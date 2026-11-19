"""
Step 2: Dimensionality reduction, clustering, and diversity picking.

Pipeline:
  fingerprints.npy  →  PCA (auto 90% variance)  →  UMAP 2D
  →  KMeans clustering  →  pick centroid-nearest molecule per cluster

The centroid-nearest strategy gives one maximally representative molecule
per cluster, ensuring uniform chemical space coverage (true max diversity).
Use --n_per_cluster > 1 to pick multiple representatives per cluster.

Input:  fingerprints.npy        (from 1_fingerprint.py)
        filtered_ligands.csv    (from 1_fingerprint.py)

Output: <out_dir>/umap_coords.npy          — 2D UMAP embedding
        <out_dir>/clustered_ligands.csv    — all filtered ligands + cluster label
        <out_dir>/picked_ligands.csv       — final diverse selection
        <out_dir>/umap_clusters.png        — scatter plot coloured by cluster
        <out_dir>/umap_picked.png          — scatter plot highlighting picked molecules

Usage:
    python3 2_cluster_pick.py \\
        --fingerprints fp_out/fingerprints.npy \\
        --ligands      fp_out/filtered_ligands.csv \\
        --out_dir      cluster_out/

    # Custom settings:
    python3 2_cluster_pick.py \\
        --fingerprints fp_out/fingerprints.npy \\
        --ligands      fp_out/filtered_ligands.csv \\
        --out_dir      cluster_out/ \\
        --n_clusters   100  --n_per_cluster 2  --seed 42
"""
import argparse
import os
from time import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')          # headless — no display needed
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans

try:
    import umap
except ImportError:
    raise ImportError("umap-learn not installed. Run: pip install umap-learn")


# ── PCA ───────────────────────────────────────────────────────────────────────

def run_pca(fp_array, variance_threshold=0.90):
    """Reduce to the number of components explaining ≥variance_threshold."""
    print(f"\n── PCA ──────────────────────────────────────────────────────")
    t0 = time()
    pca_full = PCA()
    pca_full.fit(fp_array.astype(np.float32))
    cumvar = np.cumsum(pca_full.explained_variance_ratio_)
    n_components = int(np.argmax(cumvar >= variance_threshold)) + 1
    print(f"  {n_components} components explain ≥{100*variance_threshold:.0f}% variance")

    pca = PCA(n_components=n_components)
    pca_features = pca.fit_transform(fp_array.astype(np.float32))
    print(f"  Done in {time()-t0:.0f}s  shape={pca_features.shape}")
    return pca_features


# ── UMAP ──────────────────────────────────────────────────────────────────────

def run_umap(pca_features, n_neighbors=50, min_dist=0.0, seed=42):
    """2D UMAP embedding on PCA-reduced fingerprints."""
    print(f"\n── UMAP ─────────────────────────────────────────────────────")
    print(f"  n_neighbors={n_neighbors}  min_dist={min_dist}")
    t0 = time()
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric='euclidean',    # PCA features are continuous
        compatibility_layout=True,  # keep 0.5.x layout; silence FutureWarning in umap-learn 0.6+
        n_jobs=-1,
        random_state=seed,
    )
    coords = reducer.fit_transform(pca_features)
    print(f"  Done in {time()-t0:.0f}s  shape={coords.shape}")
    return coords


# ── Clustering ────────────────────────────────────────────────────────────────

def run_kmeans(coords, n_clusters, seed=42):
    """MiniBatch KMeans on 2D UMAP coordinates."""
    print(f"\n── KMeans clustering ────────────────────────────────────────")
    print(f"  n_clusters={n_clusters}")
    t0 = time()
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        batch_size=max(5000, n_clusters * 100),
        n_init=10,
        random_state=seed,
    )
    labels = kmeans.fit_predict(coords)
    sizes  = np.bincount(labels)
    print(f"  Done in {time()-t0:.0f}s  "
          f"cluster sizes: min={sizes.min()}, median={int(np.median(sizes))}, max={sizes.max()}")
    return labels, kmeans.cluster_centers_


# ── Picking ───────────────────────────────────────────────────────────────────

def pick_representatives(coords, labels, centers, n_per_cluster=1):
    """
    For each cluster, pick the n_per_cluster molecules closest to the centroid.
    Returns a list of integer indices into the full dataset.
    """
    picked = []
    for c, center in enumerate(centers):
        idx_in_cluster = np.where(labels == c)[0]
        dists = np.linalg.norm(coords[idx_in_cluster] - center, axis=1)
        order = np.argsort(dists)
        picked.extend(idx_in_cluster[order[:n_per_cluster]].tolist())
    return sorted(picked)


# ── Plots ─────────────────────────────────────────────────────────────────────

def plot_clusters(coords, labels, out_path, n_clusters):
    cmap = 'tab10' if n_clusters <= 10 else ('tab20' if n_clusters <= 20 else 'nipy_spectral')
    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(coords[:, 0], coords[:, 1],
                    c=labels, cmap=cmap, s=1, alpha=0.4)
    plt.colorbar(sc, ax=ax, label='Cluster')
    ax.set_title(f'UMAP — {n_clusters} KMeans clusters  (n={len(coords):,})')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out_path}")


def plot_picked(coords, labels, picked_idx, out_path, n_clusters):
    cmap = 'tab10' if n_clusters <= 10 else ('tab20' if n_clusters <= 20 else 'nipy_spectral')
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(coords[:, 0], coords[:, 1],
               c=labels, cmap=cmap, s=1, alpha=0.15)
    ax.scatter(coords[picked_idx, 0], coords[picked_idx, 1],
               c='red', s=25, zorder=3, label=f'Picked ({len(picked_idx)})')
    ax.set_title(f'UMAP — picked representatives  (n={len(picked_idx)})')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out_path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main(args):
    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load ──
    print(f"\n── Loading inputs ───────────────────────────────────────────")
    fp_array = np.load(args.fingerprints)
    df       = pd.read_csv(args.ligands)
    print(f"  Fingerprints: {fp_array.shape}")
    print(f"  Ligands:      {len(df):,}")
    assert len(fp_array) == len(df), (
        f"Row count mismatch: fingerprints={len(fp_array)}, ligands={len(df)}"
    )

    # ── PCA ──
    pca_features = run_pca(fp_array, variance_threshold=args.pca_variance)

    # ── UMAP ──
    coords = run_umap(pca_features,
                      n_neighbors=args.umap_nn,
                      min_dist=args.umap_md,
                      seed=args.seed)
    umap_path = os.path.join(args.out_dir, 'umap_coords.npy')
    np.save(umap_path, coords)
    print(f"  Saved → {umap_path}")

    # ── Cluster ──
    labels, centers = run_kmeans(coords, n_clusters=args.n_clusters, seed=args.seed)

    df['cluster'] = labels
    clustered_path = os.path.join(args.out_dir, 'clustered_ligands.csv')
    df.to_csv(clustered_path, index=False)
    print(f"  Saved → {clustered_path}")

    # ── Pick ──
    print(f"\n── Picking representatives ──────────────────────────────────")
    print(f"  Strategy: centroid-nearest  |  n_per_cluster={args.n_per_cluster}")
    picked_idx = pick_representatives(coords, labels, centers, args.n_per_cluster)
    df_picked  = df.iloc[picked_idx].reset_index(drop=True)
    picked_path = os.path.join(args.out_dir, 'picked_ligands.csv')
    df_picked.to_csv(picked_path, index=False)
    print(f"  Total picked: {len(df_picked):,}  ({args.n_clusters} clusters × {args.n_per_cluster})")
    print(f"  Saved → {picked_path}")

    # ── Plots ──
    print(f"\n── Plots ────────────────────────────────────────────────────")
    plot_clusters(coords, labels,
                  os.path.join(args.out_dir, 'umap_clusters.png'),
                  args.n_clusters)
    plot_picked(coords, labels, picked_idx,
                os.path.join(args.out_dir, 'umap_picked.png'),
                args.n_clusters)

    print(f"\n── Done ─────────────────────────────────────────────────────")
    print(f"  {len(df_picked):,} diverse ligands ready for docking → {picked_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Cluster assembled ligands and pick diverse representatives.'
    )
    # Inputs
    parser.add_argument('--fingerprints', required=True,
                        help='Path to fingerprints.npy from step 1.')
    parser.add_argument('--ligands',      required=True,
                        help='Path to filtered_ligands.csv from step 1.')
    parser.add_argument('--out_dir',      default='cluster_out',
                        help='Output directory (default: cluster_out/).')
    # PCA
    parser.add_argument('--pca_variance', type=float, default=0.90,
                        help='Fraction of variance to retain with PCA (default: 0.90).')
    # UMAP
    parser.add_argument('--umap_nn',  type=int,   default=50,
                        help='UMAP n_neighbors (default: 50).')
    parser.add_argument('--umap_md',  type=float, default=0.0,
                        help='UMAP min_dist (default: 0.0 — tighter clusters).')
    # Clustering
    parser.add_argument('--n_clusters',    type=int, default=50,
                        help='Number of KMeans clusters (default: 50).')
    parser.add_argument('--n_per_cluster', type=int, default=1,
                        help='Molecules to pick per cluster (default: 1).')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for UMAP + KMeans (default: 42).')
    main(parser.parse_args())
