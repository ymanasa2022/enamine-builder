# Triazine library enumeration and diversity selection

This revises the supplied build_enamine.py around the nitrile-triazine core
N#CC1=NC(Cl)=NC(Cl)=N1. All paths, filters, route switches and diversity budgets
are editable in filters.json. The original upload is unchanged.

## Run

Requires Python 3.10+ and RDKit. Tested with RDKit 2026.09.1; pandas is no longer required.

```bash
python -m pip install -r requirements.txt
python build_enamine.py --config filters.json --validate-only
python build_enamine.py --config filters.json
```

Edit inputs.directory and the SDF patterns first. Patterns are filename globs,
e.g. "amines*.sdf", not the old command-line prefixes. Each input accepts a
string, an array of patterns, or null. Relative input directory and output paths
resolve relative to the JSON file. Amine input is required for enumeration;
aldehydes and all other pools can be null when using direct secondary amines.
A configured pattern that matches nothing is an error. Use null for unused inputs.
Choose a new output_dir for each run; existing enumeration results are protected.

## Reduce an existing 17-million-molecule CSV

The CSV must contain a smiles column. A synonyms or name column is optional.
No SDFs are read in this mode. Input CSV is left unchanged.

```bash
python build_enamine.py --config filters.json --select-existing /path/to/enamine_ligands.csv
```

This streams rows, standardizes structures, applies final-property filters,
deduplicates in SQLite, and writes enamine_ligands_selected.csv. It does not
regenerate the full CSV. Processing 17 million structures still requires substantial
CPU time and disk space; it does not load all fingerprints or all products into RAM.
SQLite stores every unique passing structure and its representative record.

## How filters work

Each stage has enabled, properties, required_smarts and excluded_smarts settings.
Supported descriptors: MW, logP, TPSA, HBD, HBA, RotBonds, RingCount,
AromaticRings, FormalCharge, FractionCSP3. Values are computed with RDKit.
logP is the RDKit Crippen estimate; it is not logD at pH 7.4. No pKa, logD,
efflux or BBB prediction is implied.

For each property, min/max bounds are inclusive. Set a bound to null to leave
that side open; set the whole descriptor to null or omit it to disable it.
Set a stage's enabled to false to bypass its property/SMARTS filter.
Unknown settings, misspelled descriptors, invalid bounds and invalid SMARTS fail validation.
required_smarts means ALL listed patterns must match. excluded_smarts rejects
ANY match. Filters use full precision, not rounded display values.

Example: to allow larger combined arms, edit filters.p23.properties.MW.max:

```json
"MW": {"min": null, "max": 375}
```

Stages:

| Stage | Object being evaluated |
|---|---|
| p1 | Free P1 secondary amine or alcohol |
| p2_precursor | Primary amine before building the combined arm |
| p3_precursor | Converted aldehyde, supplied ketone, or eligible SN2 halide |
| p23 | Complete free secondary amine, after arm construction or direct purchase |
| final | Actual assembled molecule; hard retention filter |
| preferred | Actual assembled molecule; annotation only, does not discard |

The precursor P2/P3 labels describe synthetic inputs. They do not establish the
bound S2/S3 orientation. Reversing the amine/carbonyl input roles may generate the
same arm; exact arm deduplication handles that.

The shipped settings retain final MW <=525, logP <=4.5, TPSA <=110, HBD <=2 and
RotBonds <=12. The preferred annotation uses tighter values (450, 3.5, 90, 1, 8).
These are starting heuristics for SAR exploration, not validated CNS exposure rules.
No minimum MW/logP is imposed. A strong base is not automatically rejected.
Final descriptors are recomputed after attachment.

## Diversity and duplicates

Use BOTH building-block selection and final-molecule selection.

1. Normalize SDF structures using RDKit FragmentParent and Uncharger, by default.
   Both can be switched off. Permanent ions may remain charged.
2. Deduplicate each filtered precursor pool using canonical isomeric SMILES.
   Repeat this after constructing combined P2/P3 arms.
3. Select diverse building blocks with greedy MaxMin selection on 2048-bit,
   radius-2 Morgan fingerprints, including chirality. Pool limits are budgets.
4. Assemble only selected pools. Deduplicate final canonical isomeric SMILES
   globally in SQLite, including products from swapping the two equivalent
   chloride positions and from different synthetic routes.
5. For final selection, collect a deterministic uniform hash sample of unique
   passing products, then run fingerprint MaxMin on that bounded candidate set.

Default diversity settings:

| Pool/setting | Default |
|---|---:|
| P1 amines | 40 |
| P1 alcohols | 20, route disabled by default |
| P2 primary precursors | 400 |
| P3 aldehydes/ketones | 400 |
| P3 eligible halides | 200, SN2 disabled by default |
| Complete P2/P3 arms | 1500 |
| Final candidate sample | 20000 |
| Requested final picks | 500 |
| Maximum final pairwise fingerprint similarity | 0.85 |
| Random/hash seed | 42 |

With P1 alcohols disabled, 40 x 1500 is 60,000 P1/arm pairs before multi-site
variants, symmetry deduplication and final filters. Enabling 20 alcohols adds
30,000 pairs. Budgets are not exact product counts because multiple reactive
sites can generate multiple products. Precursor combination and arm construction
are still performed before the complete-arm cap.

Set an individual pool limit to null to keep all UNIQUE blocks in that pool.
Set final_pick_count to null to skip final selection. Keep final_candidate_limit
at least as large as final_pick_count. Lower final_max_similarity to demand a
more separated selection; set it to null to disable the similarity stopping rule.
The picker can return FEWER than the requested count when the similarity limit
prevents further picks. Similarity is measured only with the stated fingerprint;
it is not a guarantee of distinct binding modes or distinct scaffolds.

Final selection is approximate when unique products exceed the candidate limit.
A 20,000-compound uniform sample of 17 million can miss rare families. For higher
coverage, reduce pools before enumeration, increase the candidate budget within
available memory/CPU, or make separate, explicitly chosen chemical-family runs.
Do not claim exhaustive clustering of the 17 million compounds.
MaxMin is diversity selection, not clustering. It avoids the need to calculate
an all-by-all 17-million-compound similarity matrix.

Exact deduplication retains explicitly specified stereoisomers. Tautomers are
NOT canonicalized into one identity. Enantiomers may still be close by the
fingerprint similarity criterion. Reusing a useful block in different final
molecules is intentional; it is not a duplicate product. Salt forms and repeated
catalog entries do not multiply the retained block pool under default normalization.

## Chemistry changes and boundaries

- Removed the interactive ipdb breakpoint.
- Corrected ester reduction and nitrile reduction products in the virtual rules.
- Added optional ketone input, allowing reductive amination to generate benzylic
  monosubstitution with both ketone branches retained.
- Preserve supplied stereochemistry. Newly generated stereocenters are left
  unspecified; stereoisomer enumeration and enantiopure availability need separate review.
- Restrict SNAr to the triazine ring; a pendant aryl chloride is preserved.
- Enumerate all matching reactive sites, rather than silently taking the first.
  Multiple unprotected amines/alcohols require synthetic route review and may
  produce alternative attachment products. Protecting groups are not removed automatically.
- Exclude amide/sulfonamide nitrogens from primary/secondary attaching amines.
- Optional SN2 uses methyl/primary sp3 chlorides, bromides and iodides, excluding fluorides.
- Optional precursor conversions are OFF by default. Enabling them represents
  a proposed multistep synthetic route, not a reagent directly ready for SNAr.
- P1 alcohol substitution is an optional exploratory branch, disabled by default;
  it is not inferred to be supported by the uploaded paper's morpholine SAR.

Virtual reaction rules establish connectivity, not yields, selectivity, stability,
or synthetic feasibility. The two chlorides in the original core are equivalent;
P1/P2/P3 labels are intended design roles, not experimentally determined binding assignments.

## Outputs

- enamine_ligands.csv: unique final products (enumeration mode only).
- enamine_ligands_full.csv: representative route, properties, preferred flag and route_count.
- enamine_ligands_selected.csv: final diversity-selected subset.
- aliases_*.csv: mappings from source/route IDs to each unique block representative,
  including blocks not selected by the diversity budget.
- props_*.csv: selected pools and their descriptors.
- products.sqlite: unique products, representative JSON record and surviving route count.
- config_used.json and summary.json: settings and counts for reproducibility.

For duplicate final products, the first representative route is retained and
route_count is incremented. It counts routes after building-block deduplication;
it is not a count of every catalog entry or every possible synthetic provenance.
Alternative block source names are in aliases CSVs. They are not all copied to
each final-product record.

Parallel processing is by precursor or P1 entry, with one shared partner pool per
worker. Results are consumed by batch and stored on disk; each worker still holds
its current task's batch and partner pools. More workers use more RAM.

## Verification

```bash
python -m unittest discover -s . -v
```

Tests cover reaction connectivity, preservation of both ketone branches,
triazine-specific substitution, salt deduplication, stereoisomer retention,
filter boundaries and validation, reproducible diversity selection, final
similarity limits, serial/parallel equivalence, existing-CSV selection and empty outputs.
No real ENAMINE SDFs or 17-million-row input were supplied, so full-scale runtime
and memory use have not been benchmarked.
