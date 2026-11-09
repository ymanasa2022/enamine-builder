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
Source-file labels locate input files; molecules are classified by their actual functional groups after conversion, so one source can supply several pools. A configured pattern that matches nothing is an error. Use null for unused inputs.
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

## Functional-group interconversion

Conversion expansion now happens BEFORE role classification and filtering.
Starting structures are retained; enabled transformations generate additional
variants that are classified into all compatible pools. Primary and secondary
alcohols can therefore supply carbonyls for arm construction, while aldehydes
and ketones can supply alcohols for the optional P1 ether branch.

| JSON route switch | Transformation |
|---|---|
| alcohol_to_aldehyde | Primary alcohol -> aldehyde |
| aldehyde_to_alcohol | Aldehyde -> primary alcohol |
| secondary_alcohol_to_ketone | Secondary alcohol -> ketone |
| ketone_to_secondary_alcohol | Ketone -> secondary alcohol |
| ester_to_aldehyde / acid_to_aldehyde | Ester or acid -> aldehyde |
| ester_to_alcohol / acid_to_alcohol | Ester or acid -> primary alcohol |
| nitrile_to_amine / imine_to_amine | Nitrile or NH-imine -> primary amine |

The four alcohol/carbonyl switches are ENABLED in the supplied JSON. Other
conversion switches remain disabled. Enabling conversion does not enable P1
alcohol SNAr automatically: set routes.p1_alcohols to true for that branch.

conversions.max_steps defaults to 2. This allows sequential changes at multiple
functional groups; set 1 for only single-step variants or 0 for no expansion.
Each unique canonical structure is visited once per original source, preventing
alcohol/carbonyl cycles from producing infinite repeats. One shortest conversion
path is retained per source/structure. Routes returning the same structure are
not all enumerated as distinct provenance records.

conversions.max_products_per_source defaults to 100 derived structures. Exceeding
it stops with an explicit error rather than silently dropping variants. This
limit and depth limit keep expansion bounded for multifunctional molecules.
conversion_provenance.csv records original entries and their retained conversion
paths BEFORE role filters and diversity selection. Pool alias files preserve
those path names through subsequent exact deduplication.

Secondary-alcohol oxidation removes stereochemistry at the oxidized carbon.
Ketone reduction generates unspecified stereochemistry at any new center.
Other supplied stereocenters are retained. Tertiary alcohols are not oxidized by
this rule; ester and amide carbonyls are not treated as ketones. The chemistry
rules do not assess protecting-group or chemoselectivity requirements.

The precursor P2/P3 labels describe synthetic inputs. They do not establish the
bound S2/S3 orientation. Reversing the amine/carbonyl input roles may generate the
same arm; exact arm deduplication handles that.

The shipped settings retain final MW <=525, logP <=5.0, TPSA <=110, HBD <=2 and
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
- Alcohol/carbonyl interconversions are ON in the supplied JSON; other precursor
  conversions are OFF. Each conversion represents a proposed synthetic step,
  not a reagent directly ready for SNAr.
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
- conversion_provenance.csv: source entries and retained conversion paths.
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
similarity limits, serial/parallel equivalence, existing-CSV selection, empty outputs,
secondary-alcohol/ketone conversion, bounded network expansion, cycle prevention,
and routing converted alcohols into carbonyl arm construction.
No real ENAMINE SDFs or 17-million-row input were supplied, so full-scale runtime
and memory use have not been benchmarked.

## SAR alignment update

The combined-arm and final hard cLogP ceilings are now 5.0 (previously 4.5).
The combined-arm RotBonds ceiling is now 8 (previously 7): free arms for the
selective cyclohexyl/cycloheptyl references 12 and 13 have eight RDKit rotatable
bonds. Preferred final settings remain annotation-only and unchanged.
The triazine counts as one aromatic ring. Basic amines are retained, with no
pKa or logD prediction. These ceilings preserve SAR exploration; they do not
validate CNS exposure. No universal structural exclusion is added for benzylic
gem-disubstitution, and no docking-based pocket assignment is enforced.

Run the editable reference audit before enumeration:

```bash
python check_sar_filters.py --config filters.json --output sar_filter_audit.csv
```

sar_references.json contains neutral structures reconstructed from the uploaded
SI names for compounds 1, 12, 13, 18 and 19, their free combined arms, morpholine,
and selected compatible reductive-amination precursors. This is a limited
reference panel, not all paper compounds or every synthetic route. Add your
own stage/name/smiles records and run again after changing filters. The report
lists exact rejection reasons, descriptors, final preferred flags and a
logP_above_4_5 flag. Exit status is 1 if any reference is rejected.
The audit does not guarantee inclusion after finite diversity selection, and
new stereocenters are still unspecified unless supplied explicitly.

All 17 supplied reference records pass the revised hard filters. Reference 13
has RDKit cLogP about 4.58 and would fail the old final cutoff. Reference 18
has TPSA about 90.2 and fails the preferred flag, but passes the hard filter.
The existing canonical isomeric SMILES and chiral Morgan fingerprints already
preserve supplied stereochemistry; no change was needed to those algorithms.
