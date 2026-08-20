# Curated tutorial dataset

A small, clean example dataset for the MicrobioLink tutorials. It is an
**inputs-only** subset of the case study in `case_study_input/`, chosen so a
seeded run produces **non-empty** results at every stage (DDIs, DMIs, TieDie
network and enrichment). No model weights, no smart-quote filenames and no
pre-downloaded FASTA/domain bulk — protein FASTA sequences and Pfam domains
are fetched live at run time.

The biology is the published case study: *Bacteroides thetaiotaomicron*
outer-membrane-vesicle (OMV) proteins interacting with human colonic
enterocyte (BEST4) proteins in Crohn's disease.

## Files

| File | What it is | Used by |
| --- | --- | --- |
| `human_proteins.csv` | Human UniProt accessions (host receptors), one per line under a `uniprot_id` header. | Membrane filter → FASTA/domain download |
| `bacterial_proteins.csv` | Bacterial (OMV) UniProt accessions, `uniprot_id` header. | FASTA/domain download |
| `degs.csv` | Endpoint differentially-expressed genes: `Gene,avg_log2FC`. | TieDie (endpoint heats) |
| `expressed_counts.csv` | Expressed-gene count matrix: `Gene,expression` (non-zero rows only). | TieDie (expressed-gene universe) |
| `background_genes.csv` | Expressed gene symbols, one per line (no header). | Enrichment background |

## How it was subset from `case_study_input`

Everything is derived deterministically (Python `random.Random(0)`), so the
selection is reproducible.

- **`human_proteins.csv` (276 accessions).** The 76 human proteins that
  survived the case study's IDR filter (guaranteed to produce forward
  domain–motif hits) plus 200 further human proteins sampled from the case
  study's full forward host–microbe interaction table. The extra 200 mostly
  drop out at the IDR / Monte-Carlo stages, so the run shows a realistic
  funnel rather than 100 % hits.
- **`bacterial_proteins.csv` (243 accessions).** The 43 bacterial proteins
  from the case study's forward interaction table plus 200 further OMV
  proteins sampled from `case_study_input/input/bacterial_protein/OMV_proteins.csv`.
- **`degs.csv` / `expressed_counts.csv` / `background_genes.csv`.** Taken
  from the case study transcriptomics
  (`Enterocytes BEST4_degs_fc05.csv` and `colon_BEST4_enterocyte_CD.csv`),
  cleaned to the two columns each stage needs.

## Notes for maintainers

- **The bacterial membrane filter is not applied to these OMV proteins.** They
  are an experimentally isolated outer-membrane-vesicle fraction and, in
  UniProt, carry essentially no subcellular-location annotation (41 of the 43
  case-study hits have none), so the location-text filter would remove the
  real hits. OMV fractionation is itself the localization step. The membrane
  filter is demonstrated on the **human** side (via OmniPath Intercell). See
  the tutorials for how this is framed.
- The tutorials run with the disorder method **`iupred`** (deterministic, CPU,
  no model download) and Monte-Carlo **`seed=0`**, so the local core
  (DMI → Monte Carlo → TieDie) reproduces its committed outputs on rerun.
- The download and enrichment stages are **live** (UniProt, Pfam, Enrichr);
  their committed snapshots may drift with upstream releases.
