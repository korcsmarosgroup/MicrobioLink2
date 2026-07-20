# Microbiolink Refactoring

## Context

The aim of this plan is to specify how the code should look and be formatted into modules 
for microbioLink. This will ensure that microbioLink is easy to understand and add new features 
in the future. This plan can be used to reformat the current code into these modules and ensure any 
excess or duplicate code is removed.
This plan should be simplistic and just decribe in plain language what each module does. There should be a seperate plan
for each module implementation after this.

Look at @README.md and https://pmc.ncbi.nlm.nih.gov/articles/PMC11787512/ for further context on what MicrobioLink is.

## Why?

Currently the MicrobioLink code is not formatted to best practices and has had development done in two places that does not align correctly. This
makes it very difficult to understand, maintain and improve. This plan highlights what functionality and modules should be within MicrobioLink, allowing
a reformatting on the code to fit this outline plan.

## Where is the code?

The ground truth microbioLink code is implemented in https://github.com/korcsmarosgroup/MicrobioLink2/tree/case-study. There may be further developments
within https://github.com/korcsmarosgroup/MicrobioLink2/tree/MicrobioLink-2.1-beta. The ground truth code in the case study branch should always be prioritised in the case
where there is duplicate code.

The implementation of the plan should all take place in the refactoring branch: https://github.com/korcsmarosgroup/MicrobioLink2/tree/refactoring

## Overall Layout

There should be core functions which contain both public and private functions. The public functions should be accessible as API calls within a package. The public functions should be enough to run microbiolink without the cli.

There should be a seperate cli folder which contains argparse and any functions which are cli specific boilerplate.

## Package Layout

This should be a top-level package 'microbiolink' with two sub-packages of /workflow which directly follows the module logic and utils which has functions that repeat across modules. For any overlapping fucntionality the flat `microbiolink/` from the case study branch takes precedent over the MicrobioLink-2.1-beta and the `microbiolink_api`. 
`pyproject.toml` with `hatchling`, console-script entry points wired through `microbiolink/cli.py`.

## Testing/CI Tooling

This is deferred entirely for a future development plan. Include only the packaging essentials (pyproject.toml, build backend, entrypoints)

## Checks against old code

Delete old code once the new code has been checked. Use `case_study_input/` → `case_study_output/` to check the outputs are the same as the old code. 
Check against the case study branch for the ground truth.
If there is new functionality display the output manually to get confirmation from me that it looks correct.

## Coding Style

All functions should be implemented in the style described in @https://github.com/TobyL98/toby_verse/blob/main/guidelines/python-coding-style.md.

## Plan for Modules

### Module 1 - Z-score Filter

- Inputs:
The input should be a gene count matrix and a user defined z-score cut off.

- Module Usage:
The module should use a z-score baed method to filter lowly expressed genes.

- Outputs:
Output is a gene count matrix with filtered genes. It should give back original count values unless it is below cut-off, in which case NaN is returned always.

### Module 2 - Membrane protein Filter

- Inputs:
For proteins, the input should be either a list of uniprot IDs, a Uniprot proteome ID, gene symbols (human only) or gene symbols from a gene count matrix (we will then need a helper function to extract the gene symbols from the gene count matrix). This should be specified as either human or microbial proteins by the user.

- Module Usage:
For human proteins, it is using OmniPath InterCell for filtering by membrane proteins.
For bacterial proteins, it is using uniprot to find secreted or membrane based proteins. These should be outer membrane or plasma membrane.

- Outputs:
Same format they provided but with an extra columns with the annotations (dictionary with uniprot ID)

- Notes
Human-side mebanre filter is from `workflow/get_human_fasta.py`. 
Bacterial-side is from the fucntion `filter_bacterial_domain_table_by_location()` in `microbiolink_api/microbiome.py` in the MicrobioLink-2.1-beta branch.


### Module 3 - Downloading Fasta

- Inputs:
For proteins, the input should be either a list of uniprot IDs, a Uniprot proteome ID, gene symbols (human only) or gene symbols from a gene count matrix (we will then need a helper function to extract the gene symbols from the gene count matrix). This should be specified as either human or microbial proteins by the user. These could be the output of module 2 if they are already filtered.
NOTE - these inputs also often come from module 2

- Module Usage:
For uniprot_ids, the fasta are downloaded from uniprot api.
For gene symbols, they are downloaded by uniprot api and the organism is always human (9606).
For uniprot proteome ID, it should get all the uniprot IDs first and then get the fastas.

- Outputs:
Stored as fasta files (one file for microbe and one for human) in a folder.


### Module 4 - Downloading the Domains

- Inputs:
Same as Downloading fasta (reuse functions where possible).

- Module Usage:
The code needs to find the Pfam domains for a uniprot ID.

- Outputs:
Pfam IDs of the domains. (Dictionary with key as Uniprot ID and value as pfam domains). 

### Module 5 - Domain-Domain Interactions (DDI)

- Inputs:
The inputs should be the outputs from module 4 for both bacteria and human domains.

- Module Usage:
This function should use the resources domine and 3did to find possible interactions between the bacterial and human domains.
The resources domine and 3did should be stored within a data folder.

- Outputs:
A pandas table of possible DDIs with the following columns: bacterial uniprot ID, bacterial pfam domain, human uniport ID,
human pfam domain.

### Module 6 - Domain-Motif Interactions (DMI)

- Inputs:
The inputs should be the domain outputs from module 4 and motif outputs from module 3. 
Ther user will need to specify whether they are doing forward DMI (bacterial domain -> host motif) ro reverse DMI (human domain - bacterial motif) or both.

- Module Usage:
Known DMIs from the Eukaryotic Linear Motif (ELM) and 3did will be used to connect bacterial domains to human motifs (forward DMI) or in reverse as explained above. The motifs are found by searching within the human or bacterial protein sequence.

- Outputs:
A table with type of DMI (reverse or forward), bacterial protein, bacterial annotation (either pfam ID or motif sequence depending on forward or reverse), human protein, human annotation (either pfam ID or motif sequence).

### Module 7 - Intrinsic Disorder Regions (IDR) Prediction

- Inputs:
The inputs should be the table output from the Module 6.

- Module Usage:
This module uses either IUPRED and AIUPRED package to filter based on the likelihood of the motif being in a disordered region, a wrapper for IUPRED and AIUPRED is found in the saez lab github here: github.com/saezlab/iupred. This is the code that shows how this done for AIUPRED: https://github.com/korcsmarosgroup/MicrobioLink2/blob/MicrobioLink-2.1-beta/microbiolink/AIUPred.py and IUPRED: github.com/korcsmarosgroup/…/blob/…/idr_prediction.py.

- Outputs:
The output should be the same table with a list filtered interactions based on the disorder/binding score cut off with three extra columns: disrordered score, binding score and the combined score. For IUPRED, the binding score may be called the anchor score.

### Module 8 - Monte Carlo Simulation

- Inputs:
The input should be the table output from module 7.

- Module Usage:
The module uses monte-carlo shuffling of motifs to try and determine if the motif at that position in the protein is likely to bind to the domain compared to any random position.

- Outputs:
The output should be the same as the outputs but with more columns that describe the score from the monte carlo shuffling simulation. These columns are likely to be: monte_carlo_hits, monte_carlo_pvalue and passes_monte_carlo 

### Module 9 - TieDie

- Inputs:
The input should be for DMIs either the table output from the monte-carlo simulation (module 8), the IDR prediction (module 7) or DMI prediction (module 7) based on the user preference. If DDIs are include it will also include module 5 output as an input. Other inputs is the differentially expressed gene list with p-values and log fold change values, this will be provided by the user.

- Module Usage:
There are three seperate steps in the TieDie process. There is a first step process which converts the inputs into the correct Inputs for the TieDie package, the second step runs the TieDie network propagation and the third step converts this into an output that the user can read and intepret. The code for the tiedie algorithm is in this package: github.com/saezlab/tiedie.

- Outputs:
The output should the output from step 3 which is the final network file and the file annotation file.

### Module 10 - Functional Enrichment Analysis

- Inputs:
The input should be the whole TieDie network from Module 9 or the list of human targets of microbial proteins

- Module Usage:
Functional enrichment can either be on the whole TieDie network or the list of human targets of microbial proteins. For TieDie network functional enrichment will be done on every node on the TieDie network.

- Outputs:
The output is a table and an enrichment plot which will be the same for each input.


## Future Work

We may add extra modules in the future. These are listed here:

- Collection of known bacterial-human interactions from data (e.g., IntAct) for additional predictions