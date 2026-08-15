# Documentation for MicrobioLink
This document is a human written plan of the basic information that is required for the documemtation of MicrobioLink2. 
For information on what is included in microbioLink look at the code base in @microbiolink/

## Zensical as the documentation generator

Zensical will be used as a modern static site generator to build and maintain the documentation. The documentation for 
Zensical can be found here: https://zensical.org/docs/get-started/ and the github repository here: https://github.com/zensical/zensical
Zensical should be installed using uv as per the project guidelines.
The Zensical site generated should be hosted on read the docs. Please see this link for documentation on how to do this: https://docs.readthedocs.com/platform/latest/intro/zensical.html

## Content of Documentation

The documentation should be split into four parts. There should be:

1. An introduction page that outlines what MicrobioLink does. This can use as an initial starting point
the overview and key features part of the current microbioLink README.md located: https://github.com/korcsmarosgroup/MicrobioLink2/blob/case-study/README.md
This will need to have new functionality that is included in @microbiolink. These new requirements should be outline in @plans/microbiolink_refactoring.md.
The introduction page should also have links to three pages: Get Started, Tutorials and API documentation.

2. This should be a "get started" page which details different ways that you can install the microbiolink package and any other information that you need to get started with the project

3. This should be an initial tutorials page that describes the tutorials that are available and who should use these. There are three tutorials that there will be links to seperate pages for
    - The first tutorial goes through how run microbioLink through all 10 of the modules in @microbiolink. This should be step by step in a Jupyter notebook which is then rendered as html by Zensical. It should explain what each of the 10 modules are, how it fits into the previous step and the possible inputs that are required. It should run both forward and reverse microbiolink but should explain what both of these are. For the modules which are extra it should explain that they are optional and why you may like to add them. The tutorial should visually show what the important outputs are of each step and what the data now looks like (where it is applicable to do so far a step). The start of the tutorial should have a link to the "Getting Started" page for those who have not already installed the microbiolink package. 
    - The second tutorial should be a basic microbioLink run in a jupyter notebook with the following pattern of the modules: membrane filter of human and bacterial proteins, download the fasta sequences, download the domains, predict the DMIs, use IUPRED to filter motifs in IDRs, run monte carlo simulations then run tiedie and enrichment analysis. The enrichment analysis should be done on the initial human proteins identified and on the TieDie final network. This should all only occur for forward microbiolink. The tutorial should visually show what the important outputs are of each step and what the data now looks like (where it is applicable to do so far a step). The start of the tutorial should have a link to the "Getting Started" page for those who have not already installed the microbiolink package.
    - The third tutorial will run through how to use the command line interface cli. The code for the cli is in @microbiolink/cli.py. This should be a step by step guide of how to use the cli for each module to run microbioLink. It should explain what each of the 10 modules are, how it fits into the previous step and the possible inputs that are required. It should run both forward and reverse microbiolink but should explain what both of these are. For the modules which are extra it should explain that they are optional and why you may like to add them. The tutorial should visually show what the important outputs are of each step and what the data now looks like (where it is applicable to do so far a step). The start of the tutorial should have a link to the "Getting Started" page for those who have not already installed the microbiolink package.
    - The fourth tutorial should be a basic microbiolink run through the cli with the following pattern of the modules: membrane filter of human and bacterial proteins, download the fasta sequences, download the domains, predict the DMIs, use IUPRED to filter motifs in IDRs, run monte carlo simulations then run tiedie and enrichment analysis. The enrichment analysis should be done on the initial human proteins identified and on the TieDie final network. This should all only occur for forward microbiolink. The tutorial should visually show what the important outputs are of each step and what the data now looks like (where it is applicable to do so far a step). The start of the tutorial should have a link to the "Getting Started" page for those who have not already installed the microbiolink package.

4. API Documentation. The API documentation should use the mkdocstrings Zensical add on to generate this from the docstrings of the public functions. This should be sense checked to ensure that the wording makes sense. It would be good to divide the API documentation up into the modules (e.g., Z-score filter, membrane filter etc) and "essential" functions (i.e., the functions you would need to use to run a module) for running the pipeline versus "helper" functions (modules that are only called by the essential functions -- the user does not need to touch)


## Future Work

In the future we should add a contribution guideline as a new page to microbioLink so that other people are able to contribute to the project