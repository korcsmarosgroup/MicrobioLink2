# Tutorials

Worked, end-to-end walkthroughs of MicrobioLink. They pick up where
[Get Started](../get-started.md) leaves off and put the ideas from
[Pipeline Concepts](../concepts.md) into practice on a real dataset.

There are four tutorials. Two are **Jupyter notebooks** that drive the
pipeline from Python; two are **command-line** walkthroughs that run the
same steps from the terminal with the `microbiolink-*` console scripts.
Within each pair, one is a **comprehensive** run through every module and
the other is a **default** first run for newcomers.

|                    | Comprehensive run                              | Default run                               |
| ------------------ | ---------------------------------------------- | ----------------------------------------- |
| **Notebook**       | [Tutorial 1](t1_api_full.md)                   | [Tutorial 2](t2_api_basic.md)             |
| **Command line**   | [Tutorial 3](t3_cli_full.md)                   | [Tutorial 4](t4_cli_basic.md)             |

## Which one is for you?

Pick **notebook or command line** by how you want to run MicrobioLink:

- **Notebook** &mdash; you are working in Python and want to call
  MicrobioLink's public functions directly, wiring the modules together and
  keeping intermediate results in memory (Tutorials 1 and 2).
- **Command line** &mdash; you want to run the pipeline from the terminal
  using the installed `microbiolink-*` console scripts, one module at a
  time, passing files between steps (Tutorials 3 and 4).

Pick **comprehensive or default** by how much you want to cover:

- **Comprehensive** &mdash; all ten modules, both forward and reverse DMI,
  including the optional stages (Z-score filtering, DDI, IDR filtering,
  network propagation and enrichment). Best when you want to understand
  every stage. Needs the `[idr]`, `[tiedie]` and `[enrichment]` extras
  (see [Get Started](../get-started.md#optional-extras)).
- **Default** &mdash; the forward-only path a first-time user is most likely
  to want, running the full depth of the pipeline (membrane filter, FASTA
  and domain download, DMI, IDR filter, Monte Carlo, TieDie and enrichment)
  but dropping the Z-score Filter and DDI modules. It still uses the
  `[idr]`, `[tiedie]` and `[enrichment]` extras. Best for a first run.

## The four tutorials

- **[Tutorial 1 &mdash; A comprehensive MicrobioLink run with optional
  modules (notebook)](t1_api_full.md)** &mdash; the full pipeline driven
  from a Jupyter notebook, for users who want every module, forward and
  reverse, with the optional stages.
- **[Tutorial 2 &mdash; A default MicrobioLink run, for first-time users
  (notebook)](t2_api_basic.md)** &mdash; the forward-only default run driven
  from a Jupyter notebook, for first-time users finding their feet in
  Python.
- **[Tutorial 3 &mdash; A comprehensive MicrobioLink run with optional
  modules (command line)](t3_cli_full.md)** &mdash; the full pipeline run
  from the `microbiolink-*` console scripts, for users who want every
  module, forward and reverse, with the optional stages.
- **[Tutorial 4 &mdash; A default MicrobioLink run, for first-time users
  (command line)](t4_cli_basic.md)** &mdash; the forward-only default run
  from the `microbiolink-*` console scripts, for first-time users working
  in the terminal.

!!! info "Content coming in Phase 2"
    The four tutorials are stubbed out so navigation and links resolve. The
    step-by-step walkthroughs are written in **Phase 2** of the
    documentation work.
