# Tutorials

Worked, end-to-end walkthroughs of MicrobioLink. They pick up where
[Get Started](../get-started.md) leaves off and put the ideas from
[Pipeline Concepts](../concepts.md) into practice on a real dataset.

The tutorials come in a 2&times;2 grid. One axis is **how you drive the
pipeline** &mdash; through the Python **API** or through the
`microbiolink-*` **command-line** scripts. The other axis is **how much
ground it covers** &mdash; a **comprehensive** run through the full set of
modules, or a **basic** run through the core happy path.

|                   | Comprehensive                                  | Basic                                     |
| ----------------- | ---------------------------------------------- | ----------------------------------------- |
| **API** (Python)  | [Tutorial 1](t1_api_full.md)                   | [Tutorial 2](t2_api_basic.md)             |
| **CLI** (scripts) | [Tutorial 3](t3_cli_full.md)                   | [Tutorial 4](t4_cli_basic.md)             |

## Which one is for you?

Pick the **row** by how you want to run MicrobioLink:

- **API** &mdash; you are working in Python (a script or a notebook) and
  want to call MicrobioLink's public functions directly, wiring the modules
  together and keeping intermediate results in memory.
- **CLI** &mdash; you want to run the pipeline from the terminal using the
  installed `microbiolink-*` console scripts, one module at a time, passing
  files between steps.

Pick the **column** by how much you want to cover:

- **Comprehensive** &mdash; the whole workflow, including the optional
  modules (for example IDR filtering, membrane filtering, network
  propagation and enrichment) and both forward and reverse DMI. Best when
  you want to understand every stage.
- **Basic** &mdash; the shortest path from inputs to a domain&ndash;motif
  interaction table, using the core modules only. Best for a first run or a
  quick result.

## The four tutorials

- **[Tutorial 1 &mdash; API, comprehensive](t1_api_full.md)** &mdash; the
  full pipeline driven from Python, for programmatic users who want every
  module and full control.
- **[Tutorial 2 &mdash; API, basic](t2_api_basic.md)** &mdash; the core
  pipeline driven from Python, for programmatic users who want the shortest
  path to a result.
- **[Tutorial 3 &mdash; CLI, comprehensive](t3_cli_full.md)** &mdash; the
  full pipeline driven from the `microbiolink-*` console scripts, for
  terminal users who want every module.
- **[Tutorial 4 &mdash; CLI, basic](t4_cli_basic.md)** &mdash; the core
  pipeline driven from the `microbiolink-*` console scripts, for terminal
  users who want a quick result.

!!! info "Content coming in Phase 2"
    The four tutorials are stubbed out so navigation and links resolve. The
    step-by-step walkthroughs are written in **Phase 2** of the
    documentation work.
