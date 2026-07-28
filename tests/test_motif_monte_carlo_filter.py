from __future__ import annotations

from pathlib import Path

import pandas as pd

from microbiolink.motif_monte_carlo_filter import run_monte_carlo_filter


def write_structure_scores(
    filename: Path,
) -> None:
    """Write a deterministic per-residue structural-score table."""

    rows = ['Human Protein\tPosition\tDisorder score\tBinding score']
    for position in range(1, 51):
        disorder = 0.9 if position in {25, 26} else 0.1
        binding = 0.9 if position in {25, 26} else 0.1
        rows.append(f'PTEST\t{position}\t{disorder}\t{binding}')

    filename.write_text('\n'.join(rows) + '\n', encoding = 'utf-8')


def test_run_monte_carlo_filter_marks_expected_motif(
    tmp_path,
) -> None:
    interaction_file = tmp_path / 'dmi_with_structure.tsv'
    structure_scores = tmp_path / 'structure_scores.tsv'
    output_file = tmp_path / 'monte_carlo.tsv'
    filtered_output = tmp_path / 'monte_carlo_kept.tsv'

    interaction_file.write_text(
        'Bacterial protein\tBacterial domain\tHuman Protein\tMotif\tStart\tEnd\tAvg IUPred score\n'
        'BAC1\tPF0001\tPTEST\tMOTIF_GOOD\t24\t26\t0.9\n'
        'BAC2\tPF0002\tPTEST\tMOTIF_BAD\t0\t2\t0.1\n',
        encoding = 'utf-8',
    )
    write_structure_scores(structure_scores)

    annotated = run_monte_carlo_filter(
        interaction_file = interaction_file,
        output_file = output_file,
        filtered_output_file = filtered_output,
        structure_scores = structure_scores,
        iterations = 2000,
        alpha = 0.05,
        min_support_fraction = 1.0,
        require_binding = False,
        seed = 7,
    )

    assert output_file.exists()
    assert filtered_output.exists()
    assert annotated.shape[0] == 2

    kept = annotated.loc[annotated['motif'] == 'MOTIF_GOOD'].iloc[0]
    rejected = annotated.loc[annotated['motif'] == 'MOTIF_BAD'].iloc[0]

    assert kept['filter_status'] == 'ok'
    assert kept['passes_monte_carlo']
    assert kept['observed_support_fraction'] == 1.0
    assert kept['monte_carlo_pvalue'] < 0.05

    assert rejected['filter_status'] == 'ok'
    assert not rejected['passes_monte_carlo']
    assert rejected['observed_support_fraction'] == 0.0
    assert rejected['monte_carlo_pvalue'] == 1.0

    filtered = pd.read_csv(filtered_output, sep = '\t')
    assert filtered['motif'].tolist() == ['MOTIF_GOOD']


def test_run_monte_carlo_filter_supports_legacy_dmi_header(
    tmp_path,
) -> None:
    interaction_file = tmp_path / 'legacy_dmi.csv'
    structure_scores = tmp_path / 'structure_scores.tsv'
    output_file = tmp_path / 'legacy_monte_carlo.tsv'

    interaction_file.write_text(
        '# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein\n'
        'PTEST;MOTIF_GOOD;24;26;PF0001;BAC1\n',
        encoding = 'utf-8',
    )
    write_structure_scores(structure_scores)

    annotated = run_monte_carlo_filter(
        interaction_file = interaction_file,
        output_file = output_file,
        structure_scores = structure_scores,
        iterations = 1000,
        alpha = 0.10,
        min_support_fraction = 1.0,
        seed = 3,
    )

    assert output_file.exists()
    assert annotated.shape[0] == 1
    assert annotated.iloc[0]['human_protein'] == 'PTEST'
    assert annotated.iloc[0]['motif'] == 'MOTIF_GOOD'
    assert annotated.iloc[0]['passes_monte_carlo']
