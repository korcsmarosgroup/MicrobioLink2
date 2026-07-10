"""Tests for microbiolink.reverse_DMI."""

import argparse

import pytest

from microbiolink.reverse_DMI import main


# ── main / output format ──────────────────────────────────────────────────────

@pytest.fixture()
def fixture_dir(tmp_path):
    """Write minimal fixture files for a reverse DMI run."""
    fasta = tmp_path / 'bacterial.fasta'
    fasta.write_text(
        '>sp|BACT001|GENE_BACT BacterialProtein\n'
        'MSTPPKRPTTLNLFPGEVASEITQSLKKLIEQFNMARDGMKSVLSGG\n'
    )

    # ELM regex: one binding motif, one CLV motif
    elm_regex = tmp_path / 'elm_classes.tsv'
    elm_regex.write_text(
        'Accession\tELMIdentifier\tFunctionalSiteName\tDescription\tRegex\tProbability\t#Instances\t#Instances_in_PDB\n'
        'ELME000001\tLIG_SH2_STAT5\tSH2 binding\tDesc\tY..[VILF]\t0.001\t5\t0\n'
        'ELME000002\tCLV_C14_Caspase3-7\tCaspase site\tDesc\t[DSTE][^P][^DEWHFYC]D[GSAN]\t0.003\t3\t0\n'
    )

    # Motif-domain interactions: LIG_SH2_STAT5 interacts with PF00017
    motif_domain = tmp_path / 'motif_domain.tsv'
    motif_domain.write_text(
        '"ELM identifier"\t"Interaction Domain Id"\t"Interaction Domain Description"\t"Interaction Domain Name"\n'
        '"LIG_SH2_STAT5"\t"PF00017"\t"SH2"\t"Src homology 2"\n'
    )

    # Human domain file: human protein has PF00017
    human_domains = tmp_path / 'human_domains.tsv'
    human_domains.write_text(
        'Entry\tPfam\tGene Names\n'
        'P12345\tPF00017\tSTAT5A\n'
    )

    return {
        'fasta': fasta,
        'elm_regex': elm_regex,
        'motif_domain': motif_domain,
        'human_domains': human_domains,
        'output': tmp_path / 'results.csv',
    }


def test_reverse_dmi_main_output_format(fixture_dir):
    args = argparse.Namespace(
        fasta_file=str(fixture_dir['fasta']),
        elm_regex_file=str(fixture_dir['elm_regex']),
        motif_domain_file=str(fixture_dir['motif_domain']),
        human_domain_file=str(fixture_dir['human_domains']),
        resource_set='default',
        output_file=str(fixture_dir['output']),
    )
    main(args)

    content = fixture_dir['output'].read_text()
    assert '# Bacterial Protein;Motif;Start;End;Human Domain;Human Protein' in content

    data_rows = [l for l in content.splitlines() if not l.startswith('#') and l]
    if data_rows:
        parts = data_rows[0].split(';')
        assert len(parts) == 6


def test_reverse_dmi_excludes_clv_motifs(fixture_dir):
    args = argparse.Namespace(
        fasta_file=str(fixture_dir['fasta']),
        elm_regex_file=str(fixture_dir['elm_regex']),
        motif_domain_file=str(fixture_dir['motif_domain']),
        human_domain_file=str(fixture_dir['human_domains']),
        resource_set='default',
        output_file=str(fixture_dir['output']),
    )
    main(args)

    content = fixture_dir['output'].read_text()
    assert 'CLV_' not in content
