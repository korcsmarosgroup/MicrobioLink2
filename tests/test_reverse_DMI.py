"""Tests for microbiolink.reverse_DMI."""

import argparse

import pytest

from microbiolink.reverse_DMI import (
    ReverseDomainMotifInteraction,
    filter_cleavage_motifs,
    main,
    predict_reverse_domain_motif_interactions_from_data,
)


# ── filter_cleavage_motifs ────────────────────────────────────────────────────

def test_filter_cleavage_motifs_removes_clv():
    elm = {
        'CLV_C14_Caspase3-7': '[DSTE]D',
        'CLV_PCSK_FUR_1': 'RXX[RK]',
        'LIG_SH2_STAT5': 'Y..[VILF]',
        'DOC_PP1_RVXF_1': '[RK].{0,1}[VI][^P][FW]',
    }
    result = filter_cleavage_motifs(elm)
    assert 'CLV_C14_Caspase3-7' not in result
    assert 'CLV_PCSK_FUR_1' not in result
    assert 'LIG_SH2_STAT5' in result
    assert 'DOC_PP1_RVXF_1' in result


def test_filter_cleavage_motifs_empty_input():
    assert filter_cleavage_motifs({}) == {}


def test_filter_cleavage_motifs_no_clv_entries():
    elm = {'LIG_FHA_1': 'pattern', 'MOD_PKA_1': 'pattern2'}
    result = filter_cleavage_motifs(elm)
    assert result == elm


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
        output_file=str(fixture_dir['output']),
    )
    main(args)

    content = fixture_dir['output'].read_text()
    assert 'CLV_' not in content


# ── predict_reverse_domain_motif_interactions_from_data ───────────────────────

def test_predict_reverse_dmi_from_data_returns_interactions():
    bacterial_sequences = {
        'sp|BACT001|GENE_BACT BacterialProtein': 'MCRAAAAAAAAAAAAAAAAAAAAAAAAAAA',
    }
    elm_regex = {'^M{0,1}(C).': '^M{0,1}(C).'}
    elm_regex = {'DEG_Nend_UBRbox_4': '^M{0,1}(C).'}
    motif_domains = {'DEG_Nend_UBRbox_4': ['PF02207']}
    human_domain_table = {'PF02207': ['P99999']}

    results = predict_reverse_domain_motif_interactions_from_data(
        bacterial_sequences=bacterial_sequences,
        elm_regex=elm_regex,
        motif_domains=motif_domains,
        human_domain_table=human_domain_table,
    )

    assert len(results) == 1
    interaction = results[0]
    assert isinstance(interaction, ReverseDomainMotifInteraction)
    assert interaction.bacterial_protein == 'BACT001'
    assert interaction.motif == 'DEG_Nend_UBRbox_4'
    assert interaction.human_domain == 'PF02207'
    assert interaction.human_protein == 'P99999'


def test_predict_reverse_dmi_from_data_clv_excluded():
    bacterial_sequences = {
        'sp|BACT002|CLV_BACT CleavageSiteProtein': 'GGGRRRKRGGGGGGGGGGGGGGGGGGGGG',
    }
    elm_regex_with_clv = {
        'CLV_PCSK_FUR_1': 'R.[RK]R.',
        'LIG_SH2_STAT5': 'Y..[VILF]',
    }
    filtered = filter_cleavage_motifs(elm_regex_with_clv)
    motif_domains = {'CLV_PCSK_FUR_1': ['PF00089'], 'LIG_SH2_STAT5': ['PF00017']}
    human_domain_table = {'PF00089': ['P11111'], 'PF00017': ['P22222']}

    results = predict_reverse_domain_motif_interactions_from_data(
        bacterial_sequences=bacterial_sequences,
        elm_regex=filtered,
        motif_domains=motif_domains,
        human_domain_table=human_domain_table,
    )

    assert all(r.motif != 'CLV_PCSK_FUR_1' for r in results)


def test_predict_reverse_dmi_from_data_no_match():
    bacterial_sequences = {
        'sp|BACT003|NOMATCH NomatchProtein': 'MCRAAAAAAAAAAAAAAAAAAAAAAAAAAA',
    }
    elm_regex = {'DEG_Nend_UBRbox_4': '^M{0,1}(C).'}
    motif_domains = {'DEG_Nend_UBRbox_4': ['PF02207']}
    human_domain_table = {}

    results = predict_reverse_domain_motif_interactions_from_data(
        bacterial_sequences=bacterial_sequences,
        elm_regex=elm_regex,
        motif_domains=motif_domains,
        human_domain_table=human_domain_table,
    )

    assert results == []
