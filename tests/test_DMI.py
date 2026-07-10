"""Tests for microbiolink.DMI."""

import argparse

import pytest

import microbiolink.DMI as DMI_module
from microbiolink.DMI import main


# ── main / output format ──────────────────────────────────────────────────────

@pytest.fixture()
def fixture_dir(tmp_path):
    """Write minimal fixture files for a forward DMI run."""
    fasta = tmp_path / 'human.fasta'
    fasta.write_text(
        '>sp|P12345|GENE_HUMAN HumanProtein\n'
        'MCRAAAAAAAAAAAAAAAAAAAAAAAAAAA\n'
    )

    elm_regex = tmp_path / 'elm_classes.tsv'
    elm_regex.write_text(
        'Accession\tELMIdentifier\tFunctionalSiteName\tDescription\tRegex\tProbability\t#Instances\t#Instances_in_PDB\n'
        'ELME000001\tDEG_Nend_UBRbox_4\tN-degron\tDesc\t^M{0,1}(C).\t0.001\t5\t0\n'
    )

    motif_domain = tmp_path / 'motif_domain.tsv'
    motif_domain.write_text(
        '"ELM identifier"\t"Interaction Domain Id"\t"Interaction Domain Description"\t"Interaction Domain Name"\n'
        '"DEG_Nend_UBRbox_4"\t"PF02207"\t"Desc"\t"Name"\n'
    )

    bacterial_domains = tmp_path / 'bacterial_domains.tsv'
    bacterial_domains.write_text(
        'Protein\tDomains\n'
        'BACPROT1\tPF02207\n'
    )

    return {
        'fasta': fasta,
        'elm_regex': elm_regex,
        'motif_domain': motif_domain,
        'bacterial_domains': bacterial_domains,
        'output': tmp_path / 'results.csv',
    }


def test_dmi_main_output_format(fixture_dir):
    args = argparse.Namespace(
        fasta_file=str(fixture_dir['fasta']),
        elm_regex_file=str(fixture_dir['elm_regex']),
        motif_domain_file=str(fixture_dir['motif_domain']),
        bacterial_domain_file=str(fixture_dir['bacterial_domains']),
        resource_set='default',
        output_file=str(fixture_dir['output']),
    )
    main(args)

    content = fixture_dir['output'].read_text()
    assert '# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein' in content

    data_rows = [l for l in content.splitlines() if not l.startswith('#') and l]
    assert data_rows == ['P12345;DEG_Nend_UBRbox_4;0;3;PF02207;BACPROT1']


def test_dmi_main_uses_resource_set_when_files_omitted(
    fixture_dir,
    monkeypatch,
):
    calls = []
    real_resolve = DMI_module.resolve_dmi_resource_bundle_by_name

    def spy_resolve(resource_set):
        calls.append(resource_set)
        return real_resolve(resource_set)

    monkeypatch.setattr(DMI_module, 'resolve_dmi_resource_bundle_by_name', spy_resolve)

    args = argparse.Namespace(
        fasta_file=str(fixture_dir['fasta']),
        elm_regex_file=None,
        motif_domain_file=None,
        bacterial_domain_file=str(fixture_dir['bacterial_domains']),
        resource_set='elm',
        output_file=str(fixture_dir['output']),
    )
    main(args)

    assert calls == ['elm']
    assert fixture_dir['output'].exists()
