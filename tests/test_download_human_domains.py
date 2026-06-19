"""Tests for microbiolink.download_human_domains."""

from unittest.mock import MagicMock, patch

import pytest

from microbiolink.download_human_domains import (
    _extract_swissprot_ids,
    _is_expressed,
    main,
    read_expressed_genes,
)


@pytest.fixture()
def zscore_csv(tmp_path):
    """A minimal z-score filtered CSV with one expressed and one filtered gene."""
    path = tmp_path / 'zscore.csv'
    path.write_text(',sample1,sample2\nGENE1,2.5,NaN\nGENE2,NaN,NaN\nGENE3,0.0,NaN\n')
    return path


def test_is_expressed_non_nan_nonzero():
    assert _is_expressed('2.5') is True


def test_is_expressed_nan():
    assert _is_expressed('NaN') is False


def test_is_expressed_zero():
    assert _is_expressed('0.0') is False


def test_is_expressed_empty():
    assert _is_expressed('') is False


def test_read_expressed_genes_filters_nan(zscore_csv):
    genes = read_expressed_genes(zscore_csv, ',')
    assert 'GENE1' in genes
    assert 'GENE2' not in genes


def test_read_expressed_genes_filters_zero(zscore_csv):
    genes = read_expressed_genes(zscore_csv, ',')
    assert 'GENE3' not in genes


def test_extract_swissprot_ids_string_value():
    translation = {'GENE1': {'Swiss-Prot': 'P12345'}}
    result = _extract_swissprot_ids(translation)
    assert result == ['P12345']


def test_extract_swissprot_ids_list_value():
    translation = {'GENE1': {'Swiss-Prot': ['P12345', 'P67890']}}
    result = _extract_swissprot_ids(translation)
    assert 'P12345' in result
    assert 'P67890' in result


def test_extract_swissprot_ids_drops_trembl_only():
    translation = {'GENE1': {'TrEMBL': 'A0A000'}}
    result = _extract_swissprot_ids(translation)
    assert result == []


def test_main_uniprot_mode(tmp_path):
    zscore = tmp_path / 'zscore.csv'
    zscore.write_text(',sample1\nP12345,2.5\nP99999,NaN\n')
    output = tmp_path / 'domains.tsv'

    tsv_response = 'Entry\tPfam\tGene Names\nP12345\tPF00001\tgeneA\n'
    mock_response = MagicMock()
    mock_response.text = tsv_response

    with patch('microbiolink.download_human_domains.download_protein_list_with_fields', return_value=tsv_response):
        result = main([
            '--gene_expression', str(zscore),
            '--id_type', 'uniprot',
            '--sep', ',',
            '--output', str(output),
        ])

    assert result == 0
    content = output.read_text()
    assert 'Entry' in content
    assert 'P12345' in content
