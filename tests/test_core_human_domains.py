"""Tests for microbiolink.core.human_domains."""

from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import requests

from microbiolink.core.human_domains import (
    _extract_swissprot_ids,
    _is_expressed,
    fetch_protein_sequences,
    get_proteins,
    read_expressed_genes,
    translate_symbol_to_uniprot,
)


@pytest.fixture()
def zscore_csv(tmp_path):
    """A minimal z-score filtered CSV with one expressed and one filtered gene."""
    path = tmp_path / 'zscore.csv'
    path.write_text(',sample1,sample2\nGENE1,2.5,NaN\nGENE2,NaN,NaN\nGENE3,0.0,NaN\n')
    return path


# ── _is_expressed ────────────────────────────────────────────────────────

def test_is_expressed_non_nan_nonzero():
    assert _is_expressed('2.5') is True


def test_is_expressed_nan():
    assert _is_expressed('NaN') is False


def test_is_expressed_zero():
    assert _is_expressed('0.0') is False


def test_is_expressed_empty():
    assert _is_expressed('') is False


# ── read_expressed_genes ─────────────────────────────────────────────────

def test_read_expressed_genes_filters_nan(zscore_csv):
    genes = read_expressed_genes(zscore_csv, ',')
    assert 'GENE1' in genes
    assert 'GENE2' not in genes


def test_read_expressed_genes_filters_zero(zscore_csv):
    genes = read_expressed_genes(zscore_csv, ',')
    assert 'GENE3' not in genes


# ── _extract_swissprot_ids ───────────────────────────────────────────────

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


# ── translate_symbol_to_uniprot ──────────────────────────────────────────

def test_translate_symbol_to_uniprot_calls_mygene():
    mock_mg = MagicMock()
    mock_mg.querymany.return_value = {
        'out': [
            {'query': 'GENE1', 'uniprot': {'Swiss-Prot': 'P12345'}},
            {'query': 'GENE2'},
        ],
    }

    with patch('mygene.MyGeneInfo', return_value = mock_mg):
        result = translate_symbol_to_uniprot(['GENE1', 'GENE2'])

    assert result == {'GENE1': {'Swiss-Prot': 'P12345'}}
    mock_mg.querymany.assert_called_once_with(
        ['GENE1', 'GENE2'],
        scopes = 'symbol',
        fields = 'uniprot',
        species = 'human',
        returnall = True,
    )


# ── fetch_protein_sequences ──────────────────────────────────────────────

def test_fetch_protein_sequences_returns_list_of_fasta():
    with patch('microbiolink.core.human_domains.fetch_fasta_sequences', return_value = '>seq\nAAA\n'):
        result = fetch_protein_sequences(['P12345'])

    assert result == ['>seq\nAAA\n']


def test_fetch_protein_sequences_handles_http_error(capsys):
    with patch(
        'microbiolink.core.human_domains.fetch_fasta_sequences',
        side_effect = requests.HTTPError('500'),
    ):
        result = fetch_protein_sequences(['P12345'])

    assert result == []
    assert 'Warning' in capsys.readouterr().err


# ── get_proteins ─────────────────────────────────────────────────────────

def test_get_proteins_uniprot_no_location_filter(tmp_path):
    gene_file = tmp_path / 'expr.csv'
    gene_file.write_text('header\nP12345,2.5\nP99999,NaN\n')

    result = get_proteins(str(gene_file), 'uniprot', ',', None, str(tmp_path))

    assert set(result) == {'P12345'}


def test_get_proteins_genesymbol_no_location_filter(tmp_path):
    gene_file = tmp_path / 'expr.csv'
    gene_file.write_text('header\nGENE1,2.5\nGENE2,NaN\n')
    translation = {'GENE1': {'Swiss-Prot': 'P12345'}}

    with patch(
        'microbiolink.core.human_domains.translate_symbol_to_uniprot',
        return_value = translation,
    ) as mock_translate:
        result = get_proteins(str(gene_file), 'genesymbol', ',', None, str(tmp_path))

    assert set(result) == {'P12345'}
    mock_translate.assert_called_once_with(['GENE1'])


def test_get_proteins_uniprot_location_filter_writes_csv(tmp_path):
    gene_file = tmp_path / 'expr.csv'
    gene_file.write_text('header\nP12345,2.5\nP99999,NaN\n')
    fake_pmtm = pd.DataFrame({'genesymbol': ['UNRELATED'], 'uniprot': ['P12345']})

    with patch('omnipath.requests.Intercell.get', return_value = fake_pmtm):
        result = get_proteins(
            str(gene_file), 'uniprot', ',', ['plasma_membrane_transmembrane'], str(tmp_path),
        )

    assert set(result) == {'P12345'}
    location_output = tmp_path / 'location_filtered_genes.csv'
    assert location_output.exists()
    assert 'P12345' in location_output.read_text()


def test_get_proteins_genesymbol_location_filter(tmp_path):
    gene_file = tmp_path / 'expr.csv'
    gene_file.write_text('header\nGENE1,2.5\n')
    fake_pmtm = pd.DataFrame({'genesymbol': ['GENE1'], 'uniprot': ['P12345']})

    with patch('omnipath.requests.Intercell.get', return_value = fake_pmtm):
        result = get_proteins(str(gene_file), 'genesymbol', ',', ['secreted'], str(tmp_path))

    assert set(result) == {'P12345'}
