"""Tests for microbiolink.core.uniprot."""

from unittest.mock import MagicMock, patch

import pytest
import requests

from microbiolink.core.uniprot import (
    DEFAULT_UNIPROT_FIELDS,
    UNIPROT_BATCH_SIZE,
    build_uniprot_accession_query,
    build_uniprot_stream_url,
    download_protein_list,
    download_protein_list_with_fields,
    download_proteome,
    download_proteome_with_fields,
    fetch_fasta_sequences,
    fetch_proteome_fasta,
    read_ids,
)


SAMPLE_FASTA = '>sp|P12345|GENE_HUMAN Gene\nMSEQUENCE\n'


# ── read_ids ───────────────────────────────────────────────────────────────

def test_read_ids_parses_column(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nP12345\nQ67890\n')

    assert read_ids(id_file, ',', 1) == ['P12345', 'Q67890']


def test_read_ids_custom_column(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('name,id\nfoo,P12345\nbar,Q67890\n')

    assert read_ids(id_file, ',', 2) == ['P12345', 'Q67890']


def test_read_ids_column_out_of_range_raises(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nP12345\n')

    with pytest.raises(ValueError):
        read_ids(id_file, ',', 2)


def test_read_ids_handles_utf8_bom(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_bytes('id\nP12345\n'.encode('utf-8-sig'))

    assert read_ids(id_file, ',', 1) == ['P12345']


# ── build_uniprot_accession_query / build_uniprot_stream_url ────────────────

def test_build_uniprot_accession_query_single():
    result = build_uniprot_accession_query(['P12345'])
    assert '%28accession%3AP12345%29' in result
    assert result.startswith('%28')
    assert result.endswith('%29')


def test_build_uniprot_accession_query_multiple():
    result = build_uniprot_accession_query(['P12345', 'Q67890'])
    assert '%28accession%3AP12345%29' in result
    assert '%28accession%3AQ67890%29' in result
    assert '+OR+' in result
    assert result.startswith('%28')
    assert result.endswith('%29')


def test_build_uniprot_stream_url_default_fields():
    url = build_uniprot_stream_url('somequery')
    for field in DEFAULT_UNIPROT_FIELDS:
        assert field in url
    assert 'format=tsv' in url
    assert 'somequery' in url


def test_build_uniprot_stream_url_custom_fields():
    url = build_uniprot_stream_url('somequery', fields = ['accession', 'length'])
    assert 'accession' in url
    assert 'length' in url
    assert 'xref_pfam' not in url


# ── download_protein_list_with_fields / download_proteome_with_fields ──────

def test_download_protein_list_with_fields_calls_correct_url():
    mock_response = MagicMock()
    mock_response.text = 'Entry\tPfam\nP12345\tPF00001'

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response) as mock_get:
        result = download_protein_list_with_fields(['P12345'])

    mock_get.assert_called_once()
    called_url = mock_get.call_args[0][0]
    assert 'accession' in called_url
    assert 'P12345' in called_url
    mock_response.raise_for_status.assert_called_once()
    assert result == mock_response.text


def test_download_proteome_with_fields_calls_correct_url():
    mock_response = MagicMock()
    mock_response.text = 'Entry\tPfam\nA0A000\tPF00002'

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response) as mock_get:
        result = download_proteome_with_fields('UP000000625')

    called_url = mock_get.call_args[0][0]
    assert 'proteome' in called_url
    assert 'UP000000625' in called_url
    mock_response.raise_for_status.assert_called_once()
    assert result == mock_response.text


def test_download_protein_list_delegates_to_with_fields():
    mock_response = MagicMock()
    mock_response.text = 'Entry\tPfam\nP12345\tPF00001'

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response):
        result = download_protein_list(['P12345'])

    assert result == mock_response.text


def test_download_proteome_delegates_to_with_fields():
    mock_response = MagicMock()
    mock_response.text = 'Entry\tPfam\nA0A000\tPF00002'

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response):
        result = download_proteome('UP000000625')

    assert result == mock_response.text


def test_uniprot_batch_size_is_1000():
    assert UNIPROT_BATCH_SIZE == 1000


# ── fetch_fasta_sequences ────────────────────────────────────────────────

def test_fetch_fasta_sequences_returns_fasta():
    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response):
        result = fetch_fasta_sequences(['P12345'])

    assert result == SAMPLE_FASTA
    mock_response.raise_for_status.assert_called_once()


def test_fetch_fasta_sequences_raises_on_http_error():
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError('500')

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response):
        with pytest.raises(requests.HTTPError):
            fetch_fasta_sequences(['P12345'])


def test_fetch_fasta_sequences_url_contains_accessions():
    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response) as mock_get:
        fetch_fasta_sequences(['P12345', 'Q67890'])

    called_url = mock_get.call_args[0][0]
    assert 'P12345' in called_url
    assert 'Q67890' in called_url
    assert 'fasta' in called_url


# ── fetch_proteome_fasta ─────────────────────────────────────────────────

def test_fetch_proteome_fasta_builds_correct_url():
    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.core.uniprot.requests.get', return_value = mock_response) as mock_get:
        fetch_proteome_fasta('UP000000625')

    called_url = mock_get.call_args[0][0]
    assert 'proteome' in called_url
    assert 'UP000000625' in called_url
    assert 'fasta' in called_url
    mock_response.raise_for_status.assert_called_once()
