"""Tests for microbiolink.get_protein_fasta."""

from unittest.mock import MagicMock, patch

import pytest
import requests

from microbiolink.get_protein_fasta import fetch_fasta_sequences, main


SAMPLE_FASTA = '>sp|P12345|GENE_HUMAN Gene\nMSEQUENCE\n'


def test_fetch_fasta_sequences_returns_fasta():
    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.get_protein_fasta.requests.get', return_value=mock_response):
        result = fetch_fasta_sequences(['P12345'])

    assert result == SAMPLE_FASTA
    mock_response.raise_for_status.assert_called_once()


def test_fetch_fasta_sequences_raises_on_http_error():
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.HTTPError('500')

    with patch('microbiolink.get_protein_fasta.requests.get', return_value=mock_response):
        with pytest.raises(requests.HTTPError):
            fetch_fasta_sequences(['P12345'])


def test_fetch_fasta_sequences_url_contains_accessions():
    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.get_protein_fasta.requests.get', return_value=mock_response) as mock_get:
        fetch_fasta_sequences(['P12345', 'Q67890'])

    called_url = mock_get.call_args[0][0]
    assert 'P12345' in called_url
    assert 'Q67890' in called_url
    assert 'fasta' in called_url


def test_main_writes_fasta_output(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nP12345\n')
    output_file = tmp_path / 'output.fasta'

    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.get_protein_fasta.requests.get', return_value=mock_response):
        result = main([
            '--id_list', str(id_file),
            '--sep', ',',
            '--id_column', '1',
            '--output', str(output_file),
        ])

    assert result == 0
    assert output_file.read_text() == SAMPLE_FASTA
