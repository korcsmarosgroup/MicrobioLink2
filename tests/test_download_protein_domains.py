"""Tests for microbiolink.download_protein_domains."""

from unittest.mock import MagicMock, patch

import pytest

from microbiolink.download_protein_domains import (
    DEFAULT_UNIPROT_FIELDS,
    UNIPROT_BATCH_SIZE,
    build_uniprot_accession_query,
    build_uniprot_stream_url,
    download_protein_list_with_fields,
    download_proteome_with_fields,
    main,
)


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
    url = build_uniprot_stream_url('somequery', fields=['accession', 'length'])
    assert 'accession' in url
    assert 'length' in url
    assert 'xref_pfam' not in url


def test_download_protein_list_with_fields_calls_correct_url():
    mock_response = MagicMock()
    mock_response.text = 'Entry\tPfam\nP12345\tPF00001'

    with patch('microbiolink.download_protein_domains.requests.get', return_value=mock_response) as mock_get:
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

    with patch('microbiolink.download_protein_domains.requests.get', return_value=mock_response) as mock_get:
        result = download_proteome_with_fields('UP000000625')

    called_url = mock_get.call_args[0][0]
    assert 'proteome' in called_url
    assert 'UP000000625' in called_url
    mock_response.raise_for_status.assert_called_once()
    assert result == mock_response.text


def test_main_writes_tsv_output(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nP12345\n')
    output_file = tmp_path / 'output.tsv'

    tsv_response = 'Entry\tPfam\tGene Names\nP12345\tPF00001\tgeneA\n'
    mock_response = MagicMock()
    mock_response.text = tsv_response

    with patch('microbiolink.download_protein_domains.requests.get', return_value=mock_response):
        result = main([
            '--id_list', str(id_file),
            '--sep', ',',
            '--id_type', 'uniprot',
            '--id_column', '1',
            '--output', str(output_file),
        ])

    assert result == 0
    content = output_file.read_text()
    assert 'Entry' in content
    assert 'P12345' in content


def test_uniprot_batch_size_is_1000():
    assert UNIPROT_BATCH_SIZE == 1000
