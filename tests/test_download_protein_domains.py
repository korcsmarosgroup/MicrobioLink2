"""Tests for microbiolink.download_protein_domains."""

from unittest.mock import MagicMock, patch

from microbiolink.download_protein_domains import main


def test_main_writes_tsv_output(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nP12345\n')
    output_file = tmp_path / 'output.tsv'

    tsv_response = 'Entry\tPfam\tGene Names\nP12345\tPF00001\tgeneA\n'
    mock_response = MagicMock()
    mock_response.text = tsv_response

    with patch('microbiolink.core.uniprot.requests.get', return_value=mock_response):
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
