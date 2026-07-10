"""Tests for microbiolink.get_protein_fasta."""

from unittest.mock import MagicMock, patch

from microbiolink.get_protein_fasta import main


SAMPLE_FASTA = '>sp|P12345|GENE_HUMAN Gene\nMSEQUENCE\n'


def test_main_writes_fasta_output(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nP12345\n')
    output_file = tmp_path / 'output.fasta'

    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.core.uniprot.requests.get', return_value=mock_response):
        result = main([
            '--id_list', str(id_file),
            '--sep', ',',
            '--id_column', '1',
            '--output', str(output_file),
        ])

    assert result == 0
    assert output_file.read_text() == SAMPLE_FASTA
