"""Tests for microbiolink.get_bacterial_fasta."""

from unittest.mock import MagicMock, patch

from microbiolink.get_bacterial_fasta import main


SAMPLE_FASTA = '>sp|A0A000|GENE_BACT Gene\nMBACTSEQ\n'


def test_uniprot_mode_batches_correctly(tmp_path):
    id_file = tmp_path / 'ids.csv'
    ids = '\n'.join([f'P{i:05d}' for i in range(150)])
    id_file.write_text('id\n' + ids + '\n')
    output_file = tmp_path / 'output.fasta'

    with patch('microbiolink.get_bacterial_fasta.fetch_fasta_sequences', return_value=SAMPLE_FASTA) as mock_fetch:
        main([
            '--id_list', str(id_file),
            '--sep', ',',
            '--id_type', 'uniprot',
            '--id_column', '1',
            '--output', str(output_file),
            '--batch_size', '100',
        ])

    assert mock_fetch.call_count == 2
    first_batch = mock_fetch.call_args_list[0][0][0]
    second_batch = mock_fetch.call_args_list[1][0][0]
    assert len(first_batch) == 100
    assert len(second_batch) == 50


def test_main_proteome_mode_writes_fasta(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('proteome_id\nUP000000625\n')
    output_file = tmp_path / 'output.fasta'

    mock_response = MagicMock()
    mock_response.text = SAMPLE_FASTA

    with patch('microbiolink.core.uniprot.requests.get', return_value=mock_response):
        result = main([
            '--id_list', str(id_file),
            '--sep', ',',
            '--id_type', 'UP',
            '--id_column', '1',
            '--output', str(output_file),
        ])

    assert result == 0
    assert SAMPLE_FASTA in output_file.read_text()


def test_main_uniprot_mode_writes_fasta(tmp_path):
    id_file = tmp_path / 'ids.csv'
    id_file.write_text('id\nA0A000\n')
    output_file = tmp_path / 'output.fasta'

    with patch('microbiolink.get_bacterial_fasta.fetch_fasta_sequences', return_value=SAMPLE_FASTA):
        result = main([
            '--id_list', str(id_file),
            '--sep', ',',
            '--id_type', 'uniprot',
            '--id_column', '1',
            '--output', str(output_file),
        ])

    assert result == 0
    assert SAMPLE_FASTA in output_file.read_text()
