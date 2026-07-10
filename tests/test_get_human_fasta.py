"""Tests for microbiolink.get_human_fasta."""

from unittest.mock import patch

from microbiolink.get_human_fasta import main


def test_main_writes_fasta_output(tmp_path):
    gene_expression = tmp_path / 'expr.csv'
    gene_expression.write_text('header\nP12345,2.5\n')

    with (
        patch('microbiolink.get_human_fasta.get_proteins', return_value = ['P12345']),
        patch(
            'microbiolink.get_human_fasta.fetch_protein_sequences',
            return_value = ['>sp|P12345|GENE_HUMAN Gene\nMSEQUENCE\n'],
        ),
    ):
        result = main([
            '-genes', str(gene_expression),
            '-id', 'uniprot',
            '-s', ',',
            '-of', str(tmp_path),
            '-oseq', 'output.fasta',
        ])

    assert result == 0
    output_file = tmp_path / 'output.fasta'
    assert output_file.read_text() == '>sp|P12345|GENE_HUMAN Gene\nMSEQUENCE\n'
