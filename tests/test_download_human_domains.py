"""Tests for microbiolink.download_human_domains."""

from unittest.mock import MagicMock, patch

from microbiolink.download_human_domains import main


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
