from __future__ import annotations

import contextlib
import io
from pathlib import Path

import pandas as pd

import microbiolink
from microbiolink import cli


FIXTURES = Path(__file__).parent / 'fixtures'


def fixture_path(*parts: str) -> Path:
    return FIXTURES.joinpath(*parts)


def test_package_exposes_metadata() -> None:
    assert microbiolink.__author__ == 'Leila Potari-Gul'
    assert microbiolink.__version__


def test_cli_dmi_end_to_end(
    monkeypatch,
    tmp_path,
) -> None:
    output_file = tmp_path / 'dmi_output.tsv'

    monkeypatch.setattr(
        cli.sys,
        'argv',
        [
            'microbiolink-dmi',
            '--fasta_file',
            str(fixture_path('dmi', 'human.fasta')),
            '--elm_regex_file',
            str(fixture_path('dmi', 'elm_regex.tsv')),
            '--motif_domain_file',
            str(fixture_path('dmi', 'motif_domain.tsv')),
            '--bacterial_domain_file',
            str(fixture_path('dmi', 'bacterial_domains.tsv')),
            '--output_file',
            str(output_file),
        ],
    )

    exit_code = cli.dmi()

    assert exit_code == 0
    assert output_file.read_text(encoding='utf-8').splitlines() == [
        '# Human Protein;Motif;Start;End;Bacterial domain;Bacteria Protein',
        'P12345;DEG_Nend_UBRbox_4;0;3;PF02207;BACPROT1',
    ]


def test_cli_reverse_dmi_end_to_end(
    monkeypatch,
    tmp_path,
) -> None:
    output_file = tmp_path / 'reverse_dmi_output.tsv'

    monkeypatch.setattr(
        cli.sys,
        'argv',
        [
            'microbiolink-reverse-dmi',
            '--fasta_file',
            str(fixture_path('dmi', 'bacterial_reverse.fasta')),
            '--elm_regex_file',
            str(fixture_path('dmi', 'elm_regex.tsv')),
            '--motif_domain_file',
            str(fixture_path('dmi', 'motif_domain.tsv')),
            '--human_domain_file',
            str(fixture_path('dmi', 'human_domains.tsv')),
            '--output_file',
            str(output_file),
        ],
    )

    exit_code = cli.reverse_dmi()

    assert exit_code == 0
    assert output_file.read_text(encoding='utf-8').splitlines() == [
        '# Bacterial Protein;Motif;Start;End;Human Domain;Human Protein',
        'BACREV1;DEG_Nend_UBRbox_4;0;3;PF02207;HUMANPROT1',
    ]


def test_cli_zscore_filter_end_to_end(
    monkeypatch,
    tmp_path,
) -> None:
    output_file = tmp_path / 'zscore_output.csv'

    monkeypatch.setattr(
        cli.sys,
        'argv',
        [
            'microbiolink-zscore-filter',
            '--input_file',
            str(fixture_path('zscore', 'input.csv')),
            '--output_file',
            str(output_file),
            '--zscore',
            '0',
        ],
    )

    exit_code = cli.z_score_filter_terminal()
    filtered = pd.read_csv(output_file, index_col=0)

    assert exit_code == 0
    assert filtered.shape == (5, 2)
    assert filtered.notna().any().all()
    assert filtered.isna().any().all()


def test_cli_download_bacterial_proteins_end_to_end(
    monkeypatch,
    tmp_path,
) -> None:
    requested_urls: list[str] = []
    output_file = tmp_path / 'download_output.tsv'
    module = cli._import_module('.download_bacterial_proteins')
    helper_module = cli._import_module('.download_protein_domains')

    class FakeResponse:
        def __init__(self, text: str) -> None:
            self.text = text

        def raise_for_status(self) -> None:
            return None

    def fake_get(url: str, timeout: int) -> FakeResponse:
        del timeout
        requested_urls.append(url)
        return FakeResponse(
            'Entry\tPfam\tGene Names\n'
            'P11111\tPF0001\tgene1\n'
            'P22222\tPF0002\tgene2\n',
        )

    del module
    monkeypatch.setattr(helper_module.requests, 'get', fake_get)
    monkeypatch.setattr(
        cli.sys,
        'argv',
        [
            'microbiolink-download-bacterial-proteins',
            '--id_list',
            str(fixture_path('download', 'id_list.tsv')),
            '--sep',
            '\t',
            '--id_type',
            'Uniprot',
            '--id_column',
            '1',
            '--output',
            str(output_file),
        ],
    )

    exit_code = cli.download_bacterial_proteins()

    assert exit_code == 0
    assert len(requested_urls) == 1
    assert 'accession%3AP11111' in requested_urls[0]
    assert 'accession%3AP22222' in requested_urls[0]
    assert output_file.read_text(encoding='utf-8').splitlines() == [
        'Entry\tPfam\tGene Names',
        'P11111\tPF0001\tgene1',
        'P22222\tPF0002\tgene2',
    ]


def test_cli_ddi_end_to_end(
    monkeypatch,
    tmp_path,
) -> None:
    output_file = tmp_path / 'ddi_output.tsv'
    ddi_resource = tmp_path / 'ddi_resource.tsv'
    ddi_resource.write_text('PF02207\tPF02207\n', encoding = 'utf-8')

    monkeypatch.setattr(
        cli.sys,
        'argv',
        [
            'microbiolink-ddi',
            '--bacterial_domain_file',
            str(fixture_path('dmi', 'bacterial_domains.tsv')),
            '--human_domain_file',
            str(fixture_path('dmi', 'human_domains.tsv')),
            '--ddi_resource_file',
            str(ddi_resource),
            '--output_file',
            str(output_file),
        ],
    )

    exit_code = cli.ddi()

    assert exit_code == 0
    assert output_file.read_text(encoding = 'utf-8').splitlines() == [
        'bacterial_protein;human_protein;bacterial_domain;human_domain;resource',
        'BACPROT1;HUMANPROT1;PF02207;PF02207;custom',
    ]


def test_cli_idr_entrypoint_reports_missing_optional_dependency(
    monkeypatch,
) -> None:
    def fake_import_module(module_name: str):
        del module_name
        raise ModuleNotFoundError("No module named 'iupred'")

    monkeypatch.setattr(cli, '_import_module', fake_import_module)

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        exit_code = cli.idr_prediction()

    assert exit_code == 2
    assert 'pip install "microbiolink[idr]"' in stderr.getvalue()


def test_cli_aiupred_entrypoint_reports_missing_biopython_dependency(
    monkeypatch,
) -> None:
    def fake_import_module(module_name: str):
        del module_name
        raise ModuleNotFoundError("No module named 'Bio'")

    monkeypatch.setattr(cli, '_import_module', fake_import_module)

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        exit_code = cli.aiupred()

    assert exit_code == 2
    assert 'pip install "microbiolink[idr]"' in stderr.getvalue()


def test_cli_aiupred_entrypoint_reports_torch_numpy_stack_error(
    monkeypatch,
) -> None:
    def fake_import_module(module_name: str):
        del module_name
        raise ImportError('AIUPred requires PyTorch, which is not installed.')

    monkeypatch.setattr(cli, '_import_module', fake_import_module)

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        exit_code = cli.aiupred()

    assert exit_code == 2
    assert 'numpy<2' in stderr.getvalue()


def test_cli_tiedie_input_processing_reports_omnipath_import_error(
    monkeypatch,
) -> None:
    def fake_import_module(module_name: str):
        del module_name
        raise KeyError('adapter')

    monkeypatch.setattr(cli, '_import_module', fake_import_module)

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        exit_code = cli.tiedie_input_processing()

    assert exit_code == 1
    assert 'microbiolink-tiedie-input-processing' in stderr.getvalue()
    assert 'Python 3.10 or 3.11 environment' in stderr.getvalue()
