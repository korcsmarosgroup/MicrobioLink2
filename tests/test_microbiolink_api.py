from __future__ import annotations

from pathlib import Path

import pandas as pd
import requests

from microbiolink_api import DMIWorkflowResult
from microbiolink_api import BidirectionalDomainMotifInteraction
from microbiolink_api import DomainMotifInteraction
from microbiolink_api import InputFormatError
from microbiolink_api import fetch_bacterial_domain_table_from_file
from microbiolink_api import filter_count_matrix_file
from microbiolink_api import filter_bacterial_domain_table_by_location
from microbiolink_api import interactions_to_dataframe
from microbiolink_api import predict_bidirectional_domain_motif_interactions
from microbiolink_api import load_default_dmi_resource_bundle
from microbiolink_api import predict_domain_motif_interactions
from microbiolink_api import predict_reverse_domain_motif_interactions
from microbiolink_api import run_dmi_workflow


FIXTURES = Path(__file__).parent / 'fixtures'


def fixture_path(*parts: str) -> Path:
    return FIXTURES.joinpath(*parts)


def test_filter_count_matrix_file_returns_dataframe(
    tmp_path,
) -> None:
    output_file = tmp_path / 'filtered_counts.csv'

    filtered_counts = filter_count_matrix_file(
        fixture_path('api', 'gene_count_matrix.csv'),
        zscore_threshold = 0,
        output_file = output_file,
    )

    assert isinstance(filtered_counts, pd.DataFrame)
    assert filtered_counts.shape == (5, 2)
    assert output_file.exists()
    assert filtered_counts.isna().any().all()


def test_predict_domain_motif_interactions_returns_records() -> None:
    default_resources = load_default_dmi_resource_bundle()

    assert 'DEG_Nend_UBRbox_4' in default_resources.elm_regex
    assert default_resources.elm_regex['DEG_Nend_UBRbox_4'] == '^M{0,1}(C).'
    assert 'DEG_Nend_UBRbox_4' in default_resources.motif_domains
    assert 'PF02207' in default_resources.motif_domains['DEG_Nend_UBRbox_4']
    assert default_resources.motif_sources['DEG_Nend_UBRbox_4'] == 'ELM'
    assert '3DID_PCNA_C_LIG_0-0' in default_resources.elm_regex
    assert 'PF02747' in default_resources.motif_domains['3DID_PCNA_C_LIG_0-0']
    assert default_resources.motif_sources['3DID_PCNA_C_LIG_0-0'] == '3did'

    interactions = predict_domain_motif_interactions(
        fasta_file = fixture_path('dmi', 'human.fasta'),
        bacterial_domain_file = fixture_path('dmi', 'bacterial_domains.tsv'),
    )

    assert interactions == [
        DomainMotifInteraction(
            human_protein = 'P12345',
            motif = 'DEG_Nend_UBRbox_4',
            start = 0,
            end = 3,
            bacterial_domain = 'PF02207',
            bacterial_protein = 'BACPROT1',
            resource = 'ELM',
        ),
    ]

    interaction_frame = interactions_to_dataframe(interactions)
    assert interaction_frame.to_dict(orient = 'records') == [
        {
            'human_protein': 'P12345',
            'motif': 'DEG_Nend_UBRbox_4',
            'start': 0,
            'end': 3,
            'bacterial_domain': 'PF02207',
            'bacterial_protein': 'BACPROT1',
            'resource': 'ELM',
        },
    ]


def test_predict_reverse_domain_motif_interactions_returns_records() -> None:
    interactions = predict_reverse_domain_motif_interactions(
        bacterial_fasta_file = fixture_path('dmi', 'bacterial_reverse.fasta'),
        human_domain_file = fixture_path('dmi', 'human_domains.tsv'),
        elm_regex_file = fixture_path('dmi', 'elm_regex.tsv'),
        motif_domain_file = fixture_path('dmi', 'motif_domain.tsv'),
    )

    assert interactions == [
        BidirectionalDomainMotifInteraction(
            host_protein = 'HUMANPROT1',
            microbial_protein = 'BACREV1',
            motif = 'DEG_Nend_UBRbox_4',
            start = 0,
            end = 3,
            domain = 'PF02207',
            motif_protein_side = 'microbe',
            domain_protein_side = 'host',
            resource = 'custom',
        ),
    ]


def test_predict_reverse_domain_motif_interactions_filters_clv_motifs(
    tmp_path,
) -> None:
    bacterial_fasta = tmp_path / 'bacterial.fasta'
    bacterial_fasta.write_text(
        '>sp|BACREV1|Bacterial reverse motif protein\n'
        'MCRAAARRAARAA\n',
        encoding = 'utf-8',
    )
    elm_regex = tmp_path / 'elm_regex.tsv'
    elm_regex.write_text(
        'Accession\tELMIdentifier\tFunctionalSiteName\tDescription\tRegex\tProbability\t#Instances\t#Instances_in_PDB\n'
        'ELME000001\tDEG_Nend_UBRbox_4\tBinding\tDesc\t^M{0,1}(C).\t0.001\t5\t0\n'
        'ELME000002\tCLV_TEST\tCleavage\tDesc\tR..R\t0.001\t5\t0\n',
        encoding = 'utf-8',
    )
    motif_domain = tmp_path / 'motif_domain.tsv'
    motif_domain.write_text(
        '"ELM identifier"\t"Interaction Domain Id"\t"Interaction Domain Description"\t"Interaction Domain Name"\n'
        '"DEG_Nend_UBRbox_4"\t"PF02207"\t"Desc"\t"Name"\n'
        '"CLV_TEST"\t"PF99999"\t"Desc"\t"Name"\n',
        encoding = 'utf-8',
    )
    human_domains = tmp_path / 'human_domains.tsv'
    human_domains.write_text(
        'Protein\tDomains\n'
        'HUMANPROT1\tPF02207\n'
        'HUMANPROT2\tPF99999\n',
        encoding = 'utf-8',
    )

    interactions = predict_reverse_domain_motif_interactions(
        bacterial_fasta_file = bacterial_fasta,
        human_domain_file = human_domains,
        elm_regex_file = elm_regex,
        motif_domain_file = motif_domain,
    )

    assert len(interactions) == 1
    assert interactions[0].motif == 'DEG_Nend_UBRbox_4'


def test_predict_bidirectional_domain_motif_interactions_combines_modes() -> None:
    interactions = predict_bidirectional_domain_motif_interactions(
        human_fasta_file = fixture_path('dmi', 'human.fasta'),
        bacterial_domain_file = fixture_path('dmi', 'bacterial_domains.tsv'),
        bacterial_fasta_file = fixture_path('dmi', 'bacterial_reverse.fasta'),
        human_domain_file = fixture_path('dmi', 'human_domains.tsv'),
        elm_regex_file = fixture_path('dmi', 'elm_regex.tsv'),
        motif_domain_file = fixture_path('dmi', 'motif_domain.tsv'),
        mode = 'both',
    )

    assert len(interactions) == 2
    assert {interaction.motif_protein_side for interaction in interactions} == {
        'host',
        'microbe',
    }


def test_run_dmi_workflow_writes_outputs(
    tmp_path,
) -> None:
    filtered_counts_output = tmp_path / 'filtered_counts.csv'
    selected_fasta_output = tmp_path / 'selected_sequences.fasta'
    dmi_output = tmp_path / 'dmi_output.tsv'

    result = run_dmi_workflow(
        count_matrix_file = fixture_path('api', 'gene_count_matrix.csv'),
        human_fasta_file = fixture_path('dmi', 'human.fasta'),
        bacterial_domain_file = fixture_path('dmi', 'bacterial_domains.tsv'),
        zscore_threshold = -3,
        filtered_counts_output = filtered_counts_output,
        selected_fasta_output = selected_fasta_output,
        dmi_output = dmi_output,
    )

    assert isinstance(result, DMIWorkflowResult)
    assert result.filtered_counts.shape == (5, 2)
    assert len(result.selected_sequences) == 2
    assert len(result.interactions) == 1
    assert filtered_counts_output.exists()
    assert selected_fasta_output.exists()
    assert dmi_output.exists()
    assert dmi_output.read_text(encoding = 'utf-8').splitlines() == [
        'human_protein;motif;start;end;bacterial_domain;bacterial_protein;resource',
        'P12345;DEG_Nend_UBRbox_4;0;3;PF02207;BACPROT1;ELM',
    ]


def test_fetch_bacterial_domain_table_from_file(
    monkeypatch,
) -> None:
    from microbiolink_api import microbiome as microbiome_module

    def fake_download_protein_list_with_fields(
        uniprot_ids: list[str],
        fields: list[str] | None = None,
    ) -> str:
        assert uniprot_ids == ['BACPROT1']
        assert fields is None
        return (
            'Entry\tPfam\tGene Names\n'
            'BACPROT1\tPF02207;PF00002\tmock_gene\n'
        )

    monkeypatch.setattr(
        microbiome_module,
        'download_protein_list_with_fields',
        fake_download_protein_list_with_fields,
    )

    domain_table = fetch_bacterial_domain_table_from_file(
        fixture_path('api', 'bacterial_proteins.tsv'),
        id_type = 'Uniprot',
        separator = '\t',
        id_column = 1,
    )

    assert domain_table.to_dict(orient = 'records') == [
        {
            'Entry': 'BACPROT1',
            'Pfam': 'PF02207;PF00002',
            'Gene Names': 'mock_gene',
        },
    ]


def test_fetch_bacterial_domain_table_from_ids_splits_failed_large_batches(
    monkeypatch,
) -> None:
    from microbiolink_api import microbiome as microbiome_module

    batch_calls: list[list[str]] = []

    def fake_download_protein_list_with_fields(
        uniprot_ids: list[str],
        fields: list[str] | None = None,
    ) -> str:
        del fields
        batch_calls.append(uniprot_ids)

        if len(uniprot_ids) > 2:
            response = requests.Response()
            response.status_code = 400
            raise requests.HTTPError('mock 400', response = response)

        header = 'Entry\tPfam\tGene Names\n'
        rows = [
            f'{uniprot_id}\tPF02207\tgene_{uniprot_id}'
            for uniprot_id in uniprot_ids
        ]
        return header + '\n'.join(rows) + '\n'

    monkeypatch.setattr(
        microbiome_module,
        'download_protein_list_with_fields',
        fake_download_protein_list_with_fields,
    )

    domain_table = microbiome_module.fetch_bacterial_domain_table_from_ids(
        ['BAC1', 'BAC2', 'BAC3', 'BAC4'],
        id_type = 'Uniprot',
    )

    assert [len(batch) for batch in batch_calls] == [4, 2, 2]
    assert domain_table['Entry'].tolist() == ['BAC1', 'BAC2', 'BAC3', 'BAC4']


def test_filter_bacterial_domain_table_by_location() -> None:
    domain_table = pd.DataFrame(
        [
            {
                'Entry': 'BACPROT1',
                'Pfam': 'PF02207;PF00002',
                'Subcellular location [CC]': 'Secreted',
            },
            {
                'Entry': 'BACPROT2',
                'Pfam': 'PF00010',
                'Subcellular location [CC]': 'Cytoplasm',
            },
        ],
    )

    filtered = filter_bacterial_domain_table_by_location(
        domain_table,
        location_filters = ['secreted'],
    )

    assert filtered.to_dict(orient = 'records') == [
        {
            'Entry': 'BACPROT1',
            'Pfam': 'PF02207;PF00002',
            'Subcellular location [CC]': 'Secreted',
        },
    ]


def test_run_dmi_workflow_accepts_bacterial_identifier_file(
    monkeypatch,
    tmp_path,
) -> None:
    from microbiolink_api import workflows as workflow_module

    def fake_fetch_bacterial_domain_table_from_file(
        id_file,
        id_type,
        separator,
        id_column,
        include_location = False,
        output_file = None,
    ) -> pd.DataFrame:
        assert Path(id_file) == fixture_path('api', 'bacterial_proteins.tsv')
        assert id_type == 'Uniprot'
        assert separator == '\t'
        assert id_column == 1
        assert include_location is False

        frame = pd.DataFrame(
            [
                {
                    'Entry': 'BACPROT1',
                    'Pfam': 'PF02207;PF00002',
                    'Gene Names': 'mock_gene',
                },
            ],
        )
        if output_file is not None:
            frame.to_csv(output_file, sep = '\t', index = False)
        return frame

    monkeypatch.setattr(
        workflow_module,
        'fetch_bacterial_domain_table_from_file',
        fake_fetch_bacterial_domain_table_from_file,
    )

    bacterial_domain_output = tmp_path / 'bacterial_domains.tsv'
    result = run_dmi_workflow(
        count_matrix_file = fixture_path('api', 'gene_count_matrix.csv'),
        human_fasta_file = fixture_path('dmi', 'human.fasta'),
        bacterial_id_file = fixture_path('api', 'bacterial_proteins.tsv'),
        bacterial_id_type = 'Uniprot',
        bacterial_domain_output = bacterial_domain_output,
    )

    assert len(result.interactions) == 1
    assert result.bacterial_domain_table.to_dict(orient = 'records') == [
        {
            'Entry': 'BACPROT1',
            'Pfam': 'PF02207;PF00002',
            'Gene Names': 'mock_gene',
        },
    ]
    assert bacterial_domain_output.exists()


def test_run_dmi_workflow_filters_bacterial_locations(
    monkeypatch,
    tmp_path,
) -> None:
    from microbiolink_api import workflows as workflow_module

    def fake_fetch_bacterial_domain_table_from_file(
        id_file,
        id_type,
        separator,
        id_column,
        include_location = False,
        output_file = None,
    ) -> pd.DataFrame:
        assert include_location is True
        frame = pd.DataFrame(
            [
                {
                    'Entry': 'BACPROT1',
                    'Pfam': 'PF02207;PF00002',
                    'Gene Names': 'mock_gene',
                    'Subcellular location [CC]': 'Secreted',
                },
                {
                    'Entry': 'BACPROT2',
                    'Pfam': 'PF02207',
                    'Gene Names': 'mock_gene_2',
                    'Subcellular location [CC]': 'Cytoplasm',
                },
            ],
        )
        if output_file is not None:
            frame.to_csv(output_file, sep = '\t', index = False)
        return frame

    monkeypatch.setattr(
        workflow_module,
        'fetch_bacterial_domain_table_from_file',
        fake_fetch_bacterial_domain_table_from_file,
    )

    bacterial_domain_output = tmp_path / 'bacterial_domains.tsv'
    result = run_dmi_workflow(
        count_matrix_file = fixture_path('api', 'gene_count_matrix.csv'),
        human_fasta_file = fixture_path('dmi', 'human.fasta'),
        bacterial_id_file = fixture_path('api', 'bacterial_proteins.tsv'),
        bacterial_id_type = 'Uniprot',
        bacterial_location_filters = ['secreted'],
        bacterial_domain_output = bacterial_domain_output,
    )

    assert len(result.interactions) == 1
    assert result.bacterial_domain_table['Entry'].tolist() == ['BACPROT1']
    assert bacterial_domain_output.exists()


def test_run_dmi_workflow_requires_microbiome_input() -> None:
    try:
        run_dmi_workflow(
            count_matrix_file = fixture_path('api', 'gene_count_matrix.csv'),
            human_fasta_file = fixture_path('dmi', 'human.fasta'),
        )
    except InputFormatError as error:
        assert 'bacterial_domain_file' in str(error)
    else:
        raise AssertionError('Expected InputFormatError to be raised.')


def test_run_dmi_workflow_does_not_drop_fasta_on_identifier_mismatch(
    monkeypatch,
    tmp_path,
) -> None:
    from microbiolink_api import workflows as workflow_module

    gene_symbol_counts = tmp_path / 'gene_symbol_counts.csv'
    gene_symbol_counts.write_text(
        'gene_symbol,sample_a,sample_b\n'
        'TLR4,100,120\n'
        'EGFR,80,90\n',
        encoding = 'utf-8',
    )

    def fake_fetch_bacterial_domain_table_from_file(
        id_file,
        id_type,
        separator,
        id_column,
        include_location = False,
        output_file = None,
    ) -> pd.DataFrame:
        del id_file, id_type, separator, id_column, include_location, output_file
        return pd.DataFrame(
            [
                {
                    'Entry': 'BACPROT1',
                    'Pfam': 'PF02207',
                    'Gene Names': 'mock_gene',
                },
            ],
        )

    monkeypatch.setattr(
        workflow_module,
        'fetch_bacterial_domain_table_from_file',
        fake_fetch_bacterial_domain_table_from_file,
    )

    result = run_dmi_workflow(
        count_matrix_file = gene_symbol_counts,
        human_fasta_file = fixture_path('dmi', 'human.fasta'),
        bacterial_id_file = fixture_path('api', 'bacterial_proteins.tsv'),
        bacterial_id_type = 'Uniprot',
    )

    assert len(result.selected_sequences) == 2
    assert len(result.interactions) == 1
