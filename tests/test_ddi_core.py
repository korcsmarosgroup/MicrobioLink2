from __future__ import annotations

from pathlib import Path

from microbiolink.core import DDIResourceBundle
from microbiolink.core import DomainDomainInteraction
from microbiolink.core import ddi_interactions_to_dataframe
from microbiolink.core import load_default_3did_ddi_resource_bundle
from microbiolink.core import load_default_ddi_resource_bundle
from microbiolink.core import load_default_domine_all_ddi_resource_bundle
from microbiolink.core import load_default_domine_hc_ddi_resource_bundle
from microbiolink.core import merge_ddi_resource_bundles
from microbiolink.core import predict_domain_domain_interactions
from microbiolink.core import predict_domain_domain_interactions_from_data
from microbiolink.core import write_domain_domain_interactions


FIXTURES = Path(__file__).parent / 'fixtures'


def fixture_path(*parts: str) -> Path:
    return FIXTURES.joinpath(*parts)


def test_load_default_ddi_resource_bundles_are_available() -> None:
    bundle_3did = load_default_3did_ddi_resource_bundle()
    bundle_domine_hc = load_default_domine_hc_ddi_resource_bundle()
    bundle_domine_all = load_default_domine_all_ddi_resource_bundle()
    default_bundle = load_default_ddi_resource_bundle()

    assert ('PF00001', 'PF00048') in bundle_3did.pfam_pairs
    assert ('PF00001', 'PF00017') in bundle_domine_hc.pfam_pairs
    assert ('PF00001', 'PF00017') in bundle_domine_all.pfam_pairs
    assert ('PF00001', 'PF00048') in default_bundle.pfam_pairs
    assert default_bundle.name == '3did+DOMINE_all'


def test_merge_ddi_resource_bundles_preserves_pair_sources() -> None:
    left_bundle = DDIResourceBundle(
        pfam_pairs = {('PF00001', 'PF00002')},
        pair_sources = {('PF00001', 'PF00002'): ('left',)},
        name = 'left',
    )
    right_bundle = DDIResourceBundle(
        pfam_pairs = {('PF00001', 'PF00002'), ('PF00002', 'PF00003')},
        pair_sources = {
            ('PF00001', 'PF00002'): ('right',),
            ('PF00002', 'PF00003'): ('right',),
        },
        name = 'right',
    )

    merged_bundle = merge_ddi_resource_bundles(left_bundle, right_bundle)

    assert merged_bundle.pair_sources[('PF00001', 'PF00002')] == (
        'left',
        'right',
    )
    assert ('PF00002', 'PF00003') in merged_bundle.pfam_pairs


def test_predict_domain_domain_interactions_from_data_returns_records() -> None:
    bundle = DDIResourceBundle(
        pfam_pairs = {('PF02207', 'PF02207')},
        pair_sources = {('PF02207', 'PF02207'): ('custom_bundle',)},
        name = 'custom_bundle',
    )

    interactions = predict_domain_domain_interactions_from_data(
        bacterial_domains = {'PF02207': ['BACPROT1']},
        human_domains = {'PF02207': ['HUMANPROT1']},
        resource_bundle = bundle,
    )

    assert interactions == [
        DomainDomainInteraction(
            bacterial_protein = 'BACPROT1',
            human_protein = 'HUMANPROT1',
            bacterial_domain = 'PF02207',
            human_domain = 'PF02207',
            resource = 'custom_bundle',
        ),
    ]


def test_predict_domain_domain_interactions_from_files_uses_custom_resource(
    tmp_path,
) -> None:
    ddi_resource = tmp_path / 'ddi_resource.tsv'
    ddi_resource.write_text('PF02207\tPF02207\n', encoding = 'utf-8')

    interactions = predict_domain_domain_interactions(
        bacterial_domain_file = fixture_path('dmi', 'bacterial_domains.tsv'),
        human_domain_file = fixture_path('dmi', 'human_domains.tsv'),
        ddi_resource_file = ddi_resource,
    )

    assert interactions == [
        DomainDomainInteraction(
            bacterial_protein = 'BACPROT1',
            human_protein = 'HUMANPROT1',
            bacterial_domain = 'PF02207',
            human_domain = 'PF02207',
            resource = 'custom',
        ),
    ]


def test_ddi_dataframe_and_writer_are_deterministic(
    tmp_path,
) -> None:
    interactions = [
        DomainDomainInteraction(
            bacterial_protein = 'BACPROT1',
            human_protein = 'HUMANPROT1',
            bacterial_domain = 'PF02207',
            human_domain = 'PF02207',
            resource = '3did',
        ),
    ]

    frame = ddi_interactions_to_dataframe(interactions)
    output_file = tmp_path / 'ddi_output.tsv'
    write_domain_domain_interactions(interactions, output_file)

    assert frame.to_dict(orient = 'records') == [
        {
            'bacterial_protein': 'BACPROT1',
            'human_protein': 'HUMANPROT1',
            'bacterial_domain': 'PF02207',
            'human_domain': 'PF02207',
            'resource': '3did',
        },
    ]
    assert output_file.read_text(encoding = 'utf-8').splitlines() == [
        'bacterial_protein;human_protein;bacterial_domain;human_domain;resource',
        'BACPROT1;HUMANPROT1;PF02207;PF02207;3did',
    ]
