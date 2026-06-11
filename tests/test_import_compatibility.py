"""Verify that names moved to download_protein_domains remain importable from their
original modules, preserving backward compatibility."""


def test_import_build_uniprot_accession_query_from_download_bacterial_proteins():
    from microbiolink.download_bacterial_proteins import build_uniprot_accession_query  # noqa: F401


def test_import_build_uniprot_stream_url_from_download_bacterial_proteins():
    from microbiolink.download_bacterial_proteins import build_uniprot_stream_url  # noqa: F401


def test_import_download_protein_list_with_fields_from_download_bacterial_proteins():
    from microbiolink.download_bacterial_proteins import download_protein_list_with_fields  # noqa: F401


def test_import_download_proteome_with_fields_from_download_bacterial_proteins():
    from microbiolink.download_bacterial_proteins import download_proteome_with_fields  # noqa: F401


def test_import_constants_from_download_bacterial_proteins():
    from microbiolink.download_bacterial_proteins import (  # noqa: F401
        DEFAULT_UNIPROT_FIELDS,
        UNIPROT_BATCH_SIZE,
    )
    from microbiolink.download_bacterial_proteins import UNIPROT_BATCH_SIZE as batch_size

    assert batch_size == 1000


def test_import_fetch_protein_sequences_from_get_human_fasta():
    from microbiolink.get_human_fasta import fetch_protein_sequences  # noqa: F401
