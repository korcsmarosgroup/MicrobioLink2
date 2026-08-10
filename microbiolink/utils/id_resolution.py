"""Resolving protein identifiers of any supported id_type to UniProt accessions."""

from . import uniprot_client


def _translate_gene_symbols_to_uniprot(gene_symbols: list[str]) -> dict[str, list[str]]:
    """Translate human gene symbols to UniProt Swiss-Prot accessions via MyGene.info.

    Args:
        gene_symbols: Gene symbols to translate.

    Returns:
        Mapping of gene symbol to its resolved Swiss-Prot UniProt
        accession(s). Symbols with no UniProt mapping are omitted.
    """

    from mygene import MyGeneInfo

    mg = MyGeneInfo()
    results = mg.querymany(
        gene_symbols,
        scopes="symbol",
        fields="uniprot",
        species="human",
        returnall=True,
    )

    translation = {}
    for entry in results["out"]:
        uniprot = entry.get("uniprot")
        if not uniprot:
            continue

        swissprot = uniprot.get("Swiss-Prot")
        if not swissprot:
            continue

        translation[entry["query"]] = (
            swissprot if isinstance(swissprot, list) else [swissprot]
        )

    return translation


def translate_uniprot_to_gene_symbols(uniprot_ids: list[str]) -> dict[str, str]:
    """Translate human UniProt accessions to gene symbols via MyGene.info.

    The reverse of _translate_gene_symbols_to_uniprot. Queries are sent in
    batches of 100 (the ground-truth batch size). Accessions with no symbol
    are omitted, so callers should fall back to the accession itself.

    Args:
        uniprot_ids: Human UniProt accessions to translate.

    Returns:
        Mapping of UniProt accession to its gene symbol.
    """

    from mygene import MyGeneInfo

    mg = MyGeneInfo()
    batch_size = 100

    translation = {}
    for start in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[start : start + batch_size]
        results = mg.querymany(
            batch,
            scopes="uniprot",
            fields="symbol",
            species="human",
            returnall=True,
        )

        for entry in results["out"]:
            symbol = entry.get("symbol")
            if symbol:
                translation[entry["query"]] = symbol

    return translation


def _resolve_proteome_uniprot_ids(proteome_ids: list[str]) -> list[str]:
    """Resolve UniProt proteome identifiers to their member accessions."""

    uniprot_ids = []
    for proteome_id in proteome_ids:
        table = uniprot_client.fetch_proteome_table(proteome_id, fields=["accession"])
        uniprot_ids.extend(table["Entry"].tolist())
    return uniprot_ids


def resolve_uniprot_ids(identifiers: list[str], id_type: str) -> list[str]:
    """Resolve protein identifiers to a flat list of UniProt accessions.

    id_type alone determines the resolution strategy: 'genesymbol' is only
    ever human (decision 3 excludes microbial gene symbols) and 'proteome'
    is only ever microbial (decision 3 excludes human proteomes), so no
    separate species argument is needed here. The CLI layer already
    restricts which id_type values are offered per species.

    Args:
        identifiers: UniProt accessions, human gene symbols, or UniProt
            proteome identifiers.
        id_type: 'uniprot', 'genesymbol', or 'proteome'.

    Returns:
        A flat list of UniProt accessions.
    """

    if id_type == "uniprot":
        return list(identifiers)

    if id_type == "genesymbol":
        translation = _translate_gene_symbols_to_uniprot(identifiers)
        return sorted(
            {uniprot_id for ids in translation.values() for uniprot_id in ids}
        )

    if id_type == "proteome":
        return _resolve_proteome_uniprot_ids(identifiers)

    raise ValueError(
        f"id_type must be 'uniprot', 'genesymbol', or 'proteome', got {id_type!r}"
    )
