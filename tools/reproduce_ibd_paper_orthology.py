#!/usr/bin/env python3

"""Reproduce the bacterial proteome download and OMA step from the IBD paper."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import date
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import time
from typing import Any
from urllib.error import HTTPError
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen


UNIPROT_PROTEOME_SEARCH_URL = 'https://rest.uniprot.org/proteomes/search'
UNIPROT_FASTA_STREAM_URL = 'https://rest.uniprot.org/uniprotkb/stream'
UNIPROT_KB_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
REQUEST_RETRIES = 4
REQUEST_BACKOFF_SECONDS = 2
FASTA_ACCESSION_PATTERN = re.compile(r'^(?:sp|tr)\|([^|]+)\|')
DOWNLOADABILITY_CHECK_LIMIT = 10


@dataclass(frozen = True)
class SpeciesRequest:
    """Describe one species request from the paper."""

    paper_species: str
    oma_species_label: str
    search_queries: tuple[str, ...]
    preferred_proteome_id: str | None = None
    download_mode_override: str | None = None
    download_query_override: str | None = None


PAPER_SPECIES: tuple[SpeciesRequest, ...] = (
    SpeciesRequest(
        paper_species = 'Ruminococcus gnavus',
        oma_species_label = 'Ruminococcus_gnavus',
        search_queries = (
            'organism_name:"Ruminococcus gnavus"',
        ),
        preferred_proteome_id = 'UP000095600',
        download_mode_override = 'ncbi_assembly_protein_faa_gz',
        download_query_override = (
            'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/406/655/'
            'GCA_001406655.1_13414_6_36/'
            'GCA_001406655.1_13414_6_36_protein.faa.gz'
        ),
    ),
    SpeciesRequest(
        paper_species = 'Clostridium symbiosum',
        oma_species_label = 'Clostridium_symbiosum',
        search_queries = (
            'organism_name:"Clostridium symbiosum"',
        ),
    ),
    SpeciesRequest(
        paper_species = 'Clostridium clostridioforme',
        oma_species_label = 'Clostridium_clostridioforme',
        search_queries = (
            'organism_name:"Clostridium clostridioforme"',
        ),
    ),
    SpeciesRequest(
        paper_species = 'Clostridium hathewayi',
        oma_species_label = 'Clostridium_hathewayi',
        search_queries = (
            'organism_name:"Clostridium hathewayi"',
            'organism_name:"Hungatella hathewayi"',
        ),
    ),
    SpeciesRequest(
        paper_species = 'Veillonella parvula',
        oma_species_label = 'Veillonella_parvula',
        search_queries = (
            'organism_name:"Veillonella parvula"',
        ),
    ),
    SpeciesRequest(
        paper_species = 'Adhesive invasive E. coli (AIEC)',
        oma_species_label = 'AIEC',
        search_queries = (
            'AIEC',
            'NRG 857C',
            'LF82',
        ),
    ),
    SpeciesRequest(
        paper_species = 'Bacteroides fragilis',
        oma_species_label = 'Bacteroides_fragilis',
        search_queries = (
            'organism_name:"Bacteroides fragilis"',
        ),
    ),
    SpeciesRequest(
        paper_species = 'Klebsiella pneumoniae',
        oma_species_label = 'Klebsiella_pneumoniae',
        search_queries = (
            'organism_name:"Klebsiella pneumoniae"',
        ),
    ),
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description = (
            'Download the eight bacterial proteomes from the Crohn\'s disease '
            'paper, stage an OMA standalone run, and optionally execute the '
            'orthology analysis with HOG inference disabled.'
        ),
    )
    parser.add_argument(
        '--output_dir',
        default = 'runs/ibd_paper_orthology',
        help = 'Directory where the reproducibility run will be written.',
    )
    parser.add_argument(
        '--oma_bundle_dir',
        help = (
            'Existing unpacked OMA standalone directory containing bin/oma and '
            'parameters.drw.'
        ),
    )
    parser.add_argument(
        '--oma_archive_path',
        help = (
            'Path to an OMA standalone .tgz archive. When set, the archive is '
            'extracted into output_dir/vendor.'
        ),
    )
    parser.add_argument(
        '--run_oma',
        action = 'store_true',
        help = 'Run OMA after staging the input directory.',
    )
    parser.add_argument(
        '--parallel_jobs',
        type = int,
        default = min(8, max(1, os.cpu_count() or 1)),
        help = 'Parallel OMA jobs for the all-vs-all phase.',
    )
    return parser.parse_args()


def ensure_directory(path: Path) -> None:
    """Create a directory if needed."""

    path.mkdir(parents = True, exist_ok = True)


def build_url(
    base_url: str,
    params: dict[str, object],
) -> str:
    """Build a URL with encoded query parameters."""

    return f'{base_url}?{urlencode(params)}'


def fetch_text(
    base_url: str,
    params: dict[str, object],
) -> str:
    """Fetch text content with small retry logic."""

    delay = REQUEST_BACKOFF_SECONDS
    url = build_url(base_url, params)

    for attempt in range(REQUEST_RETRIES):
        try:
            with urlopen(url, timeout = 120) as response:
                return response.read().decode('utf-8')
        except HTTPError as error:
            retryable = error.code in {429, 500, 502, 503, 504}

            if not retryable or attempt == REQUEST_RETRIES - 1:
                raise
        except URLError:
            if attempt == REQUEST_RETRIES - 1:
                raise

        time.sleep(delay)
        delay *= 2

    raise RuntimeError(f'Failed to fetch {url}.')


def fetch_url_bytes(url: str) -> bytes:
    """Fetch bytes from a direct URL with small retry logic."""

    delay = REQUEST_BACKOFF_SECONDS

    for attempt in range(REQUEST_RETRIES):
        try:
            with urlopen(url, timeout = 120) as response:
                return response.read()
        except HTTPError as error:
            retryable = error.code in {429, 500, 502, 503, 504}

            if not retryable or attempt == REQUEST_RETRIES - 1:
                raise
        except URLError:
            if attempt == REQUEST_RETRIES - 1:
                raise

        time.sleep(delay)
        delay *= 2

    raise RuntimeError(f'Failed to fetch {url}.')


def fetch_json(
    base_url: str,
    params: dict[str, object],
) -> dict[str, Any]:
    """Fetch a JSON response."""

    return json.loads(fetch_text(base_url, params))


def search_proteomes(query: str) -> list[dict[str, Any]]:
    """Search UniProt proteomes for one query."""

    response = fetch_json(
        UNIPROT_PROTEOME_SEARCH_URL,
        {
            'query': query,
            'format': 'json',
            'size': 50,
        },
    )
    return list(response.get('results', []))


def proteome_has_uniprotkb_records(proteome_id: str) -> bool:
    """Check whether a proteome currently resolves to UniProtKB entries."""

    return uniprotkb_query_has_records(f'(proteome:{proteome_id})')


def uniprotkb_query_has_records(query: str) -> bool:
    """Check whether a UniProtKB query returns at least one entry."""

    response = fetch_json(
        UNIPROT_KB_SEARCH_URL,
        {
            'query': query,
            'format': 'json',
            'size': 1,
        },
    )
    return bool(response.get('results'))


def build_candidate_sequence_queries(candidate: dict[str, Any]) -> list[tuple[str, str]]:
    """Build UniProtKB download queries for one proteome candidate."""

    queries: list[tuple[str, str]] = []
    proteome_id = str(candidate.get('id', '')).strip()

    if proteome_id:
        queries.append(
            (
                'proteome_id',
                f'(proteome:{proteome_id})',
            ),
        )

    taxonomy = candidate.get('taxonomy', {})
    organism_names = [
        str(taxonomy.get('scientificName', '')).strip(),
        *[
            str(name).strip()
            for name in taxonomy.get('synonyms', [])
        ],
    ]
    strain_value = str(candidate.get('strain', '')).strip()
    strain_values = [strain_value] if strain_value else []

    if strain_value and ' / ' in strain_value:
        strain_values.append(strain_value.split(' / ', maxsplit = 1)[0].strip())

    seen_queries: set[tuple[str, str]] = set()

    for organism_name in organism_names:
        if not organism_name:
            continue

        for strain in strain_values:
            if not strain:
                continue

            query = (
                f'organism_name:"{organism_name}" '
                f'AND strain:"{strain}"'
            )
            pair = ('organism_strain_query', query)

            if pair not in seen_queries:
                seen_queries.add(pair)
                queries.append(pair)

    return queries


def resolve_candidate_download_query(
    candidate: dict[str, Any],
) -> tuple[str, str] | None:
    """Resolve a UniProtKB query that returns sequences for one candidate."""

    for mode, query in build_candidate_sequence_queries(candidate):
        if uniprotkb_query_has_records(query):
            return mode, query

    return None


def select_proteome(
    species_request: SpeciesRequest,
) -> tuple[dict[str, Any], list[dict[str, Any]], str, dict[str, str], str, str, str]:
    """Resolve the proteome used for one paper species."""

    fallback_candidates: list[dict[str, Any]] = []
    fallback_query = ''
    fallback_downloadability: dict[str, str] = {}

    for query_index, query in enumerate(species_request.search_queries):
        candidates = search_proteomes(query)
        if not candidates:
            continue

        non_excluded = [
            candidate
            for candidate in candidates
            if candidate.get('proteomeType') != 'Excluded'
        ]
        ranked_candidates = non_excluded or candidates
        downloadability = {
            str(candidate.get('id', '')): 'not_checked'
            for candidate in ranked_candidates
        }

        if species_request.preferred_proteome_id:
            for candidate in ranked_candidates:
                proteome_id = str(candidate.get('id', ''))

                if proteome_id != species_request.preferred_proteome_id:
                    continue

                downloadability[proteome_id] = (
                    species_request.download_mode_override or 'proteome_id'
                )
                return (
                    candidate,
                    ranked_candidates,
                    query,
                    downloadability,
                    'Selected the preferred current UniProt proteome hit for '
                    'this paper species and used the linked assembly protein '
                    'FASTA because the full proteome is not streamable from '
                    'the current UniProtKB endpoint.',
                    species_request.download_mode_override or 'proteome_id',
                    species_request.download_query_override
                    or f'(proteome:{proteome_id})',
                )

        reference_candidates = [
            candidate
            for candidate in ranked_candidates
            if candidate.get('proteomeType') == 'Reference proteome'
        ]
        non_reference_candidates = [
            candidate
            for candidate in ranked_candidates
            if candidate.get('proteomeType') != 'Reference proteome'
        ]
        preferred_candidates = reference_candidates + non_reference_candidates

        for candidate in preferred_candidates[:DOWNLOADABILITY_CHECK_LIMIT]:
            proteome_id = str(candidate.get('id', ''))
            download_query = resolve_candidate_download_query(candidate)
            downloadability[proteome_id] = download_query[0] if download_query else 'no'

            if download_query:
                selected = candidate
                download_mode, selected_download_query = download_query

                if query_index == 0:
                    if reference_candidates and selected.get('proteomeType') == 'Reference proteome':
                        selection_reason = (
                            'Selected the first current UniProt reference '
                            'proteome hit that also had downloadable '
                            'UniProtKB sequence records.'
                        )
                    elif reference_candidates:
                        selection_reason = (
                            'Current UniProt reference proteome hits lacked '
                            'downloadable UniProtKB sequence records; '
                            'selected the first non-excluded downloadable '
                            'species-level hit instead.'
                        )
                    else:
                        selection_reason = (
                            'No reference proteome was available in the '
                            'current UniProt hits; selected the first '
                            'non-excluded downloadable species-level hit.'
                        )
                else:
                    selection_reason = (
                        'The earlier paper-species query did not yield a '
                        'downloadable UniProtKB proteome, so a fallback query '
                        f'was used: {query}'
                    )

                return (
                    selected,
                    ranked_candidates,
                    query,
                    downloadability,
                    selection_reason,
                    download_mode,
                    selected_download_query,
                )

        if not fallback_candidates:
            fallback_candidates = ranked_candidates
            fallback_query = query
            fallback_downloadability = downloadability

    if not fallback_candidates:
        raise ValueError(
            f'No UniProt proteomes found for {species_request.paper_species}.',
        )

    reference_candidates = [
        candidate
        for candidate in fallback_candidates
        if candidate.get('proteomeType') == 'Reference proteome'
    ]

    if reference_candidates:
        selected = reference_candidates[0]
        selection_reason = (
            'Selected the first current UniProt reference proteome hit, but '
            'no downloadable UniProtKB sequence records were found among the '
            'checked candidates.'
        )
    else:
        selected = fallback_candidates[0]
        selection_reason = (
            'No reference proteome was available in the current UniProt hits, '
            'and no downloadable UniProtKB sequence records were found among '
            'the checked candidates.'
        )

    return (
        selected,
        fallback_candidates,
        fallback_query,
        fallback_downloadability,
        selection_reason,
        'proteome_id',
        f'(proteome:{selected.get("id", "")})',
    )


def candidate_rows_for_species(
    species_request: SpeciesRequest,
    resolved_query: str,
    ranked_candidates: list[dict[str, Any]],
    downloadability: dict[str, str],
    selected: dict[str, Any],
    selection_reason: str,
    selected_download_mode: str,
    selected_download_query: str,
) -> list[dict[str, object]]:
    """Convert raw candidate results into a flat TSV-ready table."""

    rows: list[dict[str, object]] = []

    for rank, candidate in enumerate(ranked_candidates, start = 1):
        taxonomy = candidate.get('taxonomy', {})
        row = {
            'paper_species': species_request.paper_species,
            'oma_species_label': species_request.oma_species_label,
            'resolved_query': resolved_query,
            'rank': rank,
            'selected': candidate.get('id') == selected.get('id'),
            'downloadability_status': downloadability.get(
                str(candidate.get('id', '')),
                'not_checked',
            ),
            'selection_reason': selection_reason
            if candidate.get('id') == selected.get('id')
            else '',
            'proteome_id': candidate.get('id', ''),
            'download_mode': selected_download_mode
            if candidate.get('id') == selected.get('id')
            else '',
            'download_query': selected_download_query
            if candidate.get('id') == selected.get('id')
            else '',
            'proteome_type': candidate.get('proteomeType', ''),
            'organism': taxonomy.get('scientificName', ''),
            'organism_synonyms': ';'.join(taxonomy.get('synonyms', [])),
            'organism_id': taxonomy.get('taxonId', ''),
            'strain': candidate.get('strain', ''),
            'protein_count': candidate.get('proteinCount', ''),
            'modified': candidate.get('modified', ''),
        }
        rows.append(row)

    return rows


def download_fasta_for_query(query: str) -> str:
    """Download a FASTA file for one UniProtKB query."""

    return fetch_text(
        UNIPROT_FASTA_STREAM_URL,
        {
            'format': 'fasta',
            'query': query,
        },
    )


def download_fasta_from_locator(
    download_mode: str,
    download_query: str,
) -> str:
    """Download FASTA from a query or a direct compressed FASTA URL."""

    if download_mode == 'ncbi_assembly_protein_faa_gz':
        return gzip.decompress(
            fetch_url_bytes(download_query),
        ).decode('utf-8')

    return download_fasta_for_query(download_query)


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
) -> None:
    """Write dictionaries to TSV."""

    with open(path, 'w', newline = '', encoding = 'utf-8') as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames = fieldnames,
            delimiter = '\t',
            extrasaction = 'ignore',
        )
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def extract_archive_root(archive_path: Path) -> str:
    """Return the top-level directory name stored in a tar archive."""

    with tarfile.open(archive_path, 'r:gz') as tar_file:
        members = [member.name for member in tar_file.getmembers() if member.name]

    if not members:
        raise ValueError(f'Archive {archive_path} is empty.')

    top_level = members[0].split('/', maxsplit = 1)[0]

    if not top_level:
        raise ValueError(f'Could not determine archive root for {archive_path}.')

    return top_level


def ensure_oma_bundle(
    output_dir: Path,
    oma_bundle_dir: str | None,
    oma_archive_path: str | None,
) -> Path:
    """Resolve an unpacked OMA standalone directory."""

    if oma_bundle_dir:
        bundle_dir = Path(oma_bundle_dir).expanduser().resolve()

        if not (bundle_dir / 'bin' / 'oma').exists():
            raise FileNotFoundError(
                f'Could not find bin/oma inside {bundle_dir}.',
            )

        return bundle_dir

    if not oma_archive_path:
        raise ValueError(
            'Provide either --oma_bundle_dir or --oma_archive_path.',
        )

    archive_path = Path(oma_archive_path).expanduser().resolve()

    if not archive_path.exists():
        raise FileNotFoundError(f'OMA archive does not exist: {archive_path}')

    vendor_dir = output_dir / 'vendor'
    ensure_directory(vendor_dir)
    root_name = extract_archive_root(archive_path)
    bundle_dir = vendor_dir / root_name

    if not bundle_dir.exists():
        with tarfile.open(archive_path, 'r:gz') as tar_file:
            tar_file.extractall(vendor_dir)

    if not (bundle_dir / 'bin' / 'oma').exists():
        raise FileNotFoundError(
            f'Extracted OMA bundle does not contain bin/oma: {bundle_dir}',
        )

    return bundle_dir


def update_oma_parameters(
    template_path: Path,
    destination_path: Path,
) -> None:
    """Write a paper-oriented OMA parameters file."""

    parameters_text = template_path.read_text(encoding = 'utf-8')
    replacements = {
        "DoHierarchicalGroups := 'bottom-up';": 'DoHierarchicalGroups := false;',
        'DoGroupFunctionPrediction := true;': 'DoGroupFunctionPrediction := false;',
    }

    for source, target in replacements.items():
        parameters_text = parameters_text.replace(source, target)

    destination_path.write_text(parameters_text, encoding = 'utf-8')


def stage_oma_input(
    output_dir: Path,
    oma_bundle_dir: Path,
    selected_rows: list[dict[str, object]],
) -> tuple[Path, Path]:
    """Prepare the OMA working directory."""

    downloads_dir = output_dir / '01_selected_proteomes'
    oma_run_dir = output_dir / '02_oma_run'
    oma_db_dir = oma_run_dir / 'DB'

    ensure_directory(downloads_dir)
    ensure_directory(oma_db_dir)

    for row in selected_rows:
        species_label = str(row['oma_species_label'])
        proteome_id = str(row['proteome_id'])
        download_mode = str(row['download_mode'])
        download_query = str(row['download_query'])
        download_path = downloads_dir / f'{species_label}__{proteome_id}.fa'
        oma_species_path = oma_db_dir / f'{species_label}.fa'

        if not download_path.exists():
            fasta_text = download_fasta_from_locator(
                download_mode,
                download_query,
            )

            if not fasta_text.strip():
                raise ValueError(
                    f'No FASTA content was returned for proteome {proteome_id}.',
                )

            download_path.write_text(fasta_text, encoding = 'utf-8')

        shutil.copyfile(download_path, oma_species_path)

    update_oma_parameters(
        oma_bundle_dir / 'parameters.drw',
        oma_run_dir / 'parameters.drw',
    )

    return downloads_dir, oma_run_dir


def run_oma(
    oma_bundle_dir: Path,
    oma_run_dir: Path,
    parallel_jobs: int,
) -> None:
    """Execute OMA from the staged working directory."""

    command = [str(oma_bundle_dir / 'bin' / 'oma')]

    if parallel_jobs > 1:
        command.extend(['-n', str(parallel_jobs)])

    subprocess.run(
        command,
        cwd = oma_run_dir,
        check = True,
    )


def parse_member_accession(member_label: str) -> str:
    """Extract a stable accession from one OMA group member label."""

    member_token = member_label.strip()
    accession_match = FASTA_ACCESSION_PATTERN.match(member_token)

    if accession_match:
        return accession_match.group(1)

    return member_token.split(' ', maxsplit = 1)[0]


def parse_orthologous_groups(
    orthologous_groups_path: Path,
) -> list[dict[str, str]]:
    """Convert OMA OrthologousGroups.txt into a long table."""

    rows: list[dict[str, str]] = []

    with open(orthologous_groups_path, encoding = 'utf-8') as infile:
        for line in infile:
            stripped_line = line.strip()

            if not stripped_line or stripped_line.startswith('#'):
                continue

            fields = stripped_line.split('\t')
            oma_group = fields[0]
            members = fields[1:]

            for member in members:
                species, raw_label = member.split(':', maxsplit = 1)
                microbial_protein = parse_member_accession(raw_label)
                rows.append(
                    {
                        'species': species,
                        'oma_group': oma_group,
                        'microbial_protein': microbial_protein,
                        'oma_member_label': raw_label,
                    },
                )

    return rows


def main() -> int:
    """Run the paper orthology preparation workflow."""

    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    metadata_dir = output_dir / '00_metadata'
    processed_dir = output_dir / '03_processed'

    ensure_directory(output_dir)
    ensure_directory(metadata_dir)
    ensure_directory(processed_dir)

    oma_bundle_dir = ensure_oma_bundle(
        output_dir = output_dir,
        oma_bundle_dir = args.oma_bundle_dir,
        oma_archive_path = args.oma_archive_path,
    )

    candidate_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []

    for species_request in PAPER_SPECIES:
        (
            selected,
            ranked_candidates,
            resolved_query,
            downloadability,
            selection_reason,
            download_mode,
            download_query,
        ) = select_proteome(
            species_request,
        )

        candidate_rows.extend(
            candidate_rows_for_species(
                species_request = species_request,
                resolved_query = resolved_query,
                ranked_candidates = ranked_candidates,
                downloadability = downloadability,
                selected = selected,
                selection_reason = selection_reason,
                selected_download_mode = download_mode,
                selected_download_query = download_query,
            ),
        )

        taxonomy = selected.get('taxonomy', {})
        selected_rows.append(
            {
                'paper_species': species_request.paper_species,
                'oma_species_label': species_request.oma_species_label,
                'resolved_query': resolved_query,
                'selection_date': date.today().isoformat(),
                'selection_reason': selection_reason,
                'proteome_id': selected.get('id', ''),
                'download_mode': download_mode,
                'download_query': download_query,
                'proteome_type': selected.get('proteomeType', ''),
                'organism': taxonomy.get('scientificName', ''),
                'organism_synonyms': ';'.join(taxonomy.get('synonyms', [])),
                'organism_id': taxonomy.get('taxonId', ''),
                'strain': selected.get('strain', ''),
                'protein_count': selected.get('proteinCount', ''),
                'modified': selected.get('modified', ''),
            },
        )

    write_tsv(
        metadata_dir / 'candidate_proteomes.tsv',
        candidate_rows,
        [
            'paper_species',
            'oma_species_label',
            'resolved_query',
            'rank',
            'selected',
            'downloadability_status',
            'selection_reason',
            'proteome_id',
            'download_mode',
            'download_query',
            'proteome_type',
            'organism',
            'organism_synonyms',
            'organism_id',
            'strain',
            'protein_count',
            'modified',
        ],
    )
    write_tsv(
        metadata_dir / 'selected_proteomes.tsv',
        selected_rows,
        [
            'paper_species',
            'oma_species_label',
            'resolved_query',
            'selection_date',
            'selection_reason',
            'proteome_id',
            'download_mode',
            'download_query',
            'proteome_type',
            'organism',
            'organism_synonyms',
            'organism_id',
            'strain',
            'protein_count',
            'modified',
        ],
    )

    downloads_dir, oma_run_dir = stage_oma_input(
        output_dir = output_dir,
        oma_bundle_dir = oma_bundle_dir,
        selected_rows = selected_rows,
    )

    summary = {
        'output_dir': str(output_dir),
        'oma_bundle_dir': str(oma_bundle_dir),
        'downloads_dir': str(downloads_dir),
        'oma_run_dir': str(oma_run_dir),
        'selected_species_count': len(selected_rows),
        'run_oma': bool(args.run_oma),
        'parallel_jobs': args.parallel_jobs,
    }

    if args.run_oma:
        run_oma(
            oma_bundle_dir = oma_bundle_dir,
            oma_run_dir = oma_run_dir,
            parallel_jobs = args.parallel_jobs,
        )

    orthologous_groups_path = oma_run_dir / 'Output' / 'OrthologousGroups.txt'

    if orthologous_groups_path.exists():
        membership_rows = parse_orthologous_groups(orthologous_groups_path)
        write_tsv(
            processed_dir / 'oma_group_members.tsv',
            membership_rows,
            [
                'species',
                'oma_group',
                'microbial_protein',
                'oma_member_label',
            ],
        )
        summary['oma_group_count'] = len(
            {
                row['oma_group']
                for row in membership_rows
            },
        )
        summary['oma_group_member_count'] = len(membership_rows)

    (metadata_dir / 'run_summary.json').write_text(
        json.dumps(summary, indent = 2),
        encoding = 'utf-8',
    )

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
