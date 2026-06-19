#!/usr/bin/env python

"""Console entry points for microbiolink workflow scripts."""

from __future__ import annotations

import argparse
import importlib
import sys


def _as_exit_code(result: object) -> int:
    if isinstance(result, int):
        return result
    return 0


def _missing_optional_dependency(error: ModuleNotFoundError) -> int:
    message = str(error)

    if error.name in {'iupred', 'Bio', 'torch'}:
        print(
            'Missing optional AIUPred/IDR dependency. Install it with '
            '`pip install "microbiolink[idr]"`.',
            file=sys.stderr,
        )
        return 2

    if any(name in message for name in ["'iupred'", "'Bio'", "'torch'"]):
        print(
            'Missing optional AIUPred/IDR dependency. Install it with '
            '`pip install "microbiolink[idr]"`.',
            file=sys.stderr,
        )
        return 2

    print(message, file=sys.stderr)
    return 1


def _idr_import_error(error: ImportError) -> int:
    message = str(error)

    if any(
        pattern in message
        for pattern in [
            'PyTorch',
            'torch',
            'NumPy 1.x',
            'Numpy is not available',
            'numpy is not available',
        ]
    ):
        print(
            'The optional AIUPred/IDR stack is not usable in the current '
            'environment. Install it with `pip install "microbiolink[idr]"` '
            'in a Python 3.10 or 3.11 environment.',
            file=sys.stderr,
        )
        print(
            'If you already installed the extra, ensure the environment uses '
            '`numpy<2` so PyTorch and AIUPred can load correctly.',
            file=sys.stderr,
        )
        return 2

    print(message, file=sys.stderr)
    return 1


def _import_module(module_name: str):
    return importlib.import_module(module_name, package=__package__)


def _omnipath_import_error(command_name: str, error: Exception) -> int:
    print(
        f'Failed to load `{command_name}`. '
        'This command depends on `omnipath`, which may not work reliably '
        'with Python 3.13. Use a Python 3.10 or 3.11 environment for this '
        'step.',
        file=sys.stderr,
    )
    print(str(error), file=sys.stderr)
    return 1


def aiupred() -> int:
    try:
        aiupred_module = _import_module('.AIUPred')
    except ModuleNotFoundError as error:
        return _missing_optional_dependency(error)
    except ImportError as error:
        return _idr_import_error(error)

    return _as_exit_code(aiupred_module.main())


def dmi() -> int:
    dmi_module = _import_module('.DMI')

    parser = argparse.ArgumentParser(
        description=(
            'Predict interaction between human and microbial proteins based on '
            'domain-motif interactions.'
        ),
    )
    parser.add_argument(
        '--mode',
        choices = ['forward', 'reverse', 'both'],
        default = 'forward',
        help = 'Prediction direction. Defaults to forward.',
    )
    parser.add_argument('-fasta', '--fasta_file', help='Path to human protein FASTA file')
    parser.add_argument(
        '--bacterial_fasta_file',
        help = 'Path to bacterial protein FASTA file for reverse mode',
    )
    parser.add_argument('-motif', '--elm_regex_file', help='Path to ELM regex file')
    parser.add_argument(
        '-interaction',
        '--motif_domain_file',
        help='Path to motif-domain interaction file',
    )
    parser.add_argument(
        '-domain',
        '--bacterial_domain_file',
        help='Path to bacterial protein domain file',
    )
    parser.add_argument(
        '--human_domain_file',
        help = 'Path to human protein domain file for reverse mode',
    )
    parser.add_argument('-o', '--output_file', help='Path to output file')
    args = parser.parse_args()
    dmi_module.main(args)
    return 0


def download_bacterial_proteins() -> int:
    module = _import_module('.download_bacterial_proteins')

    return _as_exit_code(module.main(sys.argv[1:]))


def download_protein_domains() -> int:
    module = _import_module('.download_protein_domains')

    return _as_exit_code(module.main(sys.argv[1:]))


def download_human_domains() -> int:
    module = _import_module('.download_human_domains')

    return _as_exit_code(module.main(sys.argv[1:]))


def ddi() -> int:
    module = _import_module('.DDI')

    module.main(module.parse_args())
    return 0


def enrichr_id_database_ranking() -> int:
    module = _import_module('.enrichr_id_database_ranking')

    module.main()
    return 0


def get_human_fasta() -> int:
    try:
        module = _import_module('.get_human_fasta')
    except Exception as error:
        return _omnipath_import_error('microbiolink-get-human-fasta', error)

    module.main()
    return 0


def get_bacterial_fasta() -> int:
    module = _import_module('.get_bacterial_fasta')

    return _as_exit_code(module.main(sys.argv[1:]))


def get_protein_fasta() -> int:
    module = _import_module('.get_protein_fasta')

    return _as_exit_code(module.main(sys.argv[1:]))


def idr_prediction() -> int:
    try:
        module = _import_module('.idr_prediction')
    except ModuleNotFoundError as error:
        return _missing_optional_dependency(error)
    except ImportError as error:
        return _idr_import_error(error)

    return _as_exit_code(module.main(sys.argv[1:]))


def idr_prediction_score() -> int:
    try:
        module = _import_module('.idr_prediction_score')
    except ModuleNotFoundError as error:
        return _missing_optional_dependency(error)
    except ImportError as error:
        return _idr_import_error(error)

    return _as_exit_code(module.main(sys.argv[1:]))


def processing_tiedie_output() -> int:
    module = _import_module('.processing_tiedie_output')

    module.main()
    return 0


def tiedie_input_processing() -> int:
    try:
        module = _import_module('.tiedie_input_processing')
    except Exception as error:
        return _omnipath_import_error(
            'microbiolink-tiedie-input-processing',
            error,
        )

    module.main()
    return 0


def z_score_filter_terminal() -> int:
    module = _import_module('.z_score_filter_terminal')

    module.cli_main()
    return 0
