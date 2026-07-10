#!/usr/bin/env python

"""Predict host-microbe interactions from domain-domain matches."""

from __future__ import annotations

import argparse

from microbiolink.core.ddi import load_default_3did_ddi_resource_bundle
from microbiolink.core.ddi import load_default_ddi_resource_bundle
from microbiolink.core.ddi import load_default_domine_all_ddi_resource_bundle
from microbiolink.core.ddi import load_default_domine_hc_ddi_resource_bundle
from microbiolink.core.ddi import predict_domain_domain_interactions
from microbiolink.core.ddi import write_domain_domain_interactions


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for standalone execution."""

    parser = argparse.ArgumentParser(
        description = (
            'Predict interactions between human and microbial proteins based '
            'on domain-domain interactions.'
        ),
    )
    parser.add_argument(
        '--bacterial_domain_file',
        required = True,
        help = 'Path to the bacterial protein domain file.',
    )
    parser.add_argument(
        '--human_domain_file',
        required = True,
        help = 'Path to the human protein domain file.',
    )
    parser.add_argument(
        '--resource_set',
        choices = ['default', '3did', 'domine_hc', 'domine_all'],
        default = 'default',
        help = 'Packaged DDI resource set to use.',
    )
    parser.add_argument(
        '--ddi_resource_file',
        help = 'Optional custom Pfam-Pfam resource file.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required = True,
        help = 'Path to the output file.',
    )

    return parser.parse_args()


def _load_resource_bundle(resource_set: str):
    """Load the packaged DDI resource bundle selected by the CLI."""

    if resource_set == '3did':
        return load_default_3did_ddi_resource_bundle()

    if resource_set == 'domine_hc':
        return load_default_domine_hc_ddi_resource_bundle()

    if resource_set == 'domine_all':
        return load_default_domine_all_ddi_resource_bundle()

    return load_default_ddi_resource_bundle()


def main(args: argparse.Namespace) -> None:
    """Run the domain-domain interaction workflow."""

    resource_bundle = None
    if args.ddi_resource_file is None:
        resource_bundle = _load_resource_bundle(args.resource_set)

    interactions = predict_domain_domain_interactions(
        bacterial_domain_file = args.bacterial_domain_file,
        human_domain_file = args.human_domain_file,
        ddi_resource_file = args.ddi_resource_file,
        resource_bundle = resource_bundle,
    )
    write_domain_domain_interactions(interactions, args.output_file)


if __name__ == '__main__':
    main(parse_args())
