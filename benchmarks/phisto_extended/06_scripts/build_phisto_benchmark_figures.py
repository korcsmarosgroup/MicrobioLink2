#!/usr/bin/env python3

"""Build publication-style SVG figures for the PHISTO benchmark bundle."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Iterable


FONT_FAMILY = (
    "'IBM Plex Sans','Aptos','Segoe UI','Helvetica Neue',Arial,sans-serif"
)

COLORS = {
    'ink': '#15202b',
    'muted': '#5c6b73',
    'grid': '#d9e2e8',
    'panel': '#f7fafc',
    'accent_blue': '#1f5a91',
    'accent_teal': '#1f8a8a',
    'accent_gold': '#d99a2b',
    'accent_red': '#d55d52',
    'accent_green': '#3a7d44',
    'accent_purple': '#6c63a6',
    'grey': '#94a3ad',
    'light_blue': '#9ecae1',
    'light_teal': '#9dd9d2',
    'light_gold': '#f2d4a2',
    'light_red': '#f3b6b0',
}


@dataclass(frozen = True)
class ApproachRow:
    """Represent one row from the benchmark approach summary."""

    approach: str
    resource_name: str
    benchmark_unique_pairs: int
    analyzable_true_pairs: int
    resource_rows: int | None
    unique_predicted_pairs: int
    raw_prediction_rows: int
    true_positives_recovered: int
    false_positive_pairs: int
    recall_total: float
    recall_on_analyzable_pairs: float
    precision_unique_pairs: float


class SvgCanvas:
    """Small helper for writing SVG primitives."""

    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        self.elements: list[str] = []

    def rect(
        self,
        x: float,
        y: float,
        width: float,
        height: float,
        fill: str = 'none',
        stroke: str = 'none',
        stroke_width: float = 1,
        rx: float = 0,
        opacity: float | None = None,
    ) -> None:
        """Add a rectangle."""

        extra = []

        if rx:
            extra.append(f'rx="{rx}"')

        if opacity is not None:
            extra.append(f'opacity="{opacity}"')

        self.elements.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" '
            f'height="{height:.2f}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{stroke_width}" {" ".join(extra)}/>',
        )

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        stroke: str = COLORS['ink'],
        stroke_width: float = 1.5,
        dash: str | None = None,
    ) -> None:
        """Add a line."""

        dash_attr = f' stroke-dasharray="{dash}"' if dash else ''
        self.elements.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" '
            f'y2="{y2:.2f}" stroke="{stroke}" stroke-width="{stroke_width}"'
            f'{dash_attr}/>',
        )

    def circle(
        self,
        cx: float,
        cy: float,
        radius: float,
        fill: str = 'none',
        stroke: str = COLORS['ink'],
        stroke_width: float = 1.5,
        opacity: float | None = None,
    ) -> None:
        """Add a circle."""

        opacity_attr = f' opacity="{opacity}"' if opacity is not None else ''
        self.elements.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{radius:.2f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"'
            f'{opacity_attr}/>',
        )

    def text(
        self,
        x: float,
        y: float,
        value: str,
        size: float = 16,
        fill: str = COLORS['ink'],
        weight: str = '400',
        anchor: str = 'start',
        family: str = FONT_FAMILY,
        letter_spacing: float | None = None,
    ) -> None:
        """Add one line of text."""

        spacing_attr = (
            f' letter-spacing="{letter_spacing}"'
            if letter_spacing is not None
            else ''
        )
        self.elements.append(
            f'<text x="{x:.2f}" y="{y:.2f}" font-family="{family}" '
            f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
            f'text-anchor="{anchor}"{spacing_attr}>{escape(value)}</text>',
        )

    def text_block(
        self,
        x: float,
        y: float,
        lines: Iterable[str],
        size: float = 16,
        line_height: float = 1.25,
        fill: str = COLORS['ink'],
        weight: str = '400',
        anchor: str = 'start',
        family: str = FONT_FAMILY,
    ) -> None:
        """Add a multi-line text block."""

        tspans = []

        for index, line in enumerate(lines):
            dy = 0 if index == 0 else size * line_height
            tspans.append(
                f'<tspan x="{x:.2f}" dy="{dy:.2f}">{escape(line)}</tspan>',
            )

        self.elements.append(
            f'<text x="{x:.2f}" y="{y:.2f}" font-family="{family}" '
            f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
            f'text-anchor="{anchor}">{"".join(tspans)}</text>',
        )

    def save(self, path: Path) -> None:
        """Write the SVG file to disk."""

        path.write_text(
            '\n'.join(
                [
                    '<?xml version="1.0" encoding="UTF-8"?>',
                    (
                        f'<svg xmlns="http://www.w3.org/2000/svg" '
                        f'width="{self.width}" height="{self.height}" '
                        f'viewBox="0 0 {self.width} {self.height}" '
                        f'fill="none">'
                    ),
                    f'<rect width="{self.width}" height="{self.height}" '
                    f'fill="white"/>',
                    *self.elements,
                    '</svg>',
                ],
            ),
            encoding = 'utf-8',
        )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description = 'Build publication-style SVG figures for the PHISTO '
        'benchmark bundle.',
    )
    parser.add_argument(
        '--bundle_dir',
        required = True,
        help = 'Path to the documented benchmark bundle directory.',
    )
    return parser.parse_args()


def read_key_value_tsv(path: Path) -> dict[str, str]:
    """Read a two-column key-value TSV."""

    values: dict[str, str] = {}

    with open(path, newline = '', encoding = 'utf-8') as infile:
        reader = csv.DictReader(infile, delimiter = '\t')

        for row in reader:
            key = row.get('metric') or row.get('field')
            value = row.get('value')

            if key is not None and value is not None:
                values[key] = value

    return values


def read_approach_summary(path: Path) -> list[ApproachRow]:
    """Read the main benchmark approach summary."""

    rows: list[ApproachRow] = []

    with open(path, newline = '', encoding = 'utf-8') as infile:
        reader = csv.DictReader(infile, delimiter = '\t')

        for row in reader:
            resource_rows = row['resource_rows'].strip()
            rows.append(
                ApproachRow(
                    approach = row['approach'],
                    resource_name = row['resource_name'],
                    benchmark_unique_pairs = int(row['benchmark_unique_pairs']),
                    analyzable_true_pairs = int(row['analyzable_true_pairs']),
                    resource_rows = int(resource_rows) if resource_rows else None,
                    unique_predicted_pairs = int(row['unique_predicted_pairs']),
                    raw_prediction_rows = int(row['raw_prediction_rows']),
                    true_positives_recovered = int(row['true_positives_recovered']),
                    false_positive_pairs = int(
                        row['false_positive_pairs_within_panel_cross_product'],
                    ),
                    recall_total = float(row['recall_total']),
                    recall_on_analyzable_pairs = float(
                        row['recall_on_analyzable_pairs'],
                    ),
                    precision_unique_pairs = float(row['precision_unique_pairs']),
                ),
            )

    return rows


def read_count_table(path: Path, key_name: str) -> list[tuple[str, int]]:
    """Read a simple count table."""

    rows: list[tuple[str, int]] = []

    with open(path, newline = '', encoding = 'utf-8') as infile:
        reader = csv.DictReader(infile, delimiter = '\t')

        for row in reader:
            rows.append((row[key_name], int(row['count'])))

    return rows


def load_pair_set(
    path: Path,
    resource_name: str | None = None,
) -> set[tuple[str, str]]:
    """Load unique bacterial-human pairs from a detail table."""

    pairs: set[tuple[str, str]] = set()

    with open(path, newline = '', encoding = 'utf-8') as infile:
        reader = csv.DictReader(infile, delimiter = '\t')

        for row in reader:
            if resource_name is not None and row.get('resource_name') != resource_name:
                continue

            pairs.add((row['bacterial_accession'], row['human_accession']))

    return pairs


def load_dmi_resource_support(
    path: Path,
) -> dict[str, int]:
    """Count unique TP pairs by DMI support resource."""

    support: dict[str, set[tuple[str, str]]] = {}

    with open(path, newline = '', encoding = 'utf-8') as infile:
        reader = csv.DictReader(infile, delimiter = '\t')

        for row in reader:
            support.setdefault(row['resource'], set()).add(
                (row['bacterial_accession'], row['human_accession']),
            )

    return {
        resource_name: len(pairs)
        for resource_name, pairs in support.items()
    }


def format_percent(value: float, digits: int = 1) -> str:
    """Format a fraction as a percentage string."""

    return f'{value * 100:.{digits}f}%'


def panel_background(
    canvas: SvgCanvas,
    x: float,
    y: float,
    width: float,
    height: float,
    label: str,
    title: str,
    subtitle: str | None = None,
) -> None:
    """Draw a panel background and heading."""

    canvas.rect(
        x = x,
        y = y,
        width = width,
        height = height,
        fill = COLORS['panel'],
        stroke = COLORS['grid'],
        stroke_width = 1.2,
        rx = 18,
    )
    canvas.text(
        x + 24,
        y + 34,
        label,
        size = 22,
        fill = COLORS['accent_blue'],
        weight = '700',
    )
    canvas.text(
        x + 58,
        y + 34,
        title,
        size = 20,
        fill = COLORS['ink'],
        weight = '700',
    )

    if subtitle:
        canvas.text(
            x + 58,
            y + 58,
            subtitle,
            size = 13,
            fill = COLORS['muted'],
        )


def draw_vertical_bar_chart(
    canvas: SvgCanvas,
    x: float,
    y: float,
    width: float,
    height: float,
    values: list[tuple[str, float, str]],
    max_value: float,
    y_tick_step: float,
    y_label: str,
    value_label_format: str = 'count',
    x_label_size: float = 12,
) -> None:
    """Draw a publication-style vertical bar chart."""

    chart_left = x + 72
    chart_right = x + width - 20
    chart_top = y + 10
    chart_bottom = y + height - 54
    plot_width = chart_right - chart_left
    plot_height = chart_bottom - chart_top
    tick_count = int(max_value / y_tick_step)

    for tick_index in range(tick_count + 1):
        tick_value = tick_index * y_tick_step
        tick_y = chart_bottom - (tick_value / max_value) * plot_height
        canvas.line(
            chart_left,
            tick_y,
            chart_right,
            tick_y,
            stroke = COLORS['grid'],
            stroke_width = 1,
        )

        tick_label = (
            f'{tick_value:.1f}'
            if y_tick_step < 1
            else f'{int(tick_value)}'
        )
        canvas.text(
            chart_left - 10,
            tick_y + 4,
            tick_label,
            size = 11,
            fill = COLORS['muted'],
            anchor = 'end',
        )

    canvas.line(
        chart_left,
        chart_top,
        chart_left,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )
    canvas.line(
        chart_left,
        chart_bottom,
        chart_right,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )

    category_width = plot_width / len(values)
    bar_width = min(56, category_width * 0.6)

    for index, (label, value, color) in enumerate(values):
        bar_x = chart_left + category_width * index + (category_width - bar_width) / 2
        bar_height = (value / max_value) * plot_height
        bar_y = chart_bottom - bar_height

        canvas.rect(
            x = bar_x,
            y = bar_y,
            width = bar_width,
            height = bar_height,
            fill = color,
            stroke = color,
            stroke_width = 0,
            rx = 6,
        )

        if value_label_format == 'count':
            value_label = f'{int(round(value)):,}'
        else:
            value_label = f'{value:.2f}%'

        canvas.text(
            bar_x + bar_width / 2,
            bar_y - 8,
            value_label,
            size = 12,
            fill = COLORS['ink'],
            weight = '600',
            anchor = 'middle',
        )

        canvas.text_block(
            bar_x + bar_width / 2,
            chart_bottom + 20,
            label.split('\n'),
            size = x_label_size,
            fill = COLORS['ink'],
            anchor = 'middle',
        )

    canvas.text(
        x + 8,
        y + height / 2 + 8,
        y_label,
        size = 13,
        fill = COLORS['muted'],
        weight = '600',
    )


def draw_grouped_bar_chart(
    canvas: SvgCanvas,
    x: float,
    y: float,
    width: float,
    height: float,
    groups: list[tuple[str, list[tuple[str, float, str]]]],
    max_value: float,
    y_tick_step: float,
    y_label: str,
    value_as_percent: bool = True,
) -> None:
    """Draw a grouped vertical bar chart."""

    chart_left = x + 74
    chart_right = x + width - 20
    chart_top = y + 10
    chart_bottom = y + height - 64
    plot_width = chart_right - chart_left
    plot_height = chart_bottom - chart_top

    for tick_index in range(int(max_value / y_tick_step) + 1):
        tick_value = tick_index * y_tick_step
        tick_y = chart_bottom - (tick_value / max_value) * plot_height
        canvas.line(
            chart_left,
            tick_y,
            chart_right,
            tick_y,
            stroke = COLORS['grid'],
            stroke_width = 1,
        )
        canvas.text(
            chart_left - 10,
            tick_y + 4,
            f'{tick_value:.0f}',
            size = 11,
            fill = COLORS['muted'],
            anchor = 'end',
        )

    canvas.line(
        chart_left,
        chart_top,
        chart_left,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )
    canvas.line(
        chart_left,
        chart_bottom,
        chart_right,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )

    group_width = plot_width / len(groups)

    for group_index, (group_name, series) in enumerate(groups):
        inner_width = group_width * 0.72
        inner_left = chart_left + group_index * group_width + (group_width - inner_width) / 2
        bar_width = inner_width / len(series) - 10

        for series_index, (series_name, value, color) in enumerate(series):
            bar_x = inner_left + series_index * (bar_width + 10)
            bar_height = (value / max_value) * plot_height
            bar_y = chart_bottom - bar_height

            canvas.rect(
                x = bar_x,
                y = bar_y,
                width = bar_width,
                height = bar_height,
                fill = color,
                stroke = color,
                stroke_width = 0,
                rx = 6,
            )

            label = f'{value:.2f}%' if value_as_percent else f'{value:.0f}'
            canvas.text(
                bar_x + bar_width / 2,
                bar_y - 8,
                label,
                size = 11.5,
                fill = COLORS['ink'],
                weight = '600',
                anchor = 'middle',
            )

        canvas.text(
            inner_left + inner_width / 2,
            chart_bottom + 22,
            group_name,
            size = 12,
            fill = COLORS['ink'],
            weight = '600',
            anchor = 'middle',
        )

    legend_x = chart_right - 180
    legend_y = chart_top + 8

    for legend_index, (series_name, _, color) in enumerate(groups[0][1]):
        y_pos = legend_y + legend_index * 24
        canvas.rect(
            x = legend_x,
            y = y_pos - 10,
            width = 14,
            height = 14,
            fill = color,
            stroke = color,
            stroke_width = 0,
            rx = 3,
        )
        canvas.text(
            legend_x + 22,
            y_pos + 1,
            series_name,
            size = 12,
            fill = COLORS['muted'],
        )

    canvas.text(
        x + 8,
        y + height / 2 + 8,
        y_label,
        size = 13,
        fill = COLORS['muted'],
        weight = '600',
    )


def draw_horizontal_bar_chart(
    canvas: SvgCanvas,
    x: float,
    y: float,
    width: float,
    height: float,
    values: list[tuple[str, int, str]],
    max_value: int,
    x_tick_step: int,
    x_label: str,
    annotation: str | None = None,
) -> None:
    """Draw a horizontal bar chart."""

    chart_left = x + 210
    chart_right = x + width - 20
    chart_top = y + 10
    chart_bottom = y + height - 42
    plot_width = chart_right - chart_left
    row_height = (chart_bottom - chart_top) / len(values)

    for tick_value in range(0, max_value + 1, x_tick_step):
        tick_x = chart_left + (tick_value / max_value) * plot_width
        canvas.line(
            tick_x,
            chart_top,
            tick_x,
            chart_bottom,
            stroke = COLORS['grid'],
            stroke_width = 1,
        )
        canvas.text(
            tick_x,
            chart_bottom + 20,
            f'{tick_value:,}',
            size = 11,
            fill = COLORS['muted'],
            anchor = 'middle',
        )

    canvas.line(
        chart_left,
        chart_top,
        chart_left,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )
    canvas.line(
        chart_left,
        chart_bottom,
        chart_right,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )

    for index, (label, value, color) in enumerate(values):
        y_center = chart_top + row_height * index + row_height * 0.5
        bar_height = min(26, row_height * 0.6)
        bar_width = (value / max_value) * plot_width

        canvas.text_block(
            chart_left - 14,
            y_center - 3,
            label.split('\n'),
            size = 12,
            fill = COLORS['ink'],
            anchor = 'end',
        )
        canvas.rect(
            x = chart_left,
            y = y_center - bar_height / 2,
            width = bar_width,
            height = bar_height,
            fill = color,
            stroke = color,
            stroke_width = 0,
            rx = 6,
        )
        canvas.text(
            chart_left + bar_width + 8,
            y_center + 4,
            f'{value:,}',
            size = 12,
            fill = COLORS['ink'],
            weight = '600',
        )

    canvas.text(
        chart_left + plot_width / 2,
        chart_bottom + 34,
        x_label,
        size = 13,
        fill = COLORS['muted'],
        weight = '600',
        anchor = 'middle',
    )

    if annotation:
        canvas.text_block(
            chart_left + 12,
            chart_top + 18,
            annotation.split('\n'),
            size = 12,
            fill = COLORS['accent_red'],
            weight = '700',
        )


def draw_precision_recall_scatter(
    canvas: SvgCanvas,
    x: float,
    y: float,
    width: float,
    height: float,
    values: list[tuple[str, float, float, str]],
) -> None:
    """Draw a precision-recall scatter chart."""

    chart_left = x + 70
    chart_right = x + width - 26
    chart_top = y + 16
    chart_bottom = y + height - 52
    plot_width = chart_right - chart_left
    plot_height = chart_bottom - chart_top
    max_x = 9.0
    max_y = 0.7

    for tick_value in range(0, 10, 2):
        tick_x = chart_left + (tick_value / max_x) * plot_width
        canvas.line(
            tick_x,
            chart_top,
            tick_x,
            chart_bottom,
            stroke = COLORS['grid'],
            stroke_width = 1,
        )
        canvas.text(
            tick_x,
            chart_bottom + 20,
            f'{tick_value}',
            size = 11,
            fill = COLORS['muted'],
            anchor = 'middle',
        )

    for tick_index in range(0, 8):
        tick_value = tick_index * 0.1
        tick_y = chart_bottom - (tick_value / max_y) * plot_height
        canvas.line(
            chart_left,
            tick_y,
            chart_right,
            tick_y,
            stroke = COLORS['grid'],
            stroke_width = 1,
        )
        canvas.text(
            chart_left - 10,
            tick_y + 4,
            f'{tick_value:.1f}',
            size = 11,
            fill = COLORS['muted'],
            anchor = 'end',
        )

    canvas.line(
        chart_left,
        chart_top,
        chart_left,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )
    canvas.line(
        chart_left,
        chart_bottom,
        chart_right,
        chart_bottom,
        stroke = COLORS['ink'],
        stroke_width = 1.6,
    )

    for label, recall_percent, precision_percent, color in values:
        px = chart_left + (recall_percent / max_x) * plot_width
        py = chart_bottom - (precision_percent / max_y) * plot_height
        canvas.circle(
            cx = px,
            cy = py,
            radius = 7.5,
            fill = color,
            stroke = 'white',
            stroke_width = 2,
        )
        canvas.text(
            px + 10,
            py - 8,
            label,
            size = 11.5,
            fill = COLORS['ink'],
            weight = '600',
        )

    canvas.text(
        chart_left + plot_width / 2,
        chart_bottom + 36,
        'Recall on full PHISTO benchmark (%)',
        size = 13,
        fill = COLORS['muted'],
        weight = '600',
        anchor = 'middle',
    )
    canvas.text(
        x + 10,
        y + height / 2 + 8,
        'Precision among predicted panel pairs (%)',
        size = 13,
        fill = COLORS['muted'],
        weight = '600',
    )


def draw_venn(
    canvas: SvgCanvas,
    x: float,
    y: float,
    width: float,
    height: float,
    counts: dict[str, int],
) -> None:
    """Draw a simple three-set Venn diagram."""

    cx_f = x + width * 0.34
    cx_r = x + width * 0.66
    cx_d = x + width * 0.50
    cy_bottom = y + height * 0.66
    cy_top = y + height * 0.42
    radius = min(width, height) * 0.24

    canvas.circle(
        cx = cx_f,
        cy = cy_bottom,
        radius = radius,
        fill = COLORS['accent_gold'],
        stroke = COLORS['accent_gold'],
        stroke_width = 2,
        opacity = 0.20,
    )
    canvas.circle(
        cx = cx_r,
        cy = cy_bottom,
        radius = radius,
        fill = COLORS['accent_teal'],
        stroke = COLORS['accent_teal'],
        stroke_width = 2,
        opacity = 0.20,
    )
    canvas.circle(
        cx = cx_d,
        cy = cy_top,
        radius = radius,
        fill = COLORS['accent_red'],
        stroke = COLORS['accent_red'],
        stroke_width = 2,
        opacity = 0.20,
    )

    canvas.text(cx_f - radius * 0.88, cy_bottom, str(counts['f_only']), size = 24, weight = '700')
    canvas.text(cx_r + radius * 0.62, cy_bottom, str(counts['r_only']), size = 24, weight = '700')
    canvas.text(cx_d, cy_top - radius * 0.82, str(counts['d_only']), size = 24, weight = '700', anchor = 'middle')
    canvas.text(x + width * 0.50, y + height * 0.78, str(counts['f_r_only']), size = 20, weight = '700', anchor = 'middle')
    canvas.text(x + width * 0.38, y + height * 0.49, str(counts['f_d_only']), size = 20, weight = '700', anchor = 'middle')
    canvas.text(x + width * 0.62, y + height * 0.49, str(counts['r_d_only']), size = 20, weight = '700', anchor = 'middle')
    canvas.text(x + width * 0.50, y + height * 0.59, str(counts['all_three']), size = 22, weight = '700', anchor = 'middle')

    canvas.text(cx_f - radius * 0.72, cy_bottom - radius - 14, 'Forward DMI', size = 14, fill = COLORS['accent_gold'], weight = '700')
    canvas.text(cx_r + radius * 0.05, cy_bottom - radius - 14, 'Reverse DMI', size = 14, fill = COLORS['accent_teal'], weight = '700')
    canvas.text(cx_d, cy_top - radius - 18, 'Best DDI', size = 14, fill = COLORS['accent_red'], weight = '700', anchor = 'middle')
    canvas.text(
        x + width * 0.50,
        y + height - 18,
        'Union = 1,023 TP pairs (11.33% of PHISTO)',
        size = 14,
        fill = COLORS['ink'],
        weight = '700',
        anchor = 'middle',
    )


def build_main_figure(bundle_dir: Path) -> Path:
    """Build the main publication figure."""

    figure_dir = bundle_dir.joinpath('08_figures')
    figure_dir.mkdir(parents = True, exist_ok = True)

    overview = read_key_value_tsv(
        bundle_dir.joinpath('05_results', 'benchmark_overview.tsv'),
    )
    approach_rows = read_approach_summary(
        bundle_dir.joinpath('05_results', 'approach_summary.tsv'),
    )

    forward_pair_set = load_pair_set(
        bundle_dir.joinpath('05_results', 'forward_dmi_true_positive_details.tsv'),
    )
    reverse_pair_set = load_pair_set(
        bundle_dir.joinpath('05_results', 'reverse_dmi_true_positive_details.tsv'),
    )
    ddi_pair_set = load_pair_set(
        bundle_dir.joinpath('05_results', 'ddi_true_positive_details.tsv'),
        resource_name = '3did_plus_domine_v2_all',
    )

    venn_counts = {
        'f_only': len(forward_pair_set - reverse_pair_set - ddi_pair_set),
        'r_only': len(reverse_pair_set - forward_pair_set - ddi_pair_set),
        'd_only': len(ddi_pair_set - forward_pair_set - reverse_pair_set),
        'f_r_only': len((forward_pair_set & reverse_pair_set) - ddi_pair_set),
        'f_d_only': len((forward_pair_set & ddi_pair_set) - reverse_pair_set),
        'r_d_only': len((reverse_pair_set & ddi_pair_set) - forward_pair_set),
        'all_three': len(forward_pair_set & reverse_pair_set & ddi_pair_set),
    }

    main_canvas = SvgCanvas(width = 1600, height = 1120)
    main_canvas.text(
        72,
        54,
        'PHISTO Benchmarking of Extended MicrobioLink Upstream Interaction Modes',
        size = 28,
        fill = COLORS['ink'],
        weight = '700',
    )
    main_canvas.text(
        72,
        82,
        'Broad PHISTO bacteria-human panel benchmarked with forward DMI, reverse DMI, and DDI using ELM, 3did, and DOMINE resources.',
        size = 15,
        fill = COLORS['muted'],
    )

    panel_background(
        main_canvas,
        56,
        110,
        720,
        440,
        'A',
        'Benchmark Panel and Structural Coverage',
        'Coverage bottlenecks are predominantly on the bacterial side.',
    )
    panel_background(
        main_canvas,
        824,
        110,
        720,
        440,
        'B',
        'Recovery by Interaction Mode',
        'Reverse DMI recovers the largest share of the PHISTO benchmark.',
    )
    panel_background(
        main_canvas,
        56,
        592,
        720,
        440,
        'C',
        'DDI Resource Comparison',
        'Combining 3did with the full DOMINE v2 table yields the best DDI recall.',
    )
    panel_background(
        main_canvas,
        824,
        592,
        720,
        440,
        'D',
        'Overlap of Best-Performing Output Sets',
        'Union coverage remains limited despite resource extension.',
    )

    coverage_values = [
        ('Bac\ntotal', float(overview['unique_bacterial_accessions']), COLORS['grey']),
        ('Bac\nsequence', float(overview['bacterial_with_sequence']), COLORS['accent_teal']),
        ('Bac\nPfam', float(overview['bacterial_with_pfam']), COLORS['accent_blue']),
        ('Hum\ntotal', float(overview['unique_human_accessions']), COLORS['grey']),
        ('Hum\nsequence', float(overview['human_with_sequence']), COLORS['accent_teal']),
        ('Hum\nPfam', float(overview['human_with_pfam']), COLORS['accent_blue']),
    ]
    draw_vertical_bar_chart(
        canvas = main_canvas,
        x = 78,
        y = 178,
        width = 678,
        height = 330,
        values = coverage_values,
        max_value = 4000,
        y_tick_step = 1000,
        y_label = 'Protein count',
    )

    main_canvas.text_block(
        102,
        528,
        [
            'Full panel: 9,027 unique PHISTO pairs',
            'Resolved UniProt IDs: 2,715 bacterial, 3,736 human',
        ],
        size = 13,
        fill = COLORS['muted'],
    )

    summary_by_key = {
        (row.approach, row.resource_name): row
        for row in approach_rows
    }
    grouped_recovery = [
        (
            'Forward\nDMI',
            [
                (
                    'Full-panel recall',
                    summary_by_key[('forward_dmi', 'elm_plus_3did')].recall_total * 100,
                    COLORS['accent_blue'],
                ),
                (
                    'Analyzable recall',
                    summary_by_key[('forward_dmi', 'elm_plus_3did')].recall_on_analyzable_pairs * 100,
                    COLORS['accent_teal'],
                ),
            ],
        ),
        (
            'Reverse\nDMI',
            [
                (
                    'Full-panel recall',
                    summary_by_key[('reverse_dmi', 'elm_plus_3did')].recall_total * 100,
                    COLORS['accent_blue'],
                ),
                (
                    'Analyzable recall',
                    summary_by_key[('reverse_dmi', 'elm_plus_3did')].recall_on_analyzable_pairs * 100,
                    COLORS['accent_teal'],
                ),
            ],
        ),
        (
            'Best\nDDI',
            [
                (
                    'Full-panel recall',
                    summary_by_key[('ddi', '3did_plus_domine_v2_all')].recall_total * 100,
                    COLORS['accent_blue'],
                ),
                (
                    'Analyzable recall',
                    summary_by_key[('ddi', '3did_plus_domine_v2_all')].recall_on_analyzable_pairs * 100,
                    COLORS['accent_teal'],
                ),
            ],
        ),
    ]
    draw_grouped_bar_chart(
        canvas = main_canvas,
        x = 846,
        y = 178,
        width = 678,
        height = 330,
        groups = grouped_recovery,
        max_value = 20,
        y_tick_step = 5,
        y_label = 'Recall (%)',
        value_as_percent = True,
    )
    main_canvas.text_block(
        866,
        528,
        [
            'Recovered TP pairs:',
            'Forward DMI 144, Reverse DMI 756, Best DDI 292',
        ],
        size = 13,
        fill = COLORS['muted'],
    )

    ddi_chart_values = [
        ('3did', 166, COLORS['accent_blue']),
        ('DOMINE\nHC', 22, COLORS['light_red']),
        ('DOMINE\nall', 217, COLORS['accent_red']),
        ('3did +\nDOMINE HC', 179, COLORS['accent_purple']),
        ('3did +\nDOMINE all', 292, COLORS['accent_green']),
    ]
    draw_vertical_bar_chart(
        canvas = main_canvas,
        x = 78,
        y = 660,
        width = 678,
        height = 308,
        values = [(label, float(value), color) for label, value, color in ddi_chart_values],
        max_value = 320,
        y_tick_step = 80,
        y_label = 'Recovered PHISTO TP pairs',
    )
    recall_labels = {
        '3did': '1.84%',
        'DOMINE\nHC': '0.24%',
        'DOMINE\nall': '2.40%',
        '3did +\nDOMINE HC': '1.98%',
        '3did +\nDOMINE all': '3.23%',
    }
    chart_left = 78 + 72
    chart_right = 78 + 678 - 20
    plot_width = chart_right - chart_left
    category_width = plot_width / len(ddi_chart_values)
    bar_width = min(56, category_width * 0.6)
    for index, (label, value, _) in enumerate(ddi_chart_values):
        bar_x = chart_left + category_width * index + (category_width - bar_width) / 2
        main_canvas.text(
            bar_x + bar_width / 2,
            962,
            recall_labels[label],
            size = 12,
            fill = COLORS['muted'],
            anchor = 'middle',
        )
    main_canvas.text(
        392,
        1006,
        'Full-panel recall labels shown below each bar.',
        size = 12,
        fill = COLORS['muted'],
        anchor = 'middle',
    )

    draw_venn(
        canvas = main_canvas,
        x = 866,
        y = 650,
        width = 636,
        height = 340,
        counts = venn_counts,
    )

    output_path = figure_dir.joinpath('figure_1_main_benchmark_overview.svg')
    main_canvas.save(output_path)
    return output_path


def build_context_figure(bundle_dir: Path) -> Path:
    """Build a supplementary context and resource figure."""

    figure_dir = bundle_dir.joinpath('08_figures')
    figure_dir.mkdir(parents = True, exist_ok = True)

    pmid_counts = read_count_table(
        bundle_dir.joinpath('05_results', 'phisto_pmid_counts.tsv'),
        'pubmed_id',
    )[:10]
    method_counts = read_count_table(
        bundle_dir.joinpath('05_results', 'phisto_method_counts.tsv'),
        'experimental_method',
    )[:10]
    approach_rows = read_approach_summary(
        bundle_dir.joinpath('05_results', 'approach_summary.tsv'),
    )
    forward_support = load_dmi_resource_support(
        bundle_dir.joinpath('05_results', 'forward_dmi_true_positive_details.tsv'),
    )
    reverse_support = load_dmi_resource_support(
        bundle_dir.joinpath('05_results', 'reverse_dmi_true_positive_details.tsv'),
    )

    canvas = SvgCanvas(width = 1600, height = 1120)
    canvas.text(
        72,
        54,
        'PHISTO Benchmark Context and Resource Effects',
        size = 28,
        fill = COLORS['ink'],
        weight = '700',
    )
    canvas.text(
        72,
        82,
        'Supplementary figure showing benchmark skew, method composition, precision-recall tradeoffs, and DMI resource support.',
        size = 15,
        fill = COLORS['muted'],
    )

    panel_background(
        canvas,
        56,
        110,
        720,
        440,
        'A',
        'PHISTO PMID Composition',
        'The benchmark is dominated by one large high-throughput study.',
    )
    panel_background(
        canvas,
        824,
        110,
        720,
        440,
        'B',
        'Experimental Method Composition',
        'Two-hybrid pooling accounts for the overwhelming majority of unique PHISTO pairs.',
    )
    panel_background(
        canvas,
        56,
        592,
        720,
        440,
        'C',
        'Precision-Recall Positioning of Tested Configurations',
        'DDI variants improve precision; reverse DMI maximizes recall.',
    )
    panel_background(
        canvas,
        824,
        592,
        720,
        440,
        'D',
        'DMI True-Positive Resource Support',
        '3did provides the majority of DMI-supported TP recovery in both directions.',
    )

    pmid_chart_values = []

    for index, (label, count) in enumerate(pmid_counts):
        color = COLORS['accent_red'] if index == 0 else COLORS['light_blue']
        pmid_chart_values.append((label, count, color))

    draw_horizontal_bar_chart(
        canvas = canvas,
        x = 78,
        y = 178,
        width = 678,
        height = 330,
        values = pmid_chart_values,
        max_value = 9000,
        x_tick_step = 3000,
        x_label = 'Unique PHISTO pairs supported by PMID',
        annotation = 'PMID 20711500 alone supports\n8,556 / 9,027 unique pairs (94.8%).',
    )

    method_chart_values = []

    for index, (label, count) in enumerate(method_counts):
        wrapped = label.replace(' pooling approach', '\npooling approach')
        color = COLORS['accent_gold'] if index == 0 else COLORS['light_teal']
        method_chart_values.append((wrapped, count, color))

    draw_horizontal_bar_chart(
        canvas = canvas,
        x = 846,
        y = 178,
        width = 678,
        height = 330,
        values = method_chart_values,
        max_value = 9000,
        x_tick_step = 3000,
        x_label = 'Unique PHISTO pairs supported by method',
        annotation = 'The panel is dominated by\ntwo hybrid pooling approach.',
    )

    scatter_values = []

    label_map = {
        ('forward_dmi', 'elm_plus_3did'): 'F-DMI',
        ('reverse_dmi', 'elm_plus_3did'): 'R-DMI',
        ('ddi', '3did_current'): '3did',
        ('ddi', 'domine_v2_hc'): 'D-HC',
        ('ddi', 'domine_v2_all'): 'D-all',
        ('ddi', '3did_plus_domine_v2_hc'): '3+HC',
        ('ddi', '3did_plus_domine_v2_all'): '3+all',
    }
    color_map = {
        ('forward_dmi', 'elm_plus_3did'): COLORS['accent_gold'],
        ('reverse_dmi', 'elm_plus_3did'): COLORS['accent_teal'],
        ('ddi', '3did_current'): COLORS['accent_blue'],
        ('ddi', 'domine_v2_hc'): COLORS['light_red'],
        ('ddi', 'domine_v2_all'): COLORS['accent_red'],
        ('ddi', '3did_plus_domine_v2_hc'): COLORS['accent_purple'],
        ('ddi', '3did_plus_domine_v2_all'): COLORS['accent_green'],
    }

    for row in approach_rows:
        key = (row.approach, row.resource_name)
        scatter_values.append(
            (
                label_map[key],
                row.recall_total * 100,
                row.precision_unique_pairs * 100,
                color_map[key],
            ),
        )

    draw_precision_recall_scatter(
        canvas = canvas,
        x = 78,
        y = 660,
        width = 678,
        height = 308,
        values = scatter_values,
    )

    grouped_support = [
        (
            'Forward\nDMI',
            [
                ('ELM support', float(forward_support.get('ELM', 0)), COLORS['accent_blue']),
                ('3did support', float(forward_support.get('3did', 0)), COLORS['accent_red']),
            ],
        ),
        (
            'Reverse\nDMI',
            [
                ('ELM support', float(reverse_support.get('ELM', 0)), COLORS['accent_blue']),
                ('3did support', float(reverse_support.get('3did', 0)), COLORS['accent_red']),
            ],
        ),
    ]
    draw_grouped_bar_chart(
        canvas = canvas,
        x = 846,
        y = 660,
        width = 678,
        height = 308,
        groups = grouped_support,
        max_value = 800,
        y_tick_step = 200,
        y_label = 'Unique TP pairs with resource support',
        value_as_percent = False,
    )
    canvas.text_block(
        868,
        970,
        [
            'Support counts are non-exclusive:',
            'the same TP pair can have both ELM and 3did support.',
        ],
        size = 12,
        fill = COLORS['muted'],
    )

    output_path = figure_dir.joinpath(
        'figure_2_context_and_resource_effects.svg',
    )
    canvas.save(output_path)
    return output_path


def build_figures_readme(bundle_dir: Path) -> Path:
    """Write a small figure-level README."""

    figure_dir = bundle_dir.joinpath('08_figures')
    figure_dir.mkdir(parents = True, exist_ok = True)

    readme_path = figure_dir.joinpath('README.md')
    readme_path.write_text(
        '\n'.join(
            [
                '# Figure Set',
                '',
                'This folder contains publication-style vector figures generated',
                'directly from the documented PHISTO benchmark bundle.',
                '',
                '## Files',
                '',
                '- `figure_1_main_benchmark_overview.svg`',
                '  Main four-panel benchmark summary figure.',
                '- `figure_2_context_and_resource_effects.svg`',
                '  Supplementary four-panel figure showing benchmark skew,',
                '  experimental-method composition, precision-recall tradeoffs,',
                '  and DMI resource support patterns.',
                '',
                '## Format',
                '',
                '- All figures are SVG vector files.',
                '- They are suitable for editing in Inkscape, Illustrator, or',
                '  other vector editors before manuscript submission.',
                '',
                '## Suggested Usage',
                '',
                '- Use `figure_1_main_benchmark_overview.svg` as the main paper',
                '  benchmark figure.',
                '- Use `figure_2_context_and_resource_effects.svg` as a',
                '  supplementary or methods-facing figure documenting benchmark',
                '  composition and resource behavior.',
            ],
        ),
        encoding = 'utf-8',
    )
    return readme_path


def main() -> None:
    """Build all figure outputs."""

    args = parse_args()
    bundle_dir = Path(args.bundle_dir).resolve()

    figure_dir = bundle_dir.joinpath('08_figures')
    figure_dir.mkdir(parents = True, exist_ok = True)

    figure_1 = build_main_figure(bundle_dir)
    figure_2 = build_context_figure(bundle_dir)
    figure_readme = build_figures_readme(bundle_dir)

    print(
        'Built PHISTO benchmark figures\n'
        f'- {figure_1}\n'
        f'- {figure_2}\n'
        f'- {figure_readme}',
    )


if __name__ == '__main__':
    main()
