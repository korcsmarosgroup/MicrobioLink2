#!/usr/bin/env python3
"""Run tests across multiple Python versions using uv.

This script replaces tox for local multi-version testing. It creates isolated
virtual environments for each Python version, installs the package with test
dependencies, and runs pytest.

Usage:
    # Run tests on all default Python versions (3.9-3.13)
    python scripts/test.py

    # Run tests on specific versions
    python scripts/test.py 3.11 3.12

    # Run tests with additional pytest arguments
    python scripts/test.py 3.12 -- -v --tb=short

    # Run tests in parallel (requires multiple Python versions installed)
    python scripts/test.py --parallel

Requirements:
    - uv (https://docs.astral.sh/uv/)
    - Python versions you want to test must be installed or discoverable by uv
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DEFAULT_VERSIONS = ['3.9', '3.10', '3.11', '3.12', '3.13']
VENV_PREFIX = '.venv-test-'


def run_command(cmd, cwd=None, capture=False):
    """Run a shell command and return the result."""
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=capture,
        text=True,
    )
    return result


def check_uv():
    """Check if uv is installed."""
    if shutil.which('uv') is None:
        print('Error: uv is not installed.')
        print('Install it with: curl -LsSf https://astral.sh/uv/install.sh | sh')
        sys.exit(1)


def check_python_version(version):
    """Check if a Python version is available via uv."""
    result = run_command(
        ['uv', 'python', 'find', version],
        capture=True,
    )
    return result.returncode == 0


def run_tests_for_version(version, project_root, pytest_args, verbose=True):
    """Run tests for a specific Python version.

    Returns:
        tuple: (version, success, output)
    """
    venv_path = project_root / f'{VENV_PREFIX}{version}'

    if verbose:
        print(f'\n{"=" * 60}')
        print(f'Testing with Python {version}')
        print(f'{"=" * 60}')

    # Check if Python version is available
    if not check_python_version(version):
        msg = f'Python {version} not found. Install with: uv python install {version}'
        if verbose:
            print(f'Warning: {msg}')
        return (version, None, msg)

    # Create virtual environment
    if verbose:
        print(f'Creating virtual environment at {venv_path}...')

    result = run_command(
        ['uv', 'venv', str(venv_path), '--python', version],
        cwd=project_root,
        capture=not verbose,
    )
    if result.returncode != 0:
        return (version, False, f'Failed to create venv: {result.stderr}')

    # Install package with test dependencies
    if verbose:
        print('Installing package with test dependencies...')

    result = run_command(
        ['uv', 'pip', 'install', '-e', '.[tests]', '--python', str(venv_path / 'bin' / 'python')],
        cwd=project_root,
        capture=not verbose,
    )
    if result.returncode != 0:
        return (version, False, f'Failed to install: {result.stderr}')

    # Run pytest
    if verbose:
        print('Running tests...')

    python_path = venv_path / 'bin' / 'python'
    pytest_cmd = [str(python_path), '-m', 'pytest'] + list(pytest_args)

    result = run_command(
        pytest_cmd,
        cwd=project_root,
        capture=not verbose,
    )

    success = result.returncode == 0

    if verbose:
        status = 'PASSED' if success else 'FAILED'
        print(f'\nPython {version}: {status}')

    return (version, success, result.stdout if not verbose else '')


def cleanup_venvs(project_root):
    """Remove test virtual environments."""
    for venv_path in project_root.glob(f'{VENV_PREFIX}*'):
        if venv_path.is_dir():
            shutil.rmtree(venv_path)


def main():
    parser = argparse.ArgumentParser(
        description='Run tests across multiple Python versions using uv.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        'versions',
        nargs='*',
        default=DEFAULT_VERSIONS,
        help=f'Python versions to test (default: {", ".join(DEFAULT_VERSIONS)})',
    )
    parser.add_argument(
        '--parallel', '-p',
        action='store_true',
        help='Run tests in parallel',
    )
    parser.add_argument(
        '--cleanup', '-c',
        action='store_true',
        help='Remove test virtual environments after running',
    )
    parser.add_argument(
        '--cleanup-only',
        action='store_true',
        help='Only remove test virtual environments, do not run tests',
    )
    parser.add_argument(
        'pytest_args',
        nargs='*',
        default=[],
        help='Additional arguments to pass to pytest (after --)',
    )

    # Handle -- separator for pytest args
    if '--' in sys.argv:
        idx = sys.argv.index('--')
        our_args = sys.argv[1:idx]
        pytest_args = sys.argv[idx + 1:]
    else:
        our_args = sys.argv[1:]
        pytest_args = []

    args = parser.parse_args(our_args)

    # Find project root (where pyproject.toml is)
    project_root = Path.cwd()
    while project_root != project_root.parent:
        if (project_root / 'pyproject.toml').exists():
            break
        project_root = project_root.parent
    else:
        print('Error: Could not find pyproject.toml')
        sys.exit(1)

    check_uv()

    # Handle cleanup-only
    if args.cleanup_only:
        print('Cleaning up test virtual environments...')
        cleanup_venvs(project_root)
        print('Done.')
        return

    # Filter versions to only those specified
    versions = args.versions

    print(f'Testing Python versions: {", ".join(versions)}')
    print(f'Project root: {project_root}')

    if pytest_args:
        print(f'Pytest args: {" ".join(pytest_args)}')

    results = {}

    if args.parallel:
        print('\nRunning tests in parallel...')
        with ThreadPoolExecutor(max_workers=len(versions)) as executor:
            futures = {
                executor.submit(
                    run_tests_for_version,
                    version,
                    project_root,
                    pytest_args,
                    verbose=False,
                ): version
                for version in versions
            }
            for future in as_completed(futures):
                version, success, output = future.result()
                results[version] = success
                status = 'PASSED' if success else ('FAILED' if success is False else 'SKIPPED')
                print(f'Python {version}: {status}')
    else:
        for version in versions:
            version, success, output = run_tests_for_version(
                version,
                project_root,
                pytest_args,
                verbose=True,
            )
            results[version] = success

    # Summary
    print(f'\n{"=" * 60}')
    print('SUMMARY')
    print(f'{"=" * 60}')

    passed = sum(1 for s in results.values() if s is True)
    failed = sum(1 for s in results.values() if s is False)
    skipped = sum(1 for s in results.values() if s is None)

    for version in versions:
        success = results.get(version)
        if success is True:
            status = 'PASSED'
        elif success is False:
            status = 'FAILED'
        else:
            status = 'SKIPPED'
        print(f'  Python {version}: {status}')

    print(f'\nTotal: {passed} passed, {failed} failed, {skipped} skipped')

    # Cleanup if requested
    if args.cleanup:
        print('\nCleaning up test virtual environments...')
        cleanup_venvs(project_root)

    # Exit with error if any tests failed
    sys.exit(1 if failed > 0 else 0)


if __name__ == '__main__':
    main()
