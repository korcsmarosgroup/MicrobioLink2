#!/usr/bin/env python3
"""PostToolUse hook: block edits that leave a Python function over 50 lines."""

import ast
import json
import sys

MAX_FUNCTION_LINES = 50


def find_long_functions(source: str) -> list[tuple[str, int, int]]:
    tree = ast.parse(source)
    violations = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            length = node.end_lineno - node.lineno + 1
            if length > MAX_FUNCTION_LINES:
                violations.append((node.name, node.lineno, length))
    return violations


def main() -> None:
    hook_input = json.load(sys.stdin)
    file_path = hook_input.get("tool_input", {}).get("file_path", "")

    if not file_path.endswith(".py"):
        sys.exit(0)

    try:
        source = open(file_path).read()
    except OSError:
        sys.exit(0)

    try:
        violations = find_long_functions(source)
    except SyntaxError:
        sys.exit(0)

    if not violations:
        sys.exit(0)

    print(f"Functions longer than {MAX_FUNCTION_LINES} lines must be refactored:", file=sys.stderr)
    for name, lineno, length in violations:
        print(f"  {name}() at {file_path}:{lineno} — {length} lines", file=sys.stderr)
    sys.exit(2)


if __name__ == "__main__":
    main()
