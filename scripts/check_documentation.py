#!/usr/bin/env python3
"""Check Rust declaration docs, including private, test, and macro-body items.

This source-level guard complements rustdoc's public API checks. It recognizes
line/block doc comments, doc attributes, and the doc-attribute forwarding used
by this workspace's quantity macro; it is not a general Rust parser.
"""

from pathlib import Path
import re
import sys

DECLARATION = re.compile(
    r"^\s*(?:pub(?:\([^)]*\))?\s+)?"
    r"(?:(?:async|unsafe|const)\s+)*(?:fn|struct|enum|trait)\s+(\w+|\$\w+)"
)


def has_documentation(lines, index):
    """Recognize documentation attached before a declaration and its attributes."""
    preceding = []
    cursor = index - 1
    attribute_depth = 0
    while cursor >= 0:
        line = lines[cursor].strip()
        # Walk backwards over possibly multiline #[...] attributes.
        attribute_depth += line.count("]") - line.count("[")
        if (
            not line
            or line.startswith(("//", "#[", "$(#["))
            or attribute_depth > 0
            or line.endswith("]")
        ):
            preceding.append(line)
            cursor -= 1
            continue
        if line.endswith("*/"):
            while cursor >= 0:
                preceding.append(lines[cursor])
                if "/*" in lines[cursor]:
                    break
                cursor -= 1
            cursor -= 1
            continue
        break
    prefix = "\n".join(preceding)
    return any(marker in prefix for marker in ("///", "/**", "#[doc", "$(#[$type_doc])"))


def missing_documentation(root):
    """Yield paths, line numbers, and names of undocumented Rust declarations."""
    for path in sorted(root.glob("*/src/**/*.rs")):
        lines = path.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            declaration = DECLARATION.match(line)
            if declaration and not has_documentation(lines, index):
                yield path.relative_to(root), index + 1, declaration.group(1)


def main():
    """Report coverage failures and return a process status suitable for CI."""
    root = Path(__file__).resolve().parents[1]
    missing = list(missing_documentation(root))
    for path, line, name in missing:
        print(f"{path}:{line}: missing documentation for {name}")
    if missing:
        return 1
    print("All Rust function, struct, enum, and trait declarations have documentation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
