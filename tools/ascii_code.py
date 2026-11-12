"""Escape non-ASCII characters in R string literals as \\uXXXX.

Usage: python3 tools/ascii_code.py R/*.R
Comments (including roxygen) are left alone. Non-ASCII characters outside
string literals, such as identifiers, are reported, not changed.
"""
import sys

def esc(ch):
    cp = ord(ch)
    return "\\u%04X" % cp if cp <= 0xFFFF else "\\U{%X}" % cp

def fix_line(line, path, n):
    out = []
    q = None          # current quote character
    i = 0
    while i < len(line):
        c = line[i]
        if q:
            if c == "\\" and i + 1 < len(line):
                out.append(line[i:i + 2]); i += 2; continue
            if c == q:
                q = None
            out.append(esc(c) if ord(c) > 127 and q is not None or (ord(c) > 127 and c != q) else c)
            i += 1
            continue
        if c == "#":
            out.append(line[i:]); break
        if c in "\"'`":
            q = c
        elif ord(c) > 127:
            print(f"{path}:{n}: non-ASCII outside a string: {c!r}", file=sys.stderr)
        out.append(c)
        i += 1
    return "".join(out)

for path in sys.argv[1:]:
    src = open(path, encoding="utf-8").read().split("\n")
    new = [l if l.lstrip().startswith("#") else fix_line(l, path, k + 1) for k, l in enumerate(src)]
    if new != src:
        open(path, "w", encoding="utf-8").write("\n".join(new))
        print("escaped", path)
