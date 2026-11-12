#!/usr/bin/env python3
"""Rewrite decimal literals that R parses inexactly as exact hex literals.

R's parser (R_strtod) accumulates the digits of a decimal literal in long
double. Where long double is the same as double (arm64 macOS), about one in
ten 17 digit literals ends up one unit in the last place away from the
correctly rounded value that Julia's parser gives. Tableau coefficients with
large cancelling weights (the error estimates of the Verner methods, for
example) then change the step sizes. This tool finds every numeric literal
in the code of R/*.R (outside strings and comments), asks R how it parses
each one, and replaces those that differ from the correctly rounded value by
the hexadecimal literal of that value (0x1.8p+1 is 3), which R parses
exactly on every platform.

Usage: python3 tools/exact_literals.py [package root]
Run it after adding R code with long decimal constants (it needs Rscript).
"""
import glob
import json
import os
import re
import subprocess
import sys
import tempfile

NUM = re.compile(r'(?<![\w.])((?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?)(?![\w.])')


def code_spans(src):
    """Spans of src outside comments and string literals."""
    out = []
    i, n, start = 0, len(src), 0
    while i < n:
        c = src[i]
        if c == '#':
            out.append((start, i))
            j = src.find('\n', i)
            i = n if j < 0 else j
            start = i
            continue
        if c in '"\'`':
            out.append((start, i))
            q = c
            i += 1
            while i < n and src[i] != q:
                i += 2 if src[i] == '\\' else 1
            i += 1
            start = i
            continue
        i += 1
    out.append((start, n))
    return out


def literals(src):
    for a, b in code_spans(src):
        for m in NUM.finditer(src, a, b):
            if m.end() <= b:
                yield m


def r_parse(texts):
    """sprintf('%a') of R's parse of each literal."""
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, 'in.txt')
        out = os.path.join(d, 'out.txt')
        with open(inp, 'w') as fh:
            fh.write('\n'.join(texts) + '\n')
        code = ('x <- readLines(commandArgs(TRUE)[1]); '
                'v <- vapply(x, function(s) as.numeric(parse(text = s)[[1L]]), 0); '
                'writeLines(sprintf("%a", v), commandArgs(TRUE)[2])')
        subprocess.run(['Rscript', '-e', code, inp, out], check=True)
        with open(out) as fh:
            return [ln.strip() for ln in fh]


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else '.'
    files = sorted(glob.glob(os.path.join(root, 'R', '*.R')))
    srcs = {f: open(f).read() for f in files}
    texts = sorted({m.group(1) for s in srcs.values() for m in literals(s)})
    parsed = dict(zip(texts, r_parse(texts)))
    bad = {}
    for t in texts:
        h = parsed[t]
        exact = float(t)
        if h in ('Inf', '-Inf', 'NaN', 'NA') or exact in (float('inf'),):
            continue
        if float.fromhex(h) != exact:
            bad[t] = exact.hex()
    changed = 0
    for f, src in srcs.items():
        pieces, last = [], 0
        for m in literals(src):
            t = m.group(1)
            if t in bad:
                pieces.append(src[last:m.start()])
                pieces.append(bad[t])
                last = m.end()
        if pieces:
            pieces.append(src[last:])
            with open(f, 'w') as fh:
                fh.write(''.join(pieces))
            changed += 1
    print('%d literals rewritten in %d files' % (len(bad), changed))


if __name__ == '__main__':
    main()
