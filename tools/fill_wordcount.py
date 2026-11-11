#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fill_wordcount.py -- keep the "Word count (main text)" line in sync with the text.

Why this exists
---------------
Every number in the manuscript is a {{TOKEN}} that 17_fill_tokens.R replaces with a
value read from the results ledger. The word count is the one number that depends on
those replacements (a token "{{AUC_TRAIN}}" is one word, the value "0.87 (0.81-0.92)"
is three), so it cannot be typed in by hand and left there -- it goes stale the moment
the pipeline is re-run on real data.

This script recomputes it from the text itself and writes the same figure into both
copies of the manuscript (HTML template + docx deliverable), so the two can never
drift apart.

Counting rule (stated in the manuscript as "main text")
------------------------------------------------------
    section 1. Introduction  ..  everything before section 6. Declarations
    -- i.e. Introduction, Methods, Results, Discussion, Conclusions.
    Excluded: title block, abstract, keywords, declarations, references,
              tables, figure legends, supplementary list.
A word is a whitespace-separated chunk containing at least one alphanumeric
character, and each {{TOKEN}} counts as one word.

Usage
-----
    python tools/fill_wordcount.py                 # rewrite HTML + docx in place
    python tools/fill_wordcount.py --check         # report only, exit 1 if stale
    python tools/fill_wordcount.py --html a.html --docx b.docx
"""

from __future__ import annotations

import argparse
import os
import re
import sys

START_MARK = "1. Introduction"
END_MARK = "6. Declarations"
LABEL = "Word count (main text):"

TOKEN_RE = re.compile(r"\{\{[^{}]*\}\}")
TAG_RE = re.compile(r"<[^>]+>")
NUM_RE = re.compile(r"[\d,]{1,9}")


def count_words(text: str) -> int:
    """Word count under the rule stated in the docstring."""
    text = TOKEN_RE.sub("X", text)
    chunks = re.split(r"\s+", text)
    return sum(1 for c in chunks if re.search(r"[A-Za-z0-9]", c))


# --------------------------------------------------------------------------- HTML
def html_main_text(raw: str) -> str:
    # Anchor on the heading tags, not on the bare strings: "1. Introduction" and
    # "6. Declarations" also occur inside running prose ("as noted in section 6.
    # Declarations"), and find() on the plain text happily lands there, which shifts
    # the window and quietly changes the count by a few dozen words.
    m1 = re.search(r"<h[1-3][^>]*>\s*1\.\s*Introduction\s*</h[1-3]>", raw)
    m2 = re.search(r"<h[1-3][^>]*>\s*6\.\s*Declarations\s*</h[1-3]>", raw)
    if not m1 or not m2:
        raise SystemExit("HTML: cannot locate the Introduction / Declarations headings")
    seg = raw[m1.start():m2.start()]
    seg = re.sub(r"<!--.*?-->", " ", seg, flags=re.S)
    blocks = re.findall(r"<(?:p|h2|h3|li)[^>]*>(.*?)</(?:p|h2|h3|li)>", seg, flags=re.S)
    # Inline tags must be REMOVED, not replaced by a space: "<i>TGF</i>-&beta;" is one
    # word in the docx (python-docx concatenates runs without separators), and turning
    # the </i> into a space would split it into two and inflate the count by ~30 words
    # over the length of a manuscript. <br> is the one tag that really is a separator.
    txt = "\n".join(TAG_RE.sub("", re.sub(r"<br\s*/?>", " ", b)) for b in blocks)
    for ent, rep in (("&nbsp;", " "), ("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"),
                     ("&minus;", "-"), ("&ndash;", "-"), ("&mdash;", "-")):
        txt = txt.replace(ent, rep)
    return txt


def read_text(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write_text(path: str, data: str) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        fh.write(data)
    os.replace(tmp, path)


def patch_html(raw: str, n: int) -> tuple[str, str]:
    """Return (new_html, old_display). old_display is '' when the line was absent."""
    want = "%s %s" % (LABEL, format(n, ","))
    # Matches "[XXXX]", "[5,161]", "[5230]" and plain "[ to be inserted ]" style stubs.
    pat = re.compile(re.escape(LABEL) + r"\s*\[[^\]]{0,20}\]")
    m = pat.search(raw)
    if m:
        old = m.group(0)
        return raw[:m.start()] + want + raw[m.end():], old
    pat2 = re.compile(re.escape(LABEL) + r"\s*[\d,]+")
    m = pat2.search(raw)
    if m:
        return raw[:m.start()] + want + raw[m.end():], m.group(0)
    return raw, ""


# --------------------------------------------------------------------------- docx
def docx_main_text(doc) -> str:
    paras = [p.text for p in doc.paragraphs]
    ai = next((k for k, t in enumerate(paras) if t.strip().startswith(START_MARK)), -1)
    bi = next((k for k, t in enumerate(paras) if t.strip().startswith(END_MARK)), -1)
    if ai < 0 or bi < 0:
        raise SystemExit("docx: cannot locate %r / %r" % (START_MARK, END_MARK))
    return "\n".join(paras[ai:bi])


def patch_docx(path: str, check_only: bool) -> tuple[int, int, int, str]:
    """Count from the docx itself, then fix its header line.

    Counting from the docx rather than reusing the HTML figure means the two files
    are each checked against their own text: if the HTML and the docx ever drift,
    this script reports two different numbers instead of silently overwriting one
    with the other.
    """
    from docx import Document  # imported late: only needed for the docx path

    doc = Document(path)
    n = count_words(docx_main_text(doc))
    want = "%s %s" % (LABEL, format(n, ","))
    found, patched, note = -1, 0, ""

    for p in doc.paragraphs:
        if LABEL not in p.text:
            continue
        m = NUM_RE.search(p.text.split(LABEL, 1)[1])
        if m:
            found = int(m.group(0).replace(",", ""))
        if check_only:
            continue
        runs = p.runs
        if not runs:
            p.add_run(want)
            patched += 1
            continue
        # Keep the first run's formatting; collapse the rest so the number cannot be
        # split across runs (a number split in two cannot be grepped or edited later).
        head = p.text.split(LABEL, 1)[1]
        tail = head[m.end():] if m else head
        if not tail.strip():
            tail = " \u00b7 Figures: 8 \u00b7 Tables: 6 \u00b7 Supplementary files: 3"
        elif not tail.lstrip().startswith("\u00b7"):
            tail = " " + tail.lstrip()
        runs[0].text = want
        for r in runs[1:]:
            r.text = ""
        if len(runs) > 1:
            runs[1].text = tail
        else:
            runs[0].text = want + tail
        patched += 1

    if patched and not check_only:
        doc.save(path)
    if not patched and not check_only:
        note = "no paragraph containing %r -- nothing written" % LABEL
    return n, found, patched, note


# --------------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    default_html = os.path.normpath(os.path.join(root, "..", "output", "manuscript",
                                                 "formatted-manuscript.html"))
    default_docx = os.path.normpath(os.path.join(root, "..", "output", "manuscript",
                                                 "POAG_TM_manuscript_draft.docx"))
    ap.add_argument("--html", default=default_html)
    ap.add_argument("--docx", default=None,
                    help="default: the standard deliverable when --html is default, "
                         "otherwise the .docx next to --html (if any)")
    ap.add_argument("--check", action="store_true", help="report only; exit 1 if stale")
    args = ap.parse_args()

    html_path = os.path.abspath(args.html)
    if args.docx is not None:
        docx_path = os.path.abspath(args.docx)
    elif os.path.normpath(html_path) == default_html:
        docx_path = default_docx
    else:
        # Called on a working copy (e.g. results/manuscript_filled.html): look for the
        # matching .docx beside it and otherwise leave the standard deliverable alone.
        # Pairing a test HTML with the shipped docx would report a bogus mismatch and,
        # worse, rewrite the deliverable's header from a throw-away file.
        docx_path = os.path.splitext(html_path)[0] + ".docx"
    stale = False
    n_html = n_docx = None

    if os.path.exists(html_path):
        raw = read_text(html_path)
        n_html = count_words(html_main_text(raw))
        new, old = patch_html(raw, n_html)
        print("HTML : %s words | header had: %s" % (format(n_html, ","), old or "<no line>"))
        if old and old != "%s %s" % (LABEL, format(n_html, ",")):
            stale = True
        if not args.check and new != raw:
            write_text(html_path, new)
            print("       -> updated to %s" % format(n_html, ","))
    else:
        print("HTML : not found, skipped (%s)" % html_path)

    if os.path.exists(docx_path):
        n_docx, found, patched, note = patch_docx(docx_path, args.check)
        print("DOCX : %s words | header had: %s" %
              (format(n_docx, ","), format(found, ",") if found >= 0 else "<no line>"))
        if found >= 0 and found != n_docx:
            stale = True
        if note:
            print("       %s" % note)
        elif not args.check and patched:
            print("       -> updated to %s (%d paragraph(s))" % (format(n_docx, ","), patched))
    else:
        print("DOCX : not found, skipped (%s)" % docx_path)

    if n_html is not None and n_docx is not None and n_html != n_docx:
        print("\nWARNING: HTML and docx disagree by %d words -- one of them is missing text."
              % (n_html - n_docx))
        return 1

    if stale:
        print("\nSTALE: header word count does not match the text.")
        return 1 if args.check else 0
    print("\nOK: header word count matches the text.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
