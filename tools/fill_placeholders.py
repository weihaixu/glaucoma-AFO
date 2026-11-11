#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fill_placeholders.py -- fill the three fields that only the author can supply.

Everything else in the manuscript is either a {{TOKEN}} (filled from the results
ledger by 17_fill_tokens.R) or a value already written in. These three are not:

    ORCID            -> orcid.org, 0000-0000-0000-0000
    repository URL   -> https://github.com/<user>/glaucoma-AFO
    DOI              -> 10.5281/zenodo.<record>   (from Zenodo, see docs section 13)

They appear in BOTH copies of the manuscript (the HTML template and the docx that
gets submitted), so editing one by hand guarantees the two drift apart. This script
writes the same value into both and refuses malformed input, because a wrong ORCID
or DOI on a title page is the kind of error that survives to publication.

Usage
-----
    python tools/fill_placeholders.py --check
    python tools/fill_placeholders.py --orcid 0000-0002-1825-0097
    python tools/fill_placeholders.py --repo https://github.com/xuweihai/glaucoma-AFO
    python tools/fill_placeholders.py --doi 10.5281/zenodo.1234567
    python tools/fill_placeholders.py --orcid ... --repo ... --doi ...   # all at once

Nothing is written for a field you do not pass, so the three can be filled on
different days (ORCID today, DOI only after the Zenodo deposit exists).
"""

from __future__ import annotations

import argparse
import os
import re
import sys

# Longest first: "[repository URL to be inserted]" contains "[repository URL]", and
# replacing the short one first would leave "... to be inserted]" behind.
PLACEHOLDERS = {
    "orcid": "[0000-0000-0000-0000]",
    "repo": "[repository URL to be inserted]",
    "repo_short": "[repository URL]",
    "doi": "[DOI]",
}

RES = {
    "orcid": re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$"),
    "repo": re.compile(r"^https?://\S+$"),
    "doi": re.compile(r"^10\.\d{4,9}/\S+$"),
}

HINTS = {
    "orcid": "expected XXXX-XXXX-XXXX-XXXX (last character may be X)",
    "repo": "expected a full URL starting with http:// or https://",
    "doi": "expected 10.<registrant>/<suffix>, e.g. 10.5281/zenodo.1234567",
}


def validate(kind: str, value: str) -> str:
    value = value.strip()
    if kind == "repo":
        value = value.rstrip("/")
    if not RES[kind].match(value):
        raise SystemExit("%s: %r -- %s" % (kind, value, HINTS[kind]))
    return value


def patch_text(text: str, pairs) -> tuple[str, dict]:
    counts = {}
    for ph, val in pairs:
        n = text.count(ph)
        if n:
            text = text.replace(ph, val)
            counts[ph] = counts.get(ph, 0) + n
    return text, counts


def report_remaining(text: str, where: str) -> list:
    left = [(k, ph) for k, ph in PLACEHOLDERS.items() if ph in text]
    for k, ph in left:
        print("  %-9s still has %s" % (where, ph))
    return left


def main() -> int:
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    ms = os.path.normpath(os.path.join(root, "..", "output", "manuscript"))
    ap.add_argument("--html", default=os.path.join(ms, "formatted-manuscript.html"))
    ap.add_argument("--docx", default=os.path.join(ms, "POAG_TM_manuscript_draft.docx"))
    ap.add_argument("--orcid")
    ap.add_argument("--repo")
    ap.add_argument("--doi")
    ap.add_argument("--check", action="store_true", help="only list what is still missing")
    args = ap.parse_args()

    pairs = []
    if args.orcid:
        pairs.append((PLACEHOLDERS["orcid"], validate("orcid", args.orcid)))
    if args.repo:
        pairs.append((PLACEHOLDERS["repo"], validate("repo", args.repo)))
        pairs.append((PLACEHOLDERS["repo_short"], validate("repo", args.repo)))
    if args.doi:
        pairs.append((PLACEHOLDERS["doi"], validate("doi", args.doi)))

    html_path = os.path.abspath(args.html)
    docx_path = os.path.abspath(args.docx)
    html = open(html_path, encoding="utf-8").read() if os.path.exists(html_path) else None
    if html is None:
        raise SystemExit("HTML not found: %s" % html_path)

    if args.check or not pairs:
        print("Placeholders still present:")
        left_h = report_remaining(html, "HTML")
        left_d = []
        if os.path.exists(docx_path):
            from docx import Document
            left_d = report_remaining("\n".join(p.text for p in Document(docx_path).paragraphs),
                                      "DOCX")
        if not left_h and not left_d:
            print("  none -- title page, code availability and DOI are all filled")
            return 0
        return 1

    # ---- HTML ----
    new_html, counts = patch_text(html, pairs)
    if new_html != html:
        with open(html_path, "w", encoding="utf-8", newline="") as fh:
            fh.write(new_html)
        print("HTML : %d replacement(s) -> %s" % (sum(counts.values()), html_path))

    # ---- docx ----
    if os.path.exists(docx_path):
        from docx import Document
        doc = Document(docx_path)
        n_docx = 0
        for p in doc.paragraphs:
            for r in p.runs:  # every placeholder sits inside a single run (checked)
                if any(ph in r.text for ph in PLACEHOLDERS.values()):
                    r.text, c = patch_text(r.text, pairs)
                    n_docx += sum(c.values())
        if n_docx:
            doc.save(docx_path)
            print("DOCX : %d replacement(s) -> %s" % (n_docx, docx_path))
        else:
            print("DOCX : nothing to replace (already filled, or run split)")

        # A placeholder split across runs is invisible to the loop above; the paragraph
        # text would still carry it. Catch it rather than reporting success.
        joined = "\n".join(p.text for p in Document(docx_path).paragraphs)
        for k, ph in PLACEHOLDERS.items():
            if ph in joined:
                print("WARNING: %s survives in the docx -- it is split across runs; "
                      "fix that paragraph by hand" % ph)
                return 1
    else:
        print("DOCX : not found, skipped (%s)" % docx_path)

    print("\nRemaining after this run:")
    left = report_remaining(new_html, "HTML")
    return 1 if left else 0


if __name__ == "__main__":
    sys.exit(main())
