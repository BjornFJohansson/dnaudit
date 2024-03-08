#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
collect all text files
read and collect snippets
parse snippets
    make snippet instance
    collect sequences in a list of tuples [(id, seqobj),]
    simulate and compare to old result
    add snippet and sequences to digraph

check if all sequences have a unique id
check if all ids have a unique sequence

produce report

visualize graph

"""

import re
from pathlib import Path
from pydna.parsers import extract_from_text
# from pydna.editor import ape
from collections import defaultdict
from itertools import pairwise
from pydna.parsers import parse
from pydna.readers import read, read_primer

from pydna.amplify import Anneal

# sequencefilesuffixes = ["gb", "gbk", "fa", "fasta"]
# sequencefilepaths = [path for i in sequencefilesuffixes for path in Path('.').rglob(f"*.{i}")]
# sequencefiles = parse(sequencefilepaths)
# id_dict = defaultdict(list)
# for sequencefile in sequencefiles:
#     id_dict[sequencefile.id].append(sequencefile)

textfilesuffixes = ["txt", "md"]
textfilepaths = [path for i in textfilesuffixes for path in Path('.').rglob(f"*.{i}")]


class snippet:

    typedict = {"Primer": read_primer
                "Dseqrecord": }

    def __init__(self, text, pth):
        self.text = text
        self.pth = pth

    def parse(self, expect, pattern):
        sequences = extract_from_text(text)
        if len(sequences) > expect:
            print("Too many sequences.")
        elif len(sequences) < expect:
            print("Too few sequences.")
        return sequences


class fusion_pcr(snippet):
    expect = 5
    pattern = "Primer Primer *Dseqrecord Dseqrecord"
    pattern = read_primer, read_primer, parse,

    # def parse(self):
    #     sequences = super().parse(self.expect, pattern)
    #     return sequences

    def run(self):
        fp, rp, *frags, result = self.parse()

        fp = read_primer(fp)
        rp = read_primer(rp)
        frags = parse(frags)
        old_result = read(result)

        from pydna.fusionpcr import fuse_by_pcr

        results = fuse_by_pcr(frags, limit=12)

        final_products = []

        for result in results:
            ann = Anneal((fp, rp), result, limit=12)
            final_products.extend(ann.products)

        from pydna.utils import eq

        return all(eq(old_result, p) for p in final_products)


class pcr(snippet):
    pass
class ligation(snippet):
    pass
class homologous_recombination(snippet):
    pass
class crispr(snippet):
    pass
class restriction_digestion(snippet):
    pass



snippet_categories = (pcr, ligation, homologous_recombination, crispr, fusion_pcr, restriction_digestion)
snippet_headers = tuple(f"#.?{'.?'.join(f.__name__.split('_'))}" for f in snippet_categories)

for pth in textfilepaths:
    text = pth.read_text(encoding="utf8")
    # parse snippets
    snippets = re.split(f"({'|'.join(snippet_headers)})", text, re.MULTILINE|re.IGNORECASE)
    snippets = [s for s in snippets if s.strip()]
    for category, content in zip(snippets[::2], snippets[1::2]):
        node = locals()[category.strip("# ").replace(" ", "_").lower()](content, pth)
    break
