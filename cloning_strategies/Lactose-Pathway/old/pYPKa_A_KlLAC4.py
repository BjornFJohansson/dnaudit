#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from pYPKa import pYPKa
from pydna import pcr, parse, Genbank
from Bio.Restriction import ZraI, AjiI, EcoRV

enz = AjiI

primers = parse('''
>758_KlLAC4_rv (20-mer)
ttattcaaaagcgagatcaa

>757_KlLAC4_fw (21-mer)
aaatgtcttgccttattcctg
''', ds=False)

pYPKa_cut = pYPKa.cut(enz).pop()

from M84410 import M84410 as gene

gene = pcr( primers, gene)

pYPKa_enz_gene = (pYPKa_cut + gene).looped()

pYPKa_enz_gene = pYPKa_enz_gene.synced(pYPKa)

pYPKa_A_KlLAC4 = pYPKa_enz_gene
