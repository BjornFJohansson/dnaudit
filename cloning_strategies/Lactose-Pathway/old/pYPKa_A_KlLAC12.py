#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from pYPKa import pYPKa
from pydna import pcr, parse
from Bio.Restriction import ZraI, AjiI, EcoRV

enz = AjiI

from X06997 import X06997 as gene # Kluyveromyces lactis LAC12 gene for lactose permease.

primers = parse('''
>760_KlLAC12_rv (20-mer)
ttaaacagattctgcctctg

>759_KlLAC12_fw (19-mer)
aaatggcagatcattcgag
''', ds=False)

pYPKa_cut = pYPKa.cut(enz).pop()

gene = pcr( primers, gene)

pYPKa_enz_gene = (pYPKa_cut + gene).looped()

pYPKa_enz_gene = pYPKa_enz_gene.synced(pYPKa)

pYPKa_A_KlLAC12 = pYPKa_enz_gene

