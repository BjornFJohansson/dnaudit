#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
PDC1tp_KlLAC4_PGI1tp
              PGI1tp_KlLAC12_TPI1tp
'''
from pydna import pcr, parse, Genbank, Assembly
from Bio.Restriction import EcoRV, ZraI

(p167,
 p166,
 p512,
 p665) = parse(   '''
                        >167_pCAPSfw (24-mer)
                        TCCTGACGGGTAATTTTGATTTGC
                        >166_pCAPSrv (24-mer) EMPTY
                        CTGTGAAGTGGCTGAAATTTCGTA
                        >512_crp_EcoRV (29-mer)
                        ttcgccaattgattcaggtaaatccggat
                        >665_crp_ZraI (29-mer) (new)
                        agcagagtctgtgcaatgcggccgctgac''', ds=False)

from pYPK0 import pYPK0

from pYPK0_PDC1tp_KlLAC4_PGI1tp  import pYPK0_PDC1tp_KlLAC4_PGI1tp  as cas1
from pYPK0_PGI1tp_KlLAC12_TPI1tp import pYPK0_PGI1tp_KlLAC12_TPI1tp as casN

cas1  = pcr( p167, p512, cas1)
casN  = pcr( p665, p166, casN)

pYPK0_E_Z, stuffer = pYPK0.cut((EcoRV, ZraI))

a = Assembly([pYPK0_E_Z, cas1, casN], limit=61)

s = a.circular_products[0]

pYPK0_PDC1tp_KlLAC4_PGI1tp_KlLAC12_TPI1tp = s.synced("tcgcgcgtttcggtgatgacggtgaaaacctctg")
