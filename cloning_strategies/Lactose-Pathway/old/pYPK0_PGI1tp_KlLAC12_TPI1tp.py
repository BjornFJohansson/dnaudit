#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna import pcr, parse, Genbank, Assembly
from Bio.Restriction import ZraI, AjiI, EcoRV

from pYPK0 import pYPK0
from pYPKa_Z_PGI1tp    import pYPKa_Z_PGI1tp    as first
from pYPKa_A_KlLAC12   import pYPKa_A_KlLAC12   as middle
from pYPKa_E_TPI1tp    import pYPKa_E_TPI1tp    as last

p167,p166,p468,p467,p567,p568  =  parse('''
                                            >167_pCAPSfw (24-mer)
                                            TCCTGACGGGTAATTTTGATTTGC

                                            >166_pCAPSrv (24-mer) EMPTY
                                            CTGTGAAGTGGCTGAAATTTCGTA

                                            >468_pCAPs_release_fw (25-mer) 79.66
                                            gtcgaggaacgccaggttgcccact

                                            >467_pCAPs_release_re (31-mer)
                                            ATTTAAatcctgatgcgtttgtctgcacaga

                                            >567_pCAPsAjiIF (23-mer)
                                            GTcggctgcaggtcactagtgag

                                            >568_pCAPsAjiIR (22-mer)
                                            GTGCcatctgtgcagacaaacg''')


first  = pcr( p167, p567, first)
middle = pcr( p468, p467, middle)
last   = pcr( p568, p166, last)

pYPK0_E_Z, stuffer = pYPK0.cut((EcoRV, ZraI))

a = Assembly([first, middle, last, pYPK0_E_Z], limit=31)

s = a.circular_products[0]

pYPK0_PGI1tp_KlLAC12_TPI1tp = s.synced("tcgcgcgtttcggtgatgacggtgaaaacctctg")

