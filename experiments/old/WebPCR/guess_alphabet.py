#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
docstring
'''

import string
import warnings

from Bio.Alphabet       import SingleLetterAlphabet
from Bio.Alphabet       import NucleotideAlphabet
from Bio.Alphabet       import ProteinAlphabet
from Bio.Alphabet.IUPAC import extended_protein
from Bio.Alphabet.IUPAC import protein
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Alphabet.IUPAC import extended_dna
from Bio.Alphabet.IUPAC import ambiguous_rna
from Bio.Alphabet.IUPAC import unambiguous_rna

from Bio.Seq            import Seq
from Bio.SeqRecord      import SeqRecord

def guess_alphabet(sequence):
    '''
    guess_alphabet(sequence)

    '''

    if isinstance(sequence, Seq):
        sequence = sequence.tostring()

    elif isinstance(sequence, SeqRecord):
        sequence = sequence.seq.tostring()

    elif not isinstance(sequence, basestring):
        sequence = str(sequence)
        warnings.warn("Input is not string, unicode string, Seq or SeqRecord objects!")

    if len(sequence)<1:
        return SingleLetterAlphabet()

    for c in sequence:
        if c not in string.printable:
            return SingleLetterAlphabet()

    xp = set(extended_protein.letters)
    pr = set(protein.letters)

    ad = set(ambiguous_dna.letters)
    ud = set(unambiguous_dna.letters)
    ed = set(extended_dna.letters)

    ar = set(ambiguous_rna.letters)
    ur = set(unambiguous_rna.letters)

    all = xp|pr|ad|ud|ed|ar|ur

    sequence_chars = set(sequence.upper())

    if sequence_chars - all - set(string.punctuation+string.whitespace):
        return SingleLetterAlphabet()

    nucleic_count = 0

    for letter in "GATCUNgatcun":
        nucleic_count += sequence.count(letter)

    if float(nucleic_count) / float(len(sequence)) >= 0.9: # DNA or RNA
        if 'T' in sequence_chars and 'U' in sequence_chars:
            alphabet = NucleotideAlphabet()
        elif not sequence_chars-ud:
            alphabet = unambiguous_dna
        elif not sequence_chars-ad :
            alphabet = ambiguous_dna
        elif not sequence_chars-ed:
            alphabet = extended_dna
        elif not sequence_chars-ur:
            alphabet = unambiguous_rna
        elif not sequence_chars-ar:
            alphabet = extended_rna
        else:
            alphabet = NucleotideAlphabet()
    else:
        threecode = ['ALA', 'ASX', 'CYS', 'ASP','GLU', 'PHE', 'GLY', 'HIS',
                     'ILE', 'LYS', 'LEU', 'MET','ASN', 'PRO', 'GLN', 'ARG',
                     'SER', 'THR', 'VAL', 'TRP','TYR', 'GLX', 'XAA', 'TER',
                     'SEL', 'PYL', 'XLE']
        tc=set(threecode)
        three_letter_alphabet = set( [ sequence[i:i+3] for i in range(0,len(sequence),3)] )
        if not three_letter_alphabet - tc:
            alphabet = "three letter code"
        elif sequence_chars - pr:
            alphabet = protein
        elif sequence_chars - xp:
            alphabet = extended_protein
        else:
            alphabet = ProteinAlphabet()
    return alphabet


if __name__=="__main__":
    pass