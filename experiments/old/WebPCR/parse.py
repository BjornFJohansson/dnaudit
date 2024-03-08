#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''
docstring
'''
import re
import sys
import copy
import string
import textwrap
import datetime
import StringIO

from Bio                    import Alphabet
from Bio                    import SeqIO
from Bio.Seq                import Seq
from Bio.SeqRecord          import SeqRecord
from Bio.SeqUtils.CheckSum  import seguid
from formatted_record       import FormattedRecord
from guess_alphabet         import guess_alphabet

def read(data, stamp = False, filter = False):
    '''
    like parse, but only returns one sequence
    '''
    return parse(data, filter).pop()

def parse_to_dict(*args,**kwargs):
    x=parse(*args,**kwargs)
    pass

def parse(data, filter = False):
    '''
    parse(data, filter = False) --> list of FormattedRecord objects

    data is a string containing:

    1. an absolute path to a local file the file will be read in text
    mode and parsed for EMBL, FASTA and Genbank sequences

    2. a path to local a directory
    all files in the directory will be parsed as in 1.

    3. a string containing one or more
    sequences in EMBL, GENBANK, or FASTA format
    mixed formats are allowed.

    4. data can be a list or other iterable of 1 - 3

    if filter is True, sequences will be silently filtered
    for allowed characters (see docs for FormattedRecord)

    '''
    import textwrap
    import os

    frs=[]
    raw=""

    if not hasattr(data, '__iter__'):
        data = (data,)

    for item in data:

        if isinstance(item,basestring):
            item = textwrap.dedent(item)
            item = item.strip()
        else:
            continue

        if os.path.isdir(item):
            for file_ in os.listdir(item):
                with open(file_,'r') as f:
                    raw = "\n\n"+f.read()
                frs.extend( parse_string_to_formatted_records(raw) )
                raw=""

        elif os.path.isfile(os.path.join(os.getcwd(),item)):
            with open(item,'r') as f:
                raw = f.read()
            frs.extend(  parse_string_to_formatted_records(raw) )
            raw=""

        else:
            frs.extend( parse_string_to_formatted_records(item) )

    return frs

def parse_seqs(*args,**kwargs):
    '''alias for parse_string_to_formatted_records'''
    return parse_string_to_formatted_records(*args,**kwargs)

def parse_string_to_formatted_records(rawstring, stamp = False, filter = False):
    '''
    parse_seqs(rawstring, stamp = False, filter = True)

    rawstring is a string containing one or more
    sequences in EMBL, GENBANK, or FASTA format
    mixed formats are allowed.

    if stamp in True each sequence is stamped with SEGUID
    if filter is True, the input sequences will be filtered
    for allowed characters in biological alphabets.

    The function returns a list of FormattedRecord objects
    '''
    from Bio.GenBank import RecordParser
    pattern =  r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"

    rawstring = rawstring.replace( '\r\n', '\n')
    rawstring = rawstring.replace( '\r',   '\n')
    rawseqs = re.findall(pattern,textwrap.dedent(rawstring+"\n\n"),re.MULTILINE)
    sequences=[]

    while rawseqs:
        circular = False
        rawseq = rawseqs.pop(0)
        handle = StringIO.StringIO(rawseq)
        try:
            parsed = SeqIO.parse(handle, "embl").next()
            original_format = "embl"
            if "circular" in rawseq.splitlines()[0]:
                circular = True
        except StopIteration:
            handle.seek(0)
            try:
                parsed = SeqIO.parse(handle, "genbank").next()
                original_format = "genbank"
                handle.seek(0)
                parser = RecordParser()
                residue_type = parser.parse(handle).residue_type
                if "circular" in residue_type:
                    circular = True
            except StopIteration:
                handle.seek(0)
                try:
                    parsed = SeqIO.parse(handle, "fasta").next()
                    original_format = "fasta"
                    if "circular" in rawseq.splitlines()[0]:
                        circular = True
                except StopIteration:
                    continue

        sequences.append(FormattedRecord(  parsed,
                                           original_format      = original_format,
                                           raw_string           = rawseq,
                                           circular             = circular,
                                           filter               = filter))
        handle.close()

    return sequences



if __name__=="__main__":

    from Bio.Alphabet import DNAAlphabet

    seqs = parse('../tests/RefDataBjorn.fas', filter=False)

    assert len(seqs) == 771
    assert list(set([len (a) for a in seqs])) == [901]
    for i,s in enumerate(seqs):
        a = s.description
        b = a.split("|")
        c =  "|".join([b[0],b[1],b[3]])
        s.id = b[2].replace(" ","_")+"_"+str(i)
        s.description = ""
        if b[3]=="Zenion hololepis":
            s.id = b[3].replace(" ","_")+"_"+str(i)
        s.seq.alphabet = DNAAlphabet()




