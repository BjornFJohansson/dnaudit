#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import copy
import string
import textwrap
import datetime
import StringIO
import copy

from Bio                    import Alphabet
from Bio                    import SeqIO
from Bio.Seq                import Seq
from Bio.SeqRecord          import SeqRecord
from Bio.SeqUtils.CheckSum  import seguid
from guess_alphabet         import guess_alphabet

#http://docs.python.org/howto/descriptor.html#invoking-descriptors
# linear, circular should be descriptors

class FormattedRecord(SeqRecord):
    '''
    FormattedRecord( record,
                     circular = False,
                     linear   = True,
                     filter   = False,
                     *args, **kwargs):

    FormattedRecord is similar to biopython SeqRecord, but optionally
    also holds information on the topology of the sequence and some
    additional properties.

    Arguments
    ---------

    record = string, biopython Seq or SeqRecord object

    circular = True or False sets the topology of the sequence
    or
    linear   = True or False sets the topology of the sequence

    only one of the arguments circular can be set,

    if filter is True, the record is filtered for allowed characters

    ACBEDGFIHKJMLONQPSRUTWVYXZacbedgfihkjmlonqpsrutwvyxz

    properties
    ----------

    self.raw = "not set"

    Can be set by the parser function that generates the
    formatted record. Usually set to the text from which the
    formatted record was parsed.

    self.parsed_from = "not defined"

    Can be set to the sequence format that was used by the parser.

    self.guessed_alphabet = False

    if record was not supplied with an alphabet (see biopython docs
    for Seq object), the alphabet will be guessed from the sequence and
    the record object will be set with this property. self.guessed_alphabet
    will be set to true to reflect this.

    self.warnings = ""

    Set to a clear text message containing characters that were
    filtered out (if any) and which alphabet was guessed (if any).
    Otherwise this is an empty string.
    '''

    def __init__(self, record,
                         circular               = None,
                         linear                 = None,
                         filter                 = False,
                         *args, **kwargs):

        self.raw                = "not set"
        self.parsed_from        = "not defined"
        self.filtered           = None
        self._circular          = None
        self._linear            = None
        self.guessed_alphabet   = False
        self.warnings           = ""

        if isinstance(record, basestring):

            SeqRecord.__init__(self, Seq(record), *args, **kwargs)

        elif isinstance(record, Seq):

            SeqRecord.__init__(self, record, *args, **kwargs)

        elif isinstance(record, SeqRecord) or hasattr(record, "watson"):
            for key, value in record.__dict__.items():
                setattr(self, key, value )
        else:
            raise TypeError("record argument needs to be a string, a Seq object or a SeqRecord object, got {}".format(type(record)))

        if filter:
            IUPAC_single_alphabet_letters = "ACBEDGFIHKJMLONQPSRUTWVYXZacbedgfihkjmlonqpsrutwvyxz"

            filtered_out = "".join([c for c in self.seq if c not in IUPAC_single_alphabet_letters])

            if filtered_out:
                filtered = "".join([c for c in self.seq if c in IUPAC_single_alphabet_letters])
                self.seq = Seq(filtered, self.seq.alphabet)
                self.filtered = filtered_out
                self.warnings += u"{} non-permitted chars were filtered from the sequence!\n".format(", ".join(set(filtered_out)))

        if not isinstance(self.seq.alphabet, (Alphabet.ProteinAlphabet,Alphabet.DNAAlphabet,Alphabet.RNAAlphabet)):
            self.seq.alphabet = guess_alphabet(self.seq)
            self.guessed_alphabet = True
            self.warnings += str(self.seq.alphabet) + " alphabet guessed from sequence"

        if self.id in ("","."):
            self.id = self.name

        if self.description ==".":
            self.description = ""

        if not 'date' in self.annotations:
            self.annotations.update({"date": datetime.date.today().strftime("%d-%b-%Y").upper()})

        if circular == linear == self._circular == self.linear == None:
            self._circular = False
            self._linear   = True
            return

        if circular == None and linear in (True, False,):
            self._linear = linear
            self._circular = not linear
        elif linear == None and circular in (True, False,):
            self._circular = circular
            self._linear = not circular

        if self._linear != (not self._circular):
            raise(ValueError('''properties circular and/or linear not set correctly
                             circular = {}, linear    = {}'''.format(self._circular,self.linear)
                              )
                   )
    def get_linear(self):
        return self._linear
    def set_linear(self, value):
        self._linear   = value
        self._circular = not value
    def get_circular(self):
        return self._circular
    def set_circular(self, value):
        self._circular = value
        self._linear   = not value
    linear   = property(get_linear  , set_linear,   "I'm the 'linear' property.")
    circular = property(get_circular, set_circular, "I'm the 'circular' property.")



    def stamp(self):
        pattern = "(SEGUID|seguid)\s*\S{27}"
        try:
            stamp = re.search(pattern, self.description).group()
        except AttributeError:
            stamp = "SEGUID {}".format(seguid(self.seq))

        if not self.description:
            self.description = stamp
        elif not re.search(pattern, self.description):
            self.description += " "+stamp

    def verify(self):
        pattern = "(SEGUID|seguid)\s*\S{27}"
        try:
            stamp = re.search(pattern, self.description).group()
        except AttributeError:
            return False
        return seguid(self.seq) == stamp[-27:]

    def format(self,f):
        if f in ("genbank","gb"):
            s = SeqRecord.format(self,f)
            if self.circular:
                return s[:55]+"circular"+s[63:]
            else:
                return s[:55]+"linear"+s[61:]
        else:
            return SeqRecord.format(self,f)

    def __str__(self):
        return ("FormattedRecord\nCircular: {0}\nSize: {1}\n".format(self.circular, len(self))+
                SeqRecord.__str__(self))

    def __repr__(self):
        return "FormattedRecord({})".format(len(self))

    def __add__(self, other):
        answer = FormattedRecord(SeqRecord.__add__(self, other))
        answer.circular = False
        return answer

    def __radd__(self, other):
        if isinstance(other, SeqRecord):
            answer = FormattedRecord(SeqRecord.__add__(other, self))
        else:
            answer = FormattedRecord(SeqRecord.__radd__(self,other))
        answer.circular = False
        return answer

#    def _reverse_complement(self):
#        answer = FormattedRecord(self)
#        answer.seq = answer.seq.reverse_complement()
#        return answer

    def reverse_complement(self):
        answer = FormattedRecord(SeqRecord.reverse_complement(self))
        answer.circular = self.circular
        return answer

    def rc(self):
        return self.reverse_complement()

if __name__=="__main__":
    a=FormattedRecord("att", circular = True)
    print a.circular





