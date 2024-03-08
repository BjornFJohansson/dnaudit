#!/usr/bin/python
# -*- coding: utf8 -*-

import                              datetime
import                              itertools
import                              math
import                              re
import                              string
import                              sys
import                              textwrap
import                              collections
import                              warnings

from StringIO                       import StringIO
from Bio.SeqUtils.CheckSum          import seguid
from math                           import log10
from math                           import log
from Bio                            import SeqIO
from Bio.Seq                        import Seq
from Bio.Seq                        import reverse_complement
from Bio.Alphabet.IUPAC             import unambiguous_dna
from Bio.Alphabet.IUPAC             import ambiguous_dna
from Bio.SeqRecord                  import SeqRecord
from Bio.SeqUtils                   import GC
from Bio.SeqUtils.MeltingTemp       import Tm_staluc
from Bio.SeqFeature                 import SeqFeature
from Bio.SeqFeature                 import FeatureLocation
from Bio.SeqFeature                 import ExactPosition

from oligo_melting_temp             import tmbresluc
from formatted_record               import FormattedRecord
from utils                          import double

def annealing_positions(primer, template, limit):
    if len(primer)<limit:
        return []
    head = primer[-limit:]
    positions = [m.start() for m in re.finditer('(?={0})'.format(head), template, re.I)]
    if positions:
        tail = primer[:-limit]
        length = len(tail)
        revtail = tail[::-1]
        results = []
        for match_start in positions:
            tm = template[max(0,match_start-length):match_start][::-1]
            footprint = "".join(reversed([b for a,b in itertools.takewhile(lambda x: x[0].lower()==x[1].lower(),zip(revtail, tm))])) + template[match_start:match_start+limit]
            results.append((match_start+limit-1, footprint, primer[:len(primer)-len(footprint)]))
        return results
    return []


class Amplicon:
    '''
    Amplicon(forward_primer,
             reverse_primer,
             template,
             saltc=50,     saltc = monovalent cations (mM) (Na,K..)
             fc=1000,      primer concentration (nM) 1000nM = 1µM)
             rc=1000):     primer concentration (nM) 1000nM = 1µM)
    '''

    def __init__(self,
                 forward_primer,
                 reverse_primer,
                 template,
                 saltc=50,
                 fc=1000,
                 rc=1000):

        self.forward_primer             = forward_primer
        self.reverse_primer             = reverse_primer
        self.fc                         = fc
        self.rc                         = rc
        self.saltc                      = saltc
        self.template                   = template
        self.product                    = None

    def pcr_product(self):
        if self.product:
            return self.product
        begin=self.forward_primer.pos-len(self.forward_primer.footprint)+1
        end  =self.reverse_primer.pos+len(self.reverse_primer.footprint)-1
        if self.template.circular:
            tmpl=self.template+self.template
        else:
            tmpl=self.template
        if begin<0:
            begin = len(tmpl)/2+begin
        if begin >= end:
            end=end+len(tmpl)/2
        prd = (self.forward_primer.tail +
               tmpl[begin:end] +
               self.reverse_primer.tail.reverse_complement())

        # description = Genbank LOCUS max 16 chars
        prd.name = "{0}bp_PCR_prod".format(len(prd))[:16]
        prd.id = "{0}bp {1}".format( str(len(prd))[:14],seguid(prd.seq) )
        prd.description="Primers {0} {1}".format( self.forward_primer.primer.name,
                                                  self.reverse_primer.primer.name
                                                 )

        prd.features.append(SeqFeature(FeatureLocation(0,len(self.forward_primer.primer)),
                                       type ="primer_bind",strand = 1,
                                       qualifiers = {"label":self.forward_primer.primer.name}))
        prd.features.append(SeqFeature(FeatureLocation(len(prd) - len(self.reverse_primer.primer),len(prd)),
                                       type ="primer_bind",
                                       strand = -1,
                                       qualifiers = {"label":self.reverse_primer.primer.name}))

        self.product=prd

        return self.product

    def flankup(self,flankuplength=50):
        return self.template.seq[self.forward_primer.pos-flankuplength-len(self.forward_primer.footprint):self.forward_primer.pos-len(self.forward_primer.footprint)]

    def flankdn(self,flankdnlength=50):
        return self.template.seq[self.reverse_primer.pos+len(self.reverse_primer.footprint):self.reverse_primer.pos+flankdnlength+len(self.reverse_primer.footprint)]

    def _tm(self):
        # Tm calculations according to SantaLucia 1998
        self.tmf = Tm_staluc(str(self.forward_primer.footprint),dnac=50, saltc=self.saltc)
        self.tmr = Tm_staluc(str(self.reverse_primer.footprint),dnac=50, saltc=self.saltc)
        # Ta calculation for enzymes with dsDNA binding domains (dbd)
        # like Pfu-Sso7d, Phusion or Phire
        # https://www.finnzymes.fi/tm_determination.html
        self.tmf_dbd = tmbresluc(str(self.forward_primer.footprint),primerc=self.fc)
        self.tmr_dbd = tmbresluc(str(self.reverse_primer.footprint),primerc=self.rc)

    def detailed_figure(self):
        self._tm()
        f = '''
            5{fp}3
             {fap:>{fplength}} tm {tmf} (dbd) {tmf_dbd}
            {sp}5{faz}...{raz}3
            {sp}3{fzc}...{rzc}5
             {sp2}{rap} tm {tmr} (dbd) {tmr_dbd}
            {sp2}3{rp}5
            '''.format( fp       = self.forward_primer.primer.seq,
                        fap      = "|"*len(self.forward_primer.footprint),
                        fplength = len(self.forward_primer.primer.seq),
                        tmf      = round(self.tmf,1),
                        tmr      = round(self.tmr,1),
                        tmf_dbd  = round(self.tmf_dbd,1),
                        tmr_dbd  = round(self.tmr_dbd,1),
                        rp       = self.reverse_primer.primer.seq[::-1],
                        rap      = "|"*len(self.reverse_primer.footprint),
                        rplength = len(self.reverse_primer.primer.seq),
                        faz      = self.forward_primer.footprint,
                        raz      = self.reverse_primer.footprint.reverse_complement(),
                        fzc      = self.forward_primer.footprint.complement(),
                        rzc      = self.reverse_primer.footprint[::-1],
                        sp       = " "*(len(self.forward_primer.primer.seq)-len(self.forward_primer.footprint)),
                        sp2      = " "*(3+len(self.forward_primer.primer.seq))
                       )
        return textwrap.dedent(f)


    def pcr_program(self):
        self._tm()
        if not self.product:
            self.pcr_product()
        # Ta calculation according to
        # Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
        # ing temperature for DNA amplification in vitro
        # http://www.ncbi.nlm.nih.gov/pubmed/2003928
        GC_prod=GC(str(self.product.seq))
        tml = min(self.tmf,self.tmr)
        #print GC(str(self.product.seq)), self.saltc/1000.0, len(self.product)
        tmp = 81.5 + 0.41*GC(str(self.product.seq)) + 16.6*log10(self.saltc/1000.0) - 675/len(self.product)
        ta = 0.3*tml+0.7*tmp-14.9
        # Fermentas recombinant taq
        taq_extension_rate = 30  # seconds/kB PCR product length
        extension_time_taq = taq_extension_rate * len(self.product) / 1000 # seconds
        f  = textwrap.dedent(u'''
                                Taq (rate {rate} nt/s)
                                Three-step|         30 cycles     |      |SantaLucia 1998
                                94.0°C    |94.0°C                 |      |SaltC {saltc:2}mM
                                __________|_____          72.0°C  |72.0°C|
                                04min00s  |30s  \         ________|______|
                                          |      \ {ta}°C/{0:2}min{1:2}s|10min |
                                          |       \_____/         |      |
                                          |         30s           |      |4-8°C
                             '''.format(rate    = taq_extension_rate,
                                        ta      = math.ceil(ta),
                                        saltc   = self.saltc,
                                        *divmod(extension_time_taq,60)))

        PfuSso7d_extension_rate = 15 #seconds/kB PCR product
        extension_time_PfuSso7d = PfuSso7d_extension_rate * len(self.product) / 1000  # seconds

        # Ta calculation for enzymes with dsDNA binding domains like Pfu-Sso7d
        # https://www.finnzymes.fi/tm_determination.html

        length_of_f = len(self.forward_primer.footprint)
        length_of_r = len(self.reverse_primer.footprint)

        if (length_of_f>20 and length_of_r>20 and self.tmf_dbd>=69.0 and self.tmr_dbd>=69.0) or (self.tmf_dbd>=72.0 and self.tmr_dbd>=72.0):
            f+=textwrap.dedent( u'''
                                    Pfu-Sso7d (rate {rate}s/kb)
                                    Two-step|    30 cycles |      |Breslauer1986,SantaLucia1998
                                    98.0°C  |98.0C         |      |SaltC {saltc:2}mM
                                    _____ __|_____         |      |Primer1C {fc:3}µM
                                    00min30s|10s  \  72.0°C|72.0°C|Primer2C {rc:3}µM
                                            |      \_______|______|
                                            |      {0:2}min{1:2}s|10min |4-8°C
                                 '''.format(rate = PfuSso7d_extension_rate,
                                            fc = self.fc,
                                            rc = self.rc,
                                            saltc = self.saltc,
                                            *divmod(extension_time_PfuSso7d,60)))
        else:

            if (length_of_f>20 and length_of_r>20):
                ta = min(self.tmf_dbd,self.tmr_dbd)+3
            else:
                ta = min(self.tmf_dbd,self.tmr_dbd)


            f+=textwrap.dedent( u'''
                                    Pfu-Sso7d (rate {rate}s/kb)
                                    Three-step|          30 cycles   |      |Breslauer1986,SantaLucia1998
                                    98.0°C    |98.0°C                |      |SaltC {saltc:2}mM
                                    __________|_____          72.0°C |72.0°C|Primer1C {fc:3}µM
                                    00min30s  |10s  \ {ta}°C ________|______|Primer2C {rc:3}µM
                                              |      \______/{0:2}min{1:2}s|10min |
                                              |        10s           |      |4-8°C
                                 '''.format(rate = PfuSso7d_extension_rate,
                                            ta   = math.ceil(ta),
                                            fc   = self.fc/1000,
                                            rc   = self.rc/1000,
                                            saltc= self.saltc,
                                            *divmod(extension_time_PfuSso7d,60)))
        return f



class Anneal:
    '''
    Anneal(primers,
           template,
           homology_limit=13, max_product_size=15000):

    primers = iterable containing PCR primers as string, Seq, SeqRecord
    or FormattedRecord objects.

    template = string, Seq, SeqRecord or FormattedRecord object.

    returns list of Amplicon objects

    '''

    def __init__(self,
                  primers,
                  template,
                  homology_limit=13,
                  max_product_size=15000):
        self.homology_limit = homology_limit
        self.max_product_size = max_product_size

        primers=[FormattedRecord(p) for p in primers]
        primers=[p for p in primers if p.seq]

        self.template = FormattedRecord(template)

        tm = str(self.template.seq).lower()
        rc = str(self.template.rc().seq).lower()
        length = len(tm)

        primer = collections.namedtuple("primer","primer pos footprint tail")
        self.fwd_primers = []
        self.rev_primers = []
        self.amplicons   = []

        tm    = str(template.seq)
        tm_rc = str(template.seq.reverse_complement())

        if not self.template.circular:
            for p in primers:
                self.fwd_primers.extend((primer(p, pos, Seq(fp), Seq(tl)) for
                                         pos, fp, tl in annealing_positions(
                                                            p.seq.tostring(),
                                                            tm,
                                                            homology_limit)))
                self.rev_primers.extend((primer(p,len(template)-pos, Seq(fp), Seq(tl))
                                         for pos, fp, tl in annealing_positions(
                                                            p.seq.tostring(),
                                                            tm_rc,
                                                            homology_limit)))
        else:
            self.template2=double(template)
            ct = 2*tm
            ct_rc = 2*tm_rc
            for p in primers:
                ann1 = annealing_positions(p.seq.tostring(), tm, homology_limit)
                ann2 = annealing_positions(p.seq.tostring(), ct, homology_limit)
                ann  = set(ann2) - set(ann1)
                ann  = [primer(p, pos%len(template), Seq(fp), Seq(tl))
                        for (pos, fp, tl) in ann]
                self.fwd_primers.extend(ann)

                ann1 = annealing_positions(p.seq.tostring(), tm_rc, homology_limit)
                ann2 = annealing_positions(p.seq.tostring(), ct_rc, homology_limit)
                ann  = set(ann2) - set(ann1)
                ann  = [primer(p, len(template)-pos%len(template), Seq(fp), Seq(tl))
                        for (pos, fp, tl) in ann]
                self.rev_primers.extend(ann)

        for fp in self.fwd_primers:
            for rp in self.rev_primers:
                if (0 < rp.pos-fp.pos < max_product_size):
                    self.amplicons.append(Amplicon(fp,rp, self.template))
                elif self.template.circular and fp.pos>rp.pos and len(template)-fp.pos-rp.pos < max_product_size:
                    self.amplicons.append(Amplicon(fp,rp, self.template)) ### changed here 2
        self.number_of_products = len(self.amplicons)

    def report(self):
        return self.__str__()

    def __str__(self):

        mystring = "Template {name} {size} nt {top}:\n".format(name=self.template.name,
                                                               size=len(self.template),
                                                               top={True:"circular",
                                                                    False:"linear"}[self.template.circular]
                                                               )
        if self.fwd_primers:
            for primer in self.fwd_primers:
                mystring += "Primer "+primer.primer.name
                mystring += " anneals at position "
                mystring += str(primer.pos)
                mystring += "\n"
        else:
            mystring += "No forward primers anneal...\n"
        if self.rev_primers:
            for primer in self.rev_primers:
                mystring += "Primer "+primer.primer.name
                mystring += " anneals reverse at position "
                mystring += str(primer.pos)
                mystring += "\n"
        else:
             mystring += "No reverse primers anneal...\n"
        return mystring

    def products(self):
        return [a.product for a in self.amplicons]

    def featured_template(self):
        template = self.template
        for primer in self.fwd_primers:
            if primer.pos-len(primer.footprint)>0:
                start = primer.pos-len(primer.footprint)
                end = primer.pos
                template.features.append(SeqFeature(FeatureLocation(start,end+1),type ="primer_bind",strand = 1, qualifiers = {"label":fp.name,"ApEinfo_fwdcolor":"green","ApEinfo_revcolor":"red"}))
            else:
                start = len(template)-len(primer.footprint)+primer.pos
                end = start+len(primer.footprint)-len(template)

                suba = SeqFeature( FeatureLocation(start,len(template)),
                                   type ="primer_bind",
                                   strand = 1,
                                   qualifiers = {"label":primer.primer.name,
                                                 "ApEinfo_fwdcolor":"green",
                                                 "ApEinfo_revcolor":"red"} )

                subb = SeqFeature( FeatureLocation(0,end),
                                   type ="primer_bind",
                                   strand = 1,
                                   qualifiers = {"label":primer.primer.name,
                                                 "ApEinfo_fwdcolor":"green",
                                                 "ApEinfo_revcolor":"red"})

                template.features.append(SeqFeature(FeatureLocation(start,end),
                                                    type ="primer_bind",
                                                    sub_features = [suba,subb],
                                                    location_operator= "join",
                                                    strand = 1,
                                                    qualifiers = {"label":primer.primer.name}))
        for primer in self.rev_primers:
            if primer.pos+len(primer.footprint)<=len(template):
                start = primer.pos
                end = primer.pos + len(footprint)
                template.features.append(SeqFeature(FeatureLocation(start-1,end),type ="primer_bind",strand = -1, qualifiers = {"label":rp.name,"ApEinfo_fwdcolor":"green","ApEinfo_revcolor":"red"}))
            else:
                start = primer.pos
                end = primer.pos+len(primer.footprint)-len(template)
                suba = SeqFeature(FeatureLocation(start,len(template)),
                                  type ="primer_bind",
                                  strand = -1,
                                  qualifiers = {"label":primer.primer.name,
                                                "ApEinfo_fwdcolor":"green",
                                                "ApEinfo_revcolor":"red"})
                subb = SeqFeature(FeatureLocation(0,end),
                                  type ="primer_bind",
                                  strand = -1,
                                  qualifiers = {"label":primer.primer.name,
                                                "ApEinfo_fwdcolor":"green",
                                                "ApEinfo_revcolor":"red"})
                template.features.append(SeqFeature(FeatureLocation(start,end),
                                                    type ="primer_bind",
                                                    sub_features = [suba,subb],
                                                    location_operator= "join",
                                                    strand = -1,
                                                    qualifiers = {"label":primer.primer.name}))
        return template


def pcr(*args,**kwargs):
    ''' Convenience function for Anneal'''

    import itertools
    from Bio.SeqRecord import SeqRecord
    ''' flatten args '''
    output = []
    stack = []
    stack.extend(reversed(args))
    while stack:
        top = stack.pop()
        if hasattr(top, "__iter__") and not isinstance(top, SeqRecord):
            stack.extend(reversed(top))
        else:
            output.append(top)
    new=[]

    #print output[2];import sys;sys.exit()

    for s in output:
        if isinstance(s, Seq):
            s = SeqRecord(s)
        elif isinstance(s, SeqRecord):
            pass
        elif hasattr(s, "watson"):
            s=s.watson
        elif isinstance(s, basestring):
            s = SeqRecord(Seq(s))
        else:
            raise TypeError("the record property needs to be a string, a Seq object or a SeqRecord object")
        new.append(s)

    anneal_primers = Anneal(new[:-1],
                            new[-1],
                            **kwargs)
    if anneal_primers:
        if anneal_primers.number_of_products==1:
            return anneal_primers.amplicons.pop().pcr_product()
        elif anneal_primers.number_of_products==0:
            raise Exception("No PCR products! {}".format(anneal_primers.report()))
        else:
            raise Exception("PCR not specific! {}".format(anneal_primers.report()))
    else:
        raise Exception(anneal_primers.report())
    return

if __name__=="__main__":

    from parse import parse

