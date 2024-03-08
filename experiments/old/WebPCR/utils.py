#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
docstring
'''

def eq(*args,**kwargs):
    return equivalent(*args,**kwargs)

def equivalent(*args,**kwargs):
    from Bio.Seq import reverse_complement
    from Bio.SeqRecord import SeqRecord
    import itertools
    args=list(args)
    for i, arg in enumerate(args):
        if not hasattr(arg, "__iter__") or isinstance(arg,SeqRecord):
            args[i] = (arg,)
    args = list(itertools.chain.from_iterable(args))
    topology = None
    if "linear" in kwargs:
        if kwargs["linear"]==True:
            topology = "linear"
        if kwargs["linear"]==False:
            topology = "circular"
    elif "circular" in kwargs:
        if kwargs["circular"]==True:
            topology = "circular"
        if kwargs["circular"]==False:
            topology = "linear"
    else:
        # topology keyword not set, look for topology associated to each sequence
        # otherwise raise exception
        topology = set([arg.circular if hasattr(arg, "circular") else None for arg in args])

        if len(topology)!=1:
            raise Exception("sequences have different topologies")
        topology = topology.pop()
        if topology==False:
            topology = "linear"
        elif topology==True:
            topology = "circular"

    args_string_list    = [str(arg.seq).lower() if hasattr(arg,"seq") else str(arg).lower() for arg in args]
    length = set((len(s) for s in args_string_list))
    if len(length)!=1:
        return False
    same = True
    if topology == "circular":
        # force circular comparison of all given sequences
            for s1,s2 in itertools.combinations(args_string_list, 2):
                if not ( s1 in s2+s2 or reverse_complement(s1) in s2+s2):
                    same = False
    elif topology == "linear":
        # force linear comparison of all given sequences
        for s1,s2 in itertools.combinations(args_string_list, 2):
            if not ( s1==s2 or s1==reverse_complement(s2) ):
                same = False
    return same


def double(seq):
    #from pydna import ape
    from Bio.SeqFeature import SeqFeature,FeatureLocation
    from Bio.SeqRecord import SeqRecord
    double = seq+seq #SeqRecord(seq.seq)
    nf=[]
    #ape(double)
    for f in double.features:
        if f.location_operator=="join" and f.sub_features:
            for sf1,sf2 in zip(f.sub_features[:-1], f.sub_features[1:]):

                newstart = sf1.location.end+sf2.location.end
                newend   = sf1.location.start+sf2.location.start

                if  ( sf1.location.end   == len(seq),
                      sf2.location.start == 0,
                      newstart - newend  == len(f) ):

                    newstart, newend = sorted((newstart, newend))
                    nf.append( SeqFeature(FeatureLocation(newstart, newend),
                               type=f.type,
                               location_operator='',
                               strand=sf1.strand,
                               id=sf1.id,
                               qualifiers=f.qualifiers,
                               sub_features=sf1.sub_features+sf2.sub_features,
                               ref=None,
                               ref_db=None
                                              )
                                  )
                else:
                    pass

        else:
            nf.append(f)

    double.features=nf
    return double

def shift_origin(seq, shift):
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    from Bio.SeqRecord  import SeqRecord
    from pydna import ape
    length=len(seq)
    assert 0<=shift<length
    new=double(seq)[shift:shift+length]
    def wraparound(feature):
        new_start = length -(shift-feature.location.start)
        new_end   = feature.location.end-shift
        a = SeqFeature(FeatureLocation(0, new_end),
                       type=feature.type,
                       location_operator=feature.location_operator,
                       strand=feature.strand,
                       id=feature.id,
                       qualifiers=feature.qualifiers,
                       sub_features=None)
        b = SeqFeature(FeatureLocation(new_start, length),
                       type=feature.type,
                       location_operator=feature.location_operator,
                       strand=feature.strand,
                       id=feature.id,
                       qualifiers=feature.qualifiers,
                       sub_features=None)
        c = SeqFeature(FeatureLocation(new_start, new_end),
                       type=feature.type,
                       location_operator="join",
                       strand=feature.strand,
                       id=feature.id,
                       qualifiers=feature.qualifiers,
                       sub_features=[a,b])
        sub_features=[]
        for sf in feature.sub_features:
            if feature.location.end<shift:
                sub_features.append(SeqFeature(FeatureLocation(length-feature.location.start,
                                                               length-feature.location.end),
                                    type=feature.type,
                                    location_operator=feature.location_operator,
                                    strand=feature.strand,
                                    id=feature.id,
                                    qualifiers=feature.qualifiers,
                                    sub_features=None))
            elif feature.location.start>shift:
                sub_features.append(SeqFeature(FeatureLocation(feature.location.start-shift,
                                                               feature.location.end-shift),
                                    type=feature.type,
                                    location_operator=feature.location_operator,
                                    strand=feature.strand,
                                    id=feature.id,
                                    qualifiers=feature.qualifiers,
                                     sub_features=None))
            else:
                sub_features.extend(wraparound(sf))
        c.sub_features.extend(sub_features)
        return c

    for feature in seq.features:
        if shift in feature:
            new.features.append(wraparound(feature))
    return new


def sync(sequence, ref):
    import itertools
    from Bio.Seq import reverse_complement, Seq
    from Bio.SeqRecord import SeqRecord
    from pydna.find_sub_strings import find_sub_strings as fs
    from pydna import eq, dsdna
    from utils import double, shift_origin

    original = sequence

    if hasattr(sequence, "seq"):
        a               = str(sequence.seq.lower())
        a_rc            = str(sequence.seq.reverse_complement()).lower()
        sequence_rc     = sequence.reverse_complement()
        double_sequence = double(sequence)
    elif hasattr(sequence, "watson"):
        sequence = sequence.watson
        a    = str(sequence.seq).lower()
        a_rc = str(sequence.seq.reverse_complement()).lower()
        sequence_rc     = sequence.reverse_complement()
        double_sequence = double(sequence)
    else:
        a    = str(sequence).lower()
        a_rc = str(reverse_complement(sequence)).lower()
        sequence_rc = reverse_complement(sequence)
        double_sequence = a+a

    if hasattr(ref, "seq"):
        b    = str(ref.seq).lower()
    elif hasattr(ref, "watson"):
        b    = str(ref.watson.seq).lower()
    else:
        b    = str(ref.lower())

    c = fs.common_sub_strings(a+a,b)
    d = fs.common_sub_strings(a_rc+a_rc,b)

    if c:
        starta, startb, length = c.pop(0)
    else:
        starta, startb, length = 0,0,0

    if d:
        starta_rc, startb_rc, length_rc = d.pop(0)
    else:
        starta_rc, startb_rc, length_rc = 0,0,0

    if not c and not d:
        raise Exception("no overlap !")

    if length_rc>length:
        starta, startb = starta_rc, startb_rc
        sequence = sequence_rc

    if starta>startb:
        if len(a)<len(b):
            ofs = starta-startb + len(b)-len(a)
        else:
            ofs = starta-startb
    elif starta<startb:
        ofs = startb-starta + len(a)-len(b)
        ofs = len(a)-ofs
    elif starta==startb:
        ofs=0

    sequence = shift_origin(sequence, ofs)
    sequence.circular=True
    if hasattr(original,"watson"):
        assert eq(original.watson, sequence, circular =True)
        return dsdna(sequence)
    assert eq(original, sequence, circular =True)
    return sequence

def test_sync():
    from Bio.Seq import reverse_complement

    pcaps = read("/home/bjorn/Dropbox/wikidata/pCAPs.wiki")
    pcaps2 = sync(pcaps, "tctgacacatgcagctcccggagacggtcac")
    assert   str(pcaps2.seq).lower().startswith("tctgacacatgcagctcccggagacggtcac")
    pcaps2 = sync(pcaps, "agaaaccattattatcatgacattaacctataaaaa")
    assert   str(pcaps2.seq).lower().startswith("agaaaccattattatcatgacattaacctataaaaa")
    pcaps2 = sync(pcaps, "tctagacaaaccgtgggacgaattcttaag")
    pcaps3 = sync(pcaps2,"tcgcgcgtttcggtgatgacggtgaaaacc")
    assert  eq(pcaps,pcaps2,pcaps3)
    assert  eq(pcaps,pcaps3, linear=True)

    pcaps4 ="""
LOCUS       3128bp__1               3128 bp ds-DNA     linear       25-JUL-2012
DEFINITION  Primers lowgc_f lowgc_r
ACCESSION   3128bp 7pPxy/bQvs4+7CaOgiywQVzUFDc
VERSION     3128bp 7pPxy/bQvs4+7CaOgiywQVzUFDc
KEYWORDS    .
SOURCE      .
  ORGANISM  . .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     rep_origin      757..757
                     /direction="BOTH"
                     /label=rep_origin
                     /ApEinfo_fwdcolor=pink
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     gene            complement(1516..2376)
                     /gene="bla"
                     /label=bla
                     /ApEinfo_fwdcolor=pink
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     CDS             complement(1516..2376)
                     /product="beta-lactamase"
                     /codon_start="1"
                     /transl_table="11"
                     /db_xref="GI:2769263"
                     /db_xref="GOA:Q79DR3"
                     /db_xref="HSSP:P62593"
                     /db_xref="InterPro:IPR000871"
                     /db_xref="InterPro:IPR001466"
                     /db_xref="InterPro:IPR012338"
                     /db_xref="UniProtKB/TrEMBL:Q79DR3"
                     /translation="MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGY
                     I
                     ELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEY
                     S
                     PVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDR
                     W
                     EPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRS
                     A
                     LPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGA
                     S LIKHW"
                     /gene="bla"
                     /protein_id="CAA04868.1"
                     /label=beta-lactamase
                     /ApEinfo_fwdcolor=pink
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     primer_bind     1..57
                     /label=lowgc_f
                     /ApEinfo_fwdcolor=pink
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     primer_bind     complement(3083..3128)
                     /label=lowgc_r
                     /ApEinfo_fwdcolor=pink
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
        1 tttcactagt tacttgtagt cgacgtgcca tctgtgcaga caaacgcatc aggatatccg
       61 gatttacctg aatcaattgg cgaaattttt tgtacgaaat ttcagccact tcacaggcgg
      121 ttttcgcacg tacccatgcg ctacgttcct ggccctcttc aaacaggccc agttcgccaa
      181 taaaatcacc ctgattcaga taggagagga tcatttcttt accctcttcg tctttgatca
      241 gcactgccac agagccttta acgatgtagt acagcgtttc cgctttttca ccctggtgaa
      301 taagcgtgct cttggatggg tacttatgaa tgtggcaatg agacaagaac cattcgagag
      361 taggatccgt ttgaggttta ccaagtacca taagatcctt aaatttttat tatctagcta
      421 gatgataata ttatatcaag aattgtacct gaaagcaaat aaatttttta tctggcttaa
      481 ctatgcggca tcagagcaga ttgtactgag agtgcaccat atgcggtgtg aaataccgca
      541 cagatgcgta aggagaaaat accgcatcag gcgctcttcc gcttcctcgc tcactgactc
      601 gctgcgctcg gtcgttcggc tgcggcgagc ggtatcagct cactcaaagg cggtaatacg
      661 gttatccaca gaatcagggg ataacgcagg aaagaacatg tgagcaaaag gccagcaaaa
      721 ggccaggaac cgtaaaaagg ccgcgttgct ggcgtttttc cataggctcc gcccccctga
      781 cgagcatcac aaaaatcgac gctcaagtca gaggtggcga aacccgacag gactataaag
      841 ataccaggcg tttccccctg gaagctccct cgtgcgctct cctgttccga ccctgccgct
      901 taccggatac ctgtccgcct ttctcccttc gggaagcgtg gcgctttctc atagctcacg
      961 ctgtaggtat ctcagttcgg tgtaggtcgt tcgctccaag ctgggctgtg tgcacgaacc
     1021 ccccgttcag cccgaccgct gcgccttatc cggtaactat cgtcttgagt ccaacccggt
     1081 aagacacgac ttatcgccac tggcagcagc cactggtaac aggattagca gagcgaggta
     1141 tgtaggcggt gctacagagt tcttgaagtg gtggcctaac tacggctaca ctagaaggac
     1201 agtatttggt atctgcgctc tgctgaagcc agttaccttc ggaaaaagag ttggtagctc
     1261 ttgatccggc aaacaaacca ccgctggtag cggtggtttt tttgtttgca agcagcagat
     1321 tacgcgcaga aaaaaaggat ctcaagaaga tcctttgatc ttttctacgg ggtctgacgc
     1381 tcagtggaac gaaaactcac gttaagggat tttggtcatg agattatcaa aaaggatctt
     1441 cacctagatc cttttaaatt aaaaatgaag ttttaaatca atctaaagta tatatgagta
     1501 aacttggtct gacagttacc aatgcttaat cagtgaggca cctatctcag cgatctgtct
     1561 atttcgttca tccatagttg cctgactccc cgtcgtgtag ataactacga tacgggaggg
     1621 cttaccatct ggccccagtg ctgcaatgat accgcgagac ccacgctcac cggctccaga
     1681 tttatcagca ataaaccagc cagccggaag ggccgagcgc agaagtggtc ctgcaacttt
     1741 atccgcctcc atccagtcta ttaattgttg ccgggaagct agagtaagta gttcgccagt
     1801 taatagtttg cgcaacgttg ttgccattgc tacaggcatc gtggtgtcac gctcgtcgtt
     1861 tggtatggct tcattcagct ccggttccca acgatcaagg cgagttacat gatcccccat
     1921 gttgtgcaaa aaagcggtta gctccttcgg tcctccgatc gttgtcagaa gtaagttggc
     1981 cgcagtgtta tcactcatgg ttatggcagc actgcataat tctcttactg tcatgccatc
     2041 cgtaagatgc ttttctgtga ctggtgagta ctcaaccaag tcattctgag aatagtgtat
     2101 gcggcgaccg agttgctctt gcccggcgtc aatacgggat aataccgcgc cacatagcag
     2161 aactttaaaa gtgctcatca ttggaaaacg ttcttcgggg cgaaaactct caaggatctt
     2221 accgctgttg agatccagtt cgatgtaacc cactcgtgca cccaactgat cttcagcatc
     2281 ttttactttc accagcgttt ctgggtgagc aaaaacagga aggcaaaatg ccgcaaaaaa
     2341 gggaataagg gcgacacgga aatgttgaat actcatactc ttcctttttc aatattattg
     2401 aagcatttat cagggttatt gtctcatgag cggatacata tttgaatgta tttagaaaaa
     2461 taaacaaata ggggttccgc gcacatttcc ccgaaaagtg ccacctgcta agaaaccatt
     2521 attatcatga cattaaccta taaaaatagg cgtatcacga ggccctttcg tctcgcgcgt
     2581 ttcggtgatg acggtgaaaa cctctgacac atgcagctcc cggagacggt cacagcttgt
     2641 ctgtaagcgg atgccgggag cagacaagcc cgtcagggcg cgtcagcggg tgttggcggg
     2701 tgtcggggct ggcttaacta tgcggcatca gagcagattg tactgagagt gcaccataga
     2761 tcctgaggat cggggtgata aatcagtctg cgccacatcg ggggaaacaa aatggcgcga
     2821 gatctaaaaa aaaaggctcc aaaaggagcc tttcgcgcta ccaggtaacg cgccactccg
     2881 acgggattaa cgagtgccgt aaacgacgat ggttttaccg tgtgcggaga tcaggttctg
     2941 atcctcgagc atcttaagaa ttcgtcccac ggtttgtcta gagcagccga caatctggcc
     3001 aatttcctga cgggtaattt tgatttgcat gccgtccggg tgagtcatag cgtctggttg
     3061 ttttgccaga ttcagcagag tctgtgcaat gcggccgctg acgtcaagaa aagctagatt
     3121 tccaactt
//
    """


    plasmid = read(pcaps4)

    plasmid2 = sync(plasmid, pcaps)

    print plasmid2[:10].seq

    print "done!!"



def test_shift():

    a='''
LOCUS       3128bp__1               3128 bp ds-DNA     linear       24-JUL-2012
DEFINITION  Primers lowgc_f lowgc_r
ACCESSION   3128bp 7pPxy/bQvs4+7CaOgiywQVzUFDc
VERSION     3128bp 7pPxy/bQvs4+7CaOgiywQVzUFDc
KEYWORDS    .
SOURCE      .
  ORGANISM  . .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     rep_origin      757..757
                     /direction="BOTH"
                     /label=rep_origin
                     /ApEinfo_fwdcolor=pink
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     gene            complement(1516..2376)
                     /gene="bla"
                     /label=bla
                     /ApEinfo_fwdcolor=#36c04b
                     /ApEinfo_revcolor=pink
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     CDS             complement(1516..2376)
                     /product="beta-lactamase"
                     /codon_start="1"
                     /transl_table="11"
                     /db_xref="GI:2769263"
                     /db_xref="GOA:Q79DR3"
                     /db_xref="HSSP:P62593"
                     /db_xref="InterPro:IPR000871"
                     /db_xref="InterPro:IPR001466"
                     /db_xref="InterPro:IPR012338"
                     /db_xref="UniProtKB/TrEMBL:Q79DR3"
                     /translation="MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGY
                     I
                     ELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEY
                     S
                     PVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDR
                     W
                     EPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRS
                     A
                     LPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGA
                     S LIKHW"
                     /gene="bla"
                     /protein_id="CAA04868.1"
                     /label=beta-lactamase
                     /ApEinfo_fwdcolor=#3ec01b
                     /ApEinfo_revcolor=#ffe9cb
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     primer_bind     1..57
                     /label=lowgc_f
                     /ApEinfo_fwdcolor=#000eff
                     /ApEinfo_revcolor=#e3ff00
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 1}
                     width 5 offset 0
     primer_bind     complement(3083..3128)
                     /label=lowgc_r
                     /ApEinfo_fwdcolor=#000eff
                     /ApEinfo_revcolor=#e3ff00
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 1}
                     width 5 offset 0
ORIGIN
        1 tttcactagt tacttgtagt cgacgtgcca tctgtgcaga caaacgcatc aggatatccg
       61 gatttacctg aatcaattgg cgaaattttt tgtacgaaat ttcagccact tcacaggcgg
      121 ttttcgcacg tacccatgcg ctacgttcct ggccctcttc aaacaggccc agttcgccaa
      181 taaaatcacc ctgattcaga taggagagga tcatttcttt accctcttcg tctttgatca
      241 gcactgccac agagccttta acgatgtagt acagcgtttc cgctttttca ccctggtgaa
      301 taagcgtgct cttggatggg tacttatgaa tgtggcaatg agacaagaac cattcgagag
      361 taggatccgt ttgaggttta ccaagtacca taagatcctt aaatttttat tatctagcta
      421 gatgataata ttatatcaag aattgtacct gaaagcaaat aaatttttta tctggcttaa
      481 ctatgcggca tcagagcaga ttgtactgag agtgcaccat atgcggtgtg aaataccgca
      541 cagatgcgta aggagaaaat accgcatcag gcgctcttcc gcttcctcgc tcactgactc
      601 gctgcgctcg gtcgttcggc tgcggcgagc ggtatcagct cactcaaagg cggtaatacg
      661 gttatccaca gaatcagggg ataacgcagg aaagaacatg tgagcaaaag gccagcaaaa
      721 ggccaggaac cgtaaaaagg ccgcgttgct ggcgtttttc cataggctcc gcccccctga
      781 cgagcatcac aaaaatcgac gctcaagtca gaggtggcga aacccgacag gactataaag
      841 ataccaggcg tttccccctg gaagctccct cgtgcgctct cctgttccga ccctgccgct
      901 taccggatac ctgtccgcct ttctcccttc gggaagcgtg gcgctttctc atagctcacg
      961 ctgtaggtat ctcagttcgg tgtaggtcgt tcgctccaag ctgggctgtg tgcacgaacc
     1021 ccccgttcag cccgaccgct gcgccttatc cggtaactat cgtcttgagt ccaacccggt
     1081 aagacacgac ttatcgccac tggcagcagc cactggtaac aggattagca gagcgaggta
     1141 tgtaggcggt gctacagagt tcttgaagtg gtggcctaac tacggctaca ctagaaggac
     1201 agtatttggt atctgcgctc tgctgaagcc agttaccttc ggaaaaagag ttggtagctc
     1261 ttgatccggc aaacaaacca ccgctggtag cggtggtttt tttgtttgca agcagcagat
     1321 tacgcgcaga aaaaaaggat ctcaagaaga tcctttgatc ttttctacgg ggtctgacgc
     1381 tcagtggaac gaaaactcac gttaagggat tttggtcatg agattatcaa aaaggatctt
     1441 cacctagatc cttttaaatt aaaaatgaag ttttaaatca atctaaagta tatatgagta
     1501 aacttggtct gacagttacc aatgcttaat cagtgaggca cctatctcag cgatctgtct
     1561 atttcgttca tccatagttg cctgactccc cgtcgtgtag ataactacga tacgggaggg
     1621 cttaccatct ggccccagtg ctgcaatgat accgcgagac ccacgctcac cggctccaga
     1681 tttatcagca ataaaccagc cagccggaag ggccgagcgc agaagtggtc ctgcaacttt
     1741 atccgcctcc atccagtcta ttaattgttg ccgggaagct agagtaagta gttcgccagt
     1801 taatagtttg cgcaacgttg ttgccattgc tacaggcatc gtggtgtcac gctcgtcgtt
     1861 tggtatggct tcattcagct ccggttccca acgatcaagg cgagttacat gatcccccat
     1921 gttgtgcaaa aaagcggtta gctccttcgg tcctccgatc gttgtcagaa gtaagttggc
     1981 cgcagtgtta tcactcatgg ttatggcagc actgcataat tctcttactg tcatgccatc
     2041 cgtaagatgc ttttctgtga ctggtgagta ctcaaccaag tcattctgag aatagtgtat
     2101 gcggcgaccg agttgctctt gcccggcgtc aatacgggat aataccgcgc cacatagcag
     2161 aactttaaaa gtgctcatca ttggaaaacg ttcttcgggg cgaaaactct caaggatctt
     2221 accgctgttg agatccagtt cgatgtaacc cactcgtgca cccaactgat cttcagcatc
     2281 ttttactttc accagcgttt ctgggtgagc aaaaacagga aggcaaaatg ccgcaaaaaa
     2341 gggaataagg gcgacacgga aatgttgaat actcatactc ttcctttttc aatattattg
     2401 aagcatttat cagggttatt gtctcatgag cggatacata tttgaatgta tttagaaaaa
     2461 taaacaaata ggggttccgc gcacatttcc ccgaaaagtg ccacctgcta agaaaccatt
     2521 attatcatga cattaaccta taaaaatagg cgtatcacga ggccctttcg tctcgcgcgt
     2581 ttcggtgatg acggtgaaaa cctctgacac atgcagctcc cggagacggt cacagcttgt
     2641 ctgtaagcgg atgccgggag cagacaagcc cgtcagggcg cgtcagcggg tgttggcggg
     2701 tgtcggggct ggcttaacta tgcggcatca gagcagattg tactgagagt gcaccataga
     2761 tcctgaggat cggggtgata aatcagtctg cgccacatcg ggggaaacaa aatggcgcga
     2821 gatctaaaaa aaaaggctcc aaaaggagcc tttcgcgcta ccaggtaacg cgccactccg
     2881 acgggattaa cgagtgccgt aaacgacgat ggttttaccg tgtgcggaga tcaggttctg
     2941 atcctcgagc atcttaagaa ttcgtcccac ggtttgtcta gagcagccga caatctggcc
     3001 aatttcctga cgggtaattt tgatttgcat gccgtccggg tgagtcatag cgtctggttg
     3061 ttttgccaga ttcagcagag tctgtgcaat gcggccgctg acgtcaagaa aagctagatt
     3121 tccaactt
//
    '''

    a=read(a)
    ape(a)
    b = shift_origin(a,2572)
    assert eq(a,b,circular=True)
    #print b.format("gb")
    ape(b)
    import sys; sys.exit()

    a='''
LOCUS       New_DNA                   29 bp ds-DNA     linear       24-JUL-2012
DEFINITION  .
SOURCE      .
  ORGANISM  .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    join(27..29,1..7)
                     /label=MyFeature
ORIGIN
        1 tctactgagc gtatcccagc tacgtacta
//
    '''
    #a=read(a)
    #b = shift_origin(a,17)
    #print b.format("gb")

def test_eq():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    assert equivalent( "AAA" ,"TTT", linear   = True ) == True
    assert equivalent( "AAA" ,"TTT", linear   = False) == True

    assert equivalent( "ATA" ,"AAT", circular = True ) == True
    assert equivalent( "ATA" ,"AAT", circular = False) == False
    assert equivalent( "AAA" ,"AAA", linear   = True ) == True
    assert equivalent( "AAA" ,"AAA", linear   = False) == True

    assert equivalent( "ATA" ,Seq("AAT"), circular = True ) == True
    assert equivalent( "ATA" ,Seq("AAT"), circular = False) == False
    assert equivalent( "AAA" ,Seq("AAA"), linear   = True ) == True
    assert equivalent( "AAA" ,Seq("AAA"), linear   = False) == True

    assert equivalent( "ATA" ,SeqRecord("AAT"), circular = True ) == True
    assert equivalent( "ATA" ,SeqRecord("AAT"), circular = False) == False
    assert equivalent( "AAA" ,SeqRecord("AAA"), linear   = True ) == True
    assert equivalent( "AAA" ,SeqRecord("AAA"), linear   = False) == True

    assert equivalent( "ATA" ,FormattedRecord("AAT"), circular = True ) == True
    assert equivalent( "ATA" ,FormattedRecord("AAT"), circular = False) == False
    assert equivalent( "AAA" ,FormattedRecord("AAA"), linear   = True ) == True
    assert equivalent( "AAA" ,FormattedRecord("AAA"), linear   = False) == True

    assert equivalent( Seq("ATA") ,SeqRecord("AAT"), circular = True ) == True
    assert equivalent( Seq("ATA") ,SeqRecord("AAT"), circular = False) == False
    assert equivalent( Seq("AAA") ,SeqRecord("AAA"), linear   = True ) == True
    assert equivalent( Seq("AAA") ,SeqRecord("AAA"), linear   = False) == True

    assert equivalent( Seq("ATA") ,FormattedRecord("AAT"), circular = True ) == True
    assert equivalent( Seq("ATA") ,FormattedRecord("AAT"), circular = False) == False
    assert equivalent( Seq("AAA") ,FormattedRecord("AAA"), linear   = True ) == True
    assert equivalent( Seq("AAA") ,FormattedRecord("AAA"), linear   = False) == True

    assert equivalent( FormattedRecord("AAA",circular=False) ,FormattedRecord("AAA",circular=False))
    assert equivalent( FormattedRecord("AAA",circular=True)  ,FormattedRecord("AAA",circular=True))
    assert not equivalent( FormattedRecord("ATA",circular=False) ,FormattedRecord("AAT",circular=False))
    assert equivalent( FormattedRecord("ATA",circular=True)  ,FormattedRecord("AAT",circular=True))



if __name__=="__main__":
    from pydna import read, ape
    test_sync()