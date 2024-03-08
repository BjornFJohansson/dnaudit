import re
from pathlib import Path
from pydna.parsers import parse
# from pydna.editor import ape
from collections import defaultdict
from itertools import pairwise

# sequencefilesuffixes = ["gb", "gbk", "fa", "fasta"]
# sequencefilepaths = [path for i in sequencefilesuffixes for path in Path('.').rglob(f"*.{i}")]
# sequencefiles = parse(sequencefilepaths)
# id_dict = defaultdict(list)
# for sequencefile in sequencefiles:
#     id_dict[sequencefile.id].append(sequencefile)

textfilesuffixes = ["txt", "md"]
textfilepaths = [path for i in textfilesuffixes for path in Path('.').rglob(f"*.{i}")]

def parse2(text, expect):
    sequences = parse(text)
    if len(sequences) > expect:
        print(expect)
    elif len(sequences) < expect:
        print(expect)
    return sequences

def pcr(text):
    pass
def ligation(text):
    pass
def homologous_recombination(text):
    pass
def crispr(text):
    pass
def restriction_digestion(text):
    pass

def fusion_pcr(text):
    pass
    sequences = parse2(text, 5)

snippet_categories = (pcr, ligation, homologous_recombination, crispr, fusion_pcr, restriction_digestion)
snippet_headers = tuple(f"#.?{'.?'.join(f.__name__.split('_'))}" for f in snippet_categories)

for pth in textfilepaths:
    text = pth.read_text(encoding="utf8")
    # parse snippets
    snippets = re.split(f"({'|'.join(snippet_headers)})", text, re.MULTILINE|re.IGNORECASE)
    snippets = [s for s in snippets if s.strip()]
    for category, content in zip(snippets[::2], snippets[1::2]):
        node = locals()[category.strip("# ").replace(" ", "_").lower()](content)
    break
