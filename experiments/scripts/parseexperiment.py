#!/usr/bin/env python3
# -*- coding: utf-8 -*-

raw = """\
>aaa
aaa
aa

EcoRV > LOCUS
>bbb
bbb
bbb

LOCUS
xxx
//
>ccc
cccc
LOCUS
yyy
//
LOCUS
zzz
//
>dddd
ddddd"""

import re


rx = r"(^>(?:.(?!(?:^$|^>|LOCUS)))+)|(^LOCUS.+?//)|(^(?!(?:LOCUS|>))(?:.(?!(?:^$|^>|LOCUS)))+)"

rx = r"((?:^>(?:.(?!(?:^$|^>|LOCUS)))+)|(?:^LOCUS.+?//)|(?:^(?!(?:LOCUS|>))(?:.(?!(?:^$|^>|LOCUS)))+))"

re.findall(rx, raw, flags=re.MULTILINE|re.DOTALL)
