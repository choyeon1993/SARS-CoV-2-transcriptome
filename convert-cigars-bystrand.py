#!/usr/bin/env python3
#
# Copyright (c) 2020 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

import subprocess as sp
import gzip
import re
import sys

cigar_pat = re.compile('\\d+[A-Z]')

def get_switching_sites(cigar, flag, refstart):
    tokens = cigar_pat.findall(cigar)
    refpos = refstart

    for tok in tokens:
        op = tok[-1]
        length = int(tok[:-1])
        if op == 'N':
            yield flag, refpos, refpos + length
        if op in 'MDN':
            refpos += length


for line in sys.stdin:
    fields = line[:-1].split('\t')
    flag = fields[0]
    if flag == '99' or flag == '147':
        flag = '+'
    if flag == '163' or flag == '83':
        flag = '-'
    refstart = int(fields[1]) - 1
    cigar = fields[2]
    for flag, jfrom, jto in get_switching_sites(cigar, flag, refstart):
        print(flag, jfrom, jto, sep='\t')

