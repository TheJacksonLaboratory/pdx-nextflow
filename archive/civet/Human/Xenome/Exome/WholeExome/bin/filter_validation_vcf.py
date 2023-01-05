#! /usr/bin/env python
import sys
import re

# A pattern which will will identify absolute paths ending in bam or vcf, and
# leave just the filename, stripping out the path.
path_pat = '(.*?)(/[/a-zA-Z0-9-_.]*?([a-zA-Z0-9-_.]*(bam|vcf)))(.*)'
ppat = re.compile(path_pat)
date_epoch_pat = 'Date=".*,Epoch=.*,'
dpat = re.compile(date_epoch_pat)
for line in open(sys.argv[1]):
    line = line.rstrip()
    line = re.sub(ppat, r'\1\3\5', line)
    line = re.sub(dpat, 'DateEpoch,', line)
    print line
