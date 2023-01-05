#! /usr/bin/env python

# determin_sex.py
# By looking at gene KDM5D, particularly two of its exons, which
# is contained in the human Y chromosome, determine the sample's
# sex.
#
# PARAMETERS:
#    1) total reads for the sample
#    2) reads hitting in exon 1
#    3) reads hitting in exon 5
#    4) reads hitting in any exon in the gene
#    5) output directory
#
# OPTIONS:
#    -v, --version: output the program version and exit; 
#                   ignote all other parameters.
#
# Eventually, we may find that the signal is strong enough to
# make the male values totally disjoint from the female values.
# For now, however, we have a low threshhold indicating female
# and a high threshhold indicating male, and an "undetermined"
# region in between.

import sys
import os
import inspect


cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
lib_folder = os.path.join(cmd_folder, '../lib')
if lib_folder not in sys.path:
     sys.path.insert(0, lib_folder)

import cga_version as version

def usage():
    if len(sys.argv) != 6:
        print >> sys.stderr, 'USAGE:', sys.argv[0], 'total-reads exon1-reads exon5-reads gene-reads output-dir'
        sys.exit(1)

def main():
    # Check for a version request.
    version.parse_options()
    usage()

    points = [{'name':'exon 1', 'low':5, 'high':30},
              {'name':'exon 5', 'low':1, 'high':8},
              {'name':'whole gene', 'low':50, 'high':500} ]

    FEMALE = 0
    MALE = 1
    UNDETERMINED = 2
    NO_CONSENSUS = 3

    MIL = 1000000.0

    sex_text = ['Female', 'Male', 'Undetermined', 'No Consensus']

    all_reads = float(sys.argv[1])

    data = []
    data.append(int(sys.argv[2]))
    data.append(int(sys.argv[3]))
    data.append(int(sys.argv[4]))

    out_dir = sys.argv[5]
    
    ratio = []
    for d in data:
        ratio.append(d * MIL / all_reads)

    sex = []
    for n in range(len(ratio)):
        if ratio[n] <= points[n]['low']:
            sex.append(FEMALE)
        elif ratio[n] >= points[n]['high']:
            sex.append(MALE)
        else:
            sex.append(UNDETERMINED)

    if sex[0] == sex[1] and sex[1] == sex[2]:
        consensus_sex = sex[0]
    else:
        consensus_sex = NO_CONSENSUS

    of = open(os.path.join(out_dir, 'sex_determination.txt'), 'w')
    of.write('SEX: {0}\n'.format(sex_text[consensus_sex]))
    of.write('Breakdown:    Reads        RPM  Sex\n')
    #sum = open(os.path.join(out_dir, 'sex_summary.txt'), 'a')
    #sum.write('{0}\t'.format(all_reads))
    for n in range(len(sex)):
        of.write('  {0:10s} {1:6d}{2:11.3f}  {3}\n'.format(points[n]['name'], data[n], ratio[n], sex_text[sex[n]]))
        #sum.write('{0}\t{1}\t'.format(data[n], ratio[n]))
    #sum.write('{0}\n'.format(sex_text[consensus_sex]))
    of.close()
    #sum.close()

main()