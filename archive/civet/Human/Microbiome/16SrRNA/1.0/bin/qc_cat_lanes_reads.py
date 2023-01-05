#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
check for existing files matching pattern, cat all matching files into one
"""

from __future__ import print_function

import os
import sys
import argparse
import shutil
import gzip
import re
from fnmatch import fnmatchcase, filter

DEBUG=True

# default all (8) lanes and both (2) reads
# code is only safe for single digit (<=9) lanes and reads
LANES=12345678
READS=12
dest_file_ext = '.fastq'

def debug_print(str, debug=DEBUG, end='\n'):
    if debug:
        print(str, end=end, file=sys.stderr)

def str_to_list_with_prefix(in_string, prefix=''):
    str_list = []
    if len(in_string) != 0:
        for X in in_string:
            str_list.append(prefix + X)
    return str_list

def generate_lane_read_patterns(dirname='',
        lanes=LANES, reads=READS,
        lane_prefix='L00', read_prefix='R'):
    """pass number of Reads to generate filename patterns to match,
    with Lane numbers for Mi/HiSeq fastq names, or
    Sample numbers for NextSeq fastq names"""
    lanes = str_to_list_with_prefix(str(lanes), prefix=lane_prefix)
    if len(lanes):
        debug_print("lanes? " + str(lanes))
    else:
        lanes = ''
    reads = str_to_list_with_prefix(str(reads), prefix=read_prefix)
    D = re.sub('[- \.]', '_', dirname)

    template_key = '%s_%s_%s' # % (D,L,R)
    template_pattern = '*[\._]%s[\._]%s[\._]*fastq*' # % (L,R)

    #lrpd = dir,lane,read pattern dict
    lrpd = {}

    for L in lanes:
        L_norm = re.sub('[- \._*?]','', L)
        for R in reads:
            k = template_key % (D,L_norm,R)
            patt = template_pattern % (L,R)
            lrpd[k] = patt
            #debug_print('lrpd: %s = %s' % (k, patt))
    return lrpd

def match_files_by_pattern(pattern, dir):
    """checks for existing files matching lane and/or read patterns in dir"""
    debug_print('Checking for file pattern: %s' % (pattern), end='')
    match_files = filter(os.listdir(dir), pattern)
    debug_print(" -> %s found." % len(match_files))
    return match_files

def cat_files(filelist, dest_file):
    """cat all files matching src_pattern into dest_file"""
    try:
        dest = open(dest_file, 'w')
        debug_print(" -> Cat'ing files into '%s'" % (dest_file))
        for F in filelist:
            debug_print('file: '+F)
            if F.endswith('gz'):
                fh = gzip.open(F)
            else:
                fh = open(F, 'r')
            shutil.copyfileobj(fh, dest)
    except Exception as e:
        raise e

def rm_files(filename):
    """remove any matching cat'd files from dir"""
    if os.path.exists(filename):
        debug_print('Removing file: %s' % (filename))
        os.remove(filename)
    #else:
        #debug_print(" -> '" +filename+ "' file not found!")

def process_files(dir, patterns, remove):
    """check listdir for filenames matching pattern,
    cat or remove them, return number of matches"""
    num_matches = 0 # initiate zero var for future lane-less recheck
    for key,pattern in sorted(patterns.items()):
        dest_file = key + dest_file_ext
        match_files = match_files_by_pattern(pattern, dir)
        debug_print('match_files: ' + str(match_files[0:]))
        num_matches = len(match_files)
        if num_matches > 0:
            if remove == True:
                rm_files(dest_file)
            else:
                cat_files(match_files,dest_file)
    return num_matches

def parse_args(args=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', dest='dir', default='.',
                        help='directory to find source files. '+
                        '(default is current directory)')
    parser.add_argument('-l', dest='lanes', default=LANES, type=str,
                        help='list of lanes to process \'348\' (default '
                        +str(LANES)+ '). If NO lane numbers in the files,'
                        ' specify this as -l \'\'')
    parser.add_argument('-r', dest='reads', default=READS, type=str,
                        help='list of reads to process \'1\' (default '
                        +str(READS)+ ')')
    parser.add_argument('-rm', dest='remove', action='store_true',
                        default='False',
                        help='remove all matching fastq files in "-d"'
                        '(default False)')
    return parser.parse_args()

def make_it_happen(args):
    debug_print('Starting run of qc_cat_lanes_reads')
    # args....
    lanes = args.lanes
    reads = args.reads
    remove = args.remove

    dir = os.path.abspath(args.dir)
    dirname = os.path.basename(dir)
    # calldir = os.path.abspath(os.curdir)


    os.chdir(dir)
    lr_patts = generate_lane_read_patterns(dirname, lanes, reads)
    debug_print( "lr_patterns generated: " + str(len(lr_patts.keys())) )
    num_matches = process_files(dir, lr_patts, remove)
    if num_matches == 0:
        # check for reads without lanes (if '--no-lane-splitting' is specified)
        lr_patts = generate_lane_read_patterns(dirname, lanes='*', reads=reads,
                lane_prefix='S') # `S`ample number for nextseq runs
        debug_print( "lr_patterns generated (w/o lanes): " +
                    str(len(lr_patts.keys())) )
        num_matches = process_files(dir, lr_patts, remove)
        if num_matches == 0:
            # still nothing?! error out
            debug_print('Found zero files matching the lanes/reads '
                        'that were specified!')
            sys.exit(1)

if __name__ == '__main__':
    args = parse_args()
    make_it_happen(args)
