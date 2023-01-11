#! /usr/bin/env python

# concatenate_runs.py

# This script concatenates identically-named fastq files in a project that has
# been run in two runs.  Note that this means:
#  - the project directory must be the same between the two runs, 
#  - the indexes must be the same,
#  - and the lanes must be the same.


import sys
import os
import re
import glob
import argparse
import inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
lib_folder = os.path.join(cmd_folder, '../lib')
if lib_folder not in sys.path:
     sys.path.insert(0, lib_folder)

import cga_version as version

def usage():
    if len(sys.argv) != 4:
        print >> sys.stderr, 'USAGE:', sys.argv[0], 'run-directory-1 run-directory-2 project-directory'
        print >> sys.stderr, '      ', sys.argv[0], '-v | --version'
        print >> sys.stderr, 'project-ditectory must exist in both run-directory-1 and run-directory-2'
        print >> sys.stderr, 'output will be in run-directory-2/combined/project-directory'
        sys.exit(0)
    
def find_directories(dir):
    found_dirs = []
    for root, children, files in os.walk(dir):
        for f in files:
            if f.endswith('fastq'):
                if root not in found_dirs:
                    found_dirs.append(root)
                break
            
    if len(found_dirs) == 0:
        sys.stderr.write("ERROR: no fastq directory found!\n")
        sys.exit(1)
    return sorted(found_dirs)

def build_related_filenames(pd1, pd2, sample_dir, file1, od):
    files = {}
    # Get the fill filepath to a fastq file
    files['f1'] = os.path.join(sample_dir, file1)
    # now just the part past the first run's project directory
    f = files['f1'][len(pd1):]
    # use that to build the file paths for the other input file and output
    files['f2'] = pd2 + f 
    try:
        with open (files['f2']): pass
    except IOError:
        print >> sys.stderr, "ERROR: Oops, can't find paired fastq file in second run,", files['f2'], '    skipping.'
        return {}
    files['of'] = od + f
    return files

def main():
    # Check for a version request.
    version.parse_options()
    usage()
    rd1 = sys.argv[1]
    rd2 = sys.argv[2]
    pd = sys.argv[3]
    
    pd1 = os.path.join(rd1, 'Unaligned', pd)
    pd2 = os.path.join(rd2, 'Unaligned', pd)
    od = os.path.join(pd2, 'combined')
    print >> sys.stderr, "Starting..."
    if not os.path.exists(od):
        print 'Creating output directory', od
        try:
            os.mkdir(od)
        except:
            print >> sys.stderr, 'ERROR: Could not create combined output directory,', od
            print >> sys.stderr, '       Please fix this condition, and rerun.  Exiting...'
            sys.exit(1)
    if not os.path.exists(od):
        print >> sys.stderr, "ERROR:  GAK! the combined output directory still doesn't exist."
        print >> sys.stderr, "        This can't happen. Please log a software discrepancy report."
        print >> sys.stderr, "        Exiting..."
        sys.exit(1)

    print 'Searching', pd1
    # Find all the sample directories in the run1 project
    dirs = find_directories(pd1)
    for d in dirs:
        print 'Processing', d
        fqs = os.listdir(d)
        for n in range(len(fqs)-1, -1, -1):
            if not fqs[n].endswith('fastq'):
                del(fqs[n])
        for fq in fqs:
            print '    Combining', fq
            files = build_related_filenames(pd1,pd2, d, fq, od)
            if len(files) == 0:
                continue
            out_samp_dir = os.path.split(files['of'])[0]
            if not os.path.exists(out_samp_dir):
                os.mkdir(out_samp_dir)
            of1 = open(files['of'], 'w')
            for line in open(files['f1']):
                of1.write(line)
            for line in open(files['f2']):
                of1.write(line)
            of1.close()
    print 'Completed...'

main()