#! /usr/bin/env python

# concatenate_lanes.py

import sys
import os
import re
import glob
import inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
lib_folder = os.path.join(cmd_folder, '../lib')
if lib_folder not in sys.path:
     sys.path.insert(0, lib_folder)

import cga_version as version

def usage():
    if len(sys.argv) != 4:
        print >> sys.stderr, 'USAGE:', sys.argv[0], 'project-directory  lane-id1 lane-id2'
        print >> sys.stderr, '  e.g.', sys.argv[0], 'Project_CLIA_Assay 1 2'
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

def build_related_filenames(dir, lane1, lane2):
    files = {}
    ln1 = 'L00' + lane1
    ln2 = 'L00' + lane2
    key_file_str = '*' + ln1 + '*_R1_*.fastq'
    path=os.path.join(dir, key_file_str)
    key_file = glob.glob(path)[0]
    files['l1r1'] = key_file
    pat = re.compile('(.*?)_' + ln1 + '_R1_(.*fastq)')
    sub1 = r'\1_' + ln2 + r'_R1_\2'
    sub2 = r'\1_combined_R1_\2'
    sub3 = r'\1_' + ln1 + r'_R2_\2'
    sub4 = r'\1_' + ln2 + r'_R2_\2'
    sub5 = r'\1_combined_R2_\2'
    files['l2r1'] = re.sub(pat, sub1, key_file)
    files['cr1']  = re.sub(pat, sub2, key_file)
    files['l1r2'] = re.sub(pat, sub3, key_file)
    files['l2r2'] = re.sub(pat, sub4, key_file)
    files['cr2']  = re.sub(pat, sub5, key_file)
    return files

def main():
    # Check for a version request.
    version.parse_options()
    usage()
    dirs = find_directories(sys.argv[1])
    for d in dirs:
        print 'Processing', d
        files = build_related_filenames(d, sys.argv[2], sys.argv[3])
        of1 = open(files['cr1'], 'w')
        for line in open(files['l1r1']):
            of1.write(line)
        for line in open(files['l2r1']):
            of1.write(line)
        of1.close()
        of2 = open(files['cr2'], 'w')
        for line in open(files['l1r2']):
            of2.write(line)
        for line in open(files['l2r2']):
            of2.write(line)
        of2.close()
    print 'Completed...'

main()