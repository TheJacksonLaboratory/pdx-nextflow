#!/usr/bin/env python

# aggregate_CGA_stats.py OUT INP_QC  INP_HS
# A script to parse the quality and hybrid-selection statistics files
# and aggregate relevant metrics into an output file.

# Parameters:
# out = output file
# inp_qc = *stat file output by qualtool
# inp_hs = *Metricsfile.txt outoput by Picard CalculateHsMetrics 

import sys

if len(sys.argv) < 4:
    print >>sys.stderr, "Commandline arguments missing:\nFormat: aggregate_CGA_stats.py OUT INP_QC INP_DUP INP_HS\nout = output file\ninp_qc = *stat file output by qualtool\ninp_hs = *Metricsfile.txt outoput by Picard CalculateHsMetrics"
    sys.exit()
    
out = open(sys.argv[1],"w")
inp_qc  = open(sys.argv[2],"r")
inp_hs  = open(sys.argv[3],"r")

qc_out = [None, None]
read_data = False
for line in inp_qc:
    line = line.strip()
    elems = line.split("\t")
    
    if line.startswith("QC statistics"):
        read_data = True
    
    if line.startswith("Detailed QC statistics"):
        break
    
    if read_data:
        if None not in qc_out:
            break
        if line.startswith("Total number of reads"):
            try:
                elems = line.split() 
                qc_out[0] = str(int(elems[-1]))
            except Exception:
                qc_out[0] = "NA"
        if line.startswith("Total number of HQ filtered reads"):
            try:
                elems = line.split() 
                qc_out[1] = str(int(elems[-1]))
            except Exception:
                qc_out[1] = "NA"
print >>out, "Total number of reads\t%s\nTotal number of HQ filtered reads\t%s" %(qc_out[0],qc_out[1])



data_lines = []
for line in inp_hs:
    line = line.strip()
    if line.startswith("Mapped reads:"):
     print >>out ,line   
    if line.startswith("Forward strand:"):
     print >>out ,line
    if line.startswith("Reverse strand:"):
     print >>out ,line
