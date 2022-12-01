#! /usr/bin/env python
"""
For the given segments_raw.extend filename, determine PDX model and sample. 
Then look up the PDX passage number.  Finally, print the model_name,
sample_name, and passage number to stdout for use by the CNV plot module.

Rationale: The plots need to include passage number in title, and the png
filenames need to include model, sample, and passage for proper chronological
display on MTB.

Here's an example segments_raw.extend  filename:
      J000078366_J000092578_GES15_03760_SNP.segments_raw.extend.txt 

Usage example:  python get_msp.py pathname

'pathname' is either the full pathname of the above file (this is what we'll get when 
civet executes this) ..or..  just the name of the segments file (for convenience in testing).

For non-PDX contexts, just print something reasonable based on the above filename.

"""

from __future__ import print_function
import sys
import os
import requests
import json


def main():
    if len(sys.argv) != 2:
        sys.exit("ERROR: usage: eg. get_msp.py pathname" )

    try:
        pathparts = sys.argv[1].split("/")
        if len(pathparts) == 1:
            segfile = pathparts[0]
        else:
            segfile = pathparts[-1:][0]
    except:
        sys.exit("ERROR: get_msp.py failed to process the pathname" )

    try:
        nameparts = segfile.split("_")
        model_name = nameparts[0]
        sample_name = nameparts[1]
    except Exception as errmsg:
        print( sys.argv[1] )
        sys.exit( 0 )         # handle non-pdx contexts gracefully

    # hit up the pdxdata api to get passage number for this sample....
    url = "http://pdxdata.jax.org/api/inventory?model=" + model_name + "&sample=" + sample_name + "&reqitems=passage_num"
    response = requests.get(url)

    if response.status_code != 200:
        print( model_name + "_" + sample_name )   # handle non-pdx contexts gracefully
        sys.exit( 0 )         # handle non-pdx contexts gracefully

    try:
        jsondata = response.json()
    except:
        sys.exit( "ERROR: get_msp ... json decode fail." )

    if jsondata["count"] == 0:
        print( model_name + "_" + sample_name )   # handle non-pdx contexts gracefully

    else:
        passage_num = jsondata["data"][0]["passage_num"]
        if passage_num is None:
            passage_num = ""
        print( model_name + "_" + sample_name + "_" + passage_num )



if __name__ == '__main__':
    main()
