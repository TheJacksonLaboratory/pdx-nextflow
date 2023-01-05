#! /usr/bin/ python
"""
For the given PDX model, retrieve its patient gender via pdxdata API and 
print the gender (in all lowercase) to standard output.  
The API is documented here:  http://pdxdata.jax.org

Observed gender values include:  female  male  unspecified  unknown

"""
from __future__ import print_function
import sys
import requests
import argparse
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('id', help='The model ID whose gender is needed')
    args = parser.parse_args()

    url = 'http://pdxdata.jax.org/api/model?model=' + args.id + '&reqitems=gender'
    response = requests.get(url)
    if response.status_code != 200:
        print('Gender API call non-200 status code: ', response.status_code, file=sys.stderr)
        sys.exit(2)
    else:
        try:
            jsondata = response.json()
        except:
            print('Gender API call JSON decode fail.', file=sys.stderr)
            sys.exit(3)

        if jsondata['count'] == 0:
            print('Gender API call could not find pdx model {0}'.format(args.id), file=sys.stderr)
            sys.exit(1)

        gender = jsondata['data'][0]['gender']

        if gender is None:
            gender = 'unknown' 

        # print the lower-case gender to stdout to be grabbed by pipeline process
        print( gender.lower() )


if __name__ == '__main__':
    main()
