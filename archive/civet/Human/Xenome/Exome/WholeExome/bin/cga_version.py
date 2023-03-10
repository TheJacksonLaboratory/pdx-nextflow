#! /usr/bin/env python

import os
import sys
import subprocess
import inspect

def version_from_git():
    cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( 
        inspect.getouterframes(inspect.currentframe())[2][0] ))[0]))
    try:
        err = '' # Make sure err exists, if we hit the except clause
        p = subprocess.Popen('cd {0}; git describe --tags'.format(cmd_folder),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        out, err = p.communicate()
    except OSError, e:
        print >> sys.stderr, 'ERROR determining git version:\n', str(e)
        out = '(undetermined)'
    # Mask any errors, for instance not running in a git working directory.
    if err:
        print >> sys.stderr, 'ERROR misc error:\n', out, err
        out = '(undetermined)'
    return 'V' + out.strip()

def parse_options():
    #
    # "But there are libraries that do option parsing! Why do it by hand?"
    #
    # True, but this module is included in a variety of others that all have
    # their own (real) options.  This code only needs to see whether the
    # command has a single option which is either -v or --version. If
    # so, print out the version string and exit; else simply return and
    # let the rest of the program flow happen.
    #
    # If we use argparse, we have to teach it about all the possible options
    # in all the scripts that import this. Trust us.  We tried that way first.
    #
    options = ['-v', '--version']
    if len(sys.argv) == 2:
        if sys.argv[1] in options:
            print_version_string_and_exit()

def print_version_string_and_exit():
    path = os.path.abspath(sys.argv[0])
    script = os.path.split(path)[1]
    print script, version_from_git(), path
    sys.exit(0)

if __name__ == '__main__':
    parse_options()
