#! /usr/bin/env python
import os
import sys
import dendropy
import subprocess

#inputs: Alignment, Tree, query sequence
#outputs: New alignment including query seq, New tree with query seq in both EPA and starnderd format, EPA list of placement probabilities.

#
#mafft --auto schirrmeister2011_1220_clean.fas > schirrmeister2011_1220_clean
#raxml -m GTRCAT -n schirrmeister2011_1220_clean -p 1 -s schirrmeister2011_1220_clean 
#

if ("--help" in sys.argv) or ("-?" in sys.argv):
    sys.stderr.write("usage: placement.py [<query-sequence>] [<alignment-file-path>] [<newick-file-path>]\n")
    sys.exit(1)
   
if len(sys.argv) < 3:
     sys.stderr.write("Not enough arguments. Usage: placement.py [<query-sequence>] [<alignment-file-path>] [<newick-file-path>]")
else:
    query_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
    if not os.path.exists(query_fpath):
        sys.stderr.write('Not found: "%s"' % query_fpath)
    query = open(query_fpath)        
    align_fpath = os.path.expanduser(os.path.expandvars(sys.argv[2]))
    if not os.path.exists(align_fpath):
        sys.stderr.write('Not found: "%s"' % align_fpath)
    align = open(align_fpath)
    tree_fpath = os.path.expanduser(os.path.expandvars(sys.argv[3]))
    if not os.path.exists(tree_fpath):
        sys.stderr.write('Not found: "%s"' % tree_fpath)
    tree = open(tree_fpath)

outfile="test_out"
#print(query_fpath)
subprocess.call(['pagan', '-q', query_fpath, '-a', align_fpath, '-t', tree_fpath, '-o', outfile])

#print(" ".join(['raxmlHPC', '-f', 'v', '-s', '{}.fas'.format(outfile), '-t',  tree_fpath, '-m', 'GTRCAT', '-n', 'TEST']))
subprocess.call(['raxmlHPC', '-f', 'v', '-s', '{}.fas'.format(outfile), '-t',  tree_fpath, '-m', 'GTRCAT', '-n', 'TEST'])
