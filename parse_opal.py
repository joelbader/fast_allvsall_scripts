#!/usr/bin/python3

import argparse
import re

parser = argparse.ArgumentParser(description='Parses output from opal')
parser.add_argument('opal_filename', type=str, help="Input file (output from opal -x 1)")
parser.add_argument('fa_filename', type=str, help="FASTA file containing protein database (*.fa)")
parser.add_argument('out_filename', type=str, help="FASTA filename for output (*.fa)", default='')
parser.add_argument('--threshhold','-t', type=float, help="Minimum identity for calling a hit (0.75 = 75% identity)", default='')

args = parser.parse_args()
opal_in  = open(args.opal_filename, "r")
fasta_in  = open(args.fa_filename, "r")
outfile = open(args.out_filename, "w")
thresh = args.threshhold

opal_no_head = opal_in.readlines()[12:]

fa_lines = list(fasta_in)
fa_i = 0
for i,line in enumerate(opal_no_head):
    if line[0]=='#':
        #Processing FASTA with this while loop is faster than using Bio.Seq
        sequence=''
        found_comment = False
        while fa_i < len(fa_lines):
            fa_line = fa_lines[fa_i]
            if fa_line[0] != '>' and fa_line[0] != '!':
                sequence = sequence + fa_line
            else:
                if found_comment:
                    break
                else:
                    comment = fa_line
                    found_comment = True
            fa_i += 1

        #float(line.split(' (,)')[1])
        split_ls = re.split('[#: ,()\n]+',line)
        num_match = split_ls[2]
        query_len = split_ls[5]
        identity = float(num_match)/(1+float(query_len))
        if identity > thresh:
            #print fasta entry to new file
            #print(identity)
            outfile.write(comment)
            outfile.write(sequence)
    else:
        break

opal_in.close()
fasta_in.close()
outfile.close()
