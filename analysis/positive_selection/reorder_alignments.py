#!/usr/bin/env python3
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse

def run(args):
    alignments = AlignIO.read(args.input, "phylip-sequential")
    alignments = [a for a in alignments]
    first_index = 0
    for i, a in enumerate(alignments):
        if a.id == args.first:
            first_index = i
            break
    first = alignments.pop(first_index)

    alignments.insert(0, first)
    AlignIO.write(MultipleSeqAlignment(alignments), args.output, "phylip-sequential")

parser = argparse.ArgumentParser(
                    prog = 'Reorder Sequences')
parser.add_argument("-i", "--input", dest="input", required=True)
parser.add_argument("-o", "--output", dest="output", required=True)
parser.add_argument("-f", "--first", dest="first", required=True)
parser.set_defaults(func=run)
args = parser.parse_args()
args.func(args)