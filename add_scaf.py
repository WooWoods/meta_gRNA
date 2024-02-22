import os
import sys

from Bio import Seq
from Bio import SeqIO

def read_guides(gfa):
    guides = dict()
    with open(gfa) as fh:
        grna_parser = SeqIO.parse(fh, 'fasta')
        for record in grna_parser:
            guides[record.id] = record.seq
    return guides

def discard_for_dimer(guides):
    T7 = "gaatTAATACGACTCACTATA" #T7 promoter (minus the first G)
    scaff = "GTTTAAGAGCTATGCTGGAAAC"
    #scaff is the Cas9 sgRNA scaffold
    primer = "GGCATACTCTGCGACATCGT" #primer for fill-in of the oligos. The reverse complement
                    #gets added to the sgRNA template oligos.

    fh = open('meta_guides100.txt', 'w')
    for _id, sequence in list(guides.items()):
        if sequence[0] == "G":
            if sequence[1] == "G": #starts with two G
                oligo = T7 + sequence + scaff
            else: #starts with only one G
                oligo = T7 + "G" + sequence + scaff 
        else:
            oligo = T7 + "GG" + sequence + scaff

        fh.write(f'{_id}\t{oligo}\n')


if __name__ == '__main__':
    gfa = sys.argv[1]
    guides = read_guides(gfa)
    discard_for_dimer(guides)



