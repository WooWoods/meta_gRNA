import os
import sys

from collections import defaultdict
from Bio import SeqIO


def load_db(fbowtie):
    data = defaultdict(set)
    with open(fbowtie) as fh:
        for line in fh:
            arr = line.split()
            data[arr[0]].add(arr[2])
    return data

def main(off, db):
    with open(off) as fh, open('offtarget_detail.txt', 'w') as fo:
        grna_parser = SeqIO.parse(fh, 'fasta')
        for record in grna_parser:
            chroms = ','.join(db[record.id])
            fo.write(f'{record.id}\t{record.seq}\t{chroms}\n')

if __name__ == '__main__':
    off, fbowtie = sys.argv[1:]
    db = load_db(fbowtie)
    main(off, db)

