import os
import sys
sys.path.append('/public/home/wangycgroup/wuj/bin/MyScripts')

from fastx_parser import fa_parser

def rc(seq):
    bps = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([bps[b] for b in seq][::-1])

def main(fa):
    bases = set(['A', 'T', 'G', 'C'])
    seqs = set()
    with open('guides_rmdup.fa', 'w') as fh:
        for _id, seq in fa_parser(fa):
            if len(set(list(seq)) - bases) > 0:
                continue

            if seq in seqs or rc(seq) in seqs:
                continue
            fh.write(f'>{_id}\n{seq}\n')
            seqs.add(seq)


if __name__ == '__main__':
    fa = sys.argv[1]
    main(fa)
