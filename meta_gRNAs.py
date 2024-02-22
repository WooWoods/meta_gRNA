import re
import os
import sys
import csv
import random
import argparse
import subprocess
from uuid import uuid4
from collections import defaultdict, Counter

import primer3 #primer3-py package
from Bio import Seq
from Bio import SeqIO
from itertools import product


ambiguous_alph = {"N":["A","C","G","T"], "V":["A","C","G"], "H":["A","C","T"], "D":["A","G","T"],
    "B":["C","G","T"], "M":["A","C"], "K":["G","T"], "W":["A","T"], "S":["C","G"], "Y":["C","T"],
    "R":["A","G"], "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"]}

def reverse_dict(dic):
    return {v: k for k, v in dic.items()}

def read_guides(gfa):
    guides = dict()
    with open(gfa) as fh:
        grna_parser = SeqIO.parse(fh, 'fasta')
        for record in grna_parser:
            guides[record.id] = record.seq
    return guides

def check_dir(dirname):
    if os.path.exists(dirname):
        return None
    os.mkdir(dirname)

def remove_ambiguous(sequence):
    #returns a list of all possible sequences
    #given the ambiguous sequence input (e.g. with N,V, etc. bases)
    return ["".join(i) for i in product(*[ ambiguous_alph[j] for j in sequence ]) ]

def revcomp(sequence):
    #returns reverse complement of sequence
    complement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N','Y':'R',
    'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D',
    'V':'B', 'B':'V', 'X':'X', '-':'-'}
    string_list = list(sequence)
    string_list.reverse()
    revcompl = ''
    for base in string_list:
        revcompl += complement[base]
    return revcompl

def find_all(string, substr):
    #makes a list with the indexes of all occurrencies of substr in string
    listindex=[]
    start = 0
    i = string.find(substr, start)
    while i >= 0:
        listindex.append(i)
        i = string.find(substr, i + 1)
    return listindex

def GC_count(string):
    #returns %GC of string
    bases = list(string)
    G = bases.count("G")
    C = bases.count("C")
    GC_content = 100*(G+C)/len(string)
    return GC_content

def retrieve_guide(PAM,strand, sequence, guide_length=20, PAM_length=3):
    #get spacer sequence given the index of the pam sequence (PAM)
    #and the strand of the pam
    if strand == "bot":
        start = PAM - guide_length
        if start > 0:
            end = PAM
            guide = sequence[start:end]
        else: #if the spacer would start outside the rRNA gene
            guide = None
    else:
        start = PAM + PAM_length
        end = start + guide_length
        if end < len(sequence):
            guide = revcomp(sequence[start:end])
        else: #if the spacer would end outside the rRNA gene
            guide = None
    return guide

def offtarget(list_of_PAMs, genome_seqs, positions, coding_genes, scaffold, start, end, index):
    #checks if the offtarget is followed by a PAM
    #and if it's outside a rRNA gene.
    #If so, removes the spacer from guides_dict and returns 1.
    if any (genome_seqs[scaffold][start:end] == PAM for PAM in list_of_PAMs):
        if scaffold in positions:
            if not any(index in range(i,j+1) for [i, j] in positions[scaffold]) and any(index in range(i,j+1) for [i, j] in coding_genes[scaffold]):
                #if the index is not in at least one rRNA gene:
                return 'off'
            else:
                return 'on'
    return None

def primer_heterodimer(oligo, primer_revcomp):
    #checks for unwanted internal binding of the fill-in primer to oligo.
    #cutoffs are predicted melting temperature of the annealing >40 °C and 
    #of the 3' portion of 30 °C
    res_end = primer3.bindings.calcEndStability(primer_revcomp,oligo)
    res_het = primer3.bindings.calcHeterodimer(primer_revcomp,oligo,output_structure=True)
    if res_het.tm > 40 or res_end.tm > 30:
        print("\nThe end-filling primer has high predicted internal binding to the", "oligo\nThis oligo is discarded.")
        print("Heterodimer formation:")
        print(res_het)
        print("End Stability:")
        print(res_end)
        print("Predicted binding:")
        print(res_het.ascii_structure)
        return True
    return None

def build_index(genome_fa):
    check_dir('bowtie_files')
    bowtie_idx = 'bowtie_files/genome_index'
    cmd = (f"bowtie-build -r -q -f {genome_fa} bowtie_files/genome_index"),
    subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL)
    return bowtie_idx

def run_bowtie(input_fa, bowtie_idx):
    bowtie_out = 'bowtie.csv'
    cmd = (f"bowtie -v 3 -a {bowtie_idx} --suppress 6,7 \
        -f {input_fa} > {bowtie_out}") #matches with up to 3 mismatches
    subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL)
    return bowtie_out

def gen_seq_id(rRNA):
    num = random.randint(1000, 3000)
    string = str(uuid4())
    string_seg, *tmp = string.split('-')
    return f'{rRNA}_{num}{string_seg}'

def enumerate_rrna_regions(genome_fa, gtf):
    rrna_pattern = re.compile(r'transcript_biotype "rRNA"')
    regex_rrna = re.compile('(\d*S) ribosomal')
    genome_fh = open(genome_fa)
    gtf_fh = open(gtf)

    genome = SeqIO.parse(genome_fh, 'fasta')
    genome_seqs = {}
    rRNA_genes = {"5S":[], "16S":[], "23S":[]}
    positions = {} #dictionary of start end end positions of the rRNA genes
                   #scaffold:[[start1,end1],[start2,end2]]
    coding_genes = {}

    for record in genome:
        if "genome" in record.description: #if it's the main chromosome
            genome_id = re.findall('[^ ]+', record.description)[0]
            genome_seqs[genome_id] = str(record.seq).upper()
        else:
            plasmid_id = re.findall('[^ ]+', record.description)[0]
            genome_seqs["{0}".format(plasmid_id)] = str(record.seq).upper()
    genome_fh.close()

    for line in gtf_fh:
        if line[0] == "#":
            continue
        linelist = line.split("\t")
        if linelist[2] == "transcript" and rrna_pattern.search(linelist[8]):
            m = regex_rrna.search(linelist[8])
            if m is None:
                continue
            ID = m.group(1)
            start = int(linelist[3])
            end = int(linelist[4])
            strand = linelist[6]
            scaffold = linelist[0]
            positions.setdefault(scaffold,[]).append([start,end])
            if strand == "-":
                seq = revcomp(genome_seqs[scaffold][start:end])
            else:
                seq = str(genome_seqs[scaffold][start:end])
            rRNA_genes[ID].append(seq)
            continue
        if linelist[2] == 'gene':
            start = int(linelist[3])
            end = int(linelist[4])
            scaffold = linelist[0]
            coding_genes.setdefault(scaffold,[]).append([start,end])
            continue

    return rRNA_genes, positions, coding_genes, genome_seqs

def grna_miner(rRNA_genes, PAM_list, GC_low=30, GC_high=70):
    guides = dict() # gRNA_seq : gRNA_id
    for rRNA in rRNA_genes:
        for sequence in rRNA_genes[rRNA]: #loops through the sequence of each copy of all rRNAs
            PAMs_top = [] #PAM sequences in the top strand
            PAMs_bot =  [] #PAM sequences in the bottom strand
            for PAM_sequence in PAM_list:
                PAMs_bot.extend(find_all(sequence, PAM_sequence))
                PAMs_top.extend(find_all(sequence, revcomp(PAM_sequence)))

            for PAM in PAMs_bot:
                guide = retrieve_guide(PAM, "bot", sequence)
                if guide and guide not in guides:
                    if GC_low <= GC_count(guide) <= GC_high:
                        _id = gen_seq_id(rRNA)
                        guides[_id] = guide
                    else:
                        continue
            for PAM in PAMs_top:
                guide = retrieve_guide(PAM, "top", sequence)
                if guide and guide not in guides:
                    if GC_low <= GC_count(guide) <= GC_high:
                        _id = gen_seq_id(rRNA)
                        guides[_id] = guide
                    else:
                        continue
    return guides

def prepare_for_bowtie(guides):
    input_fa = 'input.fa'
    with open(input_fa, 'w+') as fh: #input file for bowtie
        for _id, seq in guides.items():
            fh.write(f'>{_id}\n{seq}\n')
    return input_fa

def discard_for_dimer(guides):
    T7 = "TTCTAATACGACTCACTATA" #T7 promoter (minus the first G)
    scaff = "GTTTAAGAGCTATGCTGGAAAC"
    #scaff is the Cas9 sgRNA scaffold
    primer = "GGCATACTCTGCGACATCGT" #primer for fill-in of the oligos. The reverse complement
                    #gets added to the sgRNA template oligos.
    primer_revcomp = revcomp(primer)

    for _id, sequence in list(guides.items()):
        if sequence[0] == "G":
            if sequence[1] == "G": #starts with two G
                oligo = T7 + sequence + scaff + primer_revcomp
                ID = primer_heterodimer(oligo[:-len(primer_revcomp)], primer_revcomp)
            else: #starts with only one G
                oligo = T7 + "G" + sequence + scaff + primer_revcomp
                ID = primer_heterodimer(oligo[:-len(primer_revcomp)], primer_revcomp)
        else:
            oligo = T7 + "GG" + sequence + scaff + primer_revcomp
            ID = primer_heterodimer(oligo[:-len(primer_revcomp)], primer_revcomp)
        
        if ID:
            try:
                del guides[_id]
            except KeyError:
                continue
    return guides

def clapse_by_density(guides, grna_candidates):
    input_fa = 'input.fa'
    with open(input_fa, 'w') as fh:
        for _id, seq in guides.items():
            fh.write(f'>{_id}\n{seq}\n')
        for _id, seq in grna_candidates.items():
            fh.write(f'>{_id}\n{seq}\n')

    bowtie_out = run_bowtie(input_fa)

    grna_positions = dict()
    with open(bowtie_out) as fh:
        for line in fh:
            arr = line.split()
            scaffold = arr[2] #chromosome or plasmid ID
            index = int(arr[3]) #starting index of off-target
            grna_positions[index] = arr[0]
    positions_sorted = sorted(grna_positions.keys())
    counter = Counter(grna_positions.values())

    discard_pos = []
    i, j = 0, 1
    while j < len(positions_sorted):
        if positions_sorted[j] - positions_sorted[i] > 50:
            i = j
            j += 1
            continue
        else:
            discarded = discard_for_dense(positions_sorted[i], positions_sorted[j], grna_positions, counter, grna_candidates)
            discard_pos.append(discarded)
            j += 1
    new_dict = dict()
    for _pos, _id in grna_positions.items():
        if _pos in discard_pos:
            continue
        try:
            new_dict[_id] = guides[_id]
        except KeyError:
            continue

    return new_dict

def discard_for_dense(i, j, grna_positions, counter, grna_candidates):
    grna1 = grna_positions[i]
    grna2 = grna_positions[j]

    if grna1 in grna_candidates and grna2 not in grna_candidates:
        return j
    elif grna1 in grna_candidates and grna2 in grna_candidates:
        if counter[grna1] > counter[2]:
            return j
        else:
            return i
    else:
        return random.choice((i, j))
    
def discard_for_offtarget(bowtie_out, genome_seqs, positions, coding_genes, guides, PAM_list, PAM_length=3, guide_length=20):
    revcomp_PAM_list = []
    for PAM in PAM_list:
        PAM = revcomp(PAM)
        revcomp_PAM_list.append(PAM)

    new_guides = {}
    offlist = set()
    onlist =set()
    with open(bowtie_out) as fh, open('offtarget.fa', 'w') as fo:
        for row in fh:
            arr = row.split()
            scaffold = arr[2] #chromosome or plasmid ID
            index = int(arr[3]) #starting index of off-target
            if arr[1] == "+": #strand
                start = index + guide_length
                end = index + guide_length + PAM_length
                discard = offtarget(PAM_list, genome_seqs, positions, coding_genes, scaffold, start, end, index)
                                #this sequence as a true
                                #off-target, the gRNA is
                                #discarded and offtargeting
                                #is increased by 1. Otherwise
                                #by 0.
            else:
                start = index - PAM_length
                end = index
                discard = offtarget(revcomp_PAM_list, genome_seqs, positions, coding_genes, scaffold, start, end, index)

            if discard == 'off':
                try:
                    offlist.add(arr[0])
                    fo.write(f'>{arr[0]}\n{guides[arr[0]]}\n')
                    del guides[arr[0]]
                except KeyError:
                    continue
            elif discard == 'on':
                onlist.add(arr[0])
            else:
                continue
    for _id in onlist - offlist:
        new_guides[_id] = guides[_id]
    return new_guides

def update_db(grna_candidates, grna_remained):
    grnas = list(grna_candidates.values())
    for _id, seq in grna_remained.items():
        if seq not in grnas:
            grna_candidates[_id] = seq

def main(args):
    if args.command == 'meta':
        GC_low = args.minGC
        GC_high = args.maxGC
        guide_length = args.length
        show_offtargets = args.offtargets
        annotation_file = args.manual_ann
        PAM = args.pam.upper()
        genome_dir = args.genomedir

        check_dir('bowtie_files')

        PAM_length = len(PAM)
        PAM_list = remove_ambiguous(PAM)

        grna_candidates = dict()

        genome_fa = ''
        gtf = ''
        for path in os.listdir(genome_dir):
            gpath = os.path.join(genome_dir, path)
            if not os.path.isdir(gpath):
                continue
            for f in os.listdir(gpath):
                if f.endswith('fna'):
                    genome_fa = os.path.join(gpath, f)
                if f.endswith('gtf'):
                    gtf = os.path.join(gpath, f)

            if not (genome_fa and gtf):
                print(gpath)
                continue
            
            rRNA_genes, positions, coding_genes, genome_seqs = enumerate_rrna_regions(genome_fa, gtf)
            build_index(genome_fa)
            guides = grna_miner(rRNA_genes, PAM_list)
            # remove off targets in candidates
            if grna_candidates:
                input_fa = prepare_for_bowtie(grna_candidates)
                bowtie_out = run_bowtie(input_fa)
                discard_for_offtarget(bowtie_out, genome_seqs, positions, coding_genes, grna_candidates, PAM_list)

            input_fa = prepare_for_bowtie(guides)
            bowtie_out = run_bowtie(input_fa)
            guides = discard_for_offtarget(bowtie_out, genome_seqs, positions, guides, PAM_list)
            guides = discard_for_dimer(guides)
            #grna_remained = clapse_by_density(guides, grna_candidates)
            update_db(grna_candidates, guides)

            genome_fa = ''
            gtf = ''
        
        with open('guides.fa', 'w') as fh:
            for _id, seq in grna_candidates.items():
                fh.write(f'>{_id}\n{seq}\n')

    if args.command == 'offtarget':
        genome_fa = args.genome
        gtf = args.gtf
        grnas = args.grnas
        PAM = args.pam.upper()
        PAM_list = remove_ambiguous(PAM)
        bowtie_idx = args.bowtie_index

        rRNA_genes, positions, coding_genes, genome_seqs = enumerate_rrna_regions(genome_fa, gtf)
        if not bowtie_idx:
            bowtie_idx = build_index(genome_fa)
        bowtie_out = run_bowtie(grnas, bowtie_idx)
        guides = read_guides(grnas)
        guides = discard_for_offtarget(bowtie_out, genome_seqs, positions, coding_genes, guides, PAM_list)

        with open('guides_rmoff2.fa', 'w') as fh:
            for _id, seq in guides.items():
                fh.write(f'>{_id}\n{seq}\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Designs a pool of sgRNAs '
            'targeting the ribosomal genes of a bacterial species.')

    AP_subparsers = parser.add_subparsers(dest='command',
        help="Sub-commands (use with -h for more info)"
        )

    # subparser for meta
    meta_parser = AP_subparsers.add_parser('meta', help="Design sgRNAs for meta-genomics")
    meta_parser.add_argument('--minGC', '-g', default=30, type=int,
                        help='Minimal accepted GC%% of a spacer (Default: 30).')

    meta_parser.add_argument('--maxGC', '-G', default=80, type=int,
                        help='Maximal accepted GC%% of a spacer (Default: 80).')

    meta_parser.add_argument('--length', '-l', default=20, type=int,
                        help='Spacer length (Default: 20).')

    meta_parser.add_argument('--offtargets', '-o', default=False, action='store_true',
                        help='Print the spacers that were discarded because of off-targeting.')

    meta_parser.add_argument('--manual_ann', '-ma', default=False, action='store',
                        help='Use if you want to provide a file with the manual '
                        'annotation (in bed format, tab-separated) of the rRNA genes (Default: False).')

    meta_parser.add_argument('--pam', '-p', default="NGG", type=str,
                        help='PAM sequence (Default: NGG).')

    meta_parser.add_argument('--genomedir', '-ref', default=False, type=str,
                        help='Directory of genomes')

    # subparser for off target detect
    off_parser = AP_subparsers.add_parser('offtarget', help='Remove off target gRNAs')
    off_parser.add_argument('--genome', '-ref', default=False, type=str,
                        help='Fasta of genomes')

    off_parser.add_argument('--gtf', '-gtf', default=False, type=str,
                        help='A gtf file to define rRNA regions')

    off_parser.add_argument('--grnas', '-grnas', default=False, type=str,
                        help='sgRNAs in fasta format')

    off_parser.add_argument('--pam', '-p', default="NGG", type=str,
                        help='PAM sequence (Default: NGG).')

    off_parser.add_argument('--bowtie_index', '-index', default=False, type=str,
                        help='Bowtie index of the genome')

    args = parser.parse_args()

    main(args)
