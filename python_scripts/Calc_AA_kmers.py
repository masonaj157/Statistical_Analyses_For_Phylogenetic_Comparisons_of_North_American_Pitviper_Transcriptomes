#!/usr/bin/env python

import argparse
import csv
from Bio import SeqIO

# Command line options
parser = argparse.ArgumentParser(description='This script will take an input fasta and spit out a tab-delimited file of kmers and their counts')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta to translate")
parser.add_argument("-k","--kmer_size",
					type=int,
                    default=20,
					help="kmer size for counting")                    
parser.add_argument("-o","--output",
                    type=str,
                    help="name of output file")
parser.add_argument("-hs","--have_stop",
					type=str,
                    default='True',
					help="Do your CDS include the final stop codon (True or False)?")
args=parser.parse_args()

fasta=args.fasta
kmer_size=args.kmer_size
output=args.output
have_stop=args.have_stop

#########################################################################################

def trim_final_stop(sequences):
    for seq in sequences:
        seq.seq = seq.seq[0:-3]
    return(sequences)


def translate_sequences(sequences):
    for seq in sequences:
        seq.seq = seq.seq.translate()
    return(sequences)


def extract_kmers(seq, kmer_size):
    if len(seq.seq) < kmer_size:
        raise Exception("Some sequences less than kmer_size")
    else:
        kmers = []
        start = 0
        end = start + kmer_size
        while end <= len(seq.seq):
            kmers.append(seq.seq[start:end])
            start +=1
            end = start + kmer_size
    return(kmers)


def collect_all_unique_kmers(sequences, kmer_size):
    all_kmers = []
    for seq in sequences:
        kmers = extract_kmers(seq, kmer_size)
        for kmer in kmers:
            all_kmers.append(kmer)
    
    unique_kmers = list(set(all_kmers))
    return(unique_kmers)


def calc_kmers_counts(sequences,unique_kmers):
    kmer_dict = {}
    for kmer in unique_kmers:
        kmer_dict[kmer] = 0
    
    for seq in sequences:
        for kmer in unique_kmers:
            if kmer in seq.seq:
                value = seq.seq.count_overlap(kmer)
                new_value = value + kmer_dict[kmer]
                kmer_dict[kmer] = new_value
    
    return(kmer_dict)

#########################################################################################

sequences = list(SeqIO.parse(fasta,"fasta"))

if have_stop in ['True','T','t','true']:
    sequences = trim_final_stop(sequences)

sequences = translate_sequences(sequences)

unique_kmers = collect_all_unique_kmers(sequences, kmer_size)

kmer_dict = calc_kmers_counts(sequences,unique_kmers)


final_counts = []
for kmer in unique_kmers:
    entry = [kmer,kmer_dict[kmer]]
    final_counts.append(entry)


with open(output,'w') as outFile:
	outfile_writer = csv.writer(outFile, delimiter = '\t')
	outfile_writer.writerow(['kmer', 'count'])
	for entry in final_counts:
		outfile_writer.writerow(entry)


outFile.close()




