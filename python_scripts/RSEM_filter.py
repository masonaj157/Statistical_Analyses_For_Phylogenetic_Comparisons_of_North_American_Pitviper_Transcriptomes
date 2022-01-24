#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description='reads RSEM output file to filter certain transcripts by TPM. Note, transcripts to be filtered must have a string to identify them')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-r","--rsem",
					type=str,
					default='',
					help="rsem output file (e.g. .genes)")
parser.add_argument("-t","--threshold",
					type=int,
					help="tpm threshold below which sequences will be filtered out")
parser.add_argument("-is1","--identifying_string1",
					type=str,
					help="string1 that identifies transcripts to be filtered")
parser.add_argument("-is2","--identifying_string2",
					type=str,
					help="string2 that identifies transcripts to be filtered")
parser.add_argument("-o","--output",
					type=str,
					default='',
					help="name of output fasta")
args=parser.parse_args()

seqfile = args.seqfile
rsem_file = args.rsem
threshold = args.threshold
idstring1 = args.identifying_string1
idstring2 = args.identifying_string2
output = args.output

#########################################################################################

class rsem:
	def __init__(self, rsem_entry):
		self.gene_id = rsem_entry[0]
		self.transcript_id = rsem_entry[1]
		self.length = rsem_entry[2]
		self.effective_length = float(rsem_entry[3])
		self.expected_count = float(rsem_entry[4])
		self.tpm = float(rsem_entry[5])
		self.fpkm = float(rsem_entry[6])



def rsem_parse(rsem_object) :
    rsem_list = []
    with open(rsem_object) as OF:
            reader = csv.reader(OF, delimiter='\t')
            for row in reader:
                rsem_list.append(row)
    
            rsem_list = [rsem(row) for row in rsem_list[1:] ]
    return(rsem_list)


def filter_by_strings_and_tpm(sequences, rsem_object, idstring1, idstring2,threshold) :
	transcripts_with_string = [ entry for entry in rsem_object if idstring1 in entry.gene_id and idstring2 in entry.gene_id]
	bad_seqs = [ entry.gene_id for entry in transcripts_with_string if entry.tpm < threshold ]
	good_seqs = [ seq for seq in sequences if seq.id not in bad_seqs ] 
	return(good_seqs)
		


#########################################################################################

rsem_object = rsem_parse(rsem_file)

sequences = list(SeqIO.parse(seqfile,"fasta"))


good_seqs = filter_by_strings_and_tpm(sequences, rsem_object, idstring1, idstring2,threshold)

	
handle = open(output, "w")
SeqIO.write(good_seqs, output, "fasta-2line")
handle.close()
