#!/usr/bin/env python

# ChimeraSorter was made as a followup to ChimeraChecker. It uses several filters to identify and sort filters into good and bad. If you have questions about this script, send a postcard with your information and question addressed to My Ass, You're Screwed Blvd., Imalone, WI 54848

# Additional software necessary to run this:
# (1) bwa 
# (2) gatk
# (3) samtools
# (4) bedtools
# (5) picard
# (6) bam-read from transrate
import time
start = time.time()

import argparse
import sys, os, shutil
import subprocess
import csv
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='Automated checking of coding seqeuences for chimeric transcripts')
parser.add_argument("-i","--input",
                    type=argparse.FileType('r'),
                    help="Fasta file of contigs to check. Use only CODING SEQUENCES.")
parser.add_argument("-r","--reads",
                    type=argparse.FileType('r'),
                    help="Fastq file of UNPAIRED reads.")
parser.add_argument("-o","--outfile",
                    type=str,
                    help="name of the output .csv file")
parser.add_argument("-a","--aligner",
                    type=str,
                    help="aligner you want to use. Either bwa or bowtie2. Aligner must be in the path or specified with a flag")
parser.add_argument("-b","--bwa",
                    nargs='?',
                    type=str,
                    default="/usr/local/bin/bwa",
                    help="Path to bwa executable. Default assumes it is in your PATH.")
parser.add_argument("-bw","--bowtie2",
                    nargs='?',
                    type=str,
                    default="/usr/local/bin/bowtie2",
                    help="Path to bowtie2 executable. Default assumes it is in your PATH.")
parser.add_argument("-s","--samtools",
                    nargs='?',
                    type=str,
                    default="samtools",
                    help="Path to samtools executable. Default assumes it is in your PATH.")
parser.add_argument("-bt","--bedtools",
                    nargs='?',
                    type=str,
                    default="bedtools",
                    help="Path to bedtools executable. Default assumes it is in your PATH.")
parser.add_argument("-m","--mismatches",
                    nargs='?',
                    type=int,
                    default=3,
                    help="Number of allowable mismatches to keep in alignment. The default is 0.")
parser.add_argument("-ts","--tooShort",
                    nargs='?',
                    type=int,
                    default=150,
                    help="Minimum length of clipped reads. Default is 150.")
parser.add_argument("-mq","--mapQuality",
                    nargs='?',
                    type=int,
                    default=0,
                    help="Minimum mapping quality. Default is 0. Note that reads with multiple mappings get assigned a quality of 0 by bwa.")
parser.add_argument("-min","--minRead",
                    nargs='?',
                    type=int,
                    default=150,
                    help="Minimum read length. Default is 150.")
parser.add_argument("-max","--maxRead",
                    nargs='?',
                    type=int,
                    default=1000,
                    help="Maximum read length. Default is 1000.")
parser.add_argument("-p","--picard",
                    nargs='?',
                    type=str,
                    default="java -jar /path/picard.jar ",
                    help="Picard command. Something like: java -jar /PATH/TO/PICARD/picard.jar")
parser.add_argument("-g","--gatk",
                    nargs='?',
                    type=str,
                    default="java -jar /path/GenomeAnalysisTK.jar",
                    help="GATK command. Some like: java -jar /PATH/TO/GATK/GenomeAnalysisTK.jar")
parser.add_argument("-c","--cov",
                    type=int,
                    default=5,
                    help="Coverage minimum for transcripts. Default is 5x.")
args = parser.parse_args()



# Check the input fasta file and output the number of detected sequences.
totalContigs = 0
for seq in SeqIO.parse(args.input,format="fasta") :
    totalContigs += 1
print("Total number of input contigs = " + str(totalContigs))
print("Maximum allowed mismatches per read retained in the alignment = " + str(args.mismatches))

print("*"*100)

print(args.aligner)


if (args.aligner != 'bwa') and (args.aligner != 'bowtie2') :
	print(" Need to specify an aligner")
	quit()



if args.aligner == 'bwa' :
	# Run the alignment
	grepNM = "grep -E 'NM:i:[0-" + str(args.mismatches) + "][[:space:]]|^@'"
	name = args.input.name.split(".")[0] 
	print("Generating bwa alignment: " + name + ".bam")
	# Build the bwa index
	command = args.bwa + " index " + args.input.name
	subprocess.call(command,shell=True)
	# Generate the initial sam alignment
	command = args.bwa + " mem -M -t 4 -R \'@RG\\tID:" + args.input.name + "\\tSM:" + args.reads.name + "' " + args.input.name + " " + args.reads.name + " | " + grepNM  + " > tmp1.sam"
	subprocess.call(command,shell=True)


if args.aligner == 'bowtie2' :
	name = args.input.name.split(".")[0]
	#Make bowtie2 index
	bowtie2_index = args.bowtie2 + '-build' 
	command = bowtie2_index + " " + args.input.name + " " + name
	subprocess.call(command, shell=True)
	#Run bowtie2 alignment
	command = args.bowtie2 + " -x " + name + " -U " + args.reads.name + " -S tmp1.sam --rg-id " + args.input.name + " --rg SM:" + args.reads.name
	subprocess.call(command, shell=True)

	
# Create a sorted bam file
command = args.picard + " SortSam INPUT=tmp1.sam OUTPUT=tmp2.bam SORT_ORDER=coordinate"
print(command) 
subprocess.call(command,shell=True)
command = args.picard + " BuildBamIndex INPUT=tmp2.bam"
print(command) 
subprocess.call(command,shell=True)
# Remove overclipped reads
command = args.picard + " CreateSequenceDictionary REFERENCE=" + args.input.name + " OUTPUT=" + args.input.name.split(".")[0] + ".dict"
print(command) 
subprocess.call(command,shell=True)
command = args.samtools + " faidx " + args.input.name
print(command) 
subprocess.call(command,shell=True)
command = args.gatk + " PrintReads -R " + args.input.name + " -I tmp2.bam  -RF OverclippedReadFilter --filter-too-short " + str(args.tooShort) + " --dont-require-soft-clips-both-ends -RF MappingQualityReadFilter --minimum-mapping-quality " + str(args.mapQuality) + " -RF ReadLengthReadFilter --min-read-length " + str(args.minRead) + " --max-read-length " + str(args.maxRead) + " -O " + name + ".bam"
print(command)
subprocess.call(command,shell=True) 
#Calculate coverage
command = args.bedtools + " genomecov -d -ibam " + name + ".bam > coverage.txt"
print(command) 
subprocess.call(command,shell=True)
subprocess.call("rm tmp*[bs]a[im]",shell=True)




# Read in the coverage data
print("*"*100)
print("Importing coverage information.")
X = pd.read_csv("coverage.txt",sep='\t',names=['transcript','pos','cov'])
print("Identified coverage data for " + str(len(set(X['transcript']))) + " transcripts.")
T = list(set(X['transcript']))


## check for sequences with less than 5X coverage in greater than 10 percent of sequences
transcript_list = []
for i in T :
	#length = len(list(X[X['transcript'] == T[0]]['cov']))
	length = len(list(X[X['transcript'] == i]['cov']))
	x = 0
	for j in range(0,args.cov):
	#for j in range(0,20):
		#thing = list(X[X['transcript'] == T[0]]['cov']).count(j)
		thing = list(X[X['transcript'] == i]['cov']).count(j)
		x = x + thing
	if x > (0.1 * length):
		entry = [i, 'absent']
		transcript_list.append(entry)
	else:
		entry = [i, 'present']
		transcript_list.append(entry)
		

name = args.outfile + '.csv'

with open(name, 'w') as csv_file:
	csv_writer = csv.writer(csv_file, delimiter = ',')
	for row in transcript_list:
		csv_writer.writerow(row)
		
csv_file.close()

 
end = time.time()
print(end - start)
                
