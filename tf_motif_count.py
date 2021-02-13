#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 12:56:07 2020
@author: jasonmoggridge

COUNT TRANSCRIPTION FACTOR BINDING SITES IN PROMOTERS OF TARGET GENES

usage: 
    $ python A3_JMoggridge.py <genes.txt> <promoters.txt> <gff3_file> <fasta_folder/>

example:
"$ python A3_JMoggridge.py zea_mays_genes.txt promoters.txt Zea_mays.B73_RefGen_v4.48.gff3.gz fasta/"

* ensure that files and fasta folder are in same working directory as this script

How it works:
    inputs:
        file of genes identifiers for promoters to search
        file of TFBS motifs to search for in promoters
        set of gff3 files to get positional info for genes' transcripts
        set of fasta files to extract promoter regions from for each gene
        
    processes:
        1. Parse genes, motifs, gff annotations
        1b - generate an equal sized set of random genes to act as decoys
        2. Get position info for 5' most transcript of each gene (mind strand sense)
        3. Extract promoter sequence for each gene
        4. Search promoters and count TFBS motif occurences
        5. Return two files with TFBS counts in (targets, decoys) for each TFBS
    
    outputs:
        'Genes_hits.txt' and 'Random_hits.txt' -> list of motif & count
"""

# Libraries  --------------------------------------------------------------- ##
import sys
import os
import gzip
import math
import re
from random import sample
from Bio import SeqIO

# Inputs  ------------------------------------------------------------------ ##

# take input names from arguments 
genes_file = sys.argv[1]
promoters_file = sys.argv[2]
gff3_file = sys.argv[3]
path_to_fastas = sys.argv[4]
#    
#genes_file = 'zea_mays_genes.txt' # test
#promoters_file = 'promoters.txt' # test
#gff3_file = 'Zea_mays.B73_RefGen_v4.48.gff3.gz' # test
#path_to_fastas = 'fasta/' # test

## Functions --------------------------------------------------------------- ##

def parse_promoters(promoters_file):
    """Parses motifs file to list"""
    with open(promoters_file) as f:
        TFBS = [p.strip() for p in f.readlines()]
    print(f'\tRead {len(TFBS)} TFBS queries from: "{promoters_file}"')  
    return(TFBS)

def parse_genes(genes_file):
    """Parses target genes to list"""
    with open(genes_file) as f:
        Genes = [g.strip() for g in f.readlines()]
        print(f'\tRead {len(Genes)} target genes from: "{genes_file}"')
    return(Genes)

def parse_gff3s(gff3_file):
    """Parses all gff3 files, returns single list with mRNA annotations"""
    # parse gff3 file to list of lists (each line tab-separated)
    with gzip.open(gff3_file, 'rt') as infile:
        gff3 = [line for line in infile.readlines() if line[0] != '#']
    gff3 = [line.strip().split('\t') for line in gff3]
    # keep only mRNA annotations
    all_mRNAs = [line for line in gff3 if line[2] == 'mRNA']

    print(f'\tRead {len(gff3)} gff3 annotation lines from: "{gff3_file}"')
    print(f'\tReturning {len(all_mRNAs)} mRNA annotations')
    return(all_mRNAs)

def get_random_decoy_genes(all_mRNAs, sample_size):
    """Get equal sized sample of random genes to analyze in parallel to targets"""
    # Create a list of all genes using parent gene from mRNA annotations
    # Avoid genes that are only found on a contig
    chromosomes = ['Mt','Pt'] + [str(x) for x in range(1,11)]
    all_Genes = []
    for mRNA in all_mRNAs:
        if mRNA[0] in chromosomes:
            parentGene = mRNA[-1].split(';')[1].strip('Parent=gene:')
            all_Genes.append(parentGene)
    # Remove any duplicate gene ids
    all_Genes = list(set(all_Genes))
    # Return random sample of gene ids equal size to target set
    return(sample(all_Genes, sample_size))

def get_TSS_dict(Genes, all_mRNAs):
    """
    For the given list of genes of interest (targets or decoys):
    Finds all mRNAs for each gene & get positional info.
    Then, keeps the 5'most mRNA, ie. smallest start for + strand
    genes, and largest end position for - sense genes
    Returns dictionary of gene: (TSS_dict)
    """
    print("\tFinding all transcription start sites for genes list")
    # a dict for all mRNAs from genes of interest
    gene_mRNAs_dict = {gene:[] for gene in Genes}
    # iterate over all mRNA annotations
    for mRNA in all_mRNAs:
        # check if annotation is for a gene of interest        
        mRNA_string =  '\t'.join(mRNA)
        for gene in Genes:
            if mRNA_string.find(gene) > 0:
                # store positional info for each mRNA from gene of interest
                chromo, start, end, strand = mRNA[0], int(mRNA[3]), int(mRNA[4]), mRNA[6]
                gene_mRNAs_dict[gene].append([chromo, start, end, strand])
    # keep only 5'-most mRNA start position
    TSS_dict = {}
    for gene in gene_mRNAs_dict:
        # if + strand, keep smallest 5' pos
        if gene_mRNAs_dict[gene][0][3] == '+':
            strand = '+'
            start = math.inf
            first = None
            for mRNA in gene_mRNAs_dict[gene]:
                if mRNA[1] < start:
                    start = mRNA[1]
                    first = mRNA
        # if - strand, keep largest 3' position
        elif gene_mRNAs_dict[gene][0][3] == '-':
            strand = '-'
            end = -math.inf
            for mRNA in gene_mRNAs_dict[gene]:
                if mRNA[2] > end:
                    end = mRNA[2]
                    first = mRNA
        TSS_dict[gene] = first
    return(TSS_dict)


def extract_promoters(TSS_dict, chromosome, seq):
    """Extracts promoter sequence for each TSS of interest on chromosome"""

    def reverse_complement(seq): 
        # get reverse complement of promoter seq on -strand
        complement = {x:y for x,y in zip('ACGTN', 'TGCAN')}
        return(''.join([complement[nt] for nt in seq[:: -1]]))
        
    def trim_DNA(dna):
        # if dna contains 100 N's in a row, returns tail of seq after run of Ns
        i = len(dna)-1
        consecutiveN = 0
        while i >= 0:
            if dna[i] == 'N':
                consecutiveN += 1
                if consecutiveN == 100:
                    return(dna[i + 100:])
            else: 
                consecutiveN = 0
            i -= 1
        return(dna)
    
    # init dict for {gene: promoter_sequence} pairs
    Promoters = {}
    # for each gene on current chromosome...
    for gene in TSS_dict:
        if (TSS_dict[gene][0] == chromosome):
            # if gene is + sense, take 500bp before TSS at start position
            if TSS_dict[gene][3] == '+':
                i = TSS_dict[gene][1]
                promoter_seq = seq[i-501: i-1]
            # if gene - sense, take 500bp after TSS at end and reverse complement
            else:
                i = TSS_dict[gene][2]
                promoter_seq = reverse_complement(seq[i-1: i+499])
            # trim promoter of any run of 100 Ns
            promoter_seq = trim_DNA(promoter_seq)
            # Store extracted promoter sequence
            Promoters[gene] = TSS_dict[gene] + [promoter_seq]
    return(Promoters)


def extract_all_promoters(fasta_path, genes_tss, random_tss):
    """
    wraps `extract_promoters()`; iterates over all fasta files; 
    returns a dicts of {gene: promoter sequence} for genes and randoms
    """
    # init dicts to capture promoter seqs for each gene of interest
    genes_promoters = {}
    random_promoters = {}
    # iterate over fasta.gz files and extract promoter regions for genes
    for file in os.listdir(fasta_path):
        if file.endswith('.fa.gz'):
            # parse fasta to label and sequence
            with gzip.open(fasta_path + file, 'rt') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    [chromosome, seq] = [record.id, str(record.seq)]
            print(f'\tfrom chromosome {chromosome}')
            # extract promoter seqs to dicts
            genes_promoters.update(extract_promoters(genes_tss, chromosome, seq))
            random_promoters.update(extract_promoters(random_tss, chromosome, seq))
    return(genes_promoters, random_promoters)


def count_tfbs(tfbsSet, genes_promoters, random_promoters):
    """Counts occurences of each TFBS in target promoters and decoy promoters.
    Returns |TFBS|*2 matrix of counts. 
    `TFBS` is the list from promoters.txt 
    `*_promoters` is a dictionary with a list for each gene, 
    where the promoter sequence is the last entry"""
    # concatenate promoter sequences for searching
    Gene_prom_str = '|'.join(x[-1] for x in genes_promoters.values())
    Rand_prom_str = '|'.join(x[-1] for x in random_promoters.values())
    # count occurrences of tfbs in each promoter
    TFBS_counts = [[tfbs, None, None] for tfbs in tfbsSet]
    for i in range(len(tfbsSet)):
        tfbs = tfbsSet[i].upper()
        TFBS_counts[i][1] = len(re.findall(tfbs, Gene_prom_str))
        TFBS_counts[i][2] = len(re.findall(tfbs, Rand_prom_str))
        
    return(TFBS_counts)


## Script ------------------------------------------------------------------ ##

# Parse all input data
print('\nParsing genes, motifs, and annotations files...')
Genes = parse_genes(genes_file)
TFBS = parse_promoters(promoters_file)
AllmRNAs = parse_gff3s(gff3_file)
del(genes_file, promoters_file, gff3_file)

# Generate a Random set of decoy genes
Random = get_random_decoy_genes(AllmRNAs, sample_size = len(Genes))

# Find TSSs for target Genes and Randoms decoys
print('\nFinding TSS for targets and decoys...')
GenesTSS = get_TSS_dict(Genes, AllmRNAs)
RandomTSS = get_TSS_dict(Random, AllmRNAs)
del(Genes, Random, AllmRNAs)

# Extract promoter seqs from fasta files
print('\nExtracting promoters from fasta files...')
(GenePromoters, RandomPromoters) = extract_all_promoters(path_to_fastas, GenesTSS, RandomTSS)
#del(path_to_fastas, GenesTSS, RandomTSS)

# Count occurences of each TFBS in promoter seqs for each set
print('Counting transcription factor binding sites')

TFBS_counts_matrix = count_tfbs(TFBS, GenePromoters, RandomPromoters)

## Output ------------------------------------------------------------------ ##

# Create output files Target_hits.txt and Decoys_hits.txt
with open('Genes_hits.txt', 'w') as outfile:
    for tfbs in TFBS_counts_matrix:
        outfile.write(tfbs[0] + '\t' + str(tfbs[1]) + '\n')
#        print(tfbs[0] + '\t' + str(tfbs[1]) + '\n')
with open('Random_hits.txt', 'w') as outfile:
    for tfbs in TFBS_counts_matrix:
        outfile.write(tfbs[0] + '\t' + str(tfbs[2]) + '\n')
#        print(tfbs[0] + '\t' + str(tfbs[2]) + '\n')

print("\nCreated output files 'Genes_hits.txt' and 'Random_hits.txt':")
print("\n------- Finished --------")

