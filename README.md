# TFmotif-counter
Python script to count the occurence of transcription factor binding sites (aka. TFBS or motifs) in the promoter regions of query genes; also provides a count for same motifs in an equal sized random sample of promoters.

## Quick overview
- Uses python3 with the libraries `sys`, `os`, `gzip`, `Bio`, `random`, `math`, and `re`. `Bio` handles io of sequence data.
- The inputs are a list of DNA patterns (motifs) to find in promoters (eg. `promoters.txt` file) and a list of genes to search (eg. `zea_mays_genes.txt`).
- As is, you need to download the gff (eg.`Zea_mays.B73_RefGen_v4.48.gff3.gz`) and fasta.gz genome assembly files (could work on other assemblies, haven't tried though). *Put inpu  files and fasta folder in same working directory as this script*.

- command line usage (unix):   
  `python tf_motif_count.py <genes.txt> <promoters.txt> <gff3_file> <fasta_folder/>`

- example:  
  `python tf_motif_count.py zea_mays_genes.txt promoters.txt Zea_mays.B73_RefGen_v4.48.gff3.gz fasta/"`

----


## How it works:
  
### inputs:
  file of genes identifiers for promoters to search
  file of TFBS motifs to search for in promoters
  set of gff3 files to get positional info for genes' transcripts
  set of fasta files to extract promoter regions from for each gene

### code:
  1. Parse genes, motifs, gff annotations
  1b - To act as decoys, generate an equal sized set of random genes from all gene names in gff3 file that will provide a control for comparing against our genes of interest.
  2. Get position info for 5' most transcript of each gene using the mRNA entry for the gene (mind strand sense: use 3' end of mRNA if gene is on the - strand)
  3. Extract promoter sequence for each gene from the fasta files using the position info from step 2.
  4. Search promoters (from step 3) and count TFBS motif occurences in the input list. Report the sum for each motif.
  5. Return two files with the TFBS counts (from step 4), one for the set of query genes and one for the set of random genes.

### outputs:
  'Genes_hits.txt' and 'Random_hits.txt' -> list of motif & count
