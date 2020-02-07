# ORthology Inference using Synteny (ORIS) Pipeline

## Example workflow for orthology inference between two genomes 1 and 2

## Step 1: Homology search

- Download [GenomeHistory 2.0](http://conantlab.org/GenomeHistory/GenomeHistory.html)
- Extract lists of gene names from fasta files. Use command such as:<br>
```grep ">" genome1.DNA.fa |egrep -o '^>(\w|\d|\.){1,20}' |egrep -o '(\w|\d|\.){1,20}' >genome1_genelist.txt```<br>
- Run GenomeHistory three times: genome 1 against itself, genome 2 against itself, and genome 1 versus genome 2.

## Step 2: 
- Run make_orthology_inf_files.pl
- Will need gff files for each genome.

## Step 3: Orthology inference
- Download map_orthology<br>
```cd map_orthology```<br>
```make```<br>
