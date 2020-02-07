# ORthology Inference using Synteny (ORIS) Pipeline

## Example workflow for orthology inference between two species A and B

## Step 1: Homology search

- Download [GenomeHistory 2.0](http://conantlab.org/GenomeHistory/GenomeHistory.html)
- Run GenomeHistory three times: genome A against itself, genome B against itself, and genome A versus genome B.

## Step 2: 
- Run make_orthology_inf_files.pl
- Will need gff files for each genome.

## Step 3: Orthology inference
- Download map_orthology<br>
```cd map_orthology```<br>
```make```<br>
