# ORthology Inference using Synteny (ORIS) Pipeline

## References

- Conant GC. 2009. [Neutral evolution on mammalian protein surfaces](https://www.cell.com/trends/genetics/fulltext/S0168-9525(09)00147-4). *Trends Genet.* **25**:377–381.<br>
- Bekaert M, Conant GC. 2011. [Transcriptional robustness and protein interactions are associated in yeast](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-5-62). *BMC Syst. Biol.* **5**:62.
- Bekaert M, Conant GC. 2014. [Gene duplication and phenotypic changes in the evolution of mammalian metabolic networks](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0087115). *PLoS One.* **9**:e87115.

## Example workflow for orthology inference between two genomes 1 and 2

### Step 1: Homology search

- Download [GenomeHistory 2.0](http://conantlab.org/GenomeHistory/GenomeHistory.html)
- Extract lists of gene names from fasta files. Use command such as:<br>
```grep ">" genome1.DNA.fa |egrep -o '^>(\w|\d|\.){1,20}' |egrep -o '(\w|\d|\.){1,20}' >genome1_genelist.txt```<br>
- Run GenomeHistory three times: genome 1 against itself, genome 2 against itself, and genome 1 versus genome 2.

### Step 2: Prepare input files
- Run make_orthology_inf_files.pl
- Will need gff files for each genome.

### Step 3: Orthology inference
- Download map_orthology<br>
```cd map_orthology```<br>
```make```<br>
