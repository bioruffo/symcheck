# README #

# What is symcheck #

This script verifies the existence of gene symbols in RefSeq and HGNC databases. It can also verify whether a gene symbol is listed as mitochondrial gene, whether the symbol is also an alias for other genes, and whether the gene is listed in any alternate/unplaced contig besides the standard chromosome scaffolds.

# How to install symcheck #

1 . Requisites:  

 * Python 3.5+ (I suggest [Anaconda Python](https://www.continuum.io/downloads))  

 
2 . Required data files:  
 
 The repository includes two files, **"NCBI_RefSeq_hg19_UCSC.tsv"** and **"HGNC_results.tsv"** that can serve as sources for RefSeq data and HGNC data, respectively. They are provided as-is and their inclusiveness may vary.  
If users wish to download new data themselves, this can be accomplished as follows:  
 * RefSeq data can be downloaded from the [UCSC Table Browser](genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=0&hgta_regionType=genome&hgta_outputType=primaryTable).  
 * HGNC official gene symbols can be downloaded from [HGNC BioMart](https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__hgnc_gene_id_1010%2Chgnc_gene__status_1010%2Chgnc_gene__approved_symbol_1010%2Chgnc_gene__approved_name_1010%2Chgnc_gene__chromosome_location_1010%2Chgnc_gene__hgnc_alias_symbol__alias_symbol_108%2Chgnc_gene__chromosome_1010%2Chgnc_gene__hgnc_alias_name__alias_name_107%2Chgnc_gene__hgnc_previous_symbol__previous_symbol_1012%2Chgnc_gene__hgnc_previous_name__previous_name_1011%2Chgnc_gene__locus_group_1010%2Chgnc_gene__locus_type_1010%2Chgnc_gene__hgnc_family__hgnc_family_id_109%2Chgnc_gene__hgnc_family__hgnc_family_name_109%2Chgnc_gene__date_approved_1010%2Chgnc_gene__date_modified_1010%2Chgnc_gene__date_symbol_changed_1010%2Chgnc_gene__date_name_changed_1010).


3 . Download or clone the source code for `symcheck`. The Python module can be ran directly from source.  

# Using symcheck #

An example use of `symcheck` is:  

~~~~
python symcheck.py  -maew -r ./NCBI_RefSeq_hg19_UCSC.tsv -n ./HGNC_results.tsv -s 'TP53 NOTCH1 COX1 MT-CO2 LIG1 MKL1'
~~~~

`symcheck` also accepts text files containing a list of gene symbols to check:

~~~~
python symcheck.py  -maew -r ./NCBI_RefSeq_hg19_UCSC.tsv -n ./HGNC_results.tsv -f ./example_genes.txt
~~~~

# Output file #
`symcheck` outputs the file `output.html`.

# LICENSE #

**symcheck is offered under the MIT License.**  
**Please read it here: https://opensource.org/licenses/MIT**  

