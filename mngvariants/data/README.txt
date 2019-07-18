In this folder you will find tabular format files (.txt), variant calling files (.vcf), reference fasta (sequences.fa), reference annotations (genes.gbk), regions of no coverage (.bed) and bam/bam index files (.bam, .bai).

If you want to look at the SNPs/small indels in a table you can open up the .txt files in something like Microsoft Excel or Google Sheets. These text files tell you the position of the mutation, whether it is in a gene and what effect it has on the gene. We actually give you two types of SNP calling, sensitive (there could be false negatives, i.e. you don't catch all of your SNPs, at least 90% of reads have this variant) and specific (could be false positives, i.e. you will catch all of your true SNPs but some of the calls might be wrong, at least 10% of the reads have this variant).

When you open the text file you will find that it describes the location of the variant, and if it is in a coding region the effect that the variant has. So for example synonymous coding means that the mutation results in the same amino acid, non synonymous mutation means that a different amino acid was formed.

If you want to look at your SNPs in a genome browser I would suggest using the following:

https://www.broadinstitute.org/software/igv/download

Select which type you want to launch (based on your computer and operating system).

When you have launched IGV select 'Genomes' at the top and then 'Load genomes from file' choose the sequences.fa file. Then select 'File' at the top and then 'Load from file'. Choose all of the .vcf files (variants), the genes.gbk file and the .bed files. You can also view how the reads align against the reference by opening up the .bam files in IGV (please note that these can be very big and use up a lot of memory).

Once you have loaded your VCFs and your annotation you will need to zoom in (top right). To move between variants (SNPs and Indels) select one of the tracks and the use Ctrl-F to move forward through the variants.
