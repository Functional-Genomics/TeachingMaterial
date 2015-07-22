```rconsole
split_counts=split(exon_counts_chr21, names(exon_counts_chr21) ) 
head(split_counts)
gene_counts_chr21=sapply( split_counts, function(x) sum(x) ) 
head(gene_counts_chr21)
```

