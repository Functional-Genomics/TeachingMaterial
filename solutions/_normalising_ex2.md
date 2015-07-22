In a typical differential expression analysis, the majority of the genes are assumed to be non-differentially expressed. For this reason, it is useful to examine boxplots of counts across samples in each dataset, both before and after normalization. An effective normalization should result in a stabilization of read count distributions across replicates.

```rconsole

# create a list to simplify the plotting step
l=list(
        split(counts, colnames(counts)),
        split(rpkm, colnames(rpkm)),
        split(deseq_ncounts, colnames(deseq_ncounts))
)
names(l)=c("raw counts", "normalised counts - RPKMs", "normalised counts - DESeq")

# visualise the data as a set of boxplots
        pdf(file="./normalising_ex2a.pdf", width=10)
        par(mfrow=c(1,3), las=2, mar=c(7,5,3,3), cex=1)
        boxplot(l[[1]], outline=F, ylab="raw counts", ylim=c(0, 100))
        boxplot(l[[2]], outline=F, ylab="RPKMs", ylim=c(0, 300))
        boxplot(l[[3]], outline=F, ylab="normalised counts - DESeq", ylim=c(0, 100))
        dev.off()
        
# since we have only the counts for the chromosome 21, 
# the number of 0 in the data is too great to be able to visualise anything.
# we are going to strip all the rows that have no expression.

counts_no0 = counts[rowSums(counts) != 0, ]
rpkm_no0 = rpkm[rowSums(rpkm) != 0, ]
deseq_ncounts_no0 = deseq_ncounts[rowSums(deseq_ncounts) != 0, ]

# we plot only the rows with expressions
l=list(
        split(counts_no0, colnames(counts)),
        split(rpkm_no0, colnames(rpkm)),
        split(deseq_ncounts_no0, colnames(deseq_ncounts))
)

# visualise the data as a set of boxplots
        pdf(file="./normalising_ex2b.pdf", width=10)
        par(mfrow=c(1,3), las=2, mar=c(7,5,3,3), cex=1)
        boxplot(l[[1]], outline=F, ylab="raw counts", ylim=c(0, 100))
        boxplot(l[[2]], outline=F, ylab="RPKMs", ylim=c(0, 300))
        boxplot(l[[3]], outline=F, ylab="normalised counts - DESeq", ylim=c(0, 100))
        dev.off()
```

