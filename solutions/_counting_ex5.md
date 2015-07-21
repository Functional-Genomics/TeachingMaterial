* Try the following accessor methods:
```rconsole
length(aln_chr21)
head(names(aln_chr21))
seqnames(aln_chr21)
strand(aln_chr21)
first(aln_chr21)
last(aln_chr21)
head(start(first(aln_chr21)))
```

* How many reads are properly paired?
```rconsole
table(isProperPair(aln_chr21))
```

* What is the percentage of reads that map to multiple locations?
```rconsole
t=table(names(aln_chr21))
head(t[t>1])
length(t[t>1])/length(t)*100

# let us inspect one of the multireads
aln_chr21[names(aln_chr21)=="ERR315494.10000261"]
```

* What information does the command `seqlevels(aln_chr21)` provide?
  It contains information on the chromosome names available in the BAM file.

