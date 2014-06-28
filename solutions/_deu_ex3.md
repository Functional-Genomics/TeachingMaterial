Similarly to Phred scores, mapping qualities reflect the probability of a wrong alignment. However, despite of the fact that they contain the same information, their range can vary depending on the mapping tool used.
In this example, and considering that we have mapped the reads using TopHat, *DEXSeq* will only consider those reads that map uniquely (i.e. mapping quality = 255; see [previous exercise](../solutions/_bam_ex3.md) on mapping qualities).

