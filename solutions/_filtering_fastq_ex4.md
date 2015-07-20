We can get the number of reads in each fastq file from the first table in the FastQC report.

We can also execute the following command if we want to automate the process:

```bash
grep '^@ERR' ERR315494_pe_1_filt1.fastq | wc -l
    # find the lines starting with "@" and count how many there are
    # 
```

Note that being too strict in some filtering steps might lead to an important loss of information. Discarding too many reads since the beginning is generally not a good option, and one can also rely on the mapping step to discard low quality reads.

[Filtering-Fastq-Ex4](https://github.com/Functional-Genomics/TeachingMaterial/blob/Cancer-Genomics-07-2015/doc/13.filtering_fastq.md#exercise-4-)
