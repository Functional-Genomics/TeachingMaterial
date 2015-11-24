* How many reads are properly paired?
```bash
samtools view -F 0x0002 untreated3.bam -bo - | samtools sort -n - untreated3_paired
samtools view untreated3_paired.bam | wc -l
# 13,212,182 / 2 = 6,606,091
```

* How would you sort the properly paired reads by name instead?
```bash
samtools sort -n untreated3_paired.bam ERR315494_paired_byname
samtools view untreated3_paired_byname.bam |head
```


* How would you sort the bam by coordinate instead? How would you index this bam?
```bash
samtools sort untreated3_paired.bam untreated3_paired_nsort
samtools index untreated3_paired_nsort.bam
```


* Which percentage of those properly paired reads map uniquely?
```bash
samtools view -q 255 untreated3_paired.bam | wc -l
# 7,417,488 / 2= 3,708,744
```

