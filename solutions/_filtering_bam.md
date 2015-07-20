* How many reads are properly paired?
```bash
samtools view -F 0x0002 ERR315494.bam -bo - | samtools sort -n - ERR315494_paired
samtools view ERR315494_paired.bam | wc -l
```

* How would you sort the properly paired reads by name instead?
```bash
samtools sort -n ERR315494_paired.bam ERR315494_paired_nsort 
samtools view ERR315494_paired_nsort.bam |head
```

* Which percentage of those properly paired reads map uniquely?
```bash
samtools view -q 50 ERR315494_paired.bam | wc -l
```

