Using the information provided in the first plot from the report (*per base sequence quality*), we might decide to trim the last nucleotides of all reads, due to the decrease in quality that we observe. 

In addition, we know that if we want to filter reads based on average quality, our threshold can be high, given the overall good quality of the data. 

You can also check the table under the *overrepresented sequences* section. It is often the case that there are suspicious sequences of polyA, which we can get rid of by setting a reasonable dust score. 

Similarly, you might find that some overrepresented sequences correspond to adapters and thus you might want to remove them. Further information on the filtering possibilities can be found in the section [Filtering FASTQ files](../doc/13.filtering_fastq.md).

