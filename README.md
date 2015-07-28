# RNA-seq data analysis practical

Note: Since this part spans only on one day, the tutorial could not been done thouroughly. 

*Do not worry*: the main points will be stressed during the course and this tutorial is enough detailed to be done by yourself after. If you have any question (on any part of the tutorial or on RNA-seq more generally) you can contact us: nf [at] ebi.ac.uk and mitra [at] ebi.ac.uk 

This tutorial will illustrate how to use standalone tools, together with R and Bioconductor for the analysis of RNA-seq data. We will also use one meta-pipeline [IRAP] (https://github.com/nunofonseca/irap). Keep in mind that this is a rapidly evolving field and that this document is not intended as a review of the many tools available to perform each step; instead, we will cover *one* of the many existing workflows to analyse this type of data.

We will be working with __a subset__ of a publicly available dataset from *Homo sapiens*, which is available:
* As raw data:
    * in the European Nucleotide Archive (ENA) [ERP003613](http://www.ebi.ac.uk/ena/data/view/ERP003613) and 
    * in ArrayExpress [E-MTAB-2836] (http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2836/) 
* but also as processed data 
    *  in the authors' website [Protein Atlas] (http://www.proteinatlas.org/humanproteome/tissue+specific) and 
    * [EBI Gene Expression Atlas] (http://www.ebi.ac.uk/gxa/experiments/E-MTAB-2836). 

For more information about this dataset please refer to the original publication ([Uhlen et al. Science (2015) - Tissue-based map of the human proteome. DOI: 10.1126/science.1260419] (http://dx.doi.org/10.1126/science.1260419))..  

The tools and R packages that we will be using during the practical are listed below (see [Software requirements](https://github.com/Functional-Genomics/TeachingMaterial#software-requirements)) and the necessary data files can be found [here](http://www.ebi.ac.uk/~mitra/courses/CG15/RNASeq.tar.gz). After dowloading and uncompressing the `tar.gz` file, you should have the following directory structure in your computer:

```
CancerGenomics
|-- DATA                      # data used for the practicals
|   |-- demultiplexing        # multiplexed data !!! optional only if you have time   
|   |-- eqtl                  # data used in the eqtl practical
|   |-- fastq                 # fastq files -> starting point
|   |-- mapped                # mapped data: BAM files
|   `-- QCreports             # precomputed QC report
|
`-- IRAP_example              # Directory setup for IRAP (raw_data +reference) + its output
    |-- data
    |   |-- contamination     # E.coli reference
    |   |-- raw_data          # fastq files
    |   |-- reference
    |   |   |--homo_sapiens   # All the reference files are in this directory
    `-- E-MTAB-2886           # output of IRAP
        | ...
```

You can also browse the files online and download only the needed material from [here](http://www.ebi.ac.uk/~mitra/courses/CancerGenomics)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

## Table of contents

1. **Dealing with raw data**
    1. [The FASTQ format](doc/11.fastq.md)
    2. [Quality assessment (QA)](doc/12.qa.md)
    3. [Filtering FASTQ files](doc/13.filtering_fastq.md)
    4. [Aligning reads to the genome](doc/15.aligning.md) (already processed - will not be run)
2. **Dealing with aligned data**
    1. [The SAM/BAM format](doc/21.bam.md)
    1. [Visualising aligned reads](doc/22.visualising.md) (optional)
    1. [Filtering BAM files](doc/23.filtering_bam.md)
    2. Gene-centric analyses:
        1. [Counting reads overlapping annotated genes](doc/24.counting.md)
            * With htseq-count
            * With R
            * Alternative approaches
        1. [Normalising counts](doc/25.normalising.md)
            * With RPKMs
            * With DESeq2
        1. [Differential gene expression](doc/26.de.md)
3. **Other topics - Not covered in the course**
     1. Dealing with raw data
           * [De-multiplexing](doc/14.demultiplexing.md)
     1. Exon-centric analyses:
           * [Differential exon usage](doc/27.deu.md)


## Software requirements
*Note: depending on the topics covered in the course some of these tools might not be used.*

* Standalone tools:
  * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [PRINSEQ](http://prinseq.sourceforge.net/)
  * [eautils](https://code.google.com/p/ea-utils/)
  * [samtools](http://sourceforge.net/projects/samtools/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)
  * [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

* Bioconductor packages:
  * [GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
  * [GenomicAlignments](http://www.bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)
  * [Rsamtools](http://www.bioconductor.org/packages/release/bioc/html/Rsamtools.html)
  * [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html)
  * [pasilla](http://www.bioconductor.org/packages/release/data/experiment/html/pasilla.html)
  * [DESeq](http://www.bioconductor.org/packages/2.13/bioc/html/DESeq.html) - only for some dependencies
  * [DESeq2](http://www.bioconductor.org/packages/2.13/bioc/html/DESeq2.html)
  * [DEXSeq](http://www.bioconductor.org/packages/2.13/bioc/html/DEXSeq.html)

* (Meta)Pipeline
  * [IRAP](https://github.com/nunofonseca/irap) - which will be used through [docker](http://docs.docker.com/) for this tutorial

## Other resources

### Course data
* [Complete course data](http://www.ebi.ac.uk/~mar/courses/CancerGenomics)

### Tutorials
* [Course materials available at the Bioconductor website](http://www.bioconductor.org/help/course-materials/)
* [Online training resources at the EBI website](http://www.ebi.ac.uk/training/online/course-list?topic%5B%5D=13&views_exposed_form_focused_field=)
* [R and Bioconductor tutorial by Thomas Girke](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual)
* Do not forget to check the documentation for the packages used in the practical!

### Cheat sheets
* [R reference card](http://cran.r-project.org/doc/contrib/Short-refcard.pdf)
* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)


## Aknowledgments
This tutorial has been inspired on material developed by Mar Gonzalez-Porta, Liliana Greger, Ângela Gonçalves, Nicolas Delhomme, Simon Anders and Martin Morgan, who we would like to thank and acknowledge. Many thanks to Claudia Calabrese  & Fatemeh Zamanzad for having participate in the elaboration of this course. Special thanks must go to Mar Gonzalez-Porta, Claudia Calabrese & Fatemeh Zamanzad with whom we have been teaching. 
