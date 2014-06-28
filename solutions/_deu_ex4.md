`dexseq-count` contains minor changes to deal with the newly generated annotation file. In addition, if we aimed at calculating gene counts by adding up the exon counts obtained with this tool, we'd realise that the numbers are slightly higher than with `htseq-count`. This is because reads that map to more than one exon of the same gene (i.e. split/spliced reads) are counted separately by `dexseq-count`.


