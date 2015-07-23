#!/bin/bash

# IMPORTANT: always leave two empty spaces at the end of each md file,
#			 otherwise pandonc doesn't distinguish the sections properly
#
# human readable command:
#pandoc *.md -o ../pdf/practical.pdf --toc 
##          --variable title:"RNA-seq data analysis practical - Cancer Genomics Workshop" 
##          --variable date:"2015/07/23"
##          --variable author:"Mitra P. Barzine and Nuno A. Fonseca, Claudia Calabrese and Fatemeh Ghavidel"
##          --variable links-as-notes 
##          --variable linkcolor:black 
##          --variable urlcolor:black 
##          --variable geometry:margin=3cm

pandoc *.md --latex-engine=xelatex -o ../pdf/practical.pdf --toc --variable title:"RNA-seq data analysis practical - Cancer Genomics Workshop" --variable date:"2015/07/23" --variable author:"Mitra P. Barzine, Nuno A. Fonseca, Claudia Calabrese and Fatemeh Ghavidel" --variable links-as-notes --variable linkcolor:black --variable urlcolor:black --variable geometry:margin=3cm

