#sh this.sh gene.bed fa.fa out
fa_file=$2
genebed=$1
outf=$3
bedtools flank -i $genebed -g $fa_file -l 2000 -r 0 -s > ${outf}_p.bed
