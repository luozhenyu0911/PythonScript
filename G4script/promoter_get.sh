#sh this.sh gene.bed fa.fa outprefix leftnum rightnum
fa_file=$2
genebed=$1
outffix=$3
left=$4
right=$5
bedtools flank -i $genebed -g $fa_file -l $left -r $right -s > ${outffix}.bed
