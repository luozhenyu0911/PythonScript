dir=/home/ycai/hm/dev/stlfr_results/
batch=$1
for i in `find /home/ycai/hm/dev/stlfr_results/${1} -name config.yaml`
do
echo "$i" >> $1.batch_type.txt
grep -w sequence_type  $i>> $1.batch_type.txt
grep read_len $i>> $1.batch_type.txt
grep -w bc_len $i>> $1.batch_type.txt
grep -w minreads $i>> $1.batch_type.txt
done
echo " see $1.batch_type.txt"
