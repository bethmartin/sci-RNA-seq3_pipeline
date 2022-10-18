for i in `ls ./nobackup/output/filtered_sam/*.sam`; do samtools view -c $i >> tmp1; done
echo "part 1 of 3 done"
wc -l tmp1
for i in `ls ./nobackup/output/rmdup_sam/*.sam`; do samtools view -c $i >> tmp2; done
echo "part 2 of 3 done"
wc -l tmp2
for i in `ls ./nobackup/output/rmdup_sam_2/*.sam`; do samtools view -c $i >> tmp3; done
echo "part 3 of 3 done"
wc -l tmp3
paste tmp1 tmp2 tmp3 > read_num.txt
#rm tmp1 tmp2 tmp3
echo "results written to read_num.txt"


