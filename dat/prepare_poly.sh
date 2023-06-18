date

st1=`date +%s.%N`

#conda activate polyfun

python ~/utils/polyfun/munge_polyfun_sumstats.py --sumstats ${1} --out poly.parquet

python ~/utils/polyfun/polyfun.py --compute-h2-L2 --output-prefix poly --sumstats poly.parquet --ref-ld-chr ~/scratch/SparsePro/POLY/baselineLF2.2.UKB/baselineLF2.2.UKB. --w-ld-chr ~/scratch/SparsePro/POLY/baselineLF2.2.UKB/weights.UKB. --allow-missing

ed1=`date +%s.%N`

echo POLYFUN taking $(echo $ed1-$st1 | bc -l) seconds

date
