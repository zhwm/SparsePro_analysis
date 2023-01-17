conda activate polyfun
cd ${1}/K${2}_W${3}
mkdir polyfun
cd polyfun

st1=`date +%s.%N`
for i in $(seq 1 22); do awk -v a=$i '{print a"."$3"."$4"."$5"\t"a"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9}' ../C$i\.fastGWA | grep -v POS; done | sed '1iSNP\tCHR\tPOS\tA1\tA2\tBETA\tSE' > poly_ss
python ~/utils/polyfun/munge_polyfun_sumstats.py --sumstats poly_ss --n 353570 --keep-hla --out poly_ss.parquet
python ~/utils/polyfun/polyfun.py --compute-h2-L2 --output-prefix test --sumstats poly_ss.parquet --ref-ld-chr ~/scratch/SparsePro/POLY/polysim/fake. --w-ld-chr ~/scratch/SparsePro/POLY/polysim/weights. --allow-missing
ed1=`date +%s.%N`
echo POLYFUN taking $(echo $ed1-$st1 | bc -l) seconds
