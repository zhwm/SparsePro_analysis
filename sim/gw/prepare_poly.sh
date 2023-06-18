conda activate polyfun
cd ${1}/ns_K${2}_W${3}
mkdir polyfun
cd polyfun

for i in $(seq 1 22); do awk -v a=$i '{print a"."$3"."$4"."$5"\t"a"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9}' ../C$i\.fastGWA | grep -v POS; done | sed '1iSNP\tCHR\tPOS\tA1\tA2\tBETA\tSE' > poly_ss
python ~/utils/polyfun/munge_polyfun_sumstats.py --sumstats poly_ss --n 353570 --keep-hla --out poly_ss.parquet
python ~/utils/polyfun/polyfun.py --compute-h2-L2 --output-prefix annopoly --sumstats poly_ss.parquet --ref-ld-chr ../../annopoly22/ --w-ld-chr ../../annopoly22/weights. --allow-missing
python ~/utils/polyfun/polyfun.py --compute-h2-L2 --output-prefix cspoly --sumstats poly_ss.parquet --ref-ld-chr ../../cspoly22/ --w-ld-chr ../../cspoly22/weights. --allow-missing
