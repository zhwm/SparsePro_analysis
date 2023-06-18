source ~/py3/bin/activate
cd ${1}/ns_K${2}_W${3}
for i in $(seq ${4} ${4}); do python ../../sample_causal.py --sdir . --snp ../${1}.prune.in --K ${2} --aw 0 0 0 0 0 0 0 0 ${3} 0 --beta 0.5 --seed $i --anno ~/scratch/UKB_Geno/${1}.anno; done
for i in $(seq ${4} ${4}); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --simu-qt --simu-causal-loci C$i\.causal --simu-hsq 0.01 --simu-rep 1 --out C$i; done
for i in $(seq ${4} ${4}); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --pheno C$i\.phen --fastGWA-lr --out C$i; done
python ../../matchss.py --rss C${4}.fastGWA --prefix C${4} --save . --idir  ~/utils/ukb/idx/ --CHR CHR --POS POS --A1 A1 --A2 A2 --BETA BETA --SE SE
for i in $(seq ${4} ${4}); do python ../../sparsepro_ukb.py --ukb ~/utils/ukb/lst/22.lst --zdir C$i\_22.z --N 353570 --save sparsepro --prefix C$i --verbose --LDdir ~/scratch/UKBBLD/ ; done
