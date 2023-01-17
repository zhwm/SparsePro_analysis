source ~/py3/bin/activate
cd ${1}_${2}_${3}
mkdir K${4}_W${5}
cd K${4}_W${5}
for i in $(seq 1 50); do python ../../sample_causal.py --sdir . --snp ../${1}_${2}_${3}.prune.in --K ${4} --aw ${5} ${5} ${5} ${5} 0 0 0 0 ${5} 0 --beta 0.5 --seed $i --anno ~/scratch/UKB_Geno/${1}.anno; done
for i in $(seq 1 50); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1}_${2}_${3} --simu-qt --simu-causal-loci C$i\.causal --simu-hsq ${4}e-4 --simu-rep 1 --out C$i; done
for i in $(seq 1 50); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1}_${2}_${3} --pheno C$i\.phen --fastGWA-lr --out C$i; done
