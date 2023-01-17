source ~/py3/bin/activate
cd ${1}
mkdir K${2}_W${3}
cd K${2}_W${3}
for i in $(seq 1 22); do python ../../sample_causal.py --sdir . --snp ../${1}.prune.in --K ${2} --aw ${3} ${3} ${3} ${3} 0 0 0 0 ${3} 0 --beta 0.5 --seed $i --anno ~/scratch/UKB_Geno/${1}.anno; done
for i in $(seq 1 22); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --simu-qt --simu-causal-loci C$i\.causal --simu-hsq 0.01 --simu-rep 1 --out C$i; done
for i in $(seq 1 22); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --pheno C$i\.phen --fastGWA-lr --out C$i; done
