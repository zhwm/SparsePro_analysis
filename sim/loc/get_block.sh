module load plink/1.9b_6.21-x86_64
source ~/py3/bin/activate
mkdir ${1}_${2}_${3}
cd ${1}_${2}_${3}
plink --bfile ~/scratch/UKB_Geno/${1} --chr ${1} --from-bp ${2} --to-bp ${3} --maf 0.001 --geno 0.1 --make-bed --out ${1}_${2}_${3}
plink --bfile ${1}_${2}_${3} --indep-pairwise 1000 80 0.1 --out ${1}_${2}_${3}
plink --bfile ${1}_${2}_${3} --matrix --r --out ${1}_${2}_${3}
python ../get_anno.py --save ${1}_${2}_${3}.anno --anno ~/scratch/UKB_Geno/${1}.anno --bim ${1}_${2}_${3}.bim
