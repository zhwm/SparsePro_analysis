module load plink/1.9b_6.21-x86_64
source ~/py3/bin/activate
mkdir ${1}
cd ${1}
plink --bfile ~/scratch/UKB_Geno/${1} --chr ${1} --maf 0.001 --geno 0.1 --make-bed --out ${1}
plink --bfile ${1} --indep-pairwise 1000 80 0.1 --out ${1}
python ../get_anno.py --save ${1}.anno --anno ~/scratch/UKB_Geno/${1}.anno --bim ${1}.bim
