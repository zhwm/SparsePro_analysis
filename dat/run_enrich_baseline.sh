source ~/py3/bin/activate
cd GWAS/${1}/${2}

python ../../../enrich_ukb.py --pip ${1}_ --anno ../../../annotation/baseline/ --annosuf baseline --prefix ${1}_baseline --save . --pthres 1e-5
