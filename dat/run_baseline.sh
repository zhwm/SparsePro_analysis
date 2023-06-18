source ~/py3/bin/activate
cd GWAS/${1}

python ../../sparsepro_ukb.py --h2 --ukb sparsepro/${1}_${2}.h2 --zdir ${1}_${2}.z --N $(cat ${1}.N) --save sparsepro_baseline --prefix ${1}_${2} --verbose --LDdir ~/scratch/UKBBLD/ --anno ../../annotation/baseline/${2}.baseline --aW sparsepro/${1}_baseline.W1e-05
