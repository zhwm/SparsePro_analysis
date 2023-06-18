source ~/py3/bin/activate
cd GWAS/${1}
python ../../sparsepro_ukb.py --ukb ~/utils/ukb/lst/${2}.lst --zdir ${1}_${2}.z --N $(cat ${1}.N) --save sparsepro --prefix ${1}_${2} --verbose --LDdir ~/scratch/UKBBLD/
