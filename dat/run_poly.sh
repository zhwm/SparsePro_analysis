source ~/py3/bin/activate
cd GWAS/${1}
python ../../sparsepro_poly.py --h2 --ukb sparsepro/${1}_${2}.h2 --zdir ${1}_${2}.z --N $(cat ${1}.N) --save sparsepro_poly --prefix ${1}_${2} --verbose --LDdir ~/scratch/UKBBLD/ --wdir polyfun/poly.${2}.snpvar_ridge_constrained.gz
