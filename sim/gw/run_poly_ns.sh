source ~/py3/bin/activate
cd ${1}/ns_K${2}_W${3}

python ../../sparsepro_poly.py --h2 --ukb sparsepro/C${4}.h2 --zdir C${4}_22.z --N 353570 --save sparsepro_annopoly --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --wdir polyfun/annopoly.22.snpvar_ridge_constrained.gz
python ../../sparsepro_poly.py --h2 --ukb sparsepro/C${4}.h2 --zdir C${4}_22.z --N 353570 --save sparsepro_cspoly --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --wdir polyfun/cspoly.22.snpvar_ridge_constrained.gz