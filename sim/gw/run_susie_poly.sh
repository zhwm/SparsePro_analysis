source ~/py3/bin/activate
module load r/4.1.2
cd ${1}/K${2}_W${3}

st2=`date +%s.%N`
python ../../susie_poly.py --ukb ../../../../dat/ukb/22.lst --N 353570 --save susie_poly --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --zdir C${4}_22.z --snpvar polyfun/test.${1}.snpvar_ridge_constrained.gz
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > C${4}_susie_poly.time
