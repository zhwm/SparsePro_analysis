source ~/py3/bin/activate
cd ${1}_${2}_${3}/K${4}_W${5}

st1=`date +%s.%N`
~/utils/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --in-files master --sss --n-causal-snps ${4}
ed1=`date +%s.%N`
echo $(echo $ed1-$st1 | bc -l) > finemap.time
