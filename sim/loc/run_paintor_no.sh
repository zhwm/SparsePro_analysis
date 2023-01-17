source ~/py3/bin/activate
cd ${1}_${2}_${3}
cd K${4}_W${5}

mkdir PAINTOR_no
st2=`date +%s.%N`
~/utils/PAINTOR_V3.0/PAINTOR -in ./ -out PAINTOR_no -input input.files -Zhead z -LDname ld -mcmc
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > PAINTOR_no.time
