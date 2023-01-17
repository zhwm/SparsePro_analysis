module load r/4.1.2
source ~/py3/bin/activate
cd ${1}_${2}_${3}
cd K${4}_W${5}

st2=`date +%s.%N`
python ../../susie_rss.py --zld zld --zdir . --N 353570 --K ${4} --save susie_rss --prefix susie_rss --verbose
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > susie_rss.time
