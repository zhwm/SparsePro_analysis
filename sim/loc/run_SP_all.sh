source ~/py3/bin/activate
cd ${1}_${2}_${3}
cd K${4}_W${5}

st2=`date +%s.%N`
python ../../sparsepro_zld.py --zld zld --zdir . --N 353570 --K ${4} --save sparsepro_no --prefix sparsepro_no --verbose --anno anno --pthres 1e-5
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > sparsepro_no.time

st2=`date +%s.%N`
python ../../sparsepro_zld.py --zld sparsepro_no/sparsepro_no.h2 --zdir . --N 353570 --K ${4} --save sparsepro_anno_all --prefix sparsepro_anno_all --verbose --anno anno --aW sparsepro_no/sparsepro_no.W1.0 --h2
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > sparsepro_anno_all.time

st2=`date +%s.%N`
python ../../sparsepro_zld.py --zld sparsepro_no/sparsepro_no.h2 --zdir . --N 353570 --K ${4} --save sparsepro_anno --prefix sparsepro_anno --verbose --anno anno --aW sparsepro_no/sparsepro_no.W1e-05 --h2
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > sparsepro_anno.time
