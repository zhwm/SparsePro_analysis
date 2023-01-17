source ~/py3/bin/activate
module load r/4.1.2
cd ${1}/K${2}_W${3}

python ../../matchss.py --rss C${4}.fastGWA --prefix C${4} --save . --idir ../../../../dat/ukb/ --CHR CHR --POS POS --A1 A1 --A2 A2 --BETA BETA --SE SE
st1=`date +%s.%N`
python ../../sparsepro_ukb.py --ukb ../../../../dat/ukb/22.lst --zdir C${4}_22.z --N 353570 --save sparsepro --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --anno ../22.idx.anno
ed1=`date +%s.%N`
echo $(echo $ed1-$st1 | bc -l) > C${4}_sparsepro.time

st2=`date +%s.%N`
python ../../susie_ukb.py --ukb ../../../../dat/ukb/22.lst --N 353570 --save susie_ukb --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --zdir C${4}_22.z
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > C${4}_susie_ukb.time

st3=`date +%s.%N`
python ../../sparsepro_ukb.py --ukb ../../../../dat/ukb/22.lst --zdir C${4}_22.z --N 353570 --save sparsepro_anno --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --anno ../22.idx.anno --aW sparsepro/C${4}.W1e-05
ed3=`date +%s.%N`
echo $(echo $ed3-$st3 | bc -l) > C${4}_sparsepro_anno.time

