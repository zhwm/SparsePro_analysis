source ~/py3/bin/activate
cd ${1}/ns_K${2}_W${3}

python ../../sparsepro_ukb.py --h2 --ukb sparsepro/C${4}.h2 --zdir C${4}_22.z --N 353570 --save sparsepro_cs_all --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --aW sparsepro/cs.W1.0 --anno ../22.idx.cs
python ../../sparsepro_ukb.py --h2 --ukb sparsepro/C${4}.h2 --zdir C${4}_22.z --N 353570 --save sparsepro_anno_all --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --aW sparsepro/anno.W1.0 --anno ../22.idx.anno
python ../../sparsepro_ukb.py --h2 --ukb sparsepro/C${4}.h2 --zdir C${4}_22.z --N 353570 --save sparsepro_anno --prefix C${4} --verbose --LDdir ~/scratch/UKBBLD/ --aW sparsepro/anno.W1e-05 --anno ../22.idx.anno
