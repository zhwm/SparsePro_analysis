source ~/py3/bin/activate
cd ${1}_${2}_${3}
cd K${4}_W${5}
python ../../get_master.py --save . --N 353570 --rep 50 --ld ../${1}_${2}_${3}.ld --anno ../${1}_${2}_${3}.anno
for i in $(seq 1 50); do sed 1d C$i\.fastGWA | awk '{print $2" "$1" "$3" "$4" "$5" "$7" "$8" "$9}' | sed '1irsid chromosome position allele1 allele2 maf beta se' > C$i\.z; done

for i in $(seq 1 50); do echo C$i; done > input.files
for i in $(seq 1 50); do sed 1d C$i\.fastGWA | awk '{print $8/$9}' | sed '1iz' > C$i; done
for i in $(seq 1 50); do cut -f 2- ../${1}_${2}_${3}.anno | tr '\t' ' ' > C$i\.annotations; done
for i in $(seq 1 50); do cp ../${1}_${2}_${3}.ld C$i\.ld; done

for i in $(seq 1 50); do sed 1d C$i\.fastGWA | awk '{print $2"\t"$8/$9}' > C$i\.txt; done
