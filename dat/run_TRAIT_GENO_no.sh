#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 300G            # memory pool per process
#SBATCH -o run_TRAIT_GENO_sparsepro.out    # STDOUT
#SBATCH -t 24:00:00            # time (D-HH:MM)

source ~/py3/bin/activate
module load r/4.1.2
st1=`date +%s.%N`
python ../sim/gw/sparsepro_ukb.py --ukb ../ukb/GENO.lst --zdir TRAIT_GENO.z --N $(cat TRAIT.N) --save sparsepro --prefix TRAIT_GENO --verbose --LDdir ~/scratch/UKBBLD/
ed1=`date +%s.%N`
echo $(echo $ed1-$st1 | bc -l) > TRAIT_GENO_sparsepro.time

st2=`date +%s.%N`
python ../sim/gw/susie_ukb.py --ukb ../ukb/GENO.lst --N $(cat TRAIT.N) --save susie --prefix TRAIT_GENO --verbose --LDdir ~/scratch/UKBBLD/ --zdir TRAIT_GENO.z
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > TRAIT_GENO_susie.time
