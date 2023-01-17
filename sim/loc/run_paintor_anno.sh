source ~/py3/bin/activate
cd ${1}_${2}_${3}
cd K${4}_W${5}

mkdir PAINTOR_anno
st1=`date +%s.%N`
~/utils/PAINTOR_V3.0/PAINTOR -in ./ -out PAINTOR_anno -input input.files -Zhead z -LDname ld -mcmc -annotations Conserved_LindbladToh,DHS_Trynka,H3K27ac_Hnisz,H3K4me3_Trynka,Transcr_Hoffman,TSS_Hoffman,UTR_3_UCSC,UTR_5_UCSC,non_synonymous,Human_Promoter_Villar_ExAC
ed1=`date +%s.%N`
echo $(echo $ed1-$st1 | bc -l) > PAINTOR_anno.time
