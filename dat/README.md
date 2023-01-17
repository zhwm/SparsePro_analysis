# Fine-mapping five biomarker of vital organs

## perform GWAS
```
./process_pheno.R
for i in eGFR FFR gammaGT glucose prate; do sed "s/TRAIT/$i/g" run_TRAIT.sh > run_$i\.sh; done
```

## perform fine-mapping with GWAS summary statistics with SparsePro- and SuSiE
```
for i in eGFR FFR gammaGT glucose prate; do for j in $(seq 1 22); do sed "s/TRAIT/$i/g" run_TRAIT_GENO_no.sh | sed "s/GENO/$j/g" > run_$i\.sh; done; done
```

## run POLYFUN based LDSC to calculate functional enrichment
```
for i in eGFR FFR gammaGT glucose prate; do ./run_POLY.sh $i\_fastgwalr.fastGWA > run_$i\.sh; done
```

## integrating GWAS summary statistics and functional annotation with SparsePro+ and SuSiE+POLY
```
for i in eGFR FFR gammaGT glucose prate; do for j in $(seq 1 22); do sed "s/TRAIT/$i/g" run_TRAIT_GENO_no.sh | sed "s/GENO/$j/g" > run_$i\.sh; done; done
```
