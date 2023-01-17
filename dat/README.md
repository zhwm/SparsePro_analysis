# Fine-mapping five biomarker of vital organs

## perform GWAS
```
./process_pheno.R
for i in eGFR FFR gammaGT glucose prate; do sed "s/TRAIT/$i/g" run_TRAIT.sh > run_$i\.sh; done
```

## perform fine-mapping with GWAS summary statistics


##
