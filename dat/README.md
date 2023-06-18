# Fine-mapping five biomarker of vital organs

### 0. Perform GWAS
```
./process_pheno.R
for i in eGFR FFR gammaGT glucose prate; do sed "s/TRAIT/$i/g" run_TRAIT.sh > run_$i\.sh; done
```

### 1. Perform fine-mapping with GWAS summary statistics with SparsePro-
```
for i in eGFR FFR gammaGT glucose prate; do for j in $(seq 1 22); do ./run_sparsepro.sh $i $j; done; done
```

### 2. Run POLYFUN based LDSC to calculate functional enrichment
```
for i in eGFR FFR gammaGT glucose prate; do ./prepare_poly.sh $i; done
```

### 3. Estimate functional enrichment with SparsePro
```
for i in eGFR FFR gammaGT glucose prate; do ./run_enrich_baseline.sh $i sparsepro; done
```

### 4. Integrate GWAS summary statistics and functional annotation with SparsePro+ and SparsePro+PolyFun
```
for i in eGFR FFR gammaGT glucose prate; do for j in $(seq 1 22); do ./run_baseline.sh $i $j; done; done
for i in eGFR FFR gammaGT glucose prate; do for j in $(seq 1 22); do ./run_poly.sh $i $j; done; done
```