# Genome-wide simulation procedures

## extract genotype and annotation information
```
./get_geno.sh CHR
```

## simulate traits and perform GWAS
```
./get_simulation.sh CHR K W
```

## run SparseRro-/SparsePro+/SuSiE
```
./run_sparsepro_susie_ukb.sh CHR K W ITE
```

## run PolyFun to obtain functional prior
```
./prepare_polyfun.sh CHR K W
```

## run PolyFun informed SuSiE
```
./run_susie_poly.sh CHR K W ITE
```
