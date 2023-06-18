# Genome-wide simulation procedures

## Extract genotype and annotation information
```
./get_geno.sh CHR
```

## Simulate traits; perform GWAS; run SparePro
```
./get_simulation_ns.sh CHR K W
```

## run SparsePro+
```
./run_anno_ns.sh CHR K W ITE
```

## run PolyFun to obtain functional prior
```
./prepare_poly.sh CHR K W
```

## run PolyFun informed SparsePro
```
./run_poly_ns.sh CHR K W ITE
```
