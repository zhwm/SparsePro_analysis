# Locus simulation procedures

## extract genotype and annotation information
```
./get_block.sh CHR ST ED
```

## simulate traits and perform GWAS
```
./get_simulation.sh CHR ST ED K W
```

## prepare summary statistics for running fine-mapping
```
./prepare_ss.sh CHR ST ED K W
```

## run FINEMAP
```
./run_finemap.sh CHR ST ED K W
```

## run PAINTOR-
```
./run_paintor_no.sh CHR ST ED K W
```

## run PAINTOR+
```
./run_paintor_anno.sh CHR ST ED K W      
```

## run SuSiE
```
./run_susie_rss.sh CHR ST ED K W
```

## run SparsePro- and SparsePro+
```
./run_SP_all.sh CHR ST ED K W
```
