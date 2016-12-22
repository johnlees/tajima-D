# tajima-D
calculates tajima's D from bcftools query input

## Generate input
Use a command such as
```
bcftools norm -r FM211187:186-1547 -m - /lustre/scratch108/bacteria/jl11/mappings/23FSpn_new/recalibrated.snps.indels.vcf.gz | bcftools query -f '[%GT,]\n' | sed 's/,$//' > cds_1.csv
```

## Timing
took 0.25s on 1144 samples, 1 gene (34 variants)

