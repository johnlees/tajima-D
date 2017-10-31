# tajima-D
calculates tajima's D from bcftools query input

## Generate input
Use a command such as
```
bcftools view -f PASS -r FM211187:186-1547 snps.indels.vcf.gz | bcftools norm -m - | bcftools view -e 'ALT[*]==\"*\"' | bcftools query -f '%POS,[%GT,]\n' | sed 's/,$//' > cds_1.csv
```

## Timing
took 0.25s on 1144 samples, 1 gene (34 variants)

