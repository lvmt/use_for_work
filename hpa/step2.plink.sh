
    for name in `less  genename`;do
        plink --vcf ${name}.filter.vcf  --recode  --make-bed  --out ${out}.filter.vcf.plink
    done 