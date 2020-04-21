
    for name in `less genename`;do
        paste sample_info_ped ${name}.filter.vcf.plink.ped.tmp > ${name}.filter.HWE.ped
        awk -F "\t" -v OFS="\t" '{if($2 != ".") print $2,$4; else print $4,$4}' > ${name}.filter.HWE.info
    done 