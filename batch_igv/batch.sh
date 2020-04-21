sam=`less sample| tr "\n" "\t"`  #样本名称
cmd=''

echo '' > get_bam.sh
echo '' > igv.batch

while read chr start end type;do
	gap=200
	newstart=`expr $start - $gap`
	newend=`expr $end + $gap`
	region=${chr}:${newstart}-${newend}
	for i in $sam;do
		echo "/NJPROJ2/DISEASE/share/Software/bin/igvtools count -w 10 --query $region /NJPROJ2/DISEASE/Proj/X101SC19102913-Z01_JiNan10li/Customzied/1.igv.PINK1.20191210/bam/${i}.${i}/${i}.final.bam ${i}.chr${chr}_${start}_${end}_${type}.wig novo37 
/NJPROJ2/DISEASE/share/Software/bin/igvtools toTDF ${i}.chr${chr}_${start}_${end}_${type}.wig  ${i}.chr${chr}_${start}_${end}_${type}.tdf novo37 
" >> get_bam.sh
		tmpstr="C:\\Users\\dell\\Desktop\\IGV\\${i}.chr${chr}_${start}_${end}_${type}.tdf"
		cmd+=','$tmpstr
	done
		
	newcmd=`echo $cmd | sed -E 's#^,##g' `
	cmd=''
	echo "new
genome hg19 
load $newcmd
snapshotDirectory C:\\Users\\dell\\Desktop\\IGV\\
maxPanelHeight 500
goto chr$chr:$newstart-$newend
region chr$chr $start $end
snapshot chr${chr}_$start-${end}-${type}.IGV.png" >>  igv.batch

done <  bedfile
		
