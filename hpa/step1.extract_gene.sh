 ##extract gene for infile
    for name in `less  genename `;do
        filter_gene -i  infile -o ${name}.filter.xls -gl ${name} 
    done