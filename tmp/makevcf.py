
def makevcf(infile, outfile):
    with open(infile, 'r') as f, open(outfile, 'w') as o:
        for line in f:
            linelist = line.strip().split('\t')
            need_str = '{chr}\t{pos}\t{ID}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{genetype}\n'.\
                format(chr=linelist[0], pos=linelist[1], ID=".", ref=linelist[3], alt=linelist[4], qual="999", filter="PASS",\
                    info=".", format="GT:PL:DP:AD", genetype="0/1:255,0,255:120:60")
            o.write(need_str)

if __name__ == "__main__":
    
    import sys
    args = sys.argv
    makevcf(args[1], args[2])

