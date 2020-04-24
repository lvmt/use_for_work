#!/usr/bin/env python
#-*- coding:utf-8 -*-

class EnsembleId(object):

    def __init__(self, args):
        self.ensemble = args['ensemble']
        self.gene = args['gene']
        self.out = args['out']
        self.ensemble_dict = self.get_ensemble_dict()
        self.gene_set = self.get_gene_set()

    def get_ensemble_dict(self):
        ensemble_dict = {}
        with open(self.ensemble, 'r') as f:
            for line in f:
                if line.startswith('Gene'):
                    continue
                linelist = line.strip().split('\t')
                ensemble_dict[linelist[1]] = linelist[0]
        return ensemble_dict

    def get_gene_set(self):
        gene_set = set()
        with open(self.gene, 'r') as f:
            for line in f:
                gene_set.add(line.strip())
        return gene_set
    
    def get_out(self):
        un_find_gene = []     #  统计在background中未找到的gene
        with open(self.out, 'w') as o:
            o.write('Gene_id\tGene_name\n')
            for gene in self.gene_set:
                if gene in self.ensemble_dict:
                    o.write("{}\n".format('\t'.join((self.ensemble_dict[gene], gene))))
                else:
                    un_find_gene.append(gene)
        print("有{}个gene在background中没有找到".format(len(un_find_gene)))
    

def main():
    en = EnsembleId(args)
    en.get_out()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('ensemble')
    parser.add_argument('gene')
    parser.add_argument('out')

    args = vars(parser.parse_args())
    main()
