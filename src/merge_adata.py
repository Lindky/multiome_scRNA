#!/usr/bin/env python

#Script to make config file for fastq sample

import os
import sys
import anndata as ad
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -r [rna data] -o [output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-o", "--out", help="output file")
    optparser.add_option("-r", "--rna", help="path of rna seq location")
    optparser.add_option("-a", "--atac", help="path of atac seq location")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.rna or not options.out:
        optparser.print_help()
        sys.exit(-1)

    options.atacpath = options.atac
    
    rna_adata = ad.read_h5ad(options.rna)
    atac_adata = ad.read_h5ad(options.atac)
  
    combined = ad.concat([rna_adata, atac_adata])
    combined.write("combined.h5ad", compression="gzip")
if __name__=='__main__':
    main()
