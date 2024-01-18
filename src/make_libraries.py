#!/usr/bin/env python

#Script to make config file for fastq sample

import os
import sys

from optparse import OptionParser

def main():
    usage = "USAGE: %prog -s [sample] -o [output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-s", "--sample", help="sample id")
    optparser.add_option("-o", "--out", help="output file")
    optparser.add_option("-r", "--rna", help="path of rna seq location")
    optparser.add_option("-a", "--atac", help="path of atac seq location")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.sample or not options.out:
        optparser.print_help()
        sys.exit(-1)

    sample_id = options.sample
    options.rnapath = options.rna
    options.atacpath = options.atac

    attr1 = [options.rnapath, sample_id, "Gene Expression"]
    attr2 = [options.atacpath, sample_id, "Chromatin Accessibility"] 
    line1 = ",".join(str(i) for i in attr1)
    line2 = ",".join(str(i) for i in attr2)
    
    out = open(options.out, "w")
    out.write("fastqs,sample,library_type\n")
    out.write(str(line1))
    out.write("\n")
    out.write(str(line2))
    out.write("\n")
            
    out.close()

if __name__=='__main__':
    main()
