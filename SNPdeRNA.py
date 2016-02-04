import csv
from sklearn import svm
import numpy as np
import sys,argparse
import os.path
from _functions import *
from _functions_SAM2SNP import SAM2SNP
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq
import re
import csv
import sys,argparse
import os.path

# Unification of the 2 remaining scripts,
# add the SAM2SNP first, move as much as possible to _functions

parser = argparse.ArgumentParser(description='Looking for SNPs in RNAseq data. This part accepts SAM or BAM input files, mapped against an available reference. In a first step it counts the occurrances of Polymorphisms against a reference. In second pass mode, only already detected SNPs are being analyzed')



parser.add_argument('-fasta',
                    dest='fasta',
                    required = True,
                    help='Input the Fasta file of one of the Parentals or the reference for the species',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-gtf',
                    dest='gtf',
                    required = True,
                    help='Input a gtf file to check for the SNPS on gene, ### Try to remove necessity for this, there is no technical need for it',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-origin',
                    dest='origin',
                    required = False,
                    default = 'SGD',
                    help='Input origin of the GTF file   valid options are SGD/CGD  ENSEMBL  ', #  Code a SNIFFER for this one....
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )


parser.add_argument('-bam1',
                    dest='bam1',
                    required = True,
                    default = 0,
                    help='Input the bam file of the Alignment.  Tophat2 output works,  NextGenMapper or STAR.',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-bam2',
                    dest='bam2',
                    required = False,
                    default = 0,
                    help='Additional Input of bam file of the Alignment.  Tophat2 output works,  NextGenMapper or STAR.',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-bam3',
                    dest='bam3',
                    required = False,
                    help='Additional Input of bam file of the Alignment.  Tophat2 output works,  NextGenMapper or STAR.',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )


parser.add_argument('-out',
                    dest='out',
                    required = False,
                    default='output.tab',
                    help='Output file,  so far in vcf-like  gene snp-position original_base snp_base count_of occurance accumulated_probability',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-feature',
                    dest='feature',
                    required = False,
                    default='gene',
                    help='give an arbitrary name to the first column, default chromosome',
                    metavar = 'string',
                    #type=argparse.FileType('w')
                    )

parser.add_argument('-secondpass',
                    dest='spass',
                    required = False,
                    default='No',
                    help='running a second pass ?',
                    metavar = 'Bool',
                    #type=argparse.FileType('w')
                    )

parser.add_argument('-minqual',
                    dest='minqual',
                    required = False,
                    default='1',
                    help='Minimum cutoff, a simple p value from a binomial test',
                    metavar = 'FLOAT',
                    #type=argparse.FileType('w')
                    )

args = parser.parse_args()

#####################################################

with open("%s" %(args.fasta), "rU") as fasta_raw, open("%s"%(args.gtf),"r") as gff_raw, open('%s' %(args.out),'w') as out_raw:
    print '[IMPORT] Loading BAM File %s' %(args.bam1)
    
    print '[IMPORT] Loading Fasta file %s' %fasta_raw
    print '[IMPORT] Loading GFF/GTF File %s ' %(gff_raw)
    print '[IMPORT] Loading Complete.'
    
    gff_file = csv.reader(gff_raw, delimiter = '\t')
    outfile = csv.writer(out_raw, delimiter = '\t')
    cutoff = 0.05
    origin = args.origin
    createIntermediate = False
    
    print '[STATUS] Initiating gff parsing'
    gffDict = parse_gff(gff_file, origin,args.feature)
    if (len(gffDict)) == 0:
        print "[STATUS] The gff/gtf file has 0 features,  check the files structure or the existence of the selected feature"
        exit 
    else:
        print "[STATUS] gff file with %s features" %(len(gffDict))


    # Starting SNPdeRNA,     First step,  open files
    SAMList = []


    # one BAM file is essential
    ## CATCH EXCEPTION for missing indexing file

    try:
        samfile1 = pysam.AlignmentFile("%s" %(args.bam1),"rb")
##    print len(samfile1)
        samDict1 = {}
        
        SAMList = [samDict1]
        samDict1 = SAM2SNP(args.feature,fasta_raw,samfile1,gffDict,outfile,cutoff,args.spass, createIntermediate)
        
    except :
        print '[ERROR] Read the replicate 1 incorrectly, please verify that the format is correct, and the file exists'
        exit

        
    
    
    try:
##        with open('%s' %(args.bam2),'r') as in_2raw:
        samfile2 = pysam.AlignmentFile("%s" %(args.bam2),"rb")
        samDict2 = {}
        samDict2 = SAM2SNP(args.feature,fasta_raw,samfile2,gffDict,outfile,cutoff,args.spass, createIntermediate)
        SAMList.append(samDict2)
    except:
        samDict2 = 0
        print '[OPEN] No replicate 2 ?  Valiadation would greatly profit from a second bam file'
        pass

    try:
##        with open('%s' %(args.bam3),'r') as in_3raw:
        samfile3 = pysam.AlignmentFile("%s" %(args.bam3),"rb")
        samDict3 = {}
        samDict3 = SAM2SNP(args.feature,fasta_raw,samfile3,gffDict,outfile,cutoff,args.spass, createIntermediate)
        SAMList.append(samDict3)

    except:
        samDict3 = 0
##        print '[OPEN] No replicate 3 ?  no problem   This should not heavily decrease predictive power'
        pass

    print len(samDict1)
    for element, value in samDict1.items():
        print element, value
