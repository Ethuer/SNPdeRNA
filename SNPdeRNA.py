import csv
##from sklearn import svm
##import numpy as np
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


import thread


def SAMresultsCount(samDict1,samDict2,samDict3,verbose = True):
    if verbose == True:
        if samDict3 > 0:
            print '[Status]  Bam files loaded with %s , %s , %s Preliminary SNPs read' %(len(samDict1), len(samDict2),len(samDict3))
            exit
        if samDict2 > 0:
            print '[Status]  Bam files loaded with %s , %s Preliminary SNPs read' %(len(samDict1), len(samDict2))
            exit
        if samDict1 > 0:
            print '[Status]  Bam files loaded with %s Preliminary SNPs read' %(len(samDict1))
            exit

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

parser.add_argument('-fastaOut',
                    dest='fastaOut',
                    required = False,
                    default = 'none',
                    help='Output a fasta file, that will have been masked for remapping',
                    metavar = 'FILE',
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
    
    # Create a conversion list between genes and chromosomes, for silencing of the fasta...
    convDict = {}

    for row in gffDict:
        print row[0],row[1]
        convDict[row[1]] = row[0]
##    convDict = 
    
    if (len(gffDict)) == 0:
        print "[STATUS] The gff/gtf file has 0 features,  check the files structure or the existence of the selected feature"
        exit 
    else:
        print "[STATUS] gff file with %s features" %(len(gffDict))


    # Starting SNPdeRNA,     First step,  open files
    SAMList = []

    # one BAM file is essential
    ## CATCH EXCEPTION for missing indexing file

    
    samfile1 = pysam.AlignmentFile("%s" %(args.bam1),"rb")
    samDict1 = {}
    SAMList = [samDict1]
    samDict1 = SAM2SNP(args.feature,fasta_raw,samfile1,gffDict,outfile,cutoff,args.spass, createIntermediate)
##    except:
##    print '[ERROR] Read the replicate 1 incorrectly, please verify that the format is correct, and the file exists'
##    exit
    
    try:
        samfile2 = pysam.AlignmentFile("%s" %(args.bam2),"rb")
        samDict2 = {}
        samDict2 = SAM2SNP(args.feature,fasta_raw,samfile2,gffDict,outfile,cutoff,args.spass, createIntermediate)
        SAMList.append(samDict2)
    except:
        samDict2 = 0
        print '[OPEN] No replicate 2 ?  Valiadation would greatly profit from a second bam file'
        pass

    try:
        samfile3 = pysam.AlignmentFile("%s" %(args.bam3),"rb")
        samDict3 = {}
        samDict3 = SAM2SNP(args.feature,fasta_raw,samfile3,gffDict,outfile,cutoff,args.spass, createIntermediate)
        SAMList.append(samDict3)
    except:
        samDict3 = 0
        pass

##    print len(samDict1)
##    for element, value in samDict1.items():
##        print element, value
    
    SAMresultsCount(samDict1,samDict2,samDict3,verbose = True)

    for element, val in samDict1.items():
        for sub_element,sub_val in val.items():
            print element, sub_element,sub_val
##        print element[0]

    # Now starts the SNP2Quant part
    # with the vcf like dictionaries loaded as samDict 1 to 3




    resDict = {}
    geneDict = {}
    genecount = 0
    posDict = {}
##    noiseDict = {}
##    outDict = {}


    writeout = True
    
##    neg_count = 0
##    pos_count = 0
##    wrong_count = 0
##    origin = 'SGD'
    
        
##    writevcfheader(writeout, outfile)


    # Populate the master dictionary
    # This creates an empty dictionary containing all SNPs 
    masterDict = populateMasterDict(samDict1,samDict2,samDict3) 

##    for karg, marg in masterDict.items():
##        print karg, karg[0]


    print '[STATUS] Creating and populating Dictionaries'
##    for repeat in range(0,20):
##        extendMasterDict(masterDict,inDict1,inDict2,inDict3)

    print '[STATUS] Master Dictionary with %s possible SNPs' %(len(masterDict))


    print '[STATUS] Initiating SNP classification and dispersion estimation'

    
    ### THe dicitonary changes with the function pass to a string, instead of List

    silenceDict = {}

    
    for silenceGene, silenceSNP in masterDict.items():
        for silencePOS, emptydata in silenceSNP.items():
            
            silenceDict[silencekey[1]] = silencekey[0]


    print len(silenceDict)
    # Start classifying SNPs as ASE, and mutate fasta

    if args.spass == 'No' and args.fastaOut != 'none' :
        print '[STATUS] Writing masked fasta file '
        print '[STATUS] Please use this as a new reference to eleminate alignment derived errors'
##        try:
        if args.fasta != 'none':
            with open('%s'%(args.fastaOut),'w') as fasta_raw:
                record_dict = SeqIO.index('%s' %(args.fasta), "fasta")
                print '[STATUS] Silencing the fasta sequence'
                MaskAFasta(silenceDict,record_dict, fasta_raw, convDict)
                print '[STATUS] Masked Fasta generated'
        if args.fasta =='none':
            print '[STATUS] No input FASTA file found'


    
    # lets get to the quantification

    # populate masterdict to randomized resampling ??

    
    

    if args.spass != 'No':
        print '[STATUS] Initiate Quantification'
        
        extendMasterDict(masterDict,samDict1,samDict2,samDict3)


        
    resDict = {}
    total = 100
    for masDictkey, masDictitems in masterDict.items():
##        print masDictitems
        resDict[masDictkey]={}
        elementcount = 0
        for element in masDictitems:

            resDict[masDictkey] = {elementcount: []}


            # Check if it is Noise
            NoiseProb = distFunNegBin(100,int(element),0.02)


            # Check if it is a full SNP
            negative = total - int(element)

            FullProb = distFunNegBin(100,negative,0.02)
                        
##            resDict[masDictkey] = {elementcount: ['FullSNP' ,FullProb]}
            # Check for Allele Specific expression
            ASEprob = pmfNegBin(total,element,0.5)
            
            resDict[masDictkey][elementcount].append(NoiseProb)
            resDict[masDictkey][elementcount].append(ASEprob)
            resDict[masDictkey][elementcount].append(FullProb)




            
            elementcount +=1

##            
##    for resDictkey, resDictitems in resDict.items():
##        print resDictkey, resDictkey[0]
        # now we run the auntification on those
        # each element is n/100 * 5
        # create a dictionary with the likelyhoods of joining, then add them up.
        # 0.014 for noise  density
        # 0.986 for Full SNPs
        # use density of inverse count  ( count against)
        # all else = ASE

        


    
##        except:
##            print '[STATUS] No masked Fasta generated'
##        

    
##    if args.spass != 'No':
##        


        # here goes also the test for full SNPs
####        
####        masterDict = likelihoodTest(element, masterDict)
####
####    count = 0
####    count_all = 0
####
####    for element, value in masterDict.items():
####        print element, value

        
##        if value[1] > 0.05:
##            print value
##            value.append('accepted')
##            count +=1
##            posDict[element] = value[0]
##            # outfile writer if firstpass,   not if secondpass
##        else :
##            value.append( 'rejected' )
##            count_all +=1
####        print element, value
##    print count, count_all
##
##

##
##    writtenlist = []
##    writecount = 0
##    for element, value in posDict.items():
##        # write the SNPs to file
##        if element in inDict1:
##            writtenlist.append(element)
####            element[0] == gene,   element[1] == position
##            outfile.writerow([element[0],element[1], inDict1[element][0],inDict1[element][1],inDict1[element][2],inDict1[element][3],inDict1[element][4],inDict1[element][5]])
##            resDict[element] = [inDict1[element][0],inDict1[element][1]]
##            writecount +=1               
##        if inDict2 != 0 and element in inDict2 and element not in writtenlist:
##            outfile.writerow([element[0],element[1], inDict2[element][0],inDict2[element][1],inDict2[element][2],inDict2[element][3],inDict2[element][4],inDict2[element][5]])
##            writtenlist.append(element)
##            resDict[element] = [inDict2[element][0],inDict2[element][1]]
##            writecount +=1
##        if inDict3 != 0 and element in inDict3 and element not in writtenlist:
##            outfile.writerow([element[0],element[1], inDict3[element][0],inDict3[element][1],inDict3[element][2],inDict3[element][3],inDict3[element][4],inDict3[element][5]])
##            writtenlist.append(element)
##            resDict[element] = [inDict3[element][0],inDict3[element][1]]
##            writecount+=1
##                
##
##
##
##    print '[STATUS] %s  SNPs written to file' %(writecount)
##
####    # now mask a fasta with this information
##    try:
##        if args.fasta != 'none':
##            with open('%s'%(args.gff),'r') as fasta_raw:
##                record_dict = SeqIO.index('%s' %(args.fasta), "fasta")
##                print '[STATUS] Silencing the fasta sequence'
##                MaskAFasta(resDict,record_dict)
##                print '[STATUS] Masked Fasta generated'
##        if args.fasta =='none':
##            print '[STATUS] No masked Fasta generated'
##    except:
##        print '[STATUS] No masked Fasta generated'
##        
####
##
##
##    

####class SNP(object):
####    """
####a Single Nucleotide polymorphism
####Has
####
####baseclass
####Attributes:
####Position  (gene, base)
####has Original reference Nucleotide
####has Alternative nucleotide
####has likelihood to be noise  (negbin(0.014)
####has likelihood to be homozygous (negbin (0.98))
####has likelihood to be Allele specific negbin (0.5)
####    """
####    def __init__(self,gene,position,ORIG,ALT,count,coverage,likelihoodNoise,likelihoodASE,likelihoodFullSNP):
####        self.gene = gene
####        self.position = position
####        self.orig = ORIG
####        self.alt = ALT
####        self.count = count
####        self.coverage = coverage
####        self.likelihoodNoise = likelihoodNoise
####        self.likelihoodASE = likelihoodASE
####        self.likelihoodFullSNP = likelihoodFullSNP
####
####
####        
####    def nucleotide(ORIG, ALT)
####        self
####    
