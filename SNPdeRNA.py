import csv
##from sklearn import svm
##import numpy as np
import sys,argparse
import os.path
import os
from _functions import *
from _functions_SAM2SNP import SAM2SNP
from _functions_SNP2Quant import FischerCombineP
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


##
##def FischerCombineP(ListOfP):
##    """
##
##    returns the probability of the combined p-values, according to the Fischer method for combining p values
##
##    """
##    summary = 0
##    for element in ListOfP:
##        arg = math.log(element)
##        arg = -2 * arg
##        summary = float(summary) + float(arg)
##
##    target =  float(summary)
##    df = len(ListOfP) -1
##
##    combined_prob = chisqprob(target,df)
##
##    return combined_prob


# Unification of the 2 remaining scripts,
# add the SAM2SNP first, move as much as possible to _functions

def writeVCFormat(masterDict,outfile,convDict):
    #write Header

    for gene, SNPs in masterDict.items():
        CHROM = convDict[gene]

        for SNP, info in SNPs.items():
            POS = SNP
            ID = '%s_%s'%(gene,POS)
            REF = info[0]
            ALT = info[1]
            QUAL = 30
            Filter = '.'
            INFO = 'NS=%s;DP=%s' %(info[5],info[2])
        
            outfile.writerrow([CHROM,POS,ID[0],REF,ALT,QUAL,Filter,INFO[0]])




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
    cutoff = float(args.minqual)
    origin = args.origin
    createIntermediate = False
    
    print '[STATUS] Initiating gff parsing'
    gffDict = parse_gff(gff_file, origin,args.feature)
    
    # Create a conversion list between genes and chromosomes, for silencing of the fasta...
    convDict = {}

    for row in gffDict:
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


    

##    try:
        
##        samfile1 = pysam.AlignmentFile("%s" %(args.bam1),"rb")
##        samDict1 = {}
##        SAMList = [samDict1]
##        samDict1 = SAM2SNP(args.feature,fasta_raw,samfile1,gffDict,outfile,cutoff,args.spass, createIntermediate)
    samDict1,SAMList = SamImport(args.bam1,SAMList,args.feature,fasta_raw,gffDict,outfile,cutoff,args.spass,createIntermediate)

##            ,feature,fasta_raw,gffDict,outfile,cutoff,spass,createIntermediate
##    except:
##        print '[ERROR] Read the replicate 1 incorrectly, please verify that the format is correct, and the file exists'
##        print '[ASSISTANCE] Most common error is missing index file from BAM input.  '
##        print '[ASSISTANCE] Attempting to create index file from here.'
##        try:
##            commandLine = 'samtools index %s' %(args.bam1)
##            os.system(commandLine)
##
##            # wait until this finishes
##            samDict1,SAMList = SamImport(args.bam1,SAMList)
##
##
##            # wait until this finishes, then continue
##
##        except:
##            print '[ERROR] Read the replicate 1 incorrectly, please verify that the format is correct, and the file exists'
##            exit
            
    
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

##    SAMresultsCount(samDict1,samDict2,samDict3,verbose = True)
##

    

    # Now starts the SNP2Quant part
    # with the vcf like dictionaries loaded as samDict 1 to 3




    resDict = {}
    geneDict = {}
    genecount = 0
    posDict = {}


    writeout = True

##def getSNPfromSAMs(gene,samDict):
##    SNPDict = {}
##    if gene in samDict:
##        for SNP, info in samDict[gene].items():
##            # info contains :
##            SNPDict[SNP] = [info[0],info[1],info[2],info[3],info[5]]
##    return SNPDict
####            ORIG,ALT,SNPcov,Basecov,DensfunProb,synonym
##
##def WriteVCFlikeoutput(masterDict,convDict,samDict):
##    for gene, SNPs in masterDict.items():
##
##        # horrible way of dealing with this :: 
##        SNPdict = getSNPfromSAMs(gene,samDict)
####        SNPdict2 = getSNPfromSAMs(gene,samDict2)
####        SNPdict3 = getSNPfromSAMs(gene,samDict3)
##        
##        # go through genes
##        CHROM = convDict[gene]
##
##        for SNP, val in SNPs.items():
##            if SNP in SNPdict:
##                val[0] =
##                val[1] =
##                val[2] =
##                val[3] =
                
            

    # Populate the master dictionary
    # This creates an empty dictionary containing all SNPs
    # masterdict here contains the a dictionary of genes. each gene is a dictionary of SNPs
    # each SNP gets a reference Seq, and alternative seq and in case of replicates a count score
    # to write this to file it needs to follow vcf format
    # compute NS (number of samples and DP combined pileup coverage
    # chrom from confDict,  pos is with the SNp, ID = gene+posingene, QUAL ??
    # e.g
##  CHROM POS    ID        REF  ALT     QUAL FILTER INFO                              FORMAT      Sample1        Sample2        Sample3
##    2      4370   rs6057    G    A       29   .      NS=2;DP=13;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:52,51 1|0:48:8:51,51 1/1:43:5:.,.
    masterDict = populateMasterDict(samDict1,samDict2,samDict3) 

    print 'masterDict with %s' %(len(masterDict))
    for element,val in masterDict.items():
        print element,val
    


    if writeout == True:
        print '[STATUS] Writing VCF file to %s' %(args.out)
        ## initiate writefunction
        writeVCFormat(masterDict,outfile,convDict)

    print '[STATUS] Creating and populating Dictionaries'

    print '[STATUS] Master Dictionary with %s possible SNPs' %(len(masterDict))

    print '[STATUS] Initiating SNP classification and dispersion estimation'

    
    silenceDict = {}

    
    for silenceGene, silenceSNP in masterDict.items():
        for silencePOS, emptydata in silenceSNP.items():
            
            silenceDict[silencePOS] = convDict[silenceGene]

    # Start classifying SNPs as ASE, and mutate fasta

    if args.spass == 'No' and args.fastaOut != 'none' :
        print '[STATUS] Writing masked fasta file '
        print '[STATUS] Please use this as a new reference to eleminate alignment derived errors'
        
        if args.fasta != 'none':
            with open('%s'%(args.fastaOut),'w') as fasta_raw:
                record_dict = SeqIO.index('%s' %(args.fasta), "fasta")
                print '[STATUS] Silencing the fasta sequence'
                MaskAFasta(silenceDict,record_dict, fasta_raw)
                print '[STATUS] Masked Fasta generated'
        if args.fasta =='none':
            print '[STATUS] No input FASTA file found, no output fasta created'
            pass


    ### INITIATE REMAPPING IF POSSIBLE
        

    writeToFile = True
    if writeToFile == True:
        print '[STATUS] writing vcf like intermediate output'

   # implement writetofile option here

    
    # lets get to the quantification

    # populate masterdict to randomized resampling ??

    if args.spass != 'No' :
     
        masterDictResample = extendMasterDict(masterDict,samDict1,samDict2,samDict3)



    # extended Masterdict for binomial testing
        resDict = {}
        total = 100
        for masDictgene, masDictSNPs in masterDictResample.items():
            for SNPpos, SNPdata in masDictSNPs.items():
                elementcount = 0
                if masDictgene not in resDict:
                    resDict[masDictgene] = {SNPpos: []}
                else:
                    resDict[masDictgene][SNPpos] = []
##            print masDictgene, SNPpos, SNPdata
##                print SNPdata

##            observations = 
            # Check if it is Noise
                elementcount = 0

            # create small dictionary, and unify the p values with Fischers 
                for element in SNPdata:
                    NoiseProb = distFunNegBin(100,int(element),0.02)


            # Check if it is a full SNP
                    negative = total - int(element)

                    FullProb = distFunNegBin(100,negative,0.02)
                        
            # Check for Allele Specific expression
                    ASEprob = pmfNegBin(total,element,0.5)
            
                    resDict[masDictgene][SNPpos].append(NoiseProb)
                    resDict[masDictgene][SNPpos].append(ASEprob)
                    resDict[masDictgene][SNPpos].append(FullProb)

                    elementcount +=1

        probDict = {}
        for gene, value in resDict.items():
            probDict[gene] = []
            for SNP, perc in value.items():
##            print gene, SNP, perc[0], perc[1],perc[3]
##            Noise = perc[0]
                probDict[gene].append(perc[1])
                ASE = perc[1]
##            Full = perc[3]


            probabil = FischerCombineP(probDict[gene])
            print gene , '  ASE probability ', probabil 

