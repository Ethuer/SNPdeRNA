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
import random

class Perceptron():
    """
    using class for perceptrons to streamline Neural network
    TemplateClass
    """
    
    
    def __init__(self):
        Theta = 10
        WeightList = []
        ObservList = []
        Classifier = 0
    
    def CalculateDotProduct(self,ObservList, WeightList):
        """
        Calculate the dot product for 2 equal length input vectors
    
        input numerical vectors
        len( Vector 1 ) = len(Vector2)
    
        Vector1 derives from sample observations, Vector 2 is a weight vector
    
        output dot prodcut
    
        """
    
        sumVector = []
        count = 0
        if len(Vector1) == len(Vector2):
        
            for element in Vector1:
                intermediate = float(element)*float(Vector2[count])
                sumVector.append(intermediate)
            
        
                count +=1
    
                dotprod = sum(sumVector)
    
        return dotprod 
    
    
        
    
    


parser = argparse.ArgumentParser(description='Looking for SNPs in RNAseq data. This part accepts SAM or BAM input files, mapped against an available reference. In a first step it counts the occurrences of Polymorphisms against a reference. In second pass mode, only already detected SNPs are being analyzed')



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

parser.add_argument('-bam',
                    dest='bam',
                    nargs='+',
                    required = True,
                    default = 0,
                    help='Input the bam file(s) of the Alignment.  Tophat2 output works,  NextGenMapper or STAR.  Bam files need to be indexed using "samtools index" ',
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
    #print '[IMPORT] Loading BAM File %s' %(args.bam1)
    
    print '[IMPORT] Loading Fasta file %s' %fasta_raw
    print '[IMPORT] Loading GFF/GTF File %s ' %(gff_raw)
    print '[IMPORT] Loading Complete.'
    
    gff_file = csv.reader(gff_raw, delimiter = '\t')
    outfile = csv.writer(out_raw, delimiter = '\t', quoting=csv.QUOTE_NONE, quotechar='')
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
    
    # list of dictionaries containing the pre-screened snips 
    dictList = []
    # one BAM file is essential
    ## CATCH EXCEPTION for missing indexing file

    # parse the fasta file before the Sam invocation 
    fastaDict = []
    for record in SeqIO.parse(fasta_raw, "fasta"):
        fastaDict.append(record)
    


    verbose = True
    errortolerance = 0.12
    createIntermediate = False


    bamCount = 0

    for bamFile in args.bam:     
        SAMList.append(bamFile)
        
        samDict,SAMList = SamImport(bamFile,SAMList,args.feature,fastaDict,gffDict,outfile,cutoff,args.spass,errortolerance )
        dictList.append(samDict)
        samDict = {}
        bamCount +=1
    
## Import functions for the indicidual bam files,  
## this should be made more flexible for multiple bam imports without necessitating enumeration by the user
##    try:
        
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
    # masterDict = populateMasterDict(samDict1,samDict2,samDict3) 

        
    # additional reference Dictionary, refDict is needed, for a universal classification approach.
    # this refDict will be repopulated later, only contains the pointers to the SNPs
    
    
    masterDict, refDict = populateMasterDictfromList(dictList)

    
#     for key, value in masterDict.items():
#         print key, value
        


#     if writeout == True:
#         print '[STATUS] Writing VCF file to %s' %(args.out)
#         ## initiate writefunction
#         writeVCFormat(masterDict,outfile,convDict)

    print '[STATUS] Creating and populating Dictionaries'

    print '[STATUS] Master Dictionary with %s possible SNPs' %(len(masterDict))

    print '[STATUS] Initiating SNP classification and dispersion estimation'

    
    silenceDict = {}

    
    for silenceGene, silenceSNP in masterDict.items():
        for silencePOS, emptydata in silenceSNP.items():
            
            silenceDict[silencePOS] = convDict[silenceGene]

    # Start classifying SNPs as ASE, and mutate fasta
    
    args.spass = 'No'
    

    if args.spass == 'No' and args.fastaOut != 'none' :
        print '[STATUS] Writing masked FASTA reference file '
        print '[STATUS] Please use this as a new reference, to eliminate alignment derived errors'
        
        if args.fasta != 'none':
            with open('%s'%(args.fastaOut),'w') as fasta_raw:
                record_dict = SeqIO.index('%s' %(args.fasta), "fasta")
                print '[STATUS] Silencing the fasta sequence'
                print '[STATUS] Masked Fasta generated'
        if args.fasta =='none':
            print '[STATUS] No input FASTA file found, no output fasta created'
            pass


    ### INITIATE REMAPPING IF POSSIBLE
        

    # option to write output so far to file, for a vcf like intermediate output 
    # This should be discouraged

    writeToFile = True
    if writeToFile == True:
        print '[STATUS] writing vcf-like intermediate output'

   
    
    # lets get to the quantification

    # populate masterdict to randomized resampling ??
    # deconstruct the master dictionary into an empty template Dict containing only SNP positions 
    
    # create template Master Dictionary
    # current masterDict is a dict of dicts, needs to unravel
   
    # need to calculate the weight matrix for SNPs according to their occurance
    
    
    # I need a weight dictionary,   derived from occurance across replicates
    # the interesting part for the weightdict is the occurance per sample
    
    # K_value calculator 
    # needs to be done in ways,  ifno replicates are given, noise is estimated from the behaviour of True negatives, that don't meet the min threshold.
    # if replicates are available, the error is derived from the variances of the distribution defined by the observations
    
        # find Variance in inDict observations
    
    
        # find Variance in minThreshold observations
        # alternativ model variance on parametric negBinDist 
        
    
    
    repeats = 10
    
    weightDictionary = PopulateWeightdict(masterDict,repeats)
#             
#     for key, value in weightDictionary.items():
#         
#         for val, content in value.items():
#             if key in masterDict:
#                 print key, val, content
#                 print masterDict[key][val]
#             
#     
    
    
    # theta is the issue  
    
    # d should be a minimum cutoff derived from the standard deviation of the data according to negbin distribution
    # second parameter is the calculation mode
    d_mode = 1
    d_val = CalculateDValue(dictList, d_mode)
    #k_val = 0.02
    
    # k_value is a preset,  modifies cost of false positive by changing the cutoff
    # represents the slope, 
    k_val = 10
    
    
    # xval is the minimum expected likelyhood, this is multiplied by the cost for false positive
    x_val = 1
    
    
    Theta = CalculateTheta(d_val,k_val,x_val)
    
    
    # implement linear function for computing theta from the expected error...
    # y = kx + d    
    # i need a robust way of measuring the error in data
    # d can be the expected error k should be

    

    
    
    if args.spass == 'No' :
       
        masterDictResample = extendMasterDictbyResampling(refDict,dictList, repeats)    
        
        # for the first perceptron Expectation is 0.1 
        masterDictPerc02 = SampleToLikelyhood(masterDictResample, 0.2 )
        masterDictPerc1 = SampleToRatio(masterDictResample)
        
        
        # calculate more useful information than just count,  lets add up the probabiliy function for > 1 %
        
        
        
        
        for key, value in masterDictPerc1.items():
            for POS, val in value.items():
                #print key, POS, val
                echo = SLPClassification(val, weightDictionary[key][POS], Theta)
                print echo, Theta
                
                
                
                
                
                #print key,POS, val
                #print key,POS, weightDictionary[key][POS]
                
                #echo = CalculateDotProduct(val, weightDictionary[key][POS])
                
                
                #print echo
        # run the individual observations into the first perceptron, 
        # 
        
        
#         for element, value in masterDictResample.items():
#             print element, value
#             for val, bins in value.items():
#                 for bin in bins:
#                     if int(bin) < 40:
#                         print bin
        
        #  check normalisation procedures   important for noise classification,
        # using random resampling might be good remove relative error from normalisation, and stabilize lack of replicates
        # implement SLP here for classification
        
        # first catch hypothesis,  expected distribution, likelyhood of being part of that distribution and where to update Theta
        
# 
# 
# 
#     # extended Masterdict for binomial testing
#         resDict = {}
#         total = 100
#         for masDictgene, masDictSNPs in masterDictResample.items():
#             for SNPpos, SNPdata in masDictSNPs.items():
#                 elementcount = 0
#                 if masDictgene not in resDict:
#                     resDict[masDictgene] = {SNPpos: []}
#                 else:
#                     resDict[masDictgene][SNPpos] = []
# 
#             # Check if it is Noise
#                 elementcount = 0
# 
#             # create small dictionary, and unify the p values with Fischers 
#             # Fischers test is insufficient for the question at hand
#             # Try implementing a Single Layer Perceptron
#                 for element in SNPdata:
#                     NoiseProb = distFunNegBin(100,int(element),0.02)
# 
# 
#             # Check if it is a full SNP
#                     negative = total - int(element)
# 
#                     FullProb = distFunNegBin(100,negative,0.02)
#                         
#             # Check for Allele Specific expression
#                     ASEprob = pmfNegBin(total,element,0.5)
#             
#                     resDict[masDictgene][SNPpos].append(NoiseProb)
#                     resDict[masDictgene][SNPpos].append(ASEprob)
#                     resDict[masDictgene][SNPpos].append(FullProb)
# 
#                     elementcount +=1
# 
#         probDict = {}
#         for gene, value in resDict.items():
#             probDict[gene] = []
#             for SNP, perc in value.items():
# ##            print gene, SNP, perc[0], perc[1],perc[3]
# ##            Noise = perc[0]
#                 probDict[gene].append(perc[1])
#                 ASE = perc[1]
# ##            Full = perc[3]
# 
# 
#             #probabil = FischerCombineP(probDict[gene])
#             print gene , '  ASE probability ', probabil 

