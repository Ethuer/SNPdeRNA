import csv
from sklearn import svm
import numpy as np
import sys,argparse
import os.path
from _functions import *

###################################################################################################################################################
#
#                                           Part 2 of the ALLELE SPECIFIC EXPRESSION PIPELINE
#
#
#
#
#                this script takes the output from the SAM2SNP file and filters the noise. The output of this one has to be vcf-like
#
#
####################################################(c) ERNST THUER 2015 ##########################################################################



# 
#CHROM POS    ID        REF  ALT     QUAL FILTER INFO                              FORMAT      Sample1        Sample2        Sample3


# lets implement another classification algorithm.
# HMM for a start

parser = argparse.ArgumentParser(description='Quantify the noise from the 3 replicates ')



parser.add_argument('-rep1',
                    dest='rep1',
                    required = True,
                    help='output of SAM2SNP for replicate 1',
                    metavar = 'FILE',
                    
                    )


parser.add_argument('-rep2',
                    dest='rep2',
                    required = False,
                    help='output of SAM2SNP for replicate 2',
                    metavar = 'FILE',
                    
                    )


parser.add_argument('-rep3',
                    dest='rep3',
                    required = False,
                    help='output of SAM2SNP for replicate 3',
                    metavar = 'FILE',
                    
                    )

parser.add_argument('-out',
                    dest='out',
                    required = True,
                    help='output in vcf format',
                    metavar = 'FILE',
                    
                    )

parser.add_argument('-gff',
                    dest='gff',
                    required = True,
                    help='reference gff file',
                    metavar = 'FILE',
                    
                    )

parser.add_argument('-nu',
                    dest='nu',
                    required = False,
                    default = '0.01',
                    help='nu value for SVM fitting,  defaults to 0.01',
                    metavar = 'FLOAT',
                    
                    )

parser.add_argument('-gamma',
                    dest='gamma',
                    required = False,
                    default = '0.0',
                    help='gamma value for SVM fitting,  defaults to 0.01',
                    metavar = 'FLOAT',
                    )

parser.add_argument('-fasta',
                    dest='fasta',
                    required = False,
                    default = 'none',
                    help='Input a fasta file, that will be masked for remapping',
                    metavar = 'FILE',
                    )
parser.add_argument('-fastaOut',
                    dest='fastaOut',
                    required = False,
                    default = 'none',
                    help='Output a fasta file, that will have been masked for remapping',
                    metavar = 'FILE',
                    )

parser.add_argument('-quiet',
                    dest='quiet',
                    required = False,
                    default = False,
                    help='Change this from False, to repress Prompt on nonessential occurances',
                    metavar = 'Bool',
                    )

args = parser.parse_args()

with open('%s' %(args.rep1),'r') as in_raw,  open('%s' %(args.out),'w') as out_raw, open('%s'%(args.gff),'r') as gff_raw:
    outfile = csv.writer(out_raw, delimiter = '\t')
    gfffile = csv.reader(gff_raw, delimiter = '\t')


    # open the replicates , create dictionaries
    infile = csv.reader(in_raw, delimiter = '\t')
    inDict1 = csv2dictIncSyn(infile)

    DictList = [inDict1]

    if len(inDict1) == 0:
        print '[ERROR] Read the replicate 1 incorrectly, please verify that the format is correct, and the file exists'
        exit
        
    
    try:
        with open('%s' %(args.rep2),'r') as in_2raw:
            infile2 = csv.reader(in_2raw, delimiter = '\t')
            inDict2 = csv2dictIncSyn(infile2)
            DictList.append(inDict2)
    except:
        inDict2 = 0
        print '[OPEN] No replicate 2 ?  not optimal  Valiadation of SNP expression is going to be limited '
        pass

    try:
        with open('%s' %(args.rep3),'r') as in_3raw:
            infile3 = csv.reader(in_3raw, delimiter = '\t')
            inDict3 = csv2dictIncSyn(infile3)
            DictList.append(inDict3)

    except:
        inDict3 = 0
        print '[OPEN] No replicate 3 ?  no problem   This should not heavily decrease predictive power'
        pass
    
    print '[OPEN] All files loaded and ready'
        
    

    # build dictionaries
    
    
##    inDict3 = csv2dict(infile3)

    resDict = {}
    gffDict = {}
    geneDict = {}
    vcfDict = {}
    genecount = 0
    posDict = {}
##    noiseDict = {}
##    outDict = {}


    writeout = True
    
    neg_count = 0
    pos_count = 0
    wrong_count = 0
    origin = 'ENSEMBL'
    
        
    writevcfheader(writeout, outfile)


    gffDict = parse_gff(gfffile,origin)


    # Populate the master dictionary
    # This creates an empty dictionary containing all SNPs 
    masterDict = populateMasterDict(inDict1,inDict2,inDict3)        

    print '[STATUS] Creating and populating Dictionaries'
##    for repeat in range(0,20):
##        extendMasterDict(masterDict,inDict1,inDict2,inDict3)

    print '[STATUS] Master Dictionary with %s possible SNPs' %(len(masterDict))


    # implement an alternative method of predicting directly from count data
    # set a minimum coverage of 100 ?? to penalize low coverage areas


    print '[STATUS] Initiating SNP classification and dispersion estimation'
    for element in DictList:
        
        masterDict = likelihoodTest(element, masterDict)

##
    # now classify according to probability,  use this as an output in pass 1 ,
    #
    count = 0
    count_all = 0

    for element, value in masterDict.items():
        
        if value[1] > 0.05:
            value.append('acceted')
            count +=1
            posDict[element] = value[0]
            # outfile writer if firstpass,   not if secondpass
        else :
            value.append( 'rejected' )
            count_all +=1
##        print element, value
    print count, count_all


    writtenlist = []
    writecount = 0
    for element, value in posDict.items():
        # write the SNPs to file
        if element in inDict1:
            writtenlist.append(element)
##            element[0] == gene,   element[1] == position
            outfile.writerow([element[0],element[1], inDict1[element][0],inDict1[element][1],inDict1[element][2],inDict1[element][3],inDict1[element][4],inDict1[element][5]])
            resDict[element] = [inDict1[element][0],inDict1[element][1]]
            writecount +=1               
        if inDict2 != 0 and element in inDict2 and element not in writtenlist:
            outfile.writerow([element[0],element[1], inDict2[element][0],inDict2[element][1],inDict2[element][2],inDict2[element][3],inDict2[element][4],inDict2[element][5]])
            writtenlist.append(element)
            resDict[element] = [inDict2[element][0],inDict2[element][1]]
            writecount +=1
        if inDict3 != 0 and element in inDict3 and element not in writtenlist:
            outfile.writerow([element[0],element[1], inDict3[element][0],inDict3[element][1],inDict3[element][2],inDict3[element][3],inDict3[element][4],inDict3[element][5]])
            writtenlist.append(element)
            resDict[element] = [inDict3[element][0],inDict3[element][1]]
            writecount+=1
                



    print '[STATUS] %s  SNPs written to file' %(writecount)

    # now mask a fasta with this information
    try:
        if args.fasta != 'none':
            with open('%s'%(args.gff),'r') as fasta_raw:
                record_dict = SeqIO.index('%s' %(args.fasta), "fasta")
                print '[STATUS] Silencing the fasta sequence'
                MaskAFasta(resDict,record_dict)
                print '[STATUS] Masked Fasta generated'
        if args.fasta =='none':
            print '[STATUS] No masked Fasta generated'
    except:
        print '[STATUS] No masked Fasta generated'
        
##

    




        
##    for key, value in masterDict.items():
####        print key, value
##        total = len(value)
##        observations = 0
##        probability = 0.2
##
##        for element in value:
##            if int(element) > 0:
##                observations += 1
####        print total, observations
##        likelyhood = pmfNegBin(total,observations,probability)
##
##
##        count_all += 1
##        if float(likelyhood) < 0.05:
##            
##            print key, likelyhood
##            count  +=1
##
##    print count, count_all
##        # test against binomial by reducing counts to true and false
##
##
##        

##
##
##
##    # Find false positives
##
##    # a nonparametric approach, if there are replicates, missing overlap with low expression
##    # a parametric approach, taking the lowest 5 10 15 % as teaching data,
##    negDict = {}
##    posDict = {}
##    if inDict2 == 0:
##        print '[STATUS] No replicates,  using parametric approach to evaluate data'
##        for element, values in inDict1.items():
##            # using the uncorrected probability value from the binomial distribution fitting
##            if float(values[4]) > 0.1 and 'NonSynon' in values[5]:
##                negDict[element] = masterDict[element]
##            else:
##                posDict[element] = masterDict[element]
##            
##
##    else:
##        print '[STATUS] Scanning data for False positives'
##        # check for lacking overlaps, < 34% replicates, immediately considered False
##        # > 34  (2 out of 3,  if probability > 0.1
##        for element,value in masterDict.items():
##
##
##            try:
##                
##                if inDict3 != 0 and element not in inDict2 and element not in inDict3:
##                    negDict[element] = value
##                    
##                elif element not in inDict1 and element not in inDict2:
##                    if element not in negDict:
##                        negDict[element] = value
##
##                # extend the teaching to nonsynonymous SNPS in only 2 out of 3 replicates
##
##                if element not in inDict1 and 'NonSynon' in value[5]:
##                    if element not in negDict:
##                        negDict[element] = value
##
##                if element not in inDict2 and 'NonSynon' in value[5]:
##                     if element not in negDict:
##                        negDict[element] = value
##                   
##                    
##            except:
##                pass
##        for element, value in masterDict.items():
##            posDict[element] = value
##
##
##    # convert noiseDict to np.array
##    # normalize to range [0,1],   easy since its subset binned count data,   divide by 50(maxcount)
##
##
##
##
##    # test the binomial pmf
##
##    
##    
##            
##    noiseArray = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
##    for element, values in negDict.items():
##        noiseArray = np.vstack([noiseArray,values])
##
##    # keep gamma in auto
##    noiseModel = svm.OneClassSVM(nu=0.1, kernel="poly",verbose=False,cache_size=2000.0)
##    
##    noiseModel.fit(noiseArray)
##    writtenlist = []
##    writecount = 0
##    outDict = {}
##    resDict = {}
##    for element, value in posDict.items():
##        prediction = noiseModel.predict(value)
##        
##        # write the SNPs to file
##        if int(prediction) > 0:
##            
##            if element in inDict1:
##                writtenlist.append(element)
##                outfile.writerow([element[0],element[1], inDict1[element][0],inDict1[element][1],inDict1[element][2],inDict1[element][3],inDict1[element][4],inDict1[element][5]])
##                resDict[element] = [inDict1[element][0],inDict1[element][1]]
##                writecount +=1               
##            if inDict2 != 0 and element in inDict2 and element not in writtenlist:
##                outfile.writerow([element[0],element[1], inDict2[element][0],inDict2[element][1],inDict2[element][2],inDict2[element][3],inDict2[element][4],inDict2[element][5]])
##                writtenlist.append(element)
##                resDict[element] = [inDict2[element][0],inDict2[element][1]]
##                writecount +=1
##            if inDict3 != 0 and element in inDict3 and element not in writtenlist:
##                outfile.writerow([element[0],element[1], inDict3[element][0],inDict3[element][1],inDict3[element][2],inDict3[element][3],inDict3[element][4],inDict3[element][5]])
##                writtenlist.append(element)
##                resDict[element] = [inDict3[element][0],inDict3[element][1]]
##                writecount+=1
##                
##
##
##
##    print '[STATUS] %s  SNPs written to file' %(writecount)
##
##    # now mask a fasta with this information
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
##
##
##    # re-map
##    # then recount the real SNPs with SAM2SNP, get the average expression over all 3 replicates. This may need a different script for the moment
##
##
##
##
##    # then give an aberration to a mean as Allele specific expression
##
