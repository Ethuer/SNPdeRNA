import csv
from _functions import csv2dict
from sklearn import svm
import numpy as np
import sys,argparse
import os.path
from _functions import writevcfheader
# this script takes the output from the BAM2SNP file and filters the noise. THe output of this one has to be vcf-like
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
                    required = True,
                    help='output of SAM2SNP for replicate 2',
                    metavar = 'FILE',
                    
                    )


parser.add_argument('-rep3',
                    dest='rep3',
                    required = True,
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



args = parser.parse_args()

with open('%s' %(args.rep1),'r') as in_raw, open('%s' %(args.rep2),'r') as in_2raw, open('%s' %(args.rep3),'r') as in_3raw, open('%s' %(args.out),'w') as out_raw, open('%s'%(args.gff),'r') as gff_raw:
    infile = csv.reader(in_raw, delimiter = '\t')
    infile2 = csv.reader(in_2raw, delimiter = '\t')
    infile3 = csv.reader(in_3raw, delimiter = '\t')
    outfile = csv.writer(out_raw, delimiter = '\t')
    gfffile = csv.reader(gff_raw, delimiter = '\t')

    # build dictionaries
    inDict = csv2dict(infile)
    inDict2 = csv2dict(infile2)
    inDict3 = csv2dict(infile3)

    resDict = {}
    gffDict = {}
    geneDict = {}
    vcfDict = {}
    genecount = 0
    noiseDict = {}
    outDict = {}


    writeout = True
    
    neg_count = 0
    pos_count = 0
    wrong_count = 0

    
##    for row in gfffile:
##        try:
##            gene = row[8].split('"')[1]
##            gffDict[gene]=row[0]
##        except:
##            pass
        
    writevcfheader(writeout, outfile)
    
    for row in gfffile:
##        print row
        # skip header
        if not '#' in row[0]:
            try:
                if row[2] == args.feature:
                    if not origin =='ENSEMBL':
                # legacy,  this has to be more flexible
                        gffDict[row[0],row[8].split('"')[1]] = [row[3],row[4],row[6],row[2]]
                    if origin == 'ENSEMBL':
                        name = re.split('=|;',row[8])[2]
                        gffDict[name] = [row[0]]
            except:
##                print 'wrong gff configuration in %s' %(row[0])
                pass
    print "gff file with %s features" %(len(gffDict))

    
##    print len(inDict), len(inDict2),len(inDict3)
##    make dictionary from csv reader objects

    # write the binning function

    for f in [1..50]:
        
    def binning(inDict):
        for element, value in inDict.items():
            
            
            
            
    

    
    count = 0    
    for element, value in inDict.items():
        
##        print element,
        if element in inDict2 and element in inDict3:
            Nucl_orig = value[1]
            pos_count += 1

            # nucleotide
            if value[1] == inDict2[element][1] and value[1] == inDict3[element][1]:
                Nucl_alt = value[1]
            elif value[1] != inDict2[element][1] or value[1] != inDict3[element][1]:
                Nucl_alt = 'N' 

            prob = float(inDict[element][4])
            prob2 = float(inDict2[element][4])
            prob3 = float(inDict3[element][4])

            snp1 = float(inDict[element][2])
            snp2 = float(inDict2[element][2])
            snp3 = float(inDict3[element][2])

            count1 = float(inDict[element][3])
            count2 = float(inDict2[element][3])
            count3 = float(inDict3[element][3])
            

            if Nucl_alt == 'N' or prob > 0.05 or prob2 > 0.05 or prob3 > 0.05 :
                
                neg_count +=1
                noiseDict[element] = [((snp1/count1)*count1),((snp2/count2)*count2),((snp3/count3)*count3)]

                
            else:
                resDict[element] = [((snp1/count1)*count1),((snp2/count2)*count2),((snp3/count3)*count3)]
                count +=1


    pos_count = 0
    neg_count = 0
    # IMPLEMENT svm   OneClassSVM

    # noiseDict contains teaching material,  resDict the testing material
    # convert noiseDict to np.array
    noiseArray = [0,0,0]
    
    for element, value in noiseDict.items():
        noiseArray= np.vstack([noiseArray,value])

    # Predict the SVM for the noise data
    print 'Predicting noise model'
##    noiseModel = svm.OneClassSVM(nu=0.01, kernel="rbf", gamma=0.0,verbose=True)
    noiseModel = svm.OneClassSVM(nu=(args.nu), kernel="rbf", gamma=(args.gamma),verbose=True)
    noiseModel.fit(noiseArray)

##    noisecount 

    ## ID =gene_name+counter
    out_count = 0
    for element, value in resDict.items():
        prediction = noiseModel.predict(value)

        
        if  int(prediction) == 1:
            neg_count +=1
        elif  int(prediction) == -1:
            pos_count +=1
            if element[0] in outDict:
                outDict[element[0]] += 1
            if not element[0] in outDict:
                outDict[element[0]] = 1


            #ID='%s_%s' %(element[0],outDict[element[0]])
	    ID = element[0]
            AF = value[0]
            DP = inDict[element][4]
            
            
            REF = inDict[element][0]
            ALT = inDict[element][1]
            
            QUAL = '30'
            FILTER = 'PASS'
            # TO DO replace number of samples with actual number'
            # TO DO replace depth with actual depth' 
            INFO = 'NS=%i;DP=%i;AF=%s;' %(3,30,AF,)

            FORMAT = '-'

##            Sample1 ='-'
            ratio1 = (float(inDict[element][2]) / float(inDict[element][3]))
##            if ratio1 > 0.8:
##                print ratio1, inDict[element][2] , float(inDict[element][3])

            
            ratio2 = (float(inDict2[element][2]) / float(inDict2[element][3]))
            ratio3 = (float(inDict3[element][2]) / float(inDict3[element][3]))
            
            Sample1 = inDict[element][2]
            Sample2 ='-'
            Sample2 = inDict2[element][2]
            Sample3 = '-'
            Sample3 = inDict3[element][2]
##            value[0],value[1],value[2]
            
            outfile.writerow([gffDict[element[0]],element[1], ID, REF, ALT, QUAL, FILTER, FORMAT,Sample1,ratio1,Sample2,ratio2,Sample3,ratio3 ])
            out_count += 1
    print '%s SNPS rejected, %s SNPS accepted ' %(neg_count, pos_count)
    print '%s SNPS written to file' %(out_count), 
##    print neg_count, pos_count, wrong_count
