import csv
##from sklearn import svm
##import numpy as np
import sys,argparse
import os.path
from _functions import *
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq
import re
import csv


def SAM2SNP(args.feature,fasta_raw,samfile,gffDict,outfile, createIntermediate = False , cutoff=0.05 , verbose = True):
        """

    wrapper function for coverage analysis of SAM files,
    returns Dictionary 'transdict' containing SNPs, their coverage
    ... and a precomputed likelihood
        """
        resultDict ={}
        SNPDict = {}
        transDict = {}
        fastadict = {}
        vcfSubDict = {}

        perc_list = []
        
        basecount = 0
        genecount = 0
        writecount = 0
        
        

        
        print "[STATUS] starting Gene-wise coverage analysis .. This may take a while"
        print "[STATUS] 0 percent of %ss analyzed " %(args.feature)

        for element in SeqIO.parse(fasta_raw, "fasta"):
        
            for key,value in gffDict.items():
                fastadict = {}
                count  = 0
            # if on chromosome
                if key[0] == element.id:
                    start = int(value[0])
                    stop = int(value[1])
                
                    gene = key[1]
                
                    fastadict = list2dict(element.seq[start:stop],start)
                             
                    genecount +=1

                
                # print gene count in percent of features  every 10 percent
                    getPercentage(genecount,len(gffDict),args.feature,perc_list)


                # flexpile function calls an additional probability function,  maybe redundant.
                    resultDict[gene]= flexpile(fastadict,samfile.pileup("%s" %(element.id),start,stop),element.id, start, stop, args.spass)
                
                # create the syn / nonsyn differentiation directly in the classification loop should be fastest

                # translate fasta of the gene to protein   POSSIBLE BUG FOOD    check for non-generic DNA motifs...
                    DNA = Seq(str(element.seq[start:stop]),generic_dna)
                    protein = 'NNN'

                
                    try:

                    # explicitly trim the sequence to contain only full codons
                        if len(DNA)%3 != 0:
                            overlap = len(DNA)%3
                            DNA = DNA[:-int(overlap)]

                        
                        protein = DNA.translate()
                        protein = list2dict(protein[0:len(protein)],0)
                    except:
                        print 'Error in translating Protein %s' %(element.id)


                
                    positionList = []

                # classify as SNPs
                    for sub_element, value in resultDict[gene].items():
                    
                        converted = classifydict(value)
                        for nucleotide, valuen in converted.items():
                            if int(valuen[0]) > 0 and float(valuen[2]) < float(args.minqual):
                                count +=1


                            # now check for synonymity,  move this to functions
                                alternativeSeq = MutableSeq(str(element.seq[start:stop]), generic_dna)
                                alternativeSeq = mutateSequence(alternativeSeq,sub_element,nucleotide,start)
                                alternativeSeq = Seq(str(alternativeSeq), generic_dna)
                                if len(alternativeSeq)%3 != 0:
                                    overlap = len(alternativeSeq)%3
                                    alternativeSeq = alternativeSeq[:-int(overlap)]
                            
                                altprot = alternativeSeq.translate()
                                altprot = list2dict(altprot[0:len(protein)],0)

                                protposition = int((sub_element-start)/3)
                                try:
                                    if protein[protposition] != altprot[protposition]:
                                        positionList.append(sub_element)
                                        synonym = 'NonSynon'
                                    if protein[protposition] == altprot[protposition]:
                                        synonym = 'Syn'
                                except:
                                    synonym = 'Unknown'

                            # correct against python specific error of starting count at 0. 

                                # Python starts counting at 0, 
                                truepos = int(sub_element)+1


                                # implement negative binomial classification for noisecorrection here :


                                observations = valuen[1]
                                total = valuen[2]

                                # miss chance is the illumina specific 1.4%    make that more flexible later
                                likelyhoodNoise = distFunNegBin(total, observations,0.014)



##                          result of a density function,   not a p value !
                                if float(likelyhoodNoise) < float(cutoff):
                                    transDict[gene,truepos] = [value[0], nucleotide, valuen[0],valuen[1],likelyhoodNoise,synonym]
                    
        # Output writing for intermediate results

        if createIntermediate == True:
            outfile.writerow('Gene','Pos','ORIG','ALT','SNPcoverage','TotalCoverage','likelyhood NegBin','Synonymity',)
            for key, value in transDict.items():
                outfile.writerow(key,value[0],value[1],value[2],value[3],value[4],value[5]],value[6])
                writecount += 1
            


        if createIntermediate == True and int(writecount) > 0:
            print '[RUN SUCCESSFUL] %s Preliminary SNPs written to file %s' %(writecount,args.out)
        elif createIntermediate == True and int(writecount) == 0:
            print '[RUN FAILED] No Preliminary SNPs written to file' %(writecount)


        print '[STATUS] %s SNPs found in samfile %s' %(len(transDict), samfile)

        return transDict
