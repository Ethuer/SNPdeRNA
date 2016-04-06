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


def SAM2SNP(feature,fasta_raw,samfile,gffDict,outfile,cutoff ,spass,errortolerance = 0.12, verbose = True, createIntermediate = False):
        """

    wrapper function for coverage analysis of SAM files,
    returns Dictionary containing SNPs, their coverage
    ... and a precomputed likelihood
    
    This is slow, but steady with the parametric error rate :
    
    Work-flow description:
    mpileup is used to obtain coverage of nucleotides that differ to the reference
    
    noise classification is computed against a negative binomial distribution.
    
    computes synonymity of the change in base composition by extracting local codon.
    
    
        """
        
        resultDict ={}
        SNPDict = {}
        transDict = {}
        fastadict = {}
        vcfSubDict = {}

        print len(fasta_raw)

        perc_list = []
        
        basecount = 0
        genecount = 0
        writecount = 0
        
        spass = 'No'

        if verbose == True:
            print "[STATUS] starting Gene-wise coverage analysis .. This may take a while"
            print "[STATUS] 0 percent of %ss analyzed " %(feature)
        

        extendedResultsDict = {}
        
        for element in fasta_raw:
            
            for key,value in gffDict.items():
                print key
                fastadict = {}
                count  = 0
            # if on chromosome
                if key[0] == element.id:
                    start = int(value[0])
                    stop = int(value[1])
                    gene = key[1]
##                    print gene
                    fastadict = list2dict(element.seq[start:stop],start)
                             
                    genecount +=1
                    
                    # load DNA as sequence into DNA variable, for synonymity checking, later
                    DNA = Seq(str(element.seq[start:stop]),generic_dna)
                    alternativeSeq = MutableSeq(str(element.seq[start:stop]), generic_dna)
                        
                    protein = 'NNN'
                    try:

                    # explicitly trim the sequence to contain only full codons
                        if len(DNA)%3 != 0:
                                overlap = len(DNA)%3
                                DNA = DNA[:-int(overlap)]

                        
                        protein = DNA.translate()
                        protein = list2dict(protein[0:len(protein)],0)
                        DNADict = list2dict(DNA[0:len(DNA)],0)
                    except:
                        print 'Error in translating Protein %s' %(element.id)

                    
                
                # print gene count in percent of features  every 10 percent
                    getPercentage(genecount,len(gffDict),feature,perc_list)


                # flexpile function calls an additional probability function,  maybe redundant.
                    SNPDict[gene]= flexpile(fastadict,samfile.pileup("%s" %(element.id),start,stop),element.id, start, stop, spass)

##                    for element

                    # implement significance test here, to keep the dictionaries sizable
                    for SNPindividual, SNPindivdata in SNPDict[gene].items():
                                    
##                        print  SNPindividual, SNPindivdata
                        # this fun returns a dictionary of the possible SNPs with given likelihood and coverage
                        # give a hightened minimum threshold for noise clearance
                        SNPdata = classifydict(SNPindivdata,errortolerance, 0.05,100)
                         
                        for possibleSNP, possibleContent in SNPdata.items():
                                position = SNPindividual
                                # python starts count at 0, not 1,
                                
                                ORIG = SNPindivdata[0]
                                ALT = possibleSNP
                                SNPcov = possibleContent[0]
                                Basecov = possibleContent[1]
                                DensfunProb = possibleContent[2]
                                # check position within gene,   position - start
                                synonym = 'Not analyzed'
                                synonym = getProt( position, protein,DNADict,ALT,int(value[0]))
##                                synonym = synonymity(protein,alternativeSeq, position, ALT,int(value[0]) )
                                position = int(position) + 1


                                # this has to change, the unique feature is the position / gene combination so this has to be a dict of dicts...

                                
##                                identity = [gene,position]

                                if not gene in transDict:
                                        transDict[gene] = {position : [ORIG,ALT,SNPcov,Basecov,DensfunProb,synonym]}
                                if gene in transDict:
                                        transDict[gene][position] = [ORIG,ALT,SNPcov,Basecov,DensfunProb,synonym]
##                                        position = '%s_alt' %(position)
##                                        identity = str([gene,position])
##                                        transDict[identity] = [ORIG,ALT,SNPcov,Basecov,DensfunProb,synonym]

        return transDict

def findCodon(NuclPos, ALT, Modulo, sequenceDict):
        """
        Could be simplified, but that is quite self explanatory.
        recreate the new codon, translate.

        Modulo 0,   position divisible by 3, so SNP is in the 0 position. Memento Python .
        Modulo 1,   position has 1 residue, so SNP is at pos 1 == 2nd posiion on codon.
        """
        codon = ['N','N','N']
        if Modulo == 0:
                codon[0] = ALT
                codon[1] = sequenceDict[(int(NuclPos)+1)]
                codon[2] = sequenceDict[(int(NuclPos)+2)]

        if Modulo == 1:
                codon[0] = sequenceDict[(int(NuclPos)-1)]
                codon[1] = ALT
                codon[2] = sequenceDict[(int(NuclPos)+1)]

        if Modulo == 2:
                codon[0] = sequenceDict[(int(NuclPos)-2)]
                codon[1] = sequenceDict[(int(NuclPos)-1)]
                codon[2] = ALT             

        return codon
        


def getProt(position, origDNADict,seqDict, nucleotide, start ):
        nucposition = int(position-start)
        protposition = ((nucposition)/3)

        Mod = nucposition%3
        Codon = findCodon(nucposition, nucleotide, Mod, seqDict)
        Codon = ''.join(Codon)
        Codon = Seq(str(Codon), generic_dna)

        if str(origDNADict[protposition]) == str(Codon.translate()):
                SYN = 'Syn'
        else:
                SYN = 'Nonsyn'
        



        

def synonymity(DNA,ALTSEQ, position,nucleotide,startofprotein):
        """
        calculate synonymity for SNPs,
        take a DNA sequence and the SNP
        mutate DNA, translate both to Proteins, and check for differences
        """
        #### improve this by making it localized,  just get the individual SNP position,  modulo to see where in the protein it is, then check against dict

                        
##                        print  SNPindividual, SNPindivdata

        # the vcf position is absolute,  so start of protein needs to be given
##        print type(alternativeSeq)
        alternativeSeq = mutateSequence(ALTSEQ,position,nucleotide,startofprotein)
        alternativeSeq = Seq(str(alternativeSeq),generic_dna)
        if len(alternativeSeq)%3 != 0:
                overlap = len(alternativeSeq)%3
                alternativeSeq = alternativeSeq[:-int(overlap)]
         
        altprotein = alternativeSeq.translate()
        altprotein = list2dict(altprotein[0:len(altprotein)],0)
        
##              now compare the two proteins               
        protposition = int((position-startofprotein)/3)

        
        try:
                if protein[protposition] != altprot[protposition]:
                        positionList.append(sub_element)
                        synonym = 'NonSynon'
                if protein[protposition] == altprot[protposition]:
                        synonym = 'Syn'
        except:
                        synonym = 'Unknown'

        return synonym
                            
##        print len(protein), len(altprotein)
##        for key, value in resultDict.items():
##                for item, content in value.items():
##                        print 'key', key, 'value', value, 'items',item,'content', content
####                extendedResultsDict[key] = SAM2SNPextension(value, createIntermediate = False)
##        return resultDict
##

##def SAM2SNPextension(resultDict, createIntermediate = False):
##                    
##        for key, value in resultDict.items():
####                            print key, value
##                key = value
##                    
##                # create the syn / nonsyn differentiation directly in the classification loop should be fastest
##
##                # translate fasta of the gene to protein   POSSIBLE BUG FOOD    check for non-generic DNA motifs...
##                    DNA = Seq(str(element.seq[start:stop]),generic_dna)
##                    protein = 'NNN'
##
##                
##                    try:
##
##                    # explicitly trim the sequence to contain only full codons
##                        if len(DNA)%3 != 0:
##                            overlap = len(DNA)%3
##                            DNA = DNA[:-int(overlap)]
##
##                        
##                        protein = DNA.translate()
##                        protein = list2dict(protein[0:len(protein)],0)
##                    except:
##                        print 'Error in translating Protein %s' %(element.id)
##
##
##                
##                    positionList = []
##
##                # classify as SNPs
##                    for sub_element, value in resultDict[gene].items():
##                        print value
##                        converted = classifydict(value)
##                        for nucleotide, valuen in converted.items():
####                            if int(valuen[0]) > 0 and float(valuen[2]) < float(cutoff):
##                                count +=1
##
##
##                            # now check for synonymity,  move this to functions
##                                alternativeSeq = MutableSeq(str(element.seq[start:stop]), generic_dna)
##                                alternativeSeq = mutateSequence(alternativeSeq,sub_element,nucleotide,start)
##                                alternativeSeq = Seq(str(alternativeSeq), generic_dna)
##                                if len(alternativeSeq)%3 != 0:
##                                    overlap = len(alternativeSeq)%3
##                                    alternativeSeq = alternativeSeq[:-int(overlap)]
##                            
##                                altprot = alternativeSeq.translate()
##                                altprot = list2dict(altprot[0:len(protein)],0)
##
##                                protposition = int((sub_element-start)/3)
##                                try:
##                                    if protein[protposition] != altprot[protposition]:
##                                        positionList.append(sub_element)
##                                        synonym = 'NonSynon'
##                                    if protein[protposition] == altprot[protposition]:
##                                        synonym = 'Syn'
##                                except:
##                                    synonym = 'Unknown'
##
##                            # correct against python specific error of starting count at 0. 
##
##                                # Python starts counting at 0, 
##                                truepos = int(sub_element)+1
##                                # implement negative binomial classification for noisecorrection here :
##                                observations = valuen[1]
##                                total = valuen[2]
##                                # miss chance is the illumina specific 1.4%    make that more flexible later
##                                likelyhoodNoise = distFunNegBin(total, observations,0.014)
##                                print likelyhoodNoise
####                          result of a density function,   not a p value !
##                                if float(likelyhoodNoise) < float(cutoff):
##                                        print cutoff
##                                        transDict[gene,truepos] = [value[0], nucleotide, valuen[0],valuen[1],likelyhoodNoise,synonym]
##                                        
##                                    
##                    
##        # Output writing for intermediate results
##
##        if createIntermediate == True:
##            outfile.writerow('Gene','Pos','ORIG','ALT','SNPcoverage','TotalCoverage','likelyhood NegBin','Synonymity',)
##            for key, value in transDict.items():
##                outfile.writerow(key,value[0],value[1],value[2],value[3],value[4],value[5],value[6])
##                writecount += 1
##            
##
##
##        if createIntermediate == True and int(writecount) > 0:
##            print '[RUN SUCCESSFUL] %s Preliminary SNPs written to file' %(writecount)
##        elif createIntermediate == True and int(writecount) == 0:
##            print '[RUN FAILED] No Preliminary SNPs written to file' %(writecount)
##
##
##        print '[STATUS] %s SNPs found in samfile %s' %(len(transDict), samfile)
##
##        return transDict
