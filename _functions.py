import pysam
from Bio import SeqIO
from scipy import stats
import operator
# Funtions for the RNAseq SNP expression 


def writevcfheader(switch, outfile):
    """
    simple header in the same format as GATKs vcf files.
    """
    if switch == True
        outfile.writerow(["##fileformat=VCFv4.0"])
        outfile.writerow(["##fileDate="])
        outfile.writerow(["##reference="])
##        outfile.writerow(["##phasing=partial"])
##        outfile.writerow(["##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"])
##        outfile.writerow(["##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"])
##        outfile.writerow(["##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">"])
##        outfile.writerow(["##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"])
##        outfile.writerow(["##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"])
##        outfile.writerow(["##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"])
##        outfile.writerow(["##FILTER=<ID=q10,Description=\"Quality below 10\">"])
##        outfile.writerow(["##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">"])
##        outfile.writerow(["##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"])
##        outfile.writerow(["##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"])
##        outfile.writerow(["##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"])
##        outfile.writerow(["##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"])
##        outfile.writerow(["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3"])
##
def getTranscripts(chromosome,gffDict):
    """
    getTranscript takes the fasta sequence of a Chromosome gtfDictionary [[[Chromosome,gene_id]= start,stop,orientation,feature  ]] file and returns the parsed location for the sampile
    """
    transcriptDict = {}
    
    for element, value in gffDict:
        if chromosome in element:
            print value[1]
    
def list2dict(String,start):
    """
    Takes a list(a String) and makes a dictionary, with the positions as the key, and the listcontent as values
    The index starts at the starting position of the transctipt in the chromosome
    """
    # switch indexing,   starting position of transcrip on the gene
    indexing = True
    Dict = {}
    
    if indexing  == True:
        count = start
    else:
        count = 0

    for piece in String:
        Dict[count] = piece
        count +=1
    return Dict


def wobble(indexSeq1,indexSeq2,position) :
    wobble = False
    """
    Wobble will catch the frameshift occurances, but also mask 25% of snps,so check the +-2 nucleotides for shifting
    """
    down_wobblespace = int(position+1)
    up_wobblespace = int(position-1)
    if nucleotide == sequenceDictionary[down_wobblespace] or nucleotide == sequenceDictionary[up_wobblespace]:
        wobble = True
        print "found wobble"

    return wobble
    

def Phred2prob(PHRED):
    """
    Takes the Phred score and translates it back to e-value probability 
    """
    probab = ((int(PHRED)*-1)/10)
    probability = '10e%i' %(float(probab))
    return float(probability)



def csv2dict(csv_reader):
    """
    allpile() takes a target transcript (name and sequence in dictionary form,  position and nucleotide)
    also takes the samfile pileup from pysam invocation
    returns a dictionary of the whole coverage of the transcript.
    """
    outDict = {}
    for row in csv_reader:
        outDict[row[0],row[1]] = [row[2],row[3],row[4],row[5],row[6]]
    return outDict



def flexpile(sequenceDictionary, genePileup, geneName, start, stop,secondpass):
    """
    allpile() takes a target transcript (name and sequence in dictionary form,  position and nucleotide)
    also takes the samfile pileup from pysam invocation
    returns a dictionary of the whole coverage of the transcript.
    """
    #print secondpass
    resDict = {}
    basecount = 1
    SNPcount = 0

    hitcount = 0
    negcount = 0
##    x_count = 0

##  can only go once through genePileup 
    for pileupcolumn in genePileup:

        if pileupcolumn.pos in sequenceDictionary:
            basecount = 0
            # nucount for G,A,T,C
            nucount = [0,0,0,0]
            bacount = [0,0,0,0]
            # element is the position on the chromosome
            element = pileupcolumn.pos
            
            # base is the corresponding nucleotide in the reference
            base = sequenceDictionary[element]


##            print secondpass
            if secondpass == 'No':
                

                for pileupread in pileupcolumn.pileups:
                    try:
                        ALT = pileupread.alignment.query_sequence[pileupread.query_position]
                    except:
                            # NoneType exception,  in case there is no alignment ? the read becomes 'X'    
                        ALT = 'X'

                    if ALT == base:
                        # DEBUG ??
##                    basecount +=1
                        bacount = nucl_compare(base, bacount)


                   # lets treat nonSNPs as snps for a moment
                    if ALT != base and ALT != 'X':
                        # each SNP has 3 possible states GATC - [ref]
                        # Dictionary will contain position{element} = ref_base,G, A, T, C,highest,??chance_for_random,poisson_exact_bonferroniholm??
                        
                        nucount = nucl_compare(ALT, nucount)

            if secondpass != 'No' :
                negcount +=1
##                print 'secondpass mode active'
                if base =='N':
                    for pileupread in pileupcolumn.pileups:
                        try:
                            ALT = pileupread.alignment.query_sequence[pileupread.query_position]
                        except:
                            # NoneType exception,  in case there is no alignment ? the read becomes 'X'    
                            ALT = 'X'

                        if ALT != 'X':
                            bacount = nucl_compare(base, bacount)
                        if ALT != base and ALT != 'X':
                            nucount = nucl_compare(ALT, nucount)
                else:
                    pass
            #print negcount
                
##            print basecount, nucount,   bacount
            basecount = max(bacount)
##            if max(bacount)==max(nucount):
##                print element, bacount, nucount
            nucount2dict(resDict,nucount,base,element,basecount)
                

    
    return resDict




def nucl_compare( ALT, nucount):
    """
counts occurance of SNP nucleotide and puts them in an ordered list
    """
    if ALT == 'G':
        nucount[0]+=1
    if ALT == 'A':
        nucount[1]+=1
    if ALT == 'T':
        nucount[2]+=1
    if ALT == 'C':
        nucount[3]+=1
    return nucount
    
def nucount2dict(resDict,nucount,base,element,basecount):
    """
write non-empty lists do dictionaries
    """
    if True in (t > 2 for t in nucount):

        resDict[element]=[base,nucount[0],nucount[1],nucount[2],nucount[3],basecount]
        
    return resDict
    



def classifydict(resDict_line):
    """
resDict, so far , contains the information of expression of different nucleotides
check if any of the nucleotides is expressed sufficiently, and write this as a SNP
Just compare reference,  most likely SNP and probability that this is at least heterozygous
    """
    snpDict = {}
    inDict = {}
    
    base = resDict_line[0]
    inDict['basecount'] = resDict_line[5]
    inDict['G'] = int(resDict_line[1])
    inDict['A'] = int(resDict_line[2])
    inDict['T'] = int(resDict_line[3])
    inDict['C'] = int(resDict_line[4])
    
##    allcount = sum(inDict.itervalues())
    allcount = (inDict['G'] + inDict['A'] + inDict['T'] + inDict['C'] + inDict['basecount'])
    
    # in case this is a secondpass, one of the nucleotides is larger than basecount
    most_likely = max(inDict.iteritems(), key=operator.itemgetter(1))[0]
    second = secondlargest(most_likely,inDict)
    

    # p is misscallchance for illumina..   lets say 1 (0.01) for the moment, can make this more flexible
    prob = stats.binom_test( inDict[most_likely],allcount, p=0.01)
    prob_2 = stats.binom_test( inDict[second],allcount, p=0.01)

    
    # build the output dictionary for the best and secondbest value
    # THIS IS TWOSIDED,   make it onesided, or its useless...
    mean = stats.binom.mean(allcount, p=0.01)
    stdev = stats.binom.std(allcount, p=0.01)
    minimum = (mean - (stdev*1.96))
    
    # build the output dictionary for the best and secondbest value
    # THIS IS TWOSIDED,   make it onesided, or its useless...
    if float(inDict[most_likely]) >= minimum and most_likely != 'basecount':
        snpDict[most_likely]=[inDict[most_likely],allcount,prob]

    if float(inDict[second]) >= minimum and second != 'basecount':
        snpDict[second]=[inDict[second],allcount,prob_2]    
    
    return snpDict



def secondlargest(largest,allDict):
    """
finding the second largest item
"""
    subDict = {}
    for element, values in allDict.items():
        if not element == largest:
            subDict[element] = values
    
    second_most_likely = max(subDict.iteritems(), key=operator.itemgetter(1))[0]
    return second_most_likely


def getPercentage(number, total, feature):
    percDict = {}
    percent = ((float(number) / float(total))*100)
    
    if percent % 10 < 0.02:
        
        if percent > 9:
            decadic_percent = round(percent,1)
            if not decadic_percent in percDict:
                print "%s percent of %ss analyzed " %(decadic_percent,feature)
                percDict[decadic_percent] = 1


### outdated functions

def allpile(sequenceDictionary, genePileup, geneName, start, stop,secondpass):
    """
    allpile() takes a target transcript (name and sequence in dictionary form,  position and nucleotide)
    also takes the samfile pileup from pysam invocation
    returns a dictionary of the whole coverage of the transcript.
    """
    resDict = {}
    basecount = 1
    SNPcount = 0
    x_count = 0

##  can only go once through genePileup 
    for pileupcolumn in genePileup:

        if pileupcolumn.pos in sequenceDictionary:

            # element is the position on the chromosome
            element = pileupcolumn.pos
            key = element
            
            # base is the corresponding nucleotide in the reference
            base = sequenceDictionary[element]

            # implement option to run on the second try, ignore all non masked nucleotides
            if secondpass == False or base == 'N':


                # ?? what ??
                if pileupcolumn.pos == element:
        
                    for pileupread in pileupcolumn.pileups:
                        
                        try:
                            ALT = pileupread.alignment.query_sequence[pileupread.query_position]
                        except:
                            # NoneType exception,  in case there is no alignment ? the read becomes 'X'
##                            ALT = base
                            ALT = 'X'
                            x_count +=1

                        if ALT !='X':
                            basecount +=1

##                print ALT, base
                        if ALT != base and ALT != 'X':
##                            basecount +=1

                            if not element in resDict:
                                resDict[element] = [base,ALT]
                                SNPcount +=1

                            # Translate PHRED to probability
                                probability = Phred2prob(int(pileupread.alignment.mapping_quality))
                                count_occurance = 1

                                resDict[element]=[sequenceDictionary[pileupcolumn.pos],ALT,count_occurance,basecount,probability]
                            
 

                # Add more information about the probability            
                            if element in resDict:          

                            # allow for multiple different SNPs per position
                                if resDict[element][1] == ALT:

                                    probability = Phred2prob(int(pileupread.alignment.mapping_quality))
                                
                                    counting = resDict[element][2]
                                    counting += 1
                                    resDict[element][2] = counting

                                    misscallChance = resDict[element][4]
                                    misscallChance = float(misscallChance) * float(probability)
                                    resDict[element][4] = misscallChance

                            
                                if resDict[element][1] != ALT:
                                    pileupcolumn_alternative= '%s_%s' %(element, ALT)
                                    if pileupcolumn_alternative in resDict:

                                # Translate PHRED to probability
                                        probability = Phred2prob(int(pileupread.alignment.mapping_quality))
                                        
                                        count_occurance = resDict[pileupcolumn_alternative][2]
                                        count_occurance += 1
                                        resDict[pileupcolumn_alternative][2] = count_occurance
                                        
                                        misscallChance = resDict[pileupcolumn_alternative][4]
                                        misscallChance = float(misscallChance) * float(probability)
##                                        
                                        
                                        SNPcount += 1
                          

                                    if not pileupcolumn_alternative in resDict:
                                        
                                        count_occurance = 1
                                        probability = Phred2prob(int(pileupread.alignment.mapping_quality))
                                        
                                        resDict[pileupcolumn_alternative]=[sequenceDictionary[key],pileupread.alignment.query_sequence[pileupread.query_position],count_occurance,basecount,probability]
                                        SNPcount +=1

                            # Translate PHRED to probability
                                        probability = Phred2prob(int(pileupread.alignment.mapping_quality))

                        if pileupcolumn.pos in resDict and ALT != 'X':
                            resDict[element][3] = basecount
                        
                            
        
    return resDict
