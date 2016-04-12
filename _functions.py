import pysam
from Bio import SeqIO
from scipy import stats
import operator
import random
# Funtions for the RNAseq SNP expression 
from Bio.Alphabet import IUPAC
import math



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
def CalculateDotProduct(Vector1, Vector2):
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
    
    

    

def findFutureProbabilities(currentProbability, decrement):
    """
    Support function for Hillclimbing to iteratively find probabilities
    
    """
    
    probability_up = float(currentProbability) + float(decrement)
    probability_down = float(currentProbability) - float(decrement)
    
    return probability_up, probability_down



def HillClimbing(total, normcount,baseexpectation,stepSize = 0.25 ):
    """
    Function implements simple hillclimbing algorithm to find optimal expectations from input data
    
    I could just compute that from the observations,
    
    10 iterations cover sufficient detail, starting in the middle, it will converge to an optimum,
    lack of flexibility of the possible end states means, this acts as a pre clustering . 
    
    
    """
    up = '+'
    down = '-'
    
    #stepSize = 0.25
    
    # starting with the expectation that the probability is heterozygous
    # 10 iterations, decreasing Stepsize   should do 
    
    
    current_prob = pmfNegBin(total,normcount,baseexpectation)
    print 'first probability = ',current_prob 
    # give it a default direction
    direction = up
    for iterations in range(1,4):
        
        # first step = decide directions:
        
        
        # decrement should reach all cutoffs within 3 iterations
        
        decrement = (float(stepSize) / float(iterations))
        #print decrement
        
        # decide to go up or down, give both choices, see what is better:
            
        upExpect, downExpect = findFutureProbabilities(baseexpectation,decrement)
            
        upprob = pmfNegBin(total,normcount,upExpect)
        downprob = pmfNegBin(total,normcount,downExpect)

        
         
        if float(upprob) > float(current_prob) or float(downprob) > float(current_prob):
            if float(upprob) > float(downprob):
                direction = '+'
                current_prob = upprob
                baseexpectation = upExpect
                
                
            elif float(upprob) < float(downprob):
                direction = '-'
                current_prob = downprob
                baseexpectation = downExpect
                
        
    return baseexpectation,  current_prob
            
        

#         for element in prob:
#             value = pmfNegBin(total,normcount,baseexpectation)
#             result.append(value)
#             
#             
        



def ClassifyPreprocess(ListofMeasurements, maxVal):
    """
    This function takes the normalized (to percentage on 1-100 scale) input count data and returns the likelyhood of the measurement following a distinct expected population distribution
    
    output will go to the SLP for classification
    
    input:
    list of values from measurements of replicates
    base expectation
    
    Here I could run an optimization to see which distribution best explains the occurance,
    Implement hillclimbing algorythm ?
    
    
    
    """
    
    baseexpect = 0.5
    
    for element in ListofMeasurements:
        HillClimbing(maxVal,ListofMeasurements, baseexpect)
    
    
    
    
    




def SLPClassification(ObservList, weightList, Theta):
    """
    simple implementation of linear classification algorithm 
    It's a simple formula that only permits linear classification, 
    In the case of noise or occurence, only one dimensional data is available.
    
    observList is a list of observations,  in this case most likely probabilities from a negative binomial distribution
    
    weightList refers to a list of weights applied to the observations.
    Theta is the threshold level, this has to be kept as a global variable to keep it dynamic
    
    TrueList is the weight corrected observationList
    
    Returns the Classification variable  Yes (1) or No (0)
    
    """
    
    TrueList = []
    count = 0
    
    
    dotProduct = CalculateDotProduct(ObservList, weightList)

    if float(dotProduct) < float(Theta):
        Classification = 0
        
    if float(dotProduct) >= float(Theta):
        Classification = 1
        
    return Classification


# add the SAM2SNP first, move as much as possible to _functions

def writeVCFormat(masterDict,outfile,convDict):
    
    """
    Simple output for vcf like format based on cvs reader / writer functions
    
    input is the original masterDictionary,
    an outfile handler to write into 
    and a conversion dictionary, this contains the conversion between gene names, and the chromosomes they are on 
    
    """
    #write Header

    for gene, SNPs in masterDict.items():
        CHROM = convDict[gene]

        for SNP, info in SNPs.items():
            POS = SNP
            ID = '%s_%i' %(gene,int(POS))
            REF = info[0]
            ALT = info[1]
            QUAL = 30
            Filter = '.'
            INFO = 'NS=%s;DP=%s' %(info[5],info[2])
            
            outrow = []
            for elmnt in CHROM,POS,ID,REF,ALT,QUAL,Filter,INFO:
                
                # elmnt = elmnt.replace('"','')
                outrow.append(elmnt)
        
            outfile.writerow(outrow)


def StaticWeightList( repeats, weight):
    """
    small support function for the weightlistcreation
    populates the SNPs with a weight vector of lenght repeat
    matches the reference Dictionary column line
    """
    
    weightListing = []
    
    for repeat in range(0,repeats):
        weightListing.append(weight)
    
    return weightListing


def PopulateWeightdict(masterDict,repeats):
    """
    Populate a dictionary to fit the resampled masterDict.
    this is needed to generate the weight vector for the Noise classification step
    in the initial SNP detection and validation pipeline.
    
    input:
    masterDictionary
    repeat of resampling strategy (same as number of columns that will be present in the reference Dictionary
    
    output:
    Dicitonary of Dictionary:  {Genes:{SNPs:Weights}}
    
    
    """
    weightDict = {}
    
    for key, value in masterDict.items():
        
        for element, content in value.items():
            
            if key in weightDict:
                weightDict[key][element]= StaticWeightList( repeats, content[5])
            
            if not key in weightDict:
                
                SNPWeightList = StaticWeightList( repeats, content[5])
                weightDict[key]={element:SNPWeightList}
                
    return weightDict




def likelihoodTest(subDictionary, masterDictionary, probability = 0.98, minCoverage = 100):
    """
    Function takes the subdictionaries provided by the individual samples, and tests them against the likelihood / pmf of the NegBin Distribution
    """
    for key, value in subDictionary.items():
        


            ### probability needs to be flexible
        probability = float(probability)
            
        total = int(value[3])
        observations = int(value[2])
            # penalize lack of coverage
        if int(total) < minCoverage:
            total = minCoverage
                
            
        likelihood = pmfNegBin(total,observations,probability, True)

        repeatedcoverage = masterDictionary[key][0] + 1
        masterDictionary[key] = [repeatedcoverage , likelihood]
        
    return masterDictionary


##def checkSynonymity(Dictionary):
##    """
##    Validate classification of genetic code to identify Synonymous SNPs or nonsynonoymous SNPs
##
##    Input is a Dictionary, containing flexpile() output of SNPs per gene for multiple genes / locations
##    """ 
##    
##    for element, value in Dictionary.items():
##        converted = classifydict(value)
##
##        for nucleotide, valuen in converted.items():
##            if int(valuen[0]) > 0 and float(valuen[2]) < float(args.minqual):
##                count +=1
##
##
##                            # now check for synonymity,  move this to functions
##                alternativeSeq = MutableSeq(str(element.seq[start:stop]), generic_dna)
####                            print nucleotide
##                alternativeSeq = mutateSequence(alternativeSeq,sub_element,nucleotide,start)
##                alternativeSeq = Seq(str(alternativeSeq), generic_dna)
##                    if len(alternativeSeq)%3 != 0:
##                        overlap = len(alternativeSeq)%3
##                        alternativeSeq = alternativeSeq[:-int(overlap)]
##                            
##                    altprot = alternativeSeq.translate()
##                    altprot = list2dict(altprot[0:len(protein)],0)
##
##                    protposition = int((sub_element-start)/3)
##                    try:
##                        if protein[protposition] != altprot[protposition]:
##                            positionList.append(sub_element)
##                            synonym = 'NonSynon'
##                        if protein[protposition] == altprot[protposition]:
##                            synonym = 'Syn'
##                    except:
##                        synonym = 'Unknown'
##
##    return 
def varNegBin(oserv_success,prob):
    """
    get the variance of the Negative Binomial Distribution 
    
    oserv_success = observed successes
    prob = probability of success
    
    """
    variance = (oserv_success*(1-prob))/(prob*prob)
    return variance




def pmfNegBin(total,observ,prob):
    """
    input
    Total       : total number of measurements
    Observ      : number of observations of bernoulli outcome
    Prob        : Expected likelyhood for observing an Observ
    
    Calculates probability Mass function for the negative binomial distribution
    A parametric approach to validate likelyhood of false claims.
    If prob is fixed on technical limitations, multiplicity correction of obtained p values may be necessary...
    """
    
    obsdf = int(observ - 1)
    totdf = int(total - 1)
    if totdf < 1:
        totdf = 1
    failures = int(total) - int(observ)

    if failures > 0:
        failuresdf = failures - 1
    elif failures == 0:
        failuresdf = 0
        

    if obsdf < 0:
        obsdf = 0


    a = math.factorial(totdf)
    b = math.factorial(obsdf)
    c = math.factorial(total-observ)
    bindiv = a // (b * c)

    term1 = math.pow(prob,(observ))
    antiprob = float(1-prob)
    term2 = math.pow(antiprob,(failures))
    try:
        arg = float(bindiv)*term1*term2
    except:
        arg = 0
##        print 'Error in evaluating %s   %s  ' %(observ, total)

    return arg


def distFunNegBin(total,observ,prob):
    """
    extend the pmf function into a distribution function
    adding a sum over all failures
    """
    failures = (total - observ)
    failures = int(failures)
    distProb = 0


    speedy = True

    
    for f in range(0,failures):
    
        newtotal = (f + observ)
        distProb = distProb + pmfNegBin(newtotal,observ,prob)
##        print distProb
        # add speedup by limiting distfunction

        # Speedy is for significance threshold determination only.
        # It will stop the loop at the threshold.

        # since this is a very time consuming step for large coverage,  speedy will decrease runtime 20 fold without major drawbacks if it is only
        # used for significance comparison

        
        if speedy == True:
            if distProb > 0.06:
                break
        
    return distProb
    


def mutateSequence(sequence,vcfposition,altnucleotide,start):
    """
    sequence contains protein sequence for the gene,
    vcfDict should contain the ALT nucleotide and the position of the SNP
    start is necessary since the vcfDict contains the absolute position, the relative position of the SNP is discerned from the starting value
    sequence is a mutable nucleotide fasta
    the identification number is record.id for the fasta seq,  so only the snps from that chromosome are considered
    the vcfDictionary contains the possible SNPs (prevalidation)
    the start position is the start of the gene on the chromosome, important for the conversion between vcfDict starting value and gene start
    """
    position = int(vcfposition) - int(start)
    sequence[position] = altnucleotide        
    return sequence


# transfor original sequence into dictionary,  count for keys.  speeds up search
def findSyn(position_in_prot,vcfposition,altnucleotide,start):
    """
    proteinDict was created with the above mentioned function,
    proteinDictAlt is a mutable string from Bio, containing the alternative Protein sequence
    """
##    count_func = 0
##    snpChromosome = {}
##    snpPositions = []
##    for element, value in proteinlistAlt.items():
##        
##        if proteinDict[count_func] != value:
####            print proteinDict[count_func] , value
##            
##            snpPositions.append(count_func)
####            print 'found change in %s' %(count_func)
##            count_func +=1

    
            # Python seems to always round to the lower integer. POSSIBLE BUG- need to improve the rounding algorythm
            
    position = (int(vcfposition) - int(start)) / 3
##    print 'looking for change in %s' %(position), 
    if position == snpPositions_in_prot:
        synity = 'NonSyn'
    if position not in snpPositions:
        synity = 'Syn'   
    
    return synity

def parse_gff(gff_file, origin,feature='gene'):
    """
    Parse a gff file,  extend this to encompass USCS /  SGD/CGD IDs
    'gene' is default
    """
    gffDict = {}
    
    for row in gff_file:
        if not '#' in row[0]:
            try:
                if row[2] == feature:
                    
                    # add more formats as options for GFF / GTF parsing

                    if origin == 'ENSEMBL':
                        name = re.split('=|;',row[8])[2]
                        gffDict[row[0],name] = [row[3],row[4],row[6],row[2]]

                    elif origin == 'SGD':
                        ident  = row[8].split(';')[0]
                        name = ident.split('=')[1]
                        
                        gffDict[row[0],name] = [row[3],row[4],row[6],row[2]]
                    else:
                # legacy,  this has to be more flexible
                        gffDict[row[0],row[8].split('"')[1]] = [row[3],row[4],row[6],row[2]]
                            
            except:
##                print 'wrong gff configuration in %s' %(row[0])
                pass
    return gffDict





def writevcfheader(switch, outfile):
    """
    simple header in the same format as GATKs vcf files.
    """
    if switch == True:
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
    convert the csv reader into a proper dictionary
    """
    outDict = {}
    for row in csv_reader:
        outDict[row[0],row[1]] = [row[2],row[3],row[4],row[5],row[6]]
    return outDict



def csv2dictIncSyn(csv_reader):
    """
    IF SYN, NONSyn information is given, that has to be added from row[7],  keep this optional, to enable vcf support  
    """
    outDict = {}
    for row in csv_reader:
        outDict[row[0],row[1]] = [row[2],row[3],row[4],row[5],row[6],row[7]]
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
    



def classifydict(resDict_line,ErrorExpect, probabilityCutoff, mincount):
    """
THe analysis, so far contains the information of expression of different nucleotides, for each nucleotide
check if any of the nucleotides is expressed sufficiently, and write this as a SNP
Just compare reference,  most likely SNP and probability that this is at least heterozygous
    """
    snpDict = {}
    inDict = {}

    ErrorExpect = float(ErrorExpect)
    probabilityCutoff = float(probabilityCutoff)
    mincount = int(mincount)

    
    base = resDict_line[0]
    inDict['basecount'] = resDict_line[5]
    inDict['G'] = int(resDict_line[1])
    inDict['A'] = int(resDict_line[2])
    inDict['T'] = int(resDict_line[3])
    inDict['C'] = int(resDict_line[4])
    
##    allcount = sum(inDict.itervalues())
    allcount = (inDict['G'] + inDict['A'] + inDict['T'] + inDict['C'] + inDict['basecount'])


    if allcount > int(mincount):
    # in case this is a secondpass, one of the nucleotides is larger than basecount
        most_likely = max(inDict.iteritems(), key=operator.itemgetter(1))[0]
        second = secondlargest(most_likely,inDict)
    

    # p is misscallchance for illumina..   1.4% for the moment, can make this more flexible
    # use negative binomial density function instead of binom.test

## legacy functionality   
####    prob = stats.binom_test( inDict[most_likely],allcount, p=0.01)
####    prob_2 = stats.binom_test( inDict[second],allcount, p=0.01)

    
        prob = distFunNegBin(allcount,inDict[most_likely],ErrorExpect)
        prob2 = distFunNegBin(allcount,inDict[second],ErrorExpect)

        if prob < probabilityCutoff and most_likely != 'basecount':
            snpDict[most_likely]=[inDict[most_likely],allcount,prob]
        if prob2 < probabilityCutoff and second != 'basecount':
        # an alternative second choice, in case both alleles mutated away from the reference,  or the organism is polyploid
            snpDict[second] = [inDict[second],allcount,prob2]

####    
####    # build the output dictionary for the best and secondbest value
####    # THIS IS TWOSIDED,   make it onesided, or its useless...
####    mean = stats.binom.mean(allcount, p=0.01)
####    stdev = stats.binom.std(allcount, p=0.01)
####    minimum = (mean - (stdev*1.96))
##    if float(inDict[most_likely]) >= minimum and most_likely != 'basecount':
##        snpDict[most_likely]=[inDict[most_likely],allcount,prob]
##
##    if float(inDict[second]) >= minimum and second != 'basecount':
##        snpDict[second]=[inDict[second],allcount,prob_2]    

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


def getPercentage(number, total, feature, percList):
    percent = ((float(number) / float(total))*100)
    decadic_percent = round(percent,0)
    if decadic_percent % 10 < 0.5:
        
        if percent > 9:
            
            if str(decadic_percent) not in percList:
                print "[STATUS] %s percent of %ss analyzed " %(decadic_percent,feature)
                percList.append(str(decadic_percent))

    return percList

### outdated functions
# 
# def allpile(sequenceDictionary, genePileup, geneName, start, stop,secondpass):
#     """
#     allpile() takes a target transcript (name and sequence in dictionary form,  position and nucleotide)
#     also takes the samfile pileup from pysam invocation
#     returns a dictionary of the whole coverage of the transcript.
#     """
#     resDict = {}
#     basecount = 1
#     SNPcount = 0
#     x_count = 0
# 
# ##  can only go once through genePileup 
#     for pileupcolumn in genePileup:
# 
#         if pileupcolumn.pos in sequenceDictionary:
# 
#             # element is the position on the chromosome
#             element = pileupcolumn.pos
#             key = element
#             
#             # base is the corresponding nucleotide in the reference
#             base = sequenceDictionary[element]
# 
#             # implement option to run on the second try, ignore all non masked nucleotides
#             if secondpass == False or base == 'N':
# 
# 
#                 # ?? what ??
#                 if pileupcolumn.pos == element:
#         
#                     for pileupread in pileupcolumn.pileups:
#                         
#                         try:
#                             ALT = pileupread.alignment.query_sequence[pileupread.query_position]
#                         except:
#                             # NoneType exception,  in case there is no alignment ? the read becomes 'X'
# ##                            ALT = base
#                             ALT = 'X'
#                             x_count +=1
# 
#                         if ALT !='X':
#                             basecount +=1
# 
# ##                print ALT, base
#                         if ALT != base and ALT != 'X':
# ##                            basecount +=1
# 
#                             if not element in resDict:
#                                 resDict[element] = [base,ALT]
#                                 SNPcount +=1
# 
#                             # Translate PHRED to probability
#                                 probability = Phred2prob(int(pileupread.alignment.mapping_quality))
#                                 count_occurance = 1
# 
#                                 resDict[element]=[sequenceDictionary[pileupcolumn.pos],ALT,count_occurance,basecount,probability]
#                             
#  
# 
#                 # Add more information about the probability            
#                             if element in resDict:          
# 
#                             # allow for multiple different SNPs per position
#                                 if resDict[element][1] == ALT:
# 
#                                     probability = Phred2prob(int(pileupread.alignment.mapping_quality))
#                                 
#                                     counting = resDict[element][2]
#                                     counting += 1
#                                     resDict[element][2] = counting
# 
#                                     misscallChance = resDict[element][4]
#                                     misscallChance = float(misscallChance) * float(probability)
#                                     resDict[element][4] = misscallChance
# 
#                             
#                                 if resDict[element][1] != ALT:
#                                     pileupcolumn_alternative= '%s_%s' %(element, ALT)
#                                     if pileupcolumn_alternative in resDict:
# 
#                                 # Translate PHRED to probability
#                                         probability = Phred2prob(int(pileupread.alignment.mapping_quality))
#                                         
#                                         count_occurance = resDict[pileupcolumn_alternative][2]
#                                         count_occurance += 1
#                                         resDict[pileupcolumn_alternative][2] = count_occurance
#                                         
#                                         misscallChance = resDict[pileupcolumn_alternative][4]
#                                         misscallChance = float(misscallChance) * float(probability)
# ##                                        
#                                         
#                                         SNPcount += 1
#                           
# 
#                                     if not pileupcolumn_alternative in resDict:
#                                         
#                                         count_occurance = 1
#                                         probability = Phred2prob(int(pileupread.alignment.mapping_quality))
#                                         
#                                         resDict[pileupcolumn_alternative]=[sequenceDictionary[key],pileupread.alignment.query_sequence[pileupread.query_position],count_occurance,basecount,probability]
#                                         SNPcount +=1
# 
#                             # Translate PHRED to probability
#                                         probability = Phred2prob(int(pileupread.alignment.mapping_quality))
# 
#                         if pileupcolumn.pos in resDict and ALT != 'X':
#                             resDict[element][3] = basecount
#                         
#                             
#         
#     return resDict


def MaskAFasta (silenceDictionary, fastaDictionary, outHandle):
    """
    silenceDictionary contains SNP results as position in the key, and the chromosomename in the .
    ## formatted as a Dict of Dicts, with the genes outside and the position inside
    ## changed this to Silencedictionary

    
    fastaDictionary contains Fasta sequences. and other parts,  convert to fastaDict where the sequences are mutable.
    the outpu handle opens a new writable file

    conversionDict is for conversion from gene names to chromosomes.
    
    No reason to keep this as an independant script, attach to the Quantification Script.
    Maskafasta takes the results from SNP quantification, a fasta file and the direction for an out file ( an open handle)
    replaces the nucleotides specified in the vcf with 'N' to enable a more accurate remapping
    """
##    vcffile = csv.reader(vcf_raw, delimiter ='\t')  # outside the function, add to SNP2QUant
##    record_dict = SeqIO.index('%s' %(args.fasta), "fasta") # same here,  add outside function

    fastaDict = {}
            
    for Identitiy, sequence in fastaDictionary.items():
        seq_mutable = sequence.seq.tomutable()

        # Identitiy are gene names,  seq_mutable is the sequence
        fastaDict[Identitiy] = seq_mutable


##    for gene, SNP in vcfDictionary.items():
    for Position,chromosome in silenceDictionary.items():
##        for Position, data in SNP.items():
            
##        print element, values
##        The VCF file contains the absolute position of the SNP on the sequence, so just replace it with an N
# value has to be the chromosome !!  I get those from the gtf 
        try:
            fastaDict[chromosome][int(Position)] = 'N'
        except:
            print '[ERROR] could not silence SNP with position %s on gene %s' %(element,gene)

                
    for element, value in fastaDict.items():
        out_seqs = SeqIO.SeqRecord((value), id = element)
        SeqIO.write(out_seqs, outHandle, "fasta")
            

def combineSNPdata(ListMasterdict,ListIndict):
    """
    Combine SNP data in the following fashion:
    gene and position have to be ident, obviously
    
    ALT coverage will be combined, 
    REF coverage will also be combined. 
    # optionally, low coverage areas should be sum(ALT) but add(REF)
    the Number of Samples will be incremented for new observations


    if the ALT bases are different, the one with the higher coverage will be used
    """
    
    outList = []
    outList_alt = ['0']
    
    
    ##            ORIG,ALT,SNPcov,Basecov,synonym
    # possible cases:   alt overlap, simple addition
    # index 0 should always overlap. since it's the same fasta reference
    if ListIndict[0] == ListMasterdict[0] and ListIndict[1] == ListMasterdict[1]:
        
        outList.append(ListIndict[0])
        outList.append(ListIndict[1])
        combCount = int(ListIndict[2]) + int(ListMasterdict[2])
        outList.append(combCount)
        combTotal = int(ListIndict[3]) + int(ListMasterdict[3])
        
        outList.append(combTotal) 
        outList.append(ListIndict[4])
        
        # increment the number of samples by the new observation
        NS = int(ListMasterdict[5]) + 1
        outList.append(NS)

    # different alternative SNP    create a list of 2x length, the bigger coverage first, and seperate at the writing function
    
    if ListIndict[1] != ListMasterdict[1]:
        if int(ListMasterdict[2]) > int(ListIndict[2]):
            outList = ListMasterdict
#             for element in ListIndict:
#                 outlist.append(element)

            
        if int(ListIndict[2]) > int(ListMasterdict[2]):
            NS = 1
            outList = ListIndict
            outlist.append(NS)
#             for element in ListMasterDict:
#                 outlist.append(element)
        

    if outList[0] != '0' :    
        #print outList
        return outList #, outList_alt
##        Basecov    
    


def add2masterDict(indict, masterdict):
    
    """
    Extend the MasterDictionary to contain all possible SNPs,  
    """

    for gene, SNP in indict.items():
        for SNPpos, details in SNP.items():
            # NS : Number of Samples set to 1 for initialisation, this is incremented later
            NS = 1
            

            if gene in masterdict:

                if SNPpos in masterdict[gene]:
                    #NS = masterdict[gene][SNPpos][5]
                    
                    if masterdict[gene][SNPpos]:
                        masterdict[gene][SNPpos] = combineSNPdata(masterdict[gene][SNPpos],indict[gene][SNPpos])
                        
                if SNPpos not in masterdict[gene]:
                    masterdict[gene][SNPpos] = [indict[gene][SNPpos][0],indict[gene][SNPpos][1],indict[gene][SNPpos][2],indict[gene][SNPpos][3],indict[gene][SNPpos][5],NS]       
            if gene not in masterdict:
                masterdict[gene] = {SNPpos:[indict[gene][SNPpos][0],indict[gene][SNPpos][1],indict[gene][SNPpos][2],indict[gene][SNPpos][3],indict[gene][SNPpos][5],NS]}            

    return masterdict
 
def add2RefDict(indict, refDict):
     
    """
    Extend the Reference Dictionary to contain location information for all possible SNPs
    refDict only needs to contain unique identifiers and positions as well as initialized empty lists
     
    """
 
    for gene, SNP in indict.items():
        for SNPpos, details in SNP.items():
             
            # NS : Number of Samples set to 1 for initialisation, this is incremented later
            NS = 1
             
 
            if gene in refDict:
 
                if SNPpos in refDict[gene]:
                    #NS = masterdict[gene][SNPpos][5]
                     
                    if refDict[gene][SNPpos]:
                        refDict[gene][SNPpos] = []
                         
                if SNPpos not in refDict[gene]:
                    refDict[gene][SNPpos] = []   
            if gene not in refDict:
                refDict[gene] = {SNPpos:[]}            
 
    return refDict


def binning(inDict, masterDict):
    """
    Appends the SNP count in the individual observations to the MasterDictionary
    """
    
    for gene, SNP in masterDict.items():
        for Pos, emptyvalue in SNP.items():


##            print inDict[gene][Pos]
            SNPcount = 0
            if gene in inDict:
                if Pos in inDict[gene] :
                
                    SNPs = int(inDict[gene][Pos][2])
                    totalCoverage = int(inDict[gene][Pos][3])
                if Pos not in inDict[gene] :
                    SNPs = 0
                    totalCoverage = 100
                
            elif gene not in inDict :
                SNPs = 0
                totalCoverage = 100


                
            if totalCoverage < 100:
                totalCoverage = 100

            count = 0
        # now create bins,  randomly draw from one sample 50 times (of whole coverage) and store how many SNPSs you get
            for repeat in range(0, 100):
                choice = random.randrange(0,int(totalCoverage))
            
                if choice < SNPs:
                    count +=1

            masterDict[gene][Pos].append(count)
##            print inDict[gene][Pos]

    return masterDict

def binningMonteCarlo(inDict, masterDict):
    """
    Appends the SNP count in the individual observations to the MasterDictionary
    This improves the parametric classification.
    
    
    takes an individual dictionary ,
    populates the emptied master dicionary with observations from 
    
    
    """
    
    for gene, SNP in masterDict.items():
        for Pos, emptyvalue in SNP.items():


##            print inDict[gene][Pos]
            SNPcount = 0
            if gene in inDict:
                if Pos in inDict[gene] :
                
                    SNPs = int(inDict[gene][Pos][2])
                    totalCoverage = int(inDict[gene][Pos][3])
                if Pos not in inDict[gene] :
                    SNPs = 0
                    totalCoverage = 100
                
            elif gene not in inDict :
                SNPs = 0
                totalCoverage = 100


                
            # penalize lack of coverage by adding a minimum of 100
            if totalCoverage < 100:
                totalCoverage = 100

            count = 0
        # now create bins,  randomly draw from one sample 100 times (of whole coverage) and store how many SNPSs you get
        # this is equal to a normalization by turning measurements into percentages of total
            for repeat in range(0, 100):
                choice = random.randrange(0,int(totalCoverage))
            
                if choice < SNPs:
                    count +=1

            masterDict[gene][Pos].append(count)
##            print inDict[gene][Pos]

    return masterDict


# def extendMasterDict(masterDict,inDict1,inDict2=0,inDict3=0):
#     """
#     Filling the master Dictionary with actual values, we take 20 observations, choosen from randomly selected datasets, of which we have up to 3...
#     This should result in a poisson like distribution for the low expressed data, depending on their mean expression
# 
#     This function selects one of the replicates, from which the observations are chosen
#     """
# 
# 
# 
#     for fill in range(0,5):
#         if inDict2 !=0 and inDict3 !=0 :
#             chosenOne = random.choice(['inDict1','inDict2','inDict3'])
#         if inDict2 !=0 :
#             chosenOne = random.choice(['inDict1','inDict2'])
#         if inDict2 ==0 and inDict3 ==0:
#             chosenOne = 'inDict1'
# 
#         
#            # take value from a random dictionary 
#         if chosenOne == 'inDict3':
#             chosenOne = inDict3
#         if chosenOne == 'inDict2' :
#             chosenOne = inDict2
#         if chosenOne == 'inDict1':
#             chosenOne = inDict1
#             
# ##        print len(chosenOne)
# 
#             
#         masterDict = binning(chosenOne, masterDict)
# 
# 
#         
#       
#     return masterDict
#     
# 
# def populateMasterDict(inDict1,inDict2=0,inDict3=0):
#     """
#     Creates the MasterDictionary, and populates it with the keys for SNPS from up to 3 samples
#     """
#     masterDict = {}
#     if inDict2 == 0:
#         masterDict = add2masterDict(inDict1,masterDict)
#     if inDict3 == 0:
#         masterDict = add2masterDict(inDict1,masterDict)
#         masterDict = add2masterDict(inDict2,masterDict)
#     else:
#         masterDict = add2masterDict(inDict1,masterDict)
#         masterDict = add2masterDict(inDict2,masterDict)
#         masterDict = add2masterDict(inDict3,masterDict)
#     
#     
#     
#     
#     
#     return masterDict

def extendMasterDictbyResampling(masterDict,ListofDictionaries, repeats ):
    """
    Filling the master Dictionary with actual values, 
    to stabilize the method a resampling approach is benefitial,
    especially if only one replicate is available
    
    the original master dictionary is taken and filled with observations from 
    randomly drawn samples 
    
    sampling takes random observations, converts to percent of total coverage of the highest expressed variable
    
    """
    for repeat in range(0,repeats):

        chosenOne = random.choice(ListofDictionaries)
    
        masterDict = binningMonteCarlo(chosenOne, masterDict)
        
        
    return masterDict





def populateMasterDictfromList(inList):
    """
    Creates the MasterDictionary from a list of Dictionaries
    """
    masterDict = {}
    refDict = {}
    for element in inList:
        
        masterDict = add2masterDict(element,masterDict)
        refDict = add2RefDict(element, refDict)
#     if inDict2 == 0:
#         masterDict = add2masterDict(inDict1,masterDict)
#     if inDict3 == 0:
#         masterDict = add2masterDict(inDict1,masterDict)
#         masterDict = add2masterDict(inDict2,masterDict)
#     else:
#         masterDict = add2masterDict(inDict1,masterDict)
#         masterDict = add2masterDict(inDict2,masterDict)
#         masterDict = add2masterDict(inDict3,masterDict)
# 
#     
    return masterDict, refDict




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
from _functions_SAM2SNP import SAM2SNP

def SamImport(bamLocation,SAMList,feature,fasta_raw,gffDict,outfile,cutoff,spass,createIntermediate):
    samfile = pysam.AlignmentFile("%s" %(bamLocation),"rb")
    
    samDict = {}
    SAMList.append(samDict)
    samDict = SAM2SNP(feature,fasta_raw,samfile,gffDict,outfile,cutoff,spass, createIntermediate)
    
    return samDict,SAMList
