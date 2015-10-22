#  test file for binning function
import csv
from _functions import csv2dict
import random


# populate master dictionary

    
def add2masterDict(indict, masterdict):
    """
    Extend the MasterDictionary to contain all possible SNPs
    """
    for element, values in indict.items():
        masterdict[element] = []
    return masterdict


def binning(inDict, masterDict):
    """
    Appends the SNP count in the individual observations to the MasterDictionary
    """
    for element, value in inDict.items():
        SNPcount = 0
        SNPs = int(value[2])
        totalCoverage = int(value[3])
##        print value[2], value[3]
        if totalCoverage < 20:
            totalCoverage = 20

        count = 0
        # now create bins,  randomly draw from one sample 50 times (of whole coverage) and store how many SNPSs you get
        for repeat in range(0, 50):
##            print totalCoverage
            choice = random.randrange(0,int(totalCoverage))
            
            if choice < SNPs:
                count +=1

        masterDict[element].append(count)

    return masterDict


def populateMasterDict(inDict1,inDict2=0,inDict3=0):
    """
    Creates the MAsterDictionary, and populates it with the keys for SNPS from up to 3 samples
    """
    masterDict = {}
    if inDict2 == 0:
        masterDict = add2masterDict(inDict1,masterDict)
    if inDict3 == 0:
        masterDict = add2masterDict(inDict1,masterDict)
        masterDict = add2masterDict(inDict2,masterDict)
    else:
        masterDict = add2masterDict(inDict1,masterDict)
        masterDict = add2masterDict(inDict2,masterDict)
        masterDict = add2masterDict(inDict3,masterDict)
    
    return masterDict


def extendMasterDict(inDict1,inDict2=0,inDict3=0):
    """
    Filling the master Dictionary with actual values, we take 20 observations, choosen from randomly selected datasets, of which we have up to 3...
    This should result in a poisson like distribution for the low expressed data, depending on their mean expression
    """
    for fill in range(0,20):
        chosenOne = random.choice(['inDict','inDict2','inDict3'])
        if chosenOne == 'inDict':
            chosenOne = inDict
        if chosenOne == 'inDict2':
            chosenOne = inDict2
        if chosenOne == 'inDict3':
            chosenOne = inDict3
            
        print len(chosenOne)
        binning(chosenOne, masterDict)




with open('Firstpass/akue_10_br3_R1','r') as in_raw, open('Firstpass/akue_10_br2_R1','r') as in2_raw:
    infile = csv.reader(in_raw, delimiter ='\t')
    in2file = csv.reader(in2_raw, delimiter ='\t')
##    for f in [1..50]:

    # let's assume there are less than 3 files, but one is obligaroty
    
    inDict = csv2dict(infile)
    try:
        inDict2 = csv2dict(in2file)
    except:
        pass

    try:
        inDict3 = csv2dict(in3file)
    except:
        pass

    
    masterDict = populateMasterDict(inDict,inDict2)
 
    
##    print len(masterDict)
##    for fill in range(0,20):
##        chosenOne = random.choice(['inDict','inDict2'])
##        if chosenOne == 'inDict':
##            chosenOne = inDict
##        if chosenOne == 'inDict2':
##            chosenOne = inDict2
##        if chosenOne == 'inDict3':
##            chosenOne = inDict3
##            
##        print len(chosenOne)
##        binning(chosenOne, masterDict)
##

    for element, value in masterDict.items():
        print value
