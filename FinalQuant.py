import csv
import functions


with open('../try2016/testoutquant.txt','r') as in_raw, open('../try2016/testoutquant.out.txt','w') as out_raw:
    infile = csv.reader(in_raw, delimiter = '\t')
    outfile = csv.writer(out_raw, delimiter = '\t')

    geneDict = {}

    for row in infile:
        # ignore the header
        if not '#' in row[0]:
            # populate gene dictionary
            if not row[0] in geneDict:
                geneDict[row[0]] = {}

        
            ORIG = row[2]
            ALT = row[3]
            OBSERV = row[4]
            TOTAL = row[5]
            PROB = row[6]
            SYN = row[7]
            
            geneDict[row[0]][row[1]] = [ORIG,ALT,OBSERV,TOTAL,PROB,SYN]

    
    print len(geneDict)
        
    for gene, snps in geneDict.items():
        obsList = [snps[2]]
        totlist = [snps[3]]
##    
    
