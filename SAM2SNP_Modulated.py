import pysam
from Bio import SeqIO
from _functions import flexpile
from _functions import list2dict
from _functions import wobble
from _functions import classifydict
from _functions import getPercentage
import re
import csv
import sys,argparse
import os.path

# arguments for commandline input and help
####################################################
parser = argparse.ArgumentParser(description='Looking for SNPs in RNAseq data. ')



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
                    help='Input a gtf file to check for the SNPS on gene ### Try to remove necessity for this',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-bam',
                    dest='bam',
                    required = True,
                    help='Input the bam file of the Alignment.  Tophat2 output worksfine,  NoxtGenMapper or STAR should be fine too.',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-out',
                    dest='out',
                    required = False,
                    default='output.tab',
                    help='Output file,  so far in undefined format  gene snp-position original_base snp_base count_of occurance accumulated_probability',
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


args = parser.parse_args()

#####################################################

with open("%s" %(args.fasta), "rU") as fasta_raw, open("%s"%(args.gtf),"r") as gff_raw, open('%s' %(args.out),'w') as out_raw:
    samfile = pysam.AlignmentFile("%s" %(args.bam),"rb")
    records = list(SeqIO.parse(fasta_raw, "fasta"))
    gff_file = csv.reader(gff_raw, delimiter = '\t')
    outfile = csv.writer(out_raw, delimiter = '\t')

    fastadict = {}
    gffDict = {}
    basecount = 0
    genecount = 0
    resultDict ={}
    SNPDict = {}
    writecount = 0


    # populate gff dictionary
    for row in gff_file:
        # skip header
        if not '#' in row[0]:
            try:
                if row[2] == args.feature:
                # legacy,  this has to be more flexible
                    gffDict[row[0],row[8].split('"')[1]] = [row[3],row[4],row[6],row[2]]
            except:
##                print 'wrong gff configuration in %s' %(row[0])
                pass
    print "gff file with %s features" %(len(gffDict))


    print " starting analysis "
    print "0 percent of %s analyzed " %(args.feature)
    for element in records:
        
        for key,value in gffDict.items():
            fastadict = {}
            
            # if on chromosome
            if key[0] == element.id:
                start = int(value[0])
                stop = int(value[1])
                
                gene = key[1]
                
                fastadict = list2dict(element.seq[start:stop],start)
                             
                genecount +=1
                # print gene count in percent of features  every 10 percent
                getPercentage(genecount,len(gffDict),args.feature)
                
                resultDict[gene]= flexpile(fastadict,samfile.pileup("%s" %(element.id),start,stop),element.id, start, stop, args.spass)
                    

    

    single_SNPDict = {}
    print ' Writing to file ' 
    for element, values in resultDict.items():
        for elementx, valuex in values.items():
            
            single_SNPDict = classifydict(valuex)
            for elementn, valuen in single_SNPDict.items():
                if int(valuen[0]) > 0:
##                print 'element = %s,  elementx = %s, valuex[0] = %s' %(element, elementx, values) 

                    outfile.writerow([element, elementx, valuex[0], elementn, valuen[0],valuen[1],valuen[2]])
                    writecount += 1
    print '%s SNPs written to file' %(writecount)
            


                                                                                         
