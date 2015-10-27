import pysam
from Bio import SeqIO
from _functions import flexpile
from _functions import list2dict
from _functions import wobble
from _functions import classifydict
from _functions import getPercentage
from _functions import parse_gff
from _functions import mutateSequence
from _functions import findSyn
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq
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
                    help='Input a gtf file to check for the SNPS on gene, ### Try to remove necessity for this, there is no technical need for it',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-bam',
                    dest='bam',
                    required = True,
                    help='Input the bam file of the Alignment.  Tophat2 output works,  NextGenMapper or STAR should be fine too.',
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
                    default='1.0',
                    help='Minimum cutoff, a simple p value from a binomial test to discern whether the coverage could have been caused by the Illumina technical error of 1-1.7%. Default is 1 (off), can be set to 0.05',
                    metavar = 'FLOAT',
                    #type=argparse.FileType('w')
                    )

args = parser.parse_args()

#####################################################

with open("%s" %(args.fasta), "rU") as fasta_raw, open("%s"%(args.gtf),"r") as gff_raw, open('%s' %(args.out),'w') as out_raw:
    print '[IMPORT] Loading BAM File'
    samfile = pysam.AlignmentFile("%s" %(args.bam),"rb")
    print '[IMPORT] Loading Fasta file %s' %fasta_raw
##    records = list(SeqIO.parse(args.fasta, "fasta"))
    print '[IMPORT] Loading GFF/GTF File'
    print '[IMPORT] Loading Complete.'
    
    gff_file = csv.reader(gff_raw, delimiter = '\t')
    outfile = csv.writer(out_raw, delimiter = '\t')

    fastadict = {}
    gffDict = {}
    basecount = 0
    genecount = 0
    resultDict ={}
    SNPDict = {}
    writecount = 0
    origin = 'OTHER'

    print '[STATUS] Initiating gff parsing'
    gffDict = parse_gff(gff_file, origin,args.feature)
    if (len(gffDict)) == 0:
        print "[STATUS] The gff/gtf file has 0 features,  check the files structure or the existence of the selected feature"
        exit 
    else:
        print "[STATUS] gff file with %s features" %(len(gffDict))
    

    
    print "[STATUS] starting Gene-wise analysis .. This may take a while"
    print "[STATUS] 0 percent of %ss analyzed " %(args.feature)
    perc_list = []
    vcfSubDict = {}
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
                
##                print '%s genes' %(genecount)
                resultDict[gene]= flexpile(fastadict,samfile.pileup("%s" %(element.id),start,stop),element.id, start, stop, args.spass)

                # create the syn / nonsyn differentiation directly in the classification loop should be fastest
                # this loop writes the output

                # translate fasta of the gene to protein
                DNA = Seq(str(element.seq[start:stop]),generic_dna)
                protein = DNA.translate()
                protein = list2dict(protein[0:len(protein)],0)
                positionList = []
                for sub_element, value in resultDict[gene].items():
                    
                    converted = classifydict(value)
                    for nucleotide, valuen in converted.items():
                        if int(valuen[0]) > 0 and float(valuen[2]) < float(args.minqual):
                            count +=1


                            # now check for synonymity,  move this to functions
                            alternativeSeq = MutableSeq(str(element.seq[start:stop]), generic_dna)
##                            print nucleotide
                            alternativeSeq = mutateSequence(alternativeSeq,sub_element,nucleotide,start)
                            alternativeSeq = Seq(str(alternativeSeq), generic_dna)
                            altprot = alternativeSeq.translate()
                            altprot = list2dict(altprot[0:len(protein)],0)

                            protposition = int((sub_element-start)/3)
                            if protein[protposition] != altprot[protposition]:
                                positionList.append(sub_element)
                                synonym = 'NonSynon'
                            if protein[protposition] == altprot[protposition]:
                                synonym = 'Syn'

                            # correct against python specific error of starting count at 0. 

                            truepos = int(sub_element)+1
                            outfile.writerow([gene, truepos, value[0], nucleotide, valuen[0],valuen[1],valuen[2],synonym])
                            writecount += 1
##                print count




##    single_SNPDict = {}
##    print '[EXPORTING]  Writing to file ' 
##    for element, values in resultDict.items():
##        for elementx, valuex in values.items():
##            
##            single_SNPDict = classifydict(valuex)
##            for elementn, valuen in single_SNPDict.items():
##                if int(valuen[0]) > 0:
####                print 'element = %s,  elementx = %s, valuex[0] = %s' %(element, elementx, values) 
##                    # this is the VCF like format,
##                    # gene, position, REF, ALT,...
##                    if float(valuen[2]) < float(args.minqual):
##                        outfile.writerow([element, elementx, valuex[0], elementn, valuen[0],valuen[1],valuen[2]])
##                        writecount += 1

    if int(writecount) > 0:
        print '[RUN SUCCESSFUL] %s Preliminary SNPs written to file %s' %(writecount,args.out)
    else:
        print '[RUN FAILED] No Preliminary SNPs written to file' %(writecount)
            


                                                                                         
