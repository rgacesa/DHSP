#!/usr/bin/env python
import sys

if len (sys.argv) < 2 or len(sys.argv) > 2:
    print '------------------- USAGE -----------------------'
    print '> python binarizeNCBItax <input file>             '
    print '> binarizes NCBI tax files for use with tax mapper'
    print '> WARNING: execute only on gi_taxid_prot.dmp'
    print '>          AND gi_taxid_nucl.dmp'
    print '-------------------------------------------------'
    exit(-1)

print "---------- TAXONOMY BINARIZER  ---------"
inFile = sys.argv[1]
#d inFile = '/home/ranko/Development/DBs/ncbi_tax/gi_taxid_prot.dmp'
outFile = inFile.replace('.dmp','.bin')
if '.dmp' not in inFile:
    print " *** ERROR: file not .dmp file! *** "
    exit(-1)
print " --> input file: ",inFile
print " --> output file:",outFile
print " -----> STARTING BINARIZATION <----- "
cntr = 0
gicntr = 0
with open(inFile,'r') as input:
    with open (outFile,'w') as output:
        for line in input:
            # parse line
            lsplit = line.split('\t')
            gi = int(lsplit[0].strip())            
            taxid = int(lsplit[1].strip())
            foundGi = False
            while not foundGi:
                ''' LINE COUNTER '''   
                cntr +=1 
                if cntr % 10000000 == 0:
                    print '   --> ',cntr,'lines binarized'                                    
                ''' GI SEARCH '''
                gicntr +=1                
                if gi == gicntr:
                    taxid = hex(int(taxid)).replace('x','')
                    taxbin = taxid.zfill(8)
                    #print 'gi',gicntr, 'tax:',taxid.zfill(8)
                    output.write(taxbin)
                    foundGi = True 
                else: 
                    output.write('0'.zfill(8))
                    #print 'gi',gicntr, 'tax:','0'.zfill(8)                    