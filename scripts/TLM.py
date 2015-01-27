#!/usr/bin/env python

'''
Created on 2 Oct 2014

@author: ranko
'''

''' DEPENDENCIES (codes)
 -> BioPython (for BLAST XML parser)
 -> RBioTools (mini library of various python scripts and functions)
 -> parseNCBItaxonomy (mini library of various python scripts and functions)
 -> binarizeNCBItax.py (used for binirizing NCBI taxonomy)
 -> prepareNCBItaxFiles.py (used to create whole-taxonomy maps)
 -> processGTxtFile.py (used to merge whole-taxonomy maps with results)
 -> optional: dot (for drawing graphs)
 -> optional: ParseHMMERTable.py, parseHHblitz (for parsing HMMER / HHblitz results)
 
    DEPENDENCIES (databases)
 -> binarized NCBI taxonomy database
 -> prepared Tax Files
 -> optional: for species mapping of sequences, NCBI NR database is required
 '''
 
''' ************* BRIEF INSTALLATION README: ******************
1) make sure BioPython is installed (at least NCBIXML module)
2) make sure required files are in same folder as taxMapper.py 
   (RBioTools.py, binarizeNCBItax.py, binarizeNCBItax.py, prepareNCBItaxFiles.py, ParseHMMERTable.py, parseHHblitz.py, processGTxtFile.py)
-- to prepare NCBI taxonomy database, do the following --
3) prepare NCBI taxonomy database for taxMapper use: 
    a) download NCBI taxonomy (from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/); 
    for this tutorial, it is assumed taxonomy is in <path_to_tax> path (location doesn't matter, but see below)
    b) binarize taxonomy by running binarizeNCBItax.py <path_to_tax>/gi_taxid_prot.dmp
       and binarizeNCBItax.py <path_to_tax>/gi_taxid_nuc.dmp (see notes/warnings below)   
4) prepare whole taxonomy maps: 
    a) run prepareNCBItaxFiles.py <path_to_tax>/names.dmp <path_to_tax>/nodes.dmp
    put produced files into <taxmapper folder>/taxFiles
    
5) (optional) to prepare taxMapper output species mapping (option which also creates text
output in addition to graph output), download NCBI NR database -- 

******************* TLM USE ***********************
run the code with -h for list of command line options
typical run involves: 
 -I <blast result file, must be in xml form>
 -O <output name>
 -X <path to NCBI taxonomy database; not required if DB is put into ./ncbi_tax>
 
if detailed species -> sequences list file is to be produced, following options should be used:
 --speciesOutBlastDBcmd <blastdbcmd full path, requires NCBI blast+>*
 --speciesOutNRDB <NR database full path>
 --doSpeciesOut Y 
 (optional) --doSpeciesOutF Y
 
*: note: this option might (or might not) work with HMMER or HHblits results file, as long as
search was done against something with standard GIs (NCBI NR, NCBI refseq or similar database)

resulting mapping graph is produced in .gv (graphviz) format and can be converted to PDF, postscript
or any other format by use of graphviz package
 
--- NOTES /  WARNINGS: ---
-> code was developed and tested under UBUNTU 12.04.5 (precise pangolin); it was not
tested under any other OS and might or might not work; it is based on python2.7 and might
not work with different versions of python
-> NCBI databases (NR and taxonomy) might have bugs such as misedited sequences, incorrect headers, 
taxonomical links leading nowhere... it can cause artifacts in produced results; also make
sure NCBI taxonomy and NR (or any other NCBI database) used are same date / version as they do 
change over time
-> code is supplied "as is" and might contain bugs
-> preparing taxonomy databases will take a while
-> code has not been tested for nucleotide blast results
-> taxonomy mapping of results not produced by search of NCBI databases is EXPERIMENTAL! 
   it might or might not work and will take a while 
-> drawing entire NCBI taxonomy to genus or species level can (and will) crash computer
'''

''' IMPORTS '''
import copy
from Bio.Blast import NCBIXML
from RBioTools import extractAllGIs
from RBioTools import determineSearchFileType
from RBioTools import argToListI
from RBioTools import loadAllFastas
from RBioTools import getFaHeader
from RBioTools import getEvFromMFaHdr
from RBioTools import isMFAHdr
from RBioTools import fetchFastaFromBlastDB
from RBioTools import parseHeader
from parseNCBItaxonomy import loadTaxID
from parseNCBItaxonomy import initTax
from parseNCBItaxonomy import getTaxFromTaxID
from parseNCBItaxonomy import extractFromIdMap
from parseHMMER import parseHMMResults
from parseHHblitz import parseHHResults
from parseHHblitz import filterHHblits
import processGTxtFile as pgtf

''' --- basic data classes and stuff --- '''

''' represents data for one NCBI Tax '''
class TaxIdData():
    def __init__(self):
        self.hitEVs = []
        self.hitGIs = []
        self.bestEV = -1.0
        self.avgEV = -1.0
        self.nrHits = 0
        self.nrHitsPerc = 0
        self.worstEV = -1.0
        self.rank = ''
    def calcStats(self, maxHits=1):
        totalEV = float(0.0)
        minEV = 999999.0
        maxEV = -1.0
        self.nrHits = len(self.hitEVs)
        for ev in self.hitEVs:
            totalEV += float(ev)
            if ev < minEV:
                minEV = ev
            if ev >= maxEV:
                maxEV = ev  
        if float(self.nrHits) > 0:
            self.avgEV = totalEV / float(self.nrHits)
        else: 
            self.avgEV = -1
        self.bestEV = float(minEV)
        self.worstEV = float(maxEV)
        self.nrHitsPerc = float(self.nrHits)/float(maxHits)
        
def getTaxDirChild (_tdHash, _tax):
    #d print _tax
    if _tax in _tdHash.keys():
        rc = _tdHash[_tax].targetTax   
        if not rc.find('\\n') == -1:
            return rc[0:rc.find('\\n')]
        return rc
    else:
        return 'N/F' 

def findTaxDirChild(_tdHash,_taxN,_taxL):
    doGo = True
    while doGo:
        c = getTaxDirChild (_tdHash, _taxN)
        #d print c
        if c == 'N/F':
            doGo = False
            return 'N/F'
        if _taxL in c:
            rc = c
            if not rc.find('\\n') == -1:
                return c[0:rc.find('\\n')]
            return rc
            doGo = False
        _taxN = c
    return 'N/F'
   
def getTaxDirAllChildren(_tdHash,_taxN):
    doGo = True
    cList = []
    cList.append(_taxN)
    while doGo:
        c = getTaxDirChild (_tdHash, _taxN)
        #print c
        if c == 'N/F':
            doGo = False
            return cList
        else: 
            cList.append(c)
            _taxN = c
    return cList

       
def grabTaxDataFHash(_tdHash,_taxN):
    t = []
    doGo = True
    tC = _taxN
    t.append(tC)
    while doGo:
        tC = getTaxDirChild(_tdHash,tC)
        #d print tC
        if tC == 'N/F':
            doGo = False
        else:                
            t.append(tC)    
    return t     
        
''' also represents data for one NCBI tax, but also 
includes link to other tax from this one (ex: canis -> lupus)'''
class TaxIdDataDirected():
    def __init__(self, _startTax,_targetTax):
        self.targetTax = _targetTax
        self.hitEVs = []
        self.bestEV = -1.0
        self.worstEV = -1.0        
        self.avgEV = -1.0
        self.nrHits = 0
        self.nrHitsPerc = 0
        self.hitGIs = []
        self.startTax = _startTax
        self.rank = ''
        self.key = ''
        
    def calcStats(self, maxHits=1):
        totalEV = float(0.0)
        minEV = 999999.0
        maxEV = -1.0
        self.nrHits = len(self.hitEVs)
        for ev in self.hitEVs:
            totalEV += float(ev)
            if ev < minEV:
                minEV = ev
            if ev >= maxEV:
                maxEV = ev  
        if float(self.nrHits) > 0:
            self.avgEV = totalEV / float(self.nrHits)
        else: 
            self.avgEV = -1
        self.bestEV = float(minEV)
        self.worstEV = float(maxEV)
        self.nrHitsPerc = float(self.nrHits)/float(maxHits)

''' function: defines coloring depending on evalue or number of hits '''
def getGColor(cType, cValue):
    col = 'white'
#   choose color depending on number of hits vs max nr of hits 
    if cType == 'nrHitsPerc':
        if cValue < 0.01:
            col = "#ee0000"
        elif cValue < 0.02:
            col = "#cd0000"
        elif cValue < 0.04:
            col = "#e31a1c"
        elif cValue < 0.07:
            col = "#f03b20"
        elif cValue < 0.10:
            col = "#fc4e2a"
        elif cValue < 0.15:
            col = "#ec7014"
        elif cValue < 0.20:
            col = "#fe9929"
        elif cValue < 0.25:
            col = "#fec44f"
        elif cValue < 0.30:
            col = "#fee391"
        elif cValue < 0.40:
            col = "#c7e9b4"
        elif cValue < 0.50:
            col = "#edf8b1"
        elif cValue < 0.60:
            col = "#7fcdbb"
        elif cValue < 0.70:
            col = "#41b6c4"
        elif cValue < 0.80:
            col = "#1d91c0"
        elif cValue < 0.85:
            col = "#a8ddb5"
        elif cValue < 0.90:
            col = "#74c476"
        elif cValue < 0.95:
            col = "#41ab5d"
        elif cValue < 0.97:
            col = "#238b45"
        elif cValue < 0.99:
            col = "#00cd00"
        elif cValue <= 10000.00:
            col = "#00ee00"                    
#   choose color depending on evalue of hit        
    elif cType == 'eValue':     
        if cValue > 1000.0:
            col = "#ee0000"
        elif cValue > 100.0:
            col = "#cd0000"
        elif cValue > 50.0:
            col = "#e31a1c"
        elif cValue > 25.0:
            col = "#f03b20"
        elif cValue > 10.0:
            col = "#fc4e2a"
        elif cValue > 5.00:
            col = "#ec7014"
        elif cValue > 2.5:
            col = "#fe9929"
        elif cValue > 1.0:
            col = "#fec44f"
        elif cValue > 0.5:
            col = "#fee391"
        elif cValue > 0.10:
            col = "#c7e9b4"
        elif cValue > 0.05:
            col = "#edf8b1"
        elif cValue > 0.01:
            col = "#7fcdbb"
        elif cValue > 1.0e-3:
            col = "#41b6c4"
        elif cValue > 1.0e-4:
            col = "#1d91c0"
        elif cValue > 1.0e-5:
            col = "#a8ddb5"
        elif cValue > 1.0e-10:
            col = "#74c476"
        elif cValue > 1.0e-15:
            col = "#41ab5d"
        elif cValue > 1.0e-20:
            col = "#238b45"
        elif cValue > 1.0e-25:
            col = "#00cd00"
        elif cValue >= 0:
            col = "#00ee00"             
        pass
    else:
        print " ERROR: getGColor: type of calculation does not exist! ["+cType+"]"    
    return col

def drawHitChart(_type, _brackets,_outFile):
    data = []
    if 'D:' in _type:
        __type = 'DOMAINS' 
    if 'K:' in _type:
        __type = 'KINGDOMS' 
    if 'P:' in _type:
        __type = 'PHYLA' 
    data.append([])
    data[0] += ['Distribution of hits by '+__type,'Nr of Hits','E-Value']
    #data[0] += cBs
    doSave = False 
    chartsB = _brackets                
    # go through all results:            
    for n in allTaxes.keys():
        newRow = []
        if _type in n:
            doSave = True
            newRow.append(n)
            for hit in allTaxes[n].hitEVs:
                hitDone = False
                for brMaxC in range(0,len(chartsB)):
                    if not hitDone and hit > chartsB[brMaxC][0]:
                        chartsB[brMaxC].append(hit)
                        hitDone = True                                                                                          
            for b in chartsB:
                newRow.append(len(b[1:]))
            data.append(newRow)                                        
    if doSave:                                      
        makeStackedBarChart(data, _outFile)                         

''' -------------------------------------------------------- '''
''' MAIN: STARTING HERE '''
''' -------------------------------------------------------- '''
'''    COMMAND LINE PARSING !
if len (sys.argv) < 2 or len(sys.argv) > 2:
    print '------------------- USAGE -----------------------'
    print '> python binarizeNCBItax <input file>             '
    print '> binarized NCBI tax files for use with tax mapper'
    print '> WARNING: execute only on gi_taxid_prot.dmp'
    print '>          AND gi_taxid_nucl.dmp'
    print '-------------------------------------------------'
    exit(-1)
'''
print "--> TaxMapper run started ! <-- "
import argparse
parser = argparse.ArgumentParser(description='Homology search result taxonomy mapper')
parser.add_argument('-I','--input',nargs='?', default='./input.xml',help='[DEF: ./input.xml] input file (BLAST or HMMER results file), if BLAST, MUST be in xml form')
parser.add_argument('-X','--tax',nargs='?',default='./ncbi_tax', help='[DEF: ./ncbi_tax] location of ncbi taxonomy FOLDER (names.dmp, nodes.dmp, gi_taxid_prot.bin and gi_taxid_nucl.bin files); Default: ./ncbi_tax; see NCBI taxonomy for details')
parser.add_argument('-O','--output',nargs='?',default='output.gv', help='[DEF: output.gv] output filename; will overwrite file if it exists')
parser.add_argument('-e','--eV',nargs='?',default=10.0,type=float, help='[DEF: 10.0] maximum eValue of result hit; all hits ABOVE this threshold will NOT be considered.')
parser.add_argument('-b','--bS',nargs='?',default=1.0, type=float, help='[DEF: 1.0] minimum bitscore of result hit; all hits BELOW this threshold will NOT be considered.')
parser.add_argument('-n','--nrHits',nargs='?',default=1, type=int, help='[DEF: 1] minimum number of hits of hits per taxa. If there are less hits, that taxa will NOT be displayed.')
parser.add_argument('-p','--percHits',nargs='?',default=0.0, type=float, help='[DEF: 0.0] minimum number of hits of hits per taxa [percenatage of max hits], decimal value [0.1 = 10 perc]')
parser.add_argument('-l','--minL',nargs='?',default=0, type=int, help='[DEF: 0] minimum alignment length of result hit; all hits BELOW this threshold will NOT be considered.')
parser.add_argument('-m','--maxL',nargs='?',default=1000000, type=int, help='[DEF: no limit] maximum alignment length of result hit; all hits ABOVE this threshold will NOT be considered.')
parser.add_argument('-R','--rank',nargs='?',default='SCPKD', help='[DEF: SCPKD] defines which taxonomical ranks to write (S = species, G = genus, F = family, O = order, C = clade, P = phylum, K = kingdom, D = domain; \nexample: KCS will draw species -> clade -> domain graph')
parser.add_argument('--excludeTaxa',nargs='?',default='D:Viruses,&D:N/D', help='[DEF: Viruses] defines excluded taxa (usually viruses due to strange problems associated with them and taxonomy; enter as "A,B,C,D..." (TAXIDS only)')
parser.add_argument('--borderC',nargs='?',type=bool, default=True, help='[DEF: True] True for output with borders around names colored to represent number of hits to appropriate taxa')
parser.add_argument('--nodeC',nargs='?',type=bool, default=True, help='[DEF: True] True for output with nodes colored to represent quality of hit(s) to appropriate taxa')
parser.add_argument('--borderEv',nargs='?',type=bool, default=False, help='[DEF: False] True for output with borders around names colored to represent quality of hit(s) to appropriate taxa')
parser.add_argument('--searchType',nargs='?', default='blastp', choices=['blastp','blastn','blastx','tblastx','psiblast','tblastn','deltablast','megablast','hmmer','jackhmmer', 'hhblits'], help='[DEF: blastp] Optional: enter blast type for automatic selection of taxonomy (nucl or prot)')
parser.add_argument('--type',nargs='?',default='auto', choices=['auto','prot','nucl'], help='[DEF: auto] prot for protein results of protein BLAST (blastp) or HMMER search, nucl for nucleotide blast (blastn)')
parser.add_argument('--autoDetect',nargs='?', default='Y', choices=['Y','N'], help='[DEF: Y] Y for automatic detection of input files (between blast/hmmer types and whether nucleotide or protein')
parser.add_argument('--doIdMapping',nargs='?', default='N', choices=['Y','N'], help='[DEF: N] Y for automatic linkage of various reference ids to gi (such as ref, pdb, gb ...) \n   NOTE: mapping can take a while (10 min+), especially for large number of ids')
parser.add_argument('--doIdMapFile',nargs='?', default='tax', help='[DEF: same as -X] (only if --doIdMapping = Y, FOLDER with uniprot id map file (idmapping_selected.tab)')
parser.add_argument('--doChartOutput',nargs='?', default='YDKP', help='[DEF: YDKP] if Y[DKP], produces charts with eValue distribution of results by Domain / Kingdom / Phylum')
parser.add_argument('--doTrimSName',nargs='?', default='Y', help='[DEF: Y] if Y, will trim first word of species name to one letter only (ex: Bacillus Subtilis -> B. Subtilis)')
parser.add_argument('--doTaxMerge',nargs='?', default='Y', help='[DEF: Y] if Y, will merge results with NCBI taxonomy')
parser.add_argument('--taxMergeFilter',nargs='?', default='D', help='[DEF: D] when merging with NCBI taxonomy, display only NCBI taxa where results are in same [D: domain, K: kingdom, P: phylum, C: class, O: order, F: family, G: genus], ex: DK: only within same domain and same kingdom as results')
parser.add_argument('--taxMergeDivisions',nargs='?', default='DKPC', help='[DEF: DKPC] which NCBI divisions are loaded for merging with NCBI taxonomy (avaliable: DK, DKP, DKPC, DKPCO')
parser.add_argument('--graphFontSizes',nargs='?', default='D16,K12,P10,C8', help='[DEF: D20,K18,P16,C14] font sizes for division type (D/K/P/C/F/G/S); if not entered, then default values')
parser.add_argument('--graphOutStats',nargs='?', default='D1,K1,P1,C1,F1,G1,S1', help='[DEF: D1,K1,P1,C1,F1,G1,S1] info written for division: 0: no # hits, best EV; 1: #hits, best EV in new line; 2 in same line; 3: #hits only, new line; 4: #hits only, same line')
parser.add_argument('--graphRankSep',nargs='?', default='6', help='[DEF: 6] distace between taxononmy ranks for output graph(s): set it to 6 for lots of results, 2-5 for less')
parser.add_argument('--HHblitsFilter',nargs='?', type = str, choices=['none','bacteria','archaea','fungi','protista','misc','metazoa','plants'], default='none', help='[DEF: no filter]: filter which taxa to load from HHblits result (and HHblits result only); following are supported: bacteria, archaea, fungi, metazoa, plants')
parser.add_argument('--HHblitsFilterTaxID',nargs='?', type = str, default='[]', help='[DEF: no filter]: filter which taxa to load from HHblits result (and HHblits result only); MUST be entered in following format: "[tax1,tax2,tax3...,taxn]", this is NCBI TAXID filter, correct taxa numbers must be entered')
parser.add_argument('--inFileTax',nargs='?', default='./taxFiles/',help='[DEF: taxFiles/]: location of prepared taxonomy graph files')
parser.add_argument('--doSpeciesOut',nargs='?', default='N',help='[DEF: N]: if Y, writes output in format: species -> hits (sequences)')
parser.add_argument('--doSpeciesOutF',nargs='?', default='N',help='[DEF: N]: if Y, writes output in format: species -> hits (sequences), full sequences')
parser.add_argument('--speciesOutBlastDBcmd',nargs='?', default='blastdbcmd',help='[DEF: blastdbcmd]: NCBI blast blastdbcmd command (required for --doSpeciesOut)')
parser.add_argument('--speciesOutNRDB',nargs='?', default='/NR_db',help='[DEF: /NR_db]: location of NCBI nr database (required for --doSpeciesOut)')


args = parser.parse_args()
   
''' ----- OPTIONS ----- '''
# SPECIES OUT PRETTY PRETTY OUTPUT PARAMETERSR
blastDBCmd = args.speciesOutBlastDBcmd
blastdbPath = args.speciesOutNRDB

# more drawing options
dofont = args.graphFontSizes
dofont = dofont.replace("'",'')
dofonts = dofont.split(',')
graphFonts = {}
# defaults
graphFonts['D'] = 16
graphFonts['K'] = 12
graphFonts['P'] = 10
graphFonts['C'] = 8
graphFonts['F'] = 6
graphFonts['G'] = 6
graphFonts['S'] = 6
graphFonts['U'] = 6
for a in dofonts:
    if 'D' in a: 
        graphFonts['D'] = int(a.replace('D',''))
    elif 'K' in a: 
        graphFonts['K'] = int(a.replace('K',''))
    elif 'P' in a: 
        graphFonts['P'] = int(a.replace('P',''))
    elif 'C' in a: 
        graphFonts['C'] = int(a.replace('C',''))
    elif 'F' in a: 
        graphFonts['F'] = int(a.replace('F',''))
    elif 'G' in a: 
        graphFonts['G'] = int(a.replace('G',''))
    elif 'S' in a: 
        graphFonts['S'] = int(a.replace('S',''))    
    elif 'U' in a: 
        graphFonts['S'] = int(a.replace('U',''))    
        
doouttype = args.graphOutStats
doouttype = doouttype.replace("'",'')
doouttypes = doouttype.split(',')
graphOutType = {}
graphOutType['D'] = 1
graphOutType['K'] = 1
graphOutType['P'] = 1
graphOutType['C'] = 1
graphOutType['F'] = 1
graphOutType['G'] = 1
graphOutType['S'] = 1
for a in doouttypes:
    if 'D' in a: 
        graphOutType['D'] = int(a.replace('D',''))
    elif 'K' in a: 
        graphOutType['K'] = int(a.replace('K',''))
    elif 'P' in a: 
        graphOutType['P'] = int(a.replace('P',''))
    elif 'C' in a: 
        graphOutType['C'] = int(a.replace('C',''))
    elif 'F' in a: 
        graphOutType['F'] = int(a.replace('F',''))
    elif 'G' in a: 
        graphOutType['G'] = int(a.replace('G',''))
    elif 'S' in a: 
        graphOutType['S'] = int(a.replace('S',''))

''' PARSING OPTIONS '''
parseEVmax = 1.0e+1        # max EV       (don't take those above it)
parseEVmax = args.eV
parseBSmin = args.bS       # min BitScore (don't take those above it)
parseMinHitNr = args.nrHits     # min nr of hits required for species to be considered
parseMinHitNrPerc = args.percHits   # min % of hits required for species to be considered

ltMin = args.minL            # min length of hit
ltMax = args.maxL      # max length of hit
minQL = 1  # min Query Length

excludedTaxa = []
et = args.excludeTaxa.replace('[','').replace(']','').replace(';','').replace('"','').replace("'",'')
if ',' not in et:
    excludedTaxa.append(et)
else: 
    for e in et.split(','):
        excludedTaxa.append(e)            

#print excludedTaxa

''' ERROR CHECKING:  
   -> INPUT error check
'''
try:
    with open(args.input):
        pass       
except IOError:
    print " --> ERROR: input file '"+args.input+"' does not exist! <-- "
    print "   --> -I or --input should point to BLAST, HMMER or HHblits results file "
    print "   --> ABORTING RUN! <--"
    exit (1)
''' DONE WITH INPUT ERROR CHECK '''

#d print 'detecting file'
''' auto detect file type'''
if args.autoDetect == 'Y':
    ft = 'unknown'
    #d print args.input
    ft = determineSearchFileType(args.input)
    #d print ft
    if not ft == 'unknown':
        args.searchType = ft 
    elif ft == 'empty':
        print ' --> WARNING: file is empty, ending run <-- '
        exit(0)
    elif ft == 'unknown':
        print ' --> ERROR: UNRECOGNIZED FILE TYPE !!! <-- '
        print '   --> make sure input is proper FASTA, BLAST, HMMER or HHblits results file'
        print '   --> skipping file (ending run!) <-- '
        exit(0)
#d print 'preparing paths'
        
''' search type if not autodetect '''    
if not args.searchType == 'none':
    if args.searchType == 'blastp': 
        isProtein = True
    elif args.searchType == 'blastn':
        isProtein = False
    elif args.searchType == 'megablast':
        isProtein = False
    elif args.searchType == 'blastx':
        isProtein = True        
    elif args.searchType == 'tblastn':
        isProtein = True
    elif args.searchType == 'tblastx':
        isProtein = True
    elif args.searchType == 'psiblast':
        isProtein = True
    elif args.searchType == 'deltablast':
        isProtein = True
    elif args.searchType == 'hmmer':
        isProtein = True        
    elif args.searchType == 'jackhmmer':
        isProtein = True    
    elif args.searchType == 'hhblits':
        isProtein = True  
    elif args.searchType == 'fastap':
        isProtein = True  
    elif args.searchType == 'fastan':
        isProtein = False                  
            
''' if forced (prot / nucl type): '''
if args.type == 'prot':
    isProtein = True
elif args.type == 'nucl': 
    isProtein = False

if args.doIdMapFile == 'tax':
    args.doIdMapFile = args.tax        
idMapFile = args.doIdMapFile+'/idmapping_selected.tab' # PATH to id_map file
if str.lower(args.doIdMapping) == 'y':
    try:
        with open(idMapFile):
            pass       
    except IOError:
        print " --> ERROR: idMap file '"+idMapFile+"' does not exist! <-- "
        print "   --> doIdMapFile should point to FOLDER of idmapping_selected.tab"
        print "   --> ABORTING RUN! " 
        exit (1)

''' PATHS '''
pathNm = args.tax+'/names.dmp'                         # PATH to names.dmp
pathNd = args.tax+'/nodes.dmp'                         # PATH to nodes.dmp
ncbiProtTaxBin = args.tax+'/gi_taxid_prot.bin'
ncbiNuclTaxBin = args.tax+'/gi_taxid_nucl.bin'

if isProtein:
    ncbiTaxPath = ncbiProtTaxBin
else:
    ncbiTaxPath = ncbiNuclTaxBin
    
''' ERROR CHECKING: ncbiTaxPath should point to proper file! '''
try:
    with open(ncbiTaxPath):
        pass       
except IOError:
    print " --> ERROR: input file '"+ncbiTaxPath+"' does not exist! <-- "
    print "   --> -X / --tax must point to FOLDER with '/gi_taxid_prot.bin' and/or '/gi_taxid_nucl.bin'"
    print "   --> ABORTING RUN! <--"
    exit (1)
''' DONE WITH ERROR CHECK '''   

# INPUT
inputFilePath = args.input #input
# OUTPUT
outFile = args.output #output

''' OUTPUT OPTIONS '''
doSpeciesOut = False
doSpeciesOutF = False
if args.doSpeciesOut == 'Y':
    doSpeciesOut = True 
if args.doSpeciesOutF == 'Y':
    doSpeciesOutF = True 

''' DRAWING OPTIONS '''
#drawRanks ='SGFOCPKD' # which tax ranks to draw? [S = species, G = genus...]
drawRanks = args.rank        # which tax ranks to draw? [S = species, G = genus...]
displayRanks = args.rank
if 'S' not in drawRanks:
    drawRanks+='S'

doColorBorders = args.borderC   # should we    print '    --> gis total:',len(allGIs.keys()) color borders?
doColorNodes = args.nodeC     # should we color nodes?
doColorBordersWithEV = args.borderEv
doTrimSName = args.doTrimSName
# parse arguments, decide on drawing parameters
drawCharts = False
chartsDomain = False
chartsKingdom = False
chartsPhylum = False
# brackets (>= then)
chartsBrackets = [[1.0e+3],[1.0e+2],[1.0e+1],[1.0e+0],[1.0e-1],[1.0e-2],[1.0e-5],[1.0e-10],[1.0e-20],[1.0e-30],[1.0e-40],[1.0e-50],[0.0]]

if 'Y' in args.doChartOutput: 
    drawCharts = True
if 'D' in args.doChartOutput:
    chartsDomain = True
if 'K' in args.doChartOutput:
    chartsKingdom = True
if 'P' in args.doChartOutput:
    chartsPhylum = True

''' ---- VARS ----------------------------------------------- '''
allGIs = {}
# TAX IDs
addedTaxIDs = {}
allTaxIDs = {}

# list of all references in ALL HITS
allHitsParsedHdrs = []
allHitsParsedEVs = []

''' PASS 1/d: PARSE FASTA RESULT (in case of plain ol' fasta input)
    no filtering can be done, so we just chop out IDs and send them forward!
'''
if 'fasta' in args.searchType:
    print " --> 1) PARSING FASTA INPUT <--"
    print " -->    opening and parsing "+inputFilePath 
    fastas = loadAllFastas(inputFilePath)
    for f in fastas:
        allHitGIs = extractAllGIs(getFaHeader(f))
        fHdr = getFaHeader(f)
        eV = 0
        if '__eV{' in fHdr and '}Ve__' in fHdr:
            eV = fHdr[fHdr.find('__eV{')+5:]
            #print eV
            eV = eV[:eV.find('}Ve__')]
            eV = float(eV)
            #print eV
            #exit()
        if eV <= parseEVmax:        
            ph = parseHeader(fHdr, True)
            allHitsParsedHdrs.append(ph)          
            for gi in allHitGIs:
                if not gi in allGIs.keys():
                    allGIs[gi] = []
                    # if we have modified fasta (produced by DHpipeline), map ev from it: 
                    if isMFAHdr(fHdr):
                        allGIs[gi].append(getEvFromMFaHdr(fHdr))     
                    else: allGIs[gi].append(0.00001)                
                else:
                    if isMFAHdr(fHdr):
                        allGIs[gi].append(getEvFromMFaHdr(fHdr)) 
                    else: allGIs[gi].append(0.00001)                
    #d print allGIs
''' PASS 1/c: PARSE HHblits RESULT (in case of HHblits)
    and filter results according to optons (see above) '''
if 'hhblits' in args.searchType:
    print " --> 1) PARSING HHBLITS INPUT <--"
    print " -->    opening and parsing "+inputFilePath   
    hhballres = parseHHResults(inputFilePath)
    print " -->    parsed "+inputFilePath
    #print args.HHblitsFilter
    nameFilter = args.HHblitsFilter
    print nameFilter
    taxFilter = argToListI(args.HHblitsFilterTaxID)
    print "    --> got ",len(hhballres.hhResults),"results"
    print " -->    filtering "
    hhballres = filterHHblits(hhballres, ncbiProtTaxBin, pathNm, pathNd, targetTaxInc=[taxFilter], targetTaxExc=[], targetTaxNames=nameFilter)
    print "    --> got ",len(hhballres.hhResults),"results after filtering"
    cntaF = 0
    for r in hhballres.hhResults:            
        qLt = r.qHMMEnd - r.qHMMStart
        if r.eV <= parseEVmax and r.rawScore >= parseBSmin and r.nrMatches <= ltMax and r.nrMatches >= ltMin \
        and qLt >= minQL:        
            cntaF += 1                
            allHitGIs = r.gis
            allHitsParsedEVs.append(r.eV)
            for gi in allHitGIs:
                if not gi in allGIs.keys():
                    ph = {}
                    ph['gi'] = r.gis
                    allHitsParsedHdrs.append(ph)
                    
                    allGIs[gi] = []
                    allGIs[gi].append(r.eV)
                else:
                    allGIs[gi].append(r.eV)  
    print "       --> got ",cntaF,"results after post-filtering [",len(allGIs.keys()),"] Taxa"
''' PASS 1/a: PARSE BLAST RESULT (in case of BLAST)
    and filter results according to optons (see above) '''
    
if 'blast' in args.searchType:     
    print " --> 1) PARSING BLAST INPUT <--"
    print " -->    opening "+inputFilePath     
    hsps_cnt = 0
    maxPsiIteration = 0
    with open(inputFilePath) as resultHandle:
            print " -->    opened "+inputFilePath
    # parse it
            blastRecords = NCBIXML.parse(resultHandle)
            cnt = 0
            psiIteration = 0
            # find max iterations:    
            if args.searchType == 'psiblast':        
                for br in blastRecords:                
                    maxPsiIteration+=1
                    
    with open(inputFilePath) as resultHandle:       
            blastRecords = NCBIXML.parse(resultHandle)                                    
            for blastRecord in blastRecords:
                if args.searchType == 'psiblast':
                    psiIteration+=1         
                for alignment in blastRecord.alignments:
                    for hsp in alignment.hsps:
                        hsps_cnt +=1
                        if hsp.expect <= parseEVmax and hsp.score >= parseBSmin and alignment.length <= ltMax and alignment.length >= ltMin\
                        and len(hsp.sbjct) >= minQL and psiIteration == maxPsiIteration:
                            cnt += 1
                            #d print '**** HIT # '+str(cnt)+'****'
                            #d print 'sequence:', alignment.title
                            #d print 'iteration:',psiIteration,'/',maxPsiIteration                            
                            
                            ph = parseHeader(alignment.title, True)
                            allHitsParsedHdrs.append(ph)                            
                            allHitsParsedEVs.append(hsp.expect)       
                            
                            allHitGIs = extractAllGIs(alignment.title)                                                        
                            for gi in allHitGIs:
                                if not gi in allGIs.keys():
                                    allGIs[gi] = []
                                    allGIs[gi].append(hsp.expect)
                                else:
                                    allGIs[gi].append(hsp.expect)                                                       
                        #d: print 'e value:', hsp.expect
                        #d: print 'bitscore:',hsp.score
    print "   --> total hits found:",hsps_cnt
    print " --> 1) DONE!         <--"
    
''' PASS 1/b: PARSE HMMER (in case of HMMER)
    and filter results according to optons (see above) '''
if 'hmmer' in args.searchType:     
    print " --> 1) PARSING HMMER INPUT <--"
    print " -->    opening "+inputFilePath         
    hmmAllResults = parseHMMResults(inputFilePath)
    if len(hmmAllResults.qry) < 1: 
        print '   --> ERROR: no HMMER results detected! <--- '
        exit (1)
    # go over all results in last iteration (only one in case of non JACK)             
    hmmQR = hmmAllResults.qry[0]
    hmmQI = hmmQR.getLastIteration()
    cnt = 0        
    for h in hmmQI.hmmTable:
        if h.fullSeqEV <= parseEVmax and h.fullSeqBS >= parseBSmin:
            cnt += 1
#            print '**** HIT # '+str(cnt)+'****'
#            print 'sequence:','eV:', h.fullSeqEV, h.seqID,
            ph = parseHeader(h.seqID+h.description, True)
            allHitsParsedHdrs.append(ph)
            allHitsParsedEVs.append(h.fullSeqEV)            
            allHitGIs = extractAllGIs(h.seqID)
            for gi in allHitGIs:
                if not gi in allGIs.keys():
                    allGIs[gi] = []
                    allGIs[gi].append(h.fullSeqEV)
                else:
                    allGIs[gi].append(h.fullSeqEV)     
    print "   --> total hits found: ",cnt       
    print " --> 1) DONE!         <--"    

''' ----- DONE WITH PARSING RESULT(S) ----- '''
#d for k in allGIs.keys():
#d    print 'gi:',k, allGIs[k]

''' determine how many GIs we found
compared to how many headers with GI there are '''
cWOGI = 0
cWGI = 0
for i in range(0,len(allHitsParsedHdrs)):
    if 'gi' not in allHitsParsedHdrs[i].keys():
        cWOGI +=1
    else: 
        cWGI +=1
    #d print allHitsParsedHdrs[i], '->',allHitsParsedEVs[i]
if cWOGI+cWGI == 0:
    print '    --> NO RESULTS FOUND AFTER FILTERING!'
    print '    --> ABORTING RUN!'
    exit (0)
percMM = float(cWOGI) / float(cWOGI+cWGI)*100 
print '  --> found',cWOGI,'/',cWOGI+cWGI,' hits without GI reference ({0:.1f}'.format(percMM),'%)'
if percMM > 50.0: 
    print '  --> --doIdMapping command should be used to map missing hits to NCBI taxonomy'
elif percMM > 10.0: 
    print '  --> consider using --doIdMapping command to map missing hits to NCBI taxonomy'

''' do IDMAPPING if required '''
if args.doIdMapping == 'y' or args.doIdMapping == 'Y':
    gisListList = extractFromIdMap(allHitsParsedHdrs,idMapFile,out='taxid')    
    for c in range(0,len(gisListList)):
        for g in gisListList[c]:
            g = int(g)
            d = TaxIdData()
            d.hitEVs.append(allHitsParsedEVs[c])        
            if not g in addedTaxIDs.keys():
                addedTaxIDs[g] = d
            else:
                for ev in d.hitEVs:
                    addedTaxIDs[g].hitEVs.append(ev)        
#d print 'filled added TaxID LIST: '
#d for k in addedTaxIDs.keys():
#d     print 'taxID:',k, addedTaxIDs[k].hitEVs                  
''' --------------------------------------------
 PASS # 2
 GO OVER ALL GIs, for each of them, extract TAXID
 and put it into dictionary: taxID -> taxIDdata
    ------------------------------------------- ''' 
print ' --> 2) LOADING taxIDs for GIs <--' 
print '    --> gis total:',len(allGIs.keys())
    
cntFnd = 0
cntNFnd = 0
print '    --> ncbi GI->Tax filename: '+ncbiTaxPath
for gi in allGIs.keys():
    taxID = loadTaxID(ncbiTaxPath, gi)
    if taxID != 0:
        cntFnd +=1
        d = TaxIdData()
        for ev in allGIs[gi]:
            d.hitEVs.append(float(ev))
            d.hitGIs.append(int(gi))                                        
        if not taxID in allTaxIDs.keys():
            allTaxIDs[taxID] = d
        else:
            for ev in d.hitEVs:
                allTaxIDs[taxID].hitEVs.append(ev)
                allTaxIDs[taxID].hitGIs.append(gi)
    else:
        cntNFnd +=1
        
print '    --> TAX ids: ','found: ',cntFnd,'not found:',cntNFnd

#d for k in allTaxIDs.keys():
#d    print 'taxID:',k, allTaxIDs[k].hitEVs       
    
# merge addedTaxIDs and allTaxIDs
for ka in addedTaxIDs.keys():
    # we have new taxID
    if ka not in allTaxIDs.keys():
        allTaxIDs[ka] = addedTaxIDs[ka]
    # we have repeated one, check for differences
    if ka in allTaxIDs.keys():
        newkaEvs = addedTaxIDs[ka].hitEVs
        for newkaEvsHit in newkaEvs:
            if newkaEvsHit not in allTaxIDs[ka].hitEVs:
                allTaxIDs[ka].hitEVs.append(newkaEvsHit)
        newkaGis = addedTaxIDs[ka].hitGIs
        for newkaGisHit in newkaGis:
            if newkaGisHit not in allTaxIDs[ka].hitGIs:
                allTaxIDs[ka].hitGIs.append(newkaGisHit)                
                
#d print 'merged tax Lists'                
#d for k in allTaxIDs.keys():
#d     print 'taxID:',k, allTaxIDs[k].hitEVs  

''' now filter this one, throw out all those without enough hits '''
maxHits = -1
for k in allTaxIDs.keys():
    l = len(allTaxIDs[k].hitEVs)
    if  l > maxHits:
        maxHits = l         

#d print 'maxhits:',maxHits
allTaxIDsClone = {}
for k in allTaxIDs.keys():    
    l = len(allTaxIDs[k].hitEVs)
    nrHitsPerc = float(float(l)/float(maxHits))
    #d print 'nrHits:',l,'nrHitsP:',nrHitsPerc,'maxH:','miHN:',parseMinHitNr,'miHNP',parseMinHitNrPerc    
    if l >= parseMinHitNr and nrHitsPerc >= parseMinHitNrPerc:
        allTaxIDsClone[int(k)] = copy.deepcopy(allTaxIDs[k])
    else:
        pass 
        #print 'filtered ',k,'out!'

allTaxIDs = {}
allTaxIDs = copy.deepcopy(allTaxIDsClone)
allTaxIDsClone = {}
#allTaxID
print ' --> 2) DONE!                  <--'

''' ------------------------------------------                
print " ----> DEBUG, LISTING ALL TAXIDS <--- "
for k in allTaxIDs.keys():
    allTaxIDs[k].calcStats()
    print 'taxid:',k, 'avg eV:',allTaxIDs[k].avgEV, 'best eV',allTaxIDs[k].bestEV, 'nr Hits:',allTaxIDs[k].nrHits,'hits:', str(allTaxIDs[k].hitEVs)
    
print " TOTAL taxids: ",len(allTaxIDs.keys())    
print " GI with existing taxonomy: ",cntFnd
print " GI with NO NCBI taxonomy: ",cntNFnd
 ---------------------------------------- '''
        
''' ------- DONE WITH PASS 2 ------ '''       

''' -------- PASS 3: TAX HASHING ---------------------------------
------------------------------------------------------------------ 
make tax stuff hash (species -> genus -> ...)
paired into tuples such as (species -> genus; genus -> family...)
clear out duplicates and populate with data such as bestEV... '''
        
print " --> 3) generating taxonomy graph <-- "
# cleanup stuff
allGIs = {}
# new hash
allTaxes = {}
print '    --> I) loading TAXONOMY database    '
initTax(pathNm,pathNd)
print '    --> I) done                         '
print '    --> II) generating hash             '
''' go over all TAX IDs, extract taxonomy for each and hash it '''
''' following block generates hash from tuples '''  
startTaxes = []
for k in allTaxIDs.keys():
    m = getTaxFromTaxID(k)    
    l = len(allTaxIDs[k].hitEVs)
#d    print doTrimSName
    rl = m.getTax9EwDummy(drawRanks,doTrimSName,excludedTaxa=excludedTaxa)
    if not rl == None:
#d    print rl.toStrList()
        h = rl.toTupleList()
    #for l in h:
    #     print l,
        taxLvlC = 0    
        for n in h:
            taxLvlC +=1
        # note: n[0] = start, n[1] = target
            if n[0] in allTaxes.keys():            
                pass
            else:
                allTaxes[n[0]] = TaxIdDataDirected(n[0],n[1])
                allTaxes[n[0]].rank = rl.records[taxLvlC-1].rank
                allTaxes[n[0]].key = str(n[0])
                if allTaxes[n[0]].startTax == allTaxes[n[0]].targetTax:
                    allTaxes[n[0]].targetTax = ''
                if taxLvlC == 1: 
                    allTaxes[n[0]].hitEVs = allTaxIDs[k].hitEVs
                    allTaxes[n[0]].hitGIs = allTaxIDs[k].hitGIs                
                    startTaxes.append(n[0])
#d print str(startTaxes)
print '    --> II) done             '
''' following block loads all hits (all hit eVS) in tree type fashion as following: 
  all species hits are added to genus of that species
  all genus hits are added to family of that genus
  all family hits are added to order of that genus
  ...
  to avoid recursion, this is done with pseudo-recursive while loop(s)
'''
print '    --> III) populating graph '
for st in startTaxes:
    doGo = True
    nS = allTaxes[st].startTax
    nSo = allTaxes[st].startTax
    nT = allTaxes[st].targetTax
    eVsToAdd = allTaxes[st].hitEVs
    GIsToAdd = allTaxes[st].hitGIs
    if not st == 'N/D' and not nT == 'N/D':            
        while doGo:
            if allTaxes[nT].startTax in allTaxes.keys():
                nS = allTaxes[nT].startTax
                if allTaxes[nT].targetTax in allTaxes.keys():
                    allTaxes[nS].hitEVs += eVsToAdd
                    allTaxes[nS].hitGIs += GIsToAdd
                    nSo = nS                 
                    nT = allTaxes[nT].targetTax                
                    #                print nS,allTaxes[nS].hitEVs,'->',nT,allTaxes[nT].hitEVs
                else: 
                    doGo = False
                    allTaxes[nT].hitEVs += eVsToAdd
                    allTaxes[nT].hitGIs += GIsToAdd
            else: 
                doGo = False
                allTaxes[nS].hitEVs += eVsToAdd     
                allTaxes[nT].hitEVs += eVsToAdd
                allTaxes[nS].hitGIs += GIsToAdd
                allTaxes[nT].hitGIs += GIsToAdd

print '    --> III) done '

''' following block calculates stats of all entries '''
print '    --> IV) calculating statistics '
maxHitsTotal = 1
for k in allTaxes.keys():
#d    print k, str(allTaxes[k].hitEVs)
    allTaxes[k].calcStats()
    if allTaxes[k].nrHits > maxHitsTotal:
        maxHitsTotal = allTaxes[k].nrHits
    if allTaxes[k].nrHits > 0:
        ndName = allTaxes[k].startTax
        # tweak output so it is pretty:
        if 'D:' in ndName: 
            oo = 'D'
        elif 'K:' in ndName: 
            oo = 'K'
        elif 'P:' in ndName: 
            oo = 'P'
        elif 'C:' in ndName: 
            oo = 'C'
        elif 'O:' in ndName: 
            oo = 'O'
        elif 'F:' in ndName: 
            oo = 'O'
        elif 'G:' in ndName: 
            oo = 'O'
        elif ':' in ndName: 
            oo = 'O'                        
        tm = []
        for k2 in allTaxes.keys():
            if allTaxes[k].startTax == allTaxes[k2].targetTax:
                tm.append(k2)
                
        allTaxes[k].startTax += '\\n [H:'+str(allTaxes[k].nrHits)+', eV:{:.1e}'.format(allTaxes[k].bestEV)+']'
        for k2m in tm:
            allTaxes[k2m].targetTax = str(allTaxes[k].startTax)
print '    --> IV) done '

for k in allTaxes.keys():
    allTaxes[k].calcStats(maxHits)
print " --> 3) DONE!                     <-- "

''' --------------- PASS 4 -------------
actually write down .gv (graphviz code (dot language)) 
for graph '''
         
print ' --> 4) GENERATING OUTPUT <--'
if '.gv' not in outFile:
    outFile = outFile + '.gv'
outFile = outFile.replace('.gv','_graph.gv')
print '   --> graph output file:',outFile
           
with open (outFile,'w') as out:
    # decide which is last taxa (used for removing orphan nodes)
    lastT = ''
    pranka = 'SGFOCPKD'
    for a in pranka:
        if a in displayRanks:
            lastT = a
            break
    #print lastT
    # graph header and options
    out.write('digraph TaxG'+'\n')
    p = '    '
    out.write(p+'{'+'\n')
    out.write(p+'root="N/D";'+'\n')
    out.write(p+'layout=twopi;'+'\n')
    out.write(p+'edge [style=dashed dir=back];'+'\n')
    out.write(p+'node [shape=plaintext];'+'\n')
    out.write(p+'ranksep='+args.graphRankSep+';'+'\n')
    out.write(p+'label = "TAXONOMY OF RESULTS";'+'\n')
    out.write(p+'center = 1;'+'\n')
    # add nodes
    out.write('# --- NODES FOLLOW --- '+'\n')    
    ''' NODES '''
    out.write(p+'"N/D"'+'[shape=circle,label="",width=0.1,height=0.1]'+'\n');
    for n in allTaxes.keys():
        # now decide on what to actually write down
        # and don't write orphan node (ones which have no proper ending)
        pt = getTaxDirAllChildren(allTaxes, n)
        #print '--> node: ',n, 'tax: ', pt                
        doWriteNode = False            
        for d in displayRanks:     
            if d+':' in n:
                doWriteNode = True
                if '&'+lastT+':' in pt[0]:
                #if '&' in pt[0] and ':' in pt[0] and lastT in pt[0]:
                    doWriteNode = False
                    
        
        if doWriteNode:                                      
            ''' color border nodes depending on # hits '''
            label = allTaxes[n].startTax
            st = allTaxes[n].startTax
            if not allTaxes[n].startTax.find('\\n') == -1:
                st = allTaxes[n].startTax[0:allTaxes[n].startTax.find('\\n')]
            borderCol = 'black'
            if doColorBorders:
                borderCol = getGColor('nrHitsPerc',allTaxes[n].nrHitsPerc)
            if doColorBordersWithEV:
                borderCol   = getGColor('eValue',allTaxes[n].bestEV)
            ''' color node depending on eV avg '''
            nodeCol = 'white'
            if doColorNodes:
                nodeCol   = getGColor('eValue',allTaxes[n].bestEV)                    
                ''' construct tooltip '''                     
            # NOTE: __ ALL HITS WILL NOT BE WRITTEN (IT CRASHES STUFF FOR LARGE NRs OF HITS         
            tool = 'HITS: '+str(allTaxes[n].nrHits)+' '+', best eV:{:.1e}'.format(allTaxes[n].bestEV)+', avg eV:{:.1e}'.format(allTaxes[n].avgEV)+', worst eV:{:.1e}'.format(allTaxes[n].worstEV)
            ''' write code for node '''
            add = '['
            # basic stuff (shape, style, colors...)
            add +='shape=box'
            add +=' width=0.01 height=0.01'
            add +=' style="rounded,filled,bold"' 
            add +=' color="'+borderCol+'" '
            add +=' fillcolor="'+nodeCol+'"'
            add +=' tooltip="'+tool+'"'
            ndName = allTaxes[n].startTax
            #        print ndName
            if 'D:' in ndName:                 
                add = add + ' shape=circle fontsize='+str(graphFonts['D'])+' '                                
            elif 'K:' in ndName:
                add = add + ' shape=ellipse fontsize='+str(graphFonts['K'])+' '                                                                                          
            elif 'P:' in ndName:
                add = add + ' fontsize='+str(graphFonts['P'])+' '
            elif 'C:' in ndName:
                add = add + ' fontsize='+str(graphFonts['C'])+' '
            elif 'O:' in ndName:
                add = add + ' fontsize='+str(graphFonts['O'])+' '
            elif 'F:' in ndName:
                add = add + ' fontsize='+str(graphFonts['F'])+' '     
            elif 'G:' in ndName:
                add = add + ' fontsize='+str(graphFonts['G'])+' '     
            elif 'S:' in ndName:
                add = add + ' fontsize='+str(graphFonts['S'])+' '           
            ## OUTPUT TYPE ##
            label = ndName
            oo = ''
            if 'D:' in ndName:                 
                oo = 'D'                            
            elif 'K:' in ndName:                 
                oo = 'K'                                        
            elif 'P:' in ndName:                 
                oo = 'P'                            
            elif 'C:' in ndName:                 
                oo = 'C'                                
            elif 'O:' in ndName:                 
                oo = 'O'                                
            elif 'F:' in ndName:                 
                oo = 'F'                                
            elif 'G:' in ndName:                 
                oo = 'G'                                
            elif 'S:' in ndName:                 
                oo = 'S'                                        
            #print ndName, oo, graphOutType[oo],add        
            if graphOutType[oo] == 1:
                add = add+' label="'+label[0:]+'"'
            elif graphOutType[oo] == 2:
                add = add+' label="'+label[0:].replace('\\n','')+'"'
            elif graphOutType[oo] == 3:
                add = add+' label="'+label[0:label.find(', eV')]+']"'
            elif graphOutType[oo] == 4:
                add = add+' label="'+label[0:label.find(', eV')].replace('\\n','')+']"'        
            else:                     
                add = add+' label="'+label[0:label.find('\\n')]+'"'            
            add +=']'                 
            #d print add       
            if '&' in allTaxes[n].startTax:            
                out.write(p+'"'+st+'"'+'[shape=circle,label="",width=0.05,height=0.05]'+'\n');
    #elif '&' in allTaxes[n].targetTax:            
    #   out.write(p+'"'+st+'"'+'[shape=circle,label="",width=0.05,height=0.05]'+'\n');
            else:
        #d print 'writing stuff...'+ p+'"'+st+'"'+add+'\n'
                out.write(p+'"'+st+'"'+add+'\n');
                                
    ''' EDGES '''
    out.write('# --- EDGES FOLLOW --- '+'\n')             
    for n in allTaxes.keys():        
        pt = getTaxDirAllChildren(allTaxes, n)
        #print '--> node: ',n, 'tax: ', pt                
        doWriteEdge = False            
        for d in displayRanks:     
            if d+':' in n:
                doWriteEdge = True
                #if '&' in pt[0] and ':' in pt[0] and lastT in pt[0]:
                if '&'+lastT+':' in pt[0]:
                    doWriteEdge = False
                        
        if not allTaxes[n].targetTax == 'N/D/H/>J' and doWriteEdge:
            st = allTaxes[n].startTax
            tt = allTaxes[n].targetTax
            if not allTaxes[n].startTax.find('\\n') == -1:
                st = allTaxes[n].startTax[0:allTaxes[n].startTax.find('\\n')]
            if not allTaxes[n].targetTax.find('\\n') == -1:
                tt = allTaxes[n].targetTax[0:allTaxes[n].targetTax.find('\\n')]
            
            add = '[' 
            if 'Eukaryota' in findTaxDirChild(allTaxes, n.strip(),'D:'):
                add = add + 'color = "#ff0000" '                                                        
            elif 'Archaea' in findTaxDirChild(allTaxes,  n.strip(),'D:'):
                add = add + 'color = "#00ff00" '                                                        
            elif 'Bacteria' in findTaxDirChild(allTaxes,  n.strip(),'D:'):
                add = add + 'color = "#0000ff" ' 
            add +=']'
            out.write(p+'"'+st+'"'+' -> '+'"'+tt+'"'+add+';'+'\n')
    # END OF MAIN GRAPH
    p = ''
    out.write(p+'}'+'\n')    
        
''' if required, write out text file with following: 
species [nr hits]
  -> taxonomy
  -> gi -> eV
  '''
if doSpeciesOutF:
    with open (outFile.replace('.gv','_out_full.txt'),'w') as out:
        out.write('<TR> ---- TAXONOMY OF RESULTS ----'+'\n')
        allTaxes2 = allTaxes.values()
        allTaxes2 = sorted(allTaxes2, key=lambda x: x.nrHits, reverse=True)                
        for n in allTaxes2:
            if 'S:' in n.startTax:
                #C:Enoplea\n [H:2, eV:1.1e-05]
                out.write('<HS> SPECIES: '+n.startTax[2:n.startTax.find('\\n')] +'\n')     
                taxWO = ''
                for i in grabTaxDataFHash(allTaxes, n.key.strip()):
                    #print i
                    if '&' not in i and 'N/D' not in i and 'N/F' not in i:
                        taxWO += i+' -> '
                taxWO=taxWO[:-4]                              
                out.write('<HT> TAXONOMY: '+taxWO+'\n')
                out.write('<HS> STATS: '+'HITS: '+str(n.nrHits)+' '+', best eV:{:.1e}'.format(n.bestEV)+', avg eV:{:.1e}'.format(n.avgEV)+', worst eV:{:.1e}'.format(n.worstEV)+'\n')
                hts = []                
                for a in range(0,len(n.hitGIs)):
                    hts.append([n.hitGIs[a],n.hitEVs[a]])
                hts = sorted(hts,key=lambda x:x[1])
                cntr = 0
                for h in hts:
                    cntr+=1
                    out.write('<HH>   HIT: '+str(cntr)+', gi:'+str(h[0])+', eV:'+str(h[1])+'\n')
                    fasta = fetchFastaFromBlastDB(h[0],blastDBCmd,blastdbPath)
                    out.write(fasta+'\n')                                        
                    
if doSpeciesOut:
    with open (outFile.replace('.gv','_out_s.txt'),'w') as out:
        out.write('<TR> ---- TAXONOMY OF RESULTS ----'+'\n')
        allTaxes2 = allTaxes.values()
        allTaxes2 = sorted(allTaxes2, key=lambda x: x.nrHits, reverse=True)                
        for n in allTaxes2:
            if 'S:' in n.startTax:
                #C:Enoplea\n [H:2, eV:1.1e-05]
                out.write('<HS> SPECIES: '+n.startTax[2:n.startTax.find('\\n')] +'\n')     
                taxWO = ''
                for i in grabTaxDataFHash(allTaxes, n.key.strip()):
                    #print i
                    if '&' not in i and 'N/D' not in i and 'N/F' not in i:
                        taxWO += i+' -> '
                taxWO=taxWO[:-4]                              
                out.write('<HT> TAXONOMY: '+taxWO+'\n')
                out.write('<HS> STATS: '+'HITS: '+str(n.nrHits)+' '+', best eV:{:.1e}'.format(n.bestEV)+', avg eV:{:.1e}'.format(n.avgEV)+', worst eV:{:.1e}'.format(n.worstEV)+'\n')
                hts = []                
                for a in range(0,len(n.hitGIs)):
                    hts.append([n.hitGIs[a],n.hitEVs[a]])
                hts = sorted(hts,key=lambda x:x[1])
                out.write('<HH> HITS: ')
                for h in hts:
                    out.write('gi:'+str(h[0])+',eV:'+str(h[1])+'\n')
                    
                #out.write('<HH> GIS: '+str(n.hitGIs)+'\n')
                #out.write('<HH> EVS: '+str(n.hitEVs)+'\n')
''' now do the merging with prepared ncbi tax gv file '''
#parser.add_argument('--doTaxMerge',nargs='?', default='Y', help='[DEF: Y] if Y, will merge results with NCBI taxonomy')
#parser.add_argument('--taxMergeFilter',nargs='?', default='D', help='[DEF: D] when merging with NCBI taxonomy, display only NCBI taxa where results are in same [D: domain, K: kingdom, P: phylum, C: class, O: order, F: family, G: genus], ex: DK: only within same domain and same kingdom as results')
#parser.add_argument('--taxMergeDivisions',nargs='?', default='DKPC', help='[DEF: DKPC] which NCBI divisions are loaded for merging with NCBI taxonomy (avaliable: DK, DKP, DKPC, DKPCO')
    
mergeNCBITax = False
if args.doTaxMerge == 'Y':
    mergeNCBITax = True
taxFilter = args.taxMergeFilter
inFileTax = args.inFileTax + 'taxOut_'+args.taxMergeDivisions+'.gtxt'

if mergeNCBITax:        
    pgtf.appendGTxTFile(inFileTax, outFile, outFile[:-3]+'_m.gv',taxFilter)    

print '--> RUN DONE! <--'
