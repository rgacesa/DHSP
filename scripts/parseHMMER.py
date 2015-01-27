'''
Created on Aug 2, 2013

@author: ranko
'''

''' 
parses HMMER results
      
Other Notes:
 -> hmmer can have multiple queries, and multiple iterations per query, 
 therefore parsed results are in format:
  qry<1+> -> iteration<1+> -> results<1+> 

MODEL
   f: parseHMMResults(filePath):    --> parses HMMER results (input: filename)
                                     --> returns: multiQueryResults (see below)
multiQueryResults
   -> qry <1 or more>                     {HMMER results for xth query sequence}
       -> hmmerType <1>                   {type of HMMER used (jack or normal) }
       -> queryFile <1>                   {file that was used to query HMM db}
       -> resultsFile <1>                 {results file for HMMER run}
       -> targetDB <1>                    {HMM db that was queried}
       -> qryNr    <1>                    {number of query, as entered}
       -> iterationResults <1 or more>    {all results for iterations}
            -> iteration                  {simply number of iteration}
                f: getLastIteration --> returns last iteration
                f: getIterationsNr  --> returns # of iterations 
            -> hmmTable                   {tabular results}
                    -> seqID
                    -> fullSeqEV
                    -> fullSeqBS
                    -> fullSeqBSBias
                    -> bestDomEV
                    -> bestDomBS
                    -> bestDomBSBias
                    -> nrDom
                    -> nrDomExp
                    -> seqID
                    -> description
                    
#NYI                -> domainResults [1 or more]
                        -> hitNR
                        -> bitScore
                        -> bitScoreBias
                        -> cEvalue
                        -> eValue
                        -> hmmFrom
                        -> hmmFromChar
                        -> hmmTo
                        -> hmmToChar
                        -> aliFrom
                        -> aliFromChar
                        -> aliTo
                        -> aliToChar
                        -> accuracy
                        -> reliable
                        -> alignmentQ
                        -> alignmentQTrim
                        -> alignmentM
                        -> alignmentMTrim
                        -> alignmentH
                        -> alignmentHTrim
                        

#TODO: write rest of parser!
'''

class HmmResultOneDomain:
    def __init__(self):
        self.bitScore = -1.0
        self.bitScoreBias = -1.0
        #self.cEvalue = -1.0
        self.eValue = -1.0
        self.hmmFrom = -1
        self.hmmFromChar = '*'
        self.hmmTo = -1
        self.hmmToChar = '*'
        self.aliFrom = -1
        self.aliFromChar = '*'
        self.aliTo = -1
        self.aliToChar = '*'
        self.accuracy = -1.0
        self.reliable = None
        self.alignmentQ = ''
        self.alignmentQTrim = ''
        self.alignmentM = ''
        self.alignmentMTrim = ''
        self.alignmentH = ''
        self.alignmentHTrim = '' 
        self.domainNR = 0       

class HmmResultsMQuery:
    def __init__(self):
        self.qry = []

class HmmTableResult:
    def __init__(self):
        self.fullSeqEV = 0.0
        self.fullSeqBS = 0.0
        self.fullSeqBSBias = 0.0
        self.bestDomEV = 0.0
        self.bestDomBS = 0.0
        self.bestDomBSBias = 0.0
        self.nrDom = 0.0
        self.nrDomExp = 0.0
        self.seqID = ""
        self.description = ""
        self.domainResults = []
        
class HmmResultsOneQuery:
    def __init__(self):
        self.iterationResults = []
        self.queryFile = ""
        self.resultsFile = ""
        self.targetDB = ""
        self.hmmerType = ""
        self.qryNr = 0  
    def getIterationsNr(self):
        return len(self.iterationResults)
    def getLastIteration(self):
        return self.iterationResults[self.getIterationsNr()-1]

class HmmResultsOneIteration:
    def __init__(self):
        self.hmmTable = []
        self.iteration = 0
    
def parseHMMResults (hmmRes):
#    print "parsing..."
#    hmmRes = '/home/ranko/workspace_e4/PyParsers/Results/HMM_RESULTS/run3_23.08/bp_nrf2_hs_bzip_nr_ref_e0.01_c_0.5_BAC.jhmr'
    with open(hmmRes) as hmmFile:
        # domain hit(s)
        hitDoms = []
        oneDomHit = HmmResultOneDomain()
        
        # file lines
        hmmFL = hmmFile.readlines()
        # all results
        hmmAllResults = HmmResultsMQuery()     
        # tmp results (one query)
        tmpHmmResultsOneQuery = HmmResultsOneQuery()    
        # tmp table (one table for one set of results)    
        tmpHmmTable = []
        # tmp one iteration
        tmpHmmerOneIteration = HmmResultsOneIteration()    
        # for decisions on hmmfullSeqEer type parsing
        hmmerType = ""    
        # PARSE HEADER    
        doParseHeader = True
        # PARSE TABLE
        doParseTable = False              
        startParseTable = False
        # PARSE DOMAINS & ALIGNMENTS
        doParseDomains = False
        startParseDomains = False
        hitNrForDoms = 0    
        tmpDomainResults = []    
        
        # PARSE QUERY    
        queryNR = 0
        newQuery = True
        newIteration = True
        # for managing jackHmmer
        jackHmmerIteration = 0
        lineCounter = -1
        for line in hmmFL:
            lineCounter+=1
            #print line
            # catch new query
            if '#' not in line and '//' in line:
                newIteration = True
                newQuery = True
            
            # catch jackhmmer iteration
            if 'jackhmmer' in hmmerType and '@@ Round:' in line:                        
                newIteration = True
        
            if 'Domain annotation for each sequence (and alignments)' in line:
                doParseDomains = True
                doParseHeader = False
                doParseTable = False
                startParseDomains = True
                grabDomGenericInfo = False
                
        
            if newIteration:
                doParseTable = False        
                startParseTable = False
                hitNrForDoms = 0 
                tmpHmmTable = []
                tmpHmmerOneIteration.iteration = jackHmmerIteration
                if jackHmmerIteration > 0:
                    tmpHmmResultsOneQuery.iterationResults.append(tmpHmmerOneIteration)            
                tmpHmmerOneIteration = HmmResultsOneIteration()
                jackHmmerIteration +=1            
                doParseHeader = True
                # PARSE TABLE
                doParseTable = False                    
                startParseTable = False        
                newIteration = False
                doParseDomains = False
                startParseDomains = False 
                grabDomGenericInfo = False
                                            
          
            if newQuery:
                # initialize                            
                doParseHeader = True
                hitNrForDoms = 0 
                # PARSE TABLE
                doParseTable = False        
                startParseTable = False
                # PARSE DOMAINS
                doParseDomains = False
                startParseDomains = False   
                grabDomGenericInfo = False      
                                       
                if queryNR >= 1:
                    hmmAllResults.qry.append(tmpHmmResultsOneQuery)
                queryNR += 1
                newQuery = False
                tmpHmmResultsOneQuery = HmmResultsOneQuery()
                tmpHmmerOneIteration = HmmResultsOneIteration()                                        
                
            if doParseHeader:
                if "# hmmsearch :: " in line:
                    hmmerType = "hmmsearch"
                    tmpHmmResultsOneQuery.hmmType = "hmmsearch"
                if "# jackhmmer :: " in line:
                    hmmerType = "jackhmmer"         
                    tmpHmmResultsOneQuery.hmmType = "jackhmmer"                    
                if "# query HMM file:" in line or '# query sequence file:' in line: 
                    tmpHmmResultsOneQuery.queryFile = line[35:].strip()
                if "# target sequence database:" in line: 
                    tmpHmmResultsOneQuery.targetDB = line[35:].strip()
                if "# output directed to file:" in line: 
                    tmpHmmResultsOneQuery.resultsFile = line[35:].strip()
                if "Scores for complete sequences" in line:
                    doParseHeader = False
                    doParseTable = True
                    doParseDomains = False
                    startParseDomains = False
                    startParseTable = False
                    grabDomGenericInfo = False                    
                
            elif doParseTable:
                doParseDomains = False
                startParseDomains = False

                line = line.strip()
                if len(line) > 0:
                    if line[0] == '-':
                        line = ' '+line[1:]
                    if line[0] == '+':
                        line = ' '+line[1:]
            # end of table
                if line.strip() == '':
                    startParseTable = False
                    doParseTable = False 
                    tmpHmmerOneIteration.hmmTable = tmpHmmTable  
                                                     
                if startParseTable:
                    line = line.strip()
                
                    if 'inclusion threshold' not in line:
                        while '  ' in line:                    
                            line = line.replace('  ',' ')
                        
                        sl = line.split(' ')
                        # RECORD RESULT
                        tabResult = HmmTableResult()
                        # FULL SEQ EV                            
                        tabResult.fullSeqEV = float(sl[0].strip())
                        # FULL SEQ SCORE
                        tabResult.fullSeqBS = float(sl[1].strip())
                        # FULL SEQ BIT SCORE BIAS
                        tabResult.fullSeqBSBias = float(sl[2].strip())
                    
                        # BEST DOMAIN EV
                        tabResult.bestDomEV = float(sl[3].strip())
                        # BEST DOMAIN BITSCORE
                        tabResult.bestDomBS = float(sl[4].strip())
                        # BEST DOMAIN BITSCORE BIAS
                        tabResult.bestDomBSBias = float(sl[5].strip())
                    
                        # NR DOMAINS
                        tabResult.nrDom = float(sl[6].strip())
                        # EXPECTED NR DOMAINS
                        tabResult.nrDomExp = float(sl[7].strip())
                
                        # SEQ ID
                        tabResult.seqID = sl[8].strip()
                
                        # DESCRIPTION (OPTIONAL)
                        tabResult.description = ""
                        if len(sl) > 9:
                            d = ''
                            for desL in range (9,len(sl)):
                                d = d + ' '+sl[desL] 
                            tabResult.description = d
#                            print tabResult.description
                
                        # SAVE IT TO TABLE
                        #print tabResult.seqID,tabResult.fullSeqEV,"added"                
                        tmpHmmTable.append(tabResult)            
                                
                if "    ------- ------ -----" in line: 
                    startParseTable = True
                    
            elif doParseDomains:                
                #print line.strip()
                if '>>' in line:
                    #print 'parsing domains for hit'
                    grabDomGenericInfo = False
                    if startParseDomains:
                        startParseDomains = False                     
                    elif not startParseDomains:
                        #print 'adding dom hits to previous hit'
                        #print 'nrhits:',hitNrForDoms, 'len of table', len(tmpHmmTable), ' ??? '
                        tmpHmmTable[hitNrForDoms].domainResults = tmpDomainResults
                        hitNrForDoms += 1
                        tmpDomainResults = []
                elif line.strip() == '':
                    grabDomGenericInfo = False
                    
                elif grabDomGenericInfo:
#                    print line.strip()
                    lineT = line.strip()
                    while '  ' in lineT:
                        lineT = lineT.replace('  ',' ')
#                    print lineT
                    lineTS = lineT.split(' ')                    
                    #print lineCounter, "|".join(lineTS)
                    oneDomGenericRes = HmmResultOneDomain()
                    oneDomGenericRes.domainNR = lineTS[0] 
                    oneDomGenericRes.reliable = lineTS[1]
                    oneDomGenericRes.bitScore = lineTS[2]
                    oneDomGenericRes.bitScoreBias = lineTS[3]
                    # we skip cValue, it sucks
                    oneDomGenericRes.eValue = lineTS[5]
                    oneDomGenericRes.hmmFrom = lineTS[6]
                    oneDomGenericRes.hmmTo = lineTS[7]                    
                    oneDomGenericRes.hmmFromChar = lineTS[8][0]
                    oneDomGenericRes.hmmToChar = lineTS[8][1]
                    oneDomGenericRes.aliFrom = lineTS[9]
                    oneDomGenericRes.aliTo = lineTS[10]
                    oneDomGenericRes.aliFromChar = lineTS[11][0]
                    oneDomGenericRes.aliToChar = lineTS[11][1]
                    oneDomGenericRes.accuracy = lineTS[15]
                    tmpDomainResults.append(oneDomGenericRes)                                         
                    
                elif 'Alignments for each domain:' in line:
                    doStop = False
                    lc = 0
                    ndlc = 0                    
                    offSet = 0       
                    domainNR = 0     
                    newDomain = True        
                    while not doStop:
                        lc += 1                                  
                        if newDomain and domainNR == 0:                             
                            newDomain = False    
                            domainNR +=1 
                            ndlc = 0
                            oneDomAlignments = HmmResultOneDomain()    
                        elif newDomain and len(tmpDomainResults) > domainNR:  
                            # add domain stuff to domain!
                            #print line
                            #print 'len tdr:',len(tmpDomainResults)
                            #print 'domainNR',domainNR                            
                            dom = tmpDomainResults[domainNR-1]
                            dom.alignmentMTrim = oneDomAlignments.alignmentMTrim
                            dom.alignmentM =     oneDomAlignments.alignmentM
                            dom.alignmentHTrim = oneDomAlignments.alignmentHTrim
                            dom.alignmentH =     oneDomAlignments.alignmentH
                            dom.alignmentQTrim = oneDomAlignments.alignmentQTrim
                            dom.alignmentQ =     oneDomAlignments.alignmentQ
                            tmpDomainResults[domainNR-1] = dom                            
                            domainNR +=1
                            newDomain = False
                            ndlc = 0
                            oneDomAlignments = HmmResultOneDomain()
                            
                        if ndlc == 1:
                            # ZLOK
                            #d print hmmFL[lineCounter+lc]
                            #d print hmmFL[lineCounter+lc].strip().split(' '), len(hmmFL[lineCounter+lc].strip().split(' '))
                            if len(hmmFL[lineCounter+lc].strip().split(' ')) < 3:                        
                                lc +=1
                                ndlc +=1
                                      
                        if hmmFL[lineCounter+lc].strip() == '':
                            newDomain = True
                        
                        #print ndlc,hmmFL[lineCounter+lc],
                        offSet = 0
                        if ndlc == 2:
                            oneDomAlignments.alignmentQ = str(hmmFL[lineCounter+lc])
                        elif ndlc == 3:
                            oneDomAlignments.alignmentM = str(hmmFL[lineCounter+lc])
                        elif ndlc == 4:
                            oneDomAlignments.alignmentH = str(hmmFL[lineCounter+lc])
                        elif ndlc == 5:
                            trimmer = str(hmmFL[lineCounter+lc])
                            for c in range(0,len(trimmer)):
                                if not trimmer[c] == ' ':
                                    offSet = c
                                    break 
                            oneDomAlignments.alignmentMTrim = str(oneDomAlignments.alignmentM[offSet:])
                            oneDomAlignments.alignmentQTrim = str(oneDomAlignments.alignmentQ[offSet:])
                            oneDomAlignments.alignmentHTrim = str(oneDomAlignments.alignmentH[offSet:])
                            oneDomAlignments.alignmentQTrim = str(oneDomAlignments.alignmentQTrim[:len(oneDomAlignments.alignmentMTrim)-1])
                            oneDomAlignments.alignmentHTrim = str(oneDomAlignments.alignmentHTrim[:len(oneDomAlignments.alignmentMTrim)-1])                                                                                                          
                        ndlc += 1 
                        #print lc,ndlc,hmmFL[lineCounter+lc].strip()                                                              
                        if '>>' in hmmFL[lineCounter+lc] or '//' in hmmFL[lineCounter+lc]: #or 'Internal pipeline' in hmmFL[lineCounter+lc]:                         
                            doStop = True
                            #print ' --- stopping at >> --- '
                        elif hmmFL[lineCounter+lc].strip() == '' and hmmFL[lineCounter+lc+1].strip() == '':
                            doStop = True
                            #print ' stopping at "" + "" ---'
                        else: 
                            pass               
                    
                elif ' ---   ------' in line:
                    grabDomGenericInfo = True                                        
    return hmmAllResults
# ------------------------ END OF PARSE MODULE ------------------
'''            
# MAIN: TEST
target = '/home/ranko/Development/workspace/DHPipeline/results/HMM_RESULTS/Nrf2_torafugu_vs_nr_bpsi_gi_e0.001_c_0.6_vs_zoophyte.hmr'
hmmAllResults = parseHMMResults(target)
print " --- ALL QUERIES --- "
cntQ = 0
print hmmAllResults.qry[0]
for hmmQR in hmmAllResults.qry:
    cntQ +=1
    print " **** QUERY ",cntQ," ****"
    print "res file: ", hmmQR.resultsFile                
    print "target DB: ", hmmQR.targetDB                
    print "query file:", hmmQR.queryFile
    print "iterations: ", len(hmmQR.iterationResults)    
    for hmmQI in hmmQR.iterationResults:        
        print " #### ITERATION: ", hmmQI.iteration," ####"        
        print " @@@@ RESULTS TABLE @@@@ "        
        cnt = 0
        for h in hmmQI.hmmTable:
            cnt +=1
            if cnt < 5000:
                if len(h.domainResults) > 1:
                    print ' ************ NEW HIT ********************** '
                    print cnt, h.seqID, "eV: ",h.fullSeqEV
                    for d in h.domainResults:
                        print ' --------- DOMAIN ------ '
                        print 'Dnr:',d.domainNR,'R:',d.reliable
                        print 'BS:',d.bitScore, '+/-',d.bitScoreBias
                        print 'EV:',d.eValue,'acc',d.accuracy
                        print 'HMM:',d.hmmFrom,d.hmmFromChar,'->',d.hmmTo,d.hmmToChar
                        print 'ALI:',d.aliFrom,d.aliFromChar,'->',d.aliTo,d.aliToChar
                        print d.alignmentQ,
                        print d.alignmentM,
                        print d.alignmentH,
                        print ' ---------- END OF DOMAIN -- '
                    print ' ******************************************* '
#                   print d.alignmentQTrim,
#                   print d.alignmentMTrim,
#                   print d.alignmentHTrim,           
print ' ----------------------------------------------------- '            
print ' ------------ LAST ITERATION ONLY -------------------- '
print ' ----------------------------------------------------- '
cntQ = 0
for hmmQR in hmmAllResults.qry:
    cntQ +=1
    print " **** QUERY ",cntQ," ****"
    print "res file: ", hmmQR.resultsFile                
    print "target DB: ", hmmQR.targetDB                
    print "query file:", hmmQR.queryFile
    print "iterations: ", len(hmmQR.iterationResults)    
    
    hmmQI = hmmQR.getLastIteration()        
    print " #### ITERATION: ", hmmQI.iteration," ####"        
    print " @@@@ RESULTS TABLE @@@@ "        
    cnt = 0
    for h in hmmQI.hmmTable:
        cnt +=1
        if cnt < 50:
            print cnt, h.seqID, "eV: ",h.fullSeqEV
'''