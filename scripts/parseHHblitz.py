'''
Created on Sep 1, 2013

@author: ranko
'''
''' 
parses HHblitz results
      
Other Notes:
  -> single query ONLY !!!

MODEL
   f: parseHMMResults(filePath):    --> parses HHblitz results (input: filename)
                                    --> returns: HHblitzResults (see below)
HHblitzResults
   -> query
   -> queryFile
   -> dbFile
   -> hhResults <list of HHResult>
  HHResult
      -> nr
      -> hdr
      -> hdrShort    
      -> prob         # chance to be good
      -> eV           # evalue (exp nr of false hits)
      -> pV           # eV/DB size (not really important)
      -> rawScore     # raw score for internal stuff
      -> ss score
      -> identities   # percentage of IDENTICAL columns in alignment
      -> similarity   # score based similarity (the bigger => the better)
      -> psipredScore # 2nd structure pred score
      -> nrMatches    # nr of matches in alignment
      -> qHMMStart    # start and end of alignment (query HMM model)
      -> qHMMEnd        
      -> hHMMStart    # start and end of alignmnet (target HMM modeL)
      -> hHMMEnd
      -> alignment    # list of strings for alignment display
'''

class HHblitzResults:
    def __init__(self):
        self.query = ''
        self.queryFile = ''
        self.dbFile = ''
        self.hhResults = []

class HHResult:
    def __init__(self):
        self.nr = 0
        self.qHdrShort = ''
        self.qHdr = ''
        self.hitHdr = ''
        self.hitHdrShort = ''
        self.prob = 0.0
        self.eV = 0.0
        self.pV = 0.0
        self.rawScore = 0.0
        self.identities = 0
        self.similarity = 0.0
        self.pripredScore = 0.0
        self.nrMatches = 0
        self.qHMMStart = -1
        self.qHMMEnd = -1
        self.hHMMStart = -1
        self.hHMMEnd = -1
        self.alignment = []
        self.gis = []  # NOTE THAT THOSE APPLY ONLY IF NR IS USED!!!
    
''' GRAB HHblits results, 
parse them, convert into object format (see above for details on 
format '''    
    
def parseHHResults(hhInputPath):
    hhResultsAll = HHblitzResults()
    # check if input seems ok:
    with open (hhInputPath) as hhF:
        c = 0
        looksFine = True
        for l in hhF:
            c+=1
            ls = l.split(' ')
            if c == 1 and not ls[0] == 'Query':
                looksFine = False
            if c == 2 and not ls[0] == 'Match_columns':
                looksFine = False
            if c == 3 and not ls[0] == 'No_of_seqs':
                looksFine = False
            if c == 4 and not ls[0] == 'Neff':
                looksFine = False
            if c == 5 and not ls[0] == 'Searched_HMMs':
                looksFine = False
            if c == 6 and not ls[0] == 'Date':
                looksFine = False
            if c == 7 and not ls[0] == 'Command':
                looksFine = False
            if c == 9 and not 'No Hit' in l:
                looksFine = False
            if c >= 10: 
                break
    if not looksFine:
        print ' does not look line HHblitz result file, aborting! '
        return     
    
    with open (hhInputPath) as hhF:
        parsingHeader = True
        parsingSummary = False
        parsingAlignments = False
        firstAlignment = True
        aNR = 0
        aR = HHResult()
        for l in hhF:
            l = l.strip();          
            if parsingSummary:
                if not l == '':
                    # grab NR
                    hhRnew = HHResult()
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]
                    hhRnew.nr = int(chunk)
                    # add query
                    hhRnew.qHdr = hhResultsAll.query
                    # grab short header of HIT
                    hhRnew.hHdrShort = l[0:30].strip()          
                    l = l[31:].strip()
                    # strip junk
                    while l.find('  ') > -1:
                        l = l.replace('  ',' ')
                    # grab prediction %
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]
                    hhRnew.prob = float(chunk)
                    # grab evalue                                        
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]                    
                    hhRnew.eV = float(chunk)
                    # grab pvalue
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]                    
                    hhRnew.eV = float(chunk)                                        
                    # grab score                                        
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]                    
                    hhRnew.rawScore = float(chunk)
                    # grab ss score                    
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]                    
                    hhRnew.pripredScore = float(chunk)
                    # grab columns
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]                    
                    hhRnew.nrMatches = int(chunk)
                    # grab query HMM
                    chunk = l[:l.find(' ')]
                    l = l[l.find(' ')+1:]                    
                    hhRnew.qHMMStart = int(chunk.split('-')[0])                                        
                    hhRnew.qHMMEnd = int(chunk.split('-')[1])
                    # grab target HMM
                    chunk = l[:l.find('(')]
                    l = l[l.find('(')+1:]                    
                    hhRnew.hHMMStart = int(chunk.split('-')[0])                                        
                    hhRnew.hHMMEnd = int(chunk.split('-')[1])                                        
                    hhResultsAll.hhResults.append(hhRnew)                
                if l == '':
                    parsingSummary = False
                    parsingAlignments = True                    
            if parsingHeader:
                #print 'hdr:', l
                if 'No Hit' in l:
                    parsingHeader = False
                    parsingSummary = True
                elif 'Query' in l:
                    hhResultsAll.query = l.strip(' ')[1]
                elif 'Command' in l:
                    while '  ' in l:
                        l = l.replace('  ',' ')
                    #chop out input file
                    iFname = l[l.find('-i')+3:]
                    iFname = iFname[:iFname.find(' ')]
                    hhResultsAll.queryFile = iFname
                    dFname = l[l.find('-d')+3:]
                    dFname = dFname[:dFname.find(' ')]
                    hhResultsAll.dbFile = dFname
                    
            # go over alignments                            
            if parsingAlignments:
                if 'No' in l[0:2]:
                    #d print l
                    if firstAlignment:
                        firstAlignment = False
                        aR = HHResult()                        
                        aR = hhResultsAll.hhResults[aNR]
                        #d print 'grabbed ',aNR
                    else: 
                        # add it
                        hhResultsAll.hhResults[aNR] = aR
                        aR = HHResult()
                        aNR +=1
                        # print aNR
                        aR = hhResultsAll.hhResults[aNR]
                        pass
                elif 'Done!' in l:
                    aR = hhResultsAll.hhResults[aNR]
                    pass                
                elif not l.strip() == '':
                    if l[0] == '>':
                        aR.hitHdr = l     
                        #print aR.hitHdr               
                        #print 'grabbing gis'
                        allGis = []
                        doGo = False
                        if ']|' in l:
                            doGo = True
                        while doGo: 
                            l = l[l.find(']|')+2:]
                            chunks = l[0:l.find(' ')]
                            for chunk in chunks.split('|'):
                                allGis.append(int(chunk))
                            if not ']|' in l: 
                                doGo = False
                        #print aR.gis						
                        aR.gis = allGis
                        #print aR.gis
                        #exit (-1)
                    elif 'Probab' in l[0:10]:
                        while '  ' in l:
                            l = l.replace('  ',' ')
                    else:                        
                        aR.alignment.append(l)
                    #add it          
    return hhResultsAll

''' ----- TEST ----- '''
def doTests():
    inFile = 'testNrf2.hhr' 
    res = parseHHResults(inFile)
    for hhRnew in res.hhResults: 
        print hhRnew.nr,' ',hhRnew.hHdrShort,'p',hhRnew.prob,'ev:',hhRnew.eV,'score:',hhRnew.rawScore,'ss:',hhRnew.pripredScore, 'm:',hhRnew.nrMatches, 'qry:',hhRnew.qHMMStart,'-',hhRnew.qHMMEnd, 'tar:',hhRnew.hHMMStart,'-',hhRnew.hHMMEnd, 'I:'
        print 'long header: ', hhRnew.hitHdr
        print 'gis: ',hhRnew.gis
        #    print ' ----- ALIGNMENT FOLLOWS ---------------- '
        #    for line in hhRnew.alignment:
        #        print line
        #    print ' ---------------------------------------- '
    _nrTaxBinarized = '/home/ranko/Development/DBs/ncbi_tax/gi_taxid_prot.bin'
    _names = '/home/ranko/Development/DBs/ncbi_tax/names.dmp'
    _nodes = '/home/ranko/Development/DBs/ncbi_tax/nodes.dmp'

    resF = filterHHblits(res, _nrTaxBinarized, _names, _nodes, targetTaxInc=[2])

    for hhRnew in resF.hhResults: 
        print hhRnew.nr,' ',hhRnew.hHdrShort,'p',hhRnew.prob,'ev:',hhRnew.eV,'score:',hhRnew.rawScore,'ss:',hhRnew.pripredScore, 'm:',hhRnew.nrMatches, 'qry:',hhRnew.qHMMStart,'-',hhRnew.qHMMEnd, 'tar:',hhRnew.hHMMStart,'-',hhRnew.hHMMEnd, 'I:'
        #print 'long header: ', hhRnew.hitHdr
        #print 'gis: ',hhRnew.gis
    
    resF = filterHHblits(res, _nrTaxBinarized, _names, _nodes, targetTaxInc=[9606])
    for hhRnew in resF.hhResults: 
        print hhRnew.nr,' ',hhRnew.hHdrShort,'p',hhRnew.prob,'ev:',hhRnew.eV,'score:',hhRnew.rawScore,'ss:',hhRnew.pripredScore, 'm:',hhRnew.nrMatches, 'qry:',hhRnew.qHMMStart,'-',hhRnew.qHMMEnd, 'tar:',hhRnew.hHMMStart,'-',hhRnew.hHMMEnd, 'I:'
        #print 'long header: ', hhRnew.hitHdr
        #print 'gis: ',hhRnew.gis    

    print ' RUN 3, filtering no filter (should show all hits)'    
    resF = filterHHblits(res, _nrTaxBinarized, _names, _nodes)
    for hhRnew in resF.hhResults: 
        print hhRnew.nr,' ',hhRnew.hHdrShort,'p',hhRnew.prob,'ev:',hhRnew.eV,'score:',hhRnew.rawScore,'ss:',hhRnew.pripredScore, 'm:',hhRnew.nrMatches, 'qry:',hhRnew.qHMMStart,'-',hhRnew.qHMMEnd, 'tar:',hhRnew.hHMMStart,'-',hhRnew.hHMMEnd, 'I:'
        #print 'long header: ', hhRnew.hitHdr
        #print 'gis: ',hhRnew.gis    
    
    #    print ' ----- ALIGNMENT FOLLOWS ---------------- '
    #    for line in hhRnew.alignment:
    #        print line
    #    print ' ---------------------------------------- '

    print ' RUN 4, filtering by "bacteria")'    
    resF = filterHHblits(res, _nrTaxBinarized, _names, _nodes, targetTaxNames=['bacteria'])
    for hhRnew in resF.hhResults: 
        print hhRnew.nr,' ',hhRnew.hHdrShort,'p',hhRnew.prob,'ev:',hhRnew.eV,'score:',hhRnew.rawScore,'ss:',hhRnew.pripredScore, 'm:',hhRnew.nrMatches, 'qry:',hhRnew.qHMMStart,'-',hhRnew.qHMMEnd, 'tar:',hhRnew.hHMMStart,'-',hhRnew.hHMMEnd, 'I:'
        #print 'long header: ', hhRnew.hitHdr
        #print 'gis: ',hhRnew.gis    


'''
Filters HHblits PARSED results by NCBI taxonomy,
keep following in mind:
 -> requires taxID as input
 -> taxIDs of interest: 
  Bacteria: 2
  Fungi: 4751
  Archaea: 2157
  Metazoa: 33208
  Viridiplantae: 33090
  human: 9606
  
  PROTISTA (are evil):  Alveolata: 336030;   Amoebazoa: 554915;   Apusozoa: 554296;   Breviatea: 1401294;   Centroheliozoa: 193537;
    Cryptophyta: 3027;  Euglenozoa: 33682;  Fornicata: 207245     Glaucocystophyceae: 38254     Haptophyceae: 2830  Heterolobosea: 5752
  Jakobida:556282   Katablepharidophyta: 339960   Malawimonadidae: 136087   Choanoflagellida: 28009   Nucleariidae and Fonticula group: 1001604
  Opisthokonta incertae sedis: 42461   Oxymonadida: 66288   Parabasalid: 5719   Rhizaria: 543769
  
  OTHER: unclassified eukaryotes: 42452; Rhodophyta: 2763; Stramenopiles: 33634; Chlorophyta: 3041

 -> can also accept 'named' variant, which simply translates
 string name into (list) of taxIDs
 '''
 
def filterHHblits (inputR, nrTaxBinarized, names, nodes, targetTaxInc = [], targetTaxExc = [], targetTaxNames = []):
    fn = 'none'
    outputR = HHblitzResults()
    corig = len(inputR.hhResults)
    cfiltered = 0
    gisfiltered = 0
    from debinarizeNCBItax import loadTaxID
    import parseNCBItaxonomy as pnt
    import copy

# following can be handled: 'none','bacteria','archaea','fungi','protista','misc','metazoa','plants'
    if 'none' in targetTaxNames:
        return inputR
    if 'bacteria' in targetTaxNames:
        if not 2 in targetTaxInc:
            targetTaxInc.append(2)
    if 'archaea' in targetTaxNames:
        if not 2157 in targetTaxInc:
            targetTaxInc.append(2157)
    if 'fungi' in targetTaxNames:
        if not 4751 in targetTaxInc:
            targetTaxInc.append(4751)
    if 'metazoa' in targetTaxNames:
        if not 33208 in targetTaxInc:
            targetTaxInc.append(33208)
    if 'viridiplantae' in targetTaxNames or 'plants' in targetTaxNames:
        if not 33090 in targetTaxInc:
            targetTaxInc.append(33090)
    if 'protista' in targetTaxNames:
        protaxid = [336030,554915,554296,1401294,193537,3027,33682,207245,38254,2830,5752,556282,339960,136087,28009,1001604,42461,66288,5719,543769]
        for i in protaxid:
            targetTaxInc.append(i)
    if 'misc' in targetTaxNames:
        targetTaxInc.append(33634)
        targetTaxInc.append(2763)
    
    
    #d print ' ---> FILTERING HHblits results <--- '
    #d print '     --> grabbing taxa: '+str(targetTaxInc)
    #d print '     --> excluding taxa: '+str(targetTaxExc)  
    #d print '   --> initing Taxonomy DB <-- '    
    pnt.initTax(names, nodes)
    #for a in pnt.taxMap.keys():
    #    print a
    #d print '   --> processing input HHblits results <-- '
    for oneRes in inputR.hhResults:
        newOneResult = copy.deepcopy(oneRes)   
        newOneResult.gis = []              
        #d print oneRes.hitHdr
        #d print oneRes.gis
        #d print newOneResult.gis
        for oneGI in oneRes.gis:
            taxID = loadTaxID(nrTaxBinarized, oneGI)
            if not taxID == 0: 
                #d print taxID
                doInclude = False
                allMyTaxa = pnt.getNumericTax(taxID)
                if len(targetTaxInc) > 0:
                    for incTax in targetTaxInc:
                        if incTax in allMyTaxa:
                            doInclude = True
                else: 
                    doInclude = True                
                for excTax in targetTaxExc:
                    if excTax in allMyTaxa:
                        doInclude = False            
                if doInclude:
                    #d print 'including','gi',oneGI,'taxid:',taxID
                    newOneResult.gis.append(oneGI)
                    gisfiltered +=1
                else:
                    pass
        #d print 'after filter: ',newOneResult.gis        
        if len(newOneResult.gis) > 0:
            newOneResult.hitHdr = newOneResult.hitHdr[0:newOneResult.hitHdr.find('. ')+2]
            for gi in newOneResult.gis:
                newOneResult.hitHdr += '|'+str(gi)
            newOneResult.hitHdr += '|'
            outputR.hhResults.append(newOneResult)
            cfiltered +=1
    #d print 'filtered results from',corig,'to',cfiltered,'hits',gisfiltered,'GIs' 
    return outputR
