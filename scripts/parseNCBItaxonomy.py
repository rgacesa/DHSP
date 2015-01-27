'''
Created on Aug 29, 2014
'''
taxIDNames = {}
taxMap = {}
''' tax record class 
    contains complete taxonomy record
    for single organism (taxID) '''
class TaxRecord:
    def __init__(self):
        self.records = []
        
    def add(self, rec):
        self.records.append(rec)
        
    def toStrList(self):
        ret = ''
        for l in self.records:
            ret = ret + l.toStr() + '\n'
        return ret

    def toStrRow(self):
        ret = ''
        for l in self.records:
            ret = ret + l.toStrB() + ' -> '
        return ret[0:-4]
    
    def toStrDot(self):
        ret = ''
        for l in self.records:
            if not l.name == 'N/D':
                ret = ret + '"'+l.name+'"' + ' -> '
        return ret[0:-4]+";"    
    
    def toHash(self):
        ret = {}
        for l in range(0,len(self.records)-1):
            # if not self.records[l].name == 'N/D':
            ret[self.records[l].name] = self.records[l+1].name
        return ret
         
    def toTupleList(self):
        ret = []
        for l in range(0,len(self.records)-1):
            # if not self.records[l].name == 'N/D':
            t = self.records[l].name, self.records[l+1].name
            ret.append(t)
        return ret
         
    def toStrListMR(self):
        ret = ''
        for l in self.records:
            if l.rank != 'no rank':
                ret = ret + l.toStr() + '\n'
        return ret
    
    def getRank(self, rank):
        ret = TaxRecordOneLine()
        ret.rank = rank
        ret.taxID = 0        
        self.name = 'N/D'        
        
        ret
        for l in self.records:
            if l.rank != 'no rank' and l.rank == rank:
                ret = l
        return ret        
    
    def getTax7(self):
        newRec = TaxRecord()
        newRec.add(self.getRank('species'))
        newRec.add(self.getRank('genus'))
        newRec.add(self.getRank('family'))
        newRec.add(self.getRank('order'))
        newRec.add(self.getRank('class'))
        newRec.add(self.getRank('phylum'))
        if self.getRank('kingdom').name == 'N/D' or self.getRank('kingdom').name == '':
            newRec.add(self.getRank('superkingdom'))
        else:
            newRec.add(self.getRank('kingdom'))            
        return newRec
    
    def getTax8(self):                        
        newRec = TaxRecord()
        if not self.getRank('species').name == 'N/D':
            newRec.add(self.getRank('species'))
        if not self.getRank('genus').name == 'N/D':            
            newRec.add(self.getRank('genus'))
        if not self.getRank('family').name == 'N/D':
            newRec.add(self.getRank('family'))
        if not self.getRank('order').name == 'N/D':            
            newRec.add(self.getRank('order'))
        if not self.getRank('class').name == 'N/D':            
            newRec.add(self.getRank('class'))
        if not self.getRank('phylum').name == 'N/D':            
            newRec.add(self.getRank('phylum'))
        if not self.getRank('kingdom').name == 'N/D':
            newRec.add(self.getRank('kingdom'))
        if not self.getRank('superkingdom').name == 'N/D':            
            newRec.add(self.getRank('superkingdom'))            
        return newRec
    
    def getTax9E(self, taxF='SGFOCPKD',trimSpecies='Y'):                
        newRec = TaxRecord()
        if not self.getRank('species').name == 'N/D' and 'S' in taxF:
            sR = self.getRank('species')
#            print sR.name
            if trimSpecies == 'Y' and ' ' in sR.name:
#                print 'trimming name'
                sRNS = sR.name.split(' ')
                sR.name = sRNS[0][0:3] + '. '+sRNS[1]
#                print 'trimmed: ',sR.name                               
            newRec.add(sR)                        
        if not self.getRank('genus').name == 'N/D' and 'G' in taxF:            
            newRec.add(self.getRank('genus'))
        if not self.getRank('family').name == 'N/D' and 'F' in taxF:
            newRec.add(self.getRank('family'))
        if not self.getRank('order').name == 'N/D' and 'O' in taxF:            
            newRec.add(self.getRank('order'))
        if not self.getRank('class').name == 'N/D' and 'C' in taxF:            
            newRec.add(self.getRank('class'))
        if not self.getRank('phylum').name == 'N/D' and 'P' in taxF:            
            newRec.add(self.getRank('phylum'))
        if not self.getRank('kingdom').name == 'N/D' and 'K' in taxF:
            newRec.add(self.getRank('kingdom'))
        if not self.getRank('superkingdom').name == 'N/D' and 'D' in taxF:            
            newRec.add(self.getRank('superkingdom'))
            
        newRec.add(self.getRank('end'))                    
        return newRec    

    ''' as 9E, but also keep track of 'blind spots' 
    and adds them as &X: notation 
    ideally by 10 points (0 for first, 10 for second...)
    '''
    def getTax9EwDummy(self, taxF='SGFOCPKD',trimSpecies='Y',excludedTaxa=[]):
        # order by level: 
        # S -> G -> F -> O -> K -> P -> K -> D
        isStarterDefined = False
        if 'D' in taxF:
            starget = 'superkingdom'        
        if 'K' in taxF:
            starget = 'kingdom'
        if 'P' in taxF:
            starget = 'phylum'                                                            
        if 'C' in taxF:
            starget = 'class'
        if 'O' in taxF:            
            starget = 'order'        
        if 'G' in taxF:
            starget = 'genus'
        if 'F' in taxF:
            starget = 'family'                    
        if 'S' in taxF:
            starget = 'species'        
        #d print starget
        #d print self.getRank(starget).name        
        if not self.getRank(starget).name == 'N/D':
            isStarterDefined = True            
                        
        newRec = TaxRecord()
        if not self.getRank('species').name == 'N/D' and 'S' in taxF:
            sR = self.getRank('species')            
#            print sR.name
            if trimSpecies == 'Y' and ' ' in sR.name:
#                print 'trimming name'
                sRNS = sR.name.split(' ')
                sR.name = sRNS[0][0:3] + '. '+sRNS[1]
#                print 'trimmed: ',sR.name                               
            newRec.add(sR)         
        elif 'S' in taxF and isStarterDefined: 
            tmpTR = '&S:'+self.getRank('species')
            newRec.add(tmpTR) 
                                   
        if not self.getRank('genus').name == 'N/D' and 'G' in taxF:            
            newRec.add(self.getRank('genus'))
            isStarterDefined = True
        elif 'G' in taxF and isStarterDefined:                     
            tmpTR = self.getRank('genus')
            tmpTR.name = '&G:'+str(self.getRank('species').name)
            newRec.add(tmpTR)   
                   
        if not self.getRank('family').name == 'N/D' and 'F' in taxF:
            newRec.add(self.getRank('family'))
            isStarterDefined = True
        elif 'F' in taxF and isStarterDefined: 
            tmpTR = self.getRank('family')
            tmpTR.name = '&F:'+str(self.getRank('genus').name)
            newRec.add(tmpTR)    
                    
        if not self.getRank('order').name == 'N/D' and 'O' in taxF:            
            newRec.add(self.getRank('order'))
            isStarterDefined = True
        elif 'O' in taxF and isStarterDefined: 
            tmpTR = self.getRank('order')
            tmpTR.name = '&O:'+str(self.getRank('family').name)
            newRec.add(tmpTR)     
                    
        if not self.getRank('class').name == 'N/D' and 'C' in taxF:            
            newRec.add(self.getRank('class'))
            isStarterDefined = True
        elif 'C' in taxF and isStarterDefined: 
            tmpTR = self.getRank('class')
            tmpTR.name = '&C:'+str(self.getRank('order').name)
            newRec.add(tmpTR) 
                 
        if not self.getRank('phylum').name == 'N/D' and 'P' in taxF:
            newRec.add(self.getRank('phylum'))
            isStarterDefined = True
        elif 'P' in taxF and isStarterDefined: 
            tmpTR = self.getRank('phylum')
            tmpTR.name = '&P:'+str(self.getRank('class').name)
            newRec.add(tmpTR)                                 
                           
        if not self.getRank('kingdom').name == 'N/D' and 'K' in taxF:
            newRec.add(self.getRank('kingdom'))
            isStarterDefined = True
        elif 'K' in taxF and isStarterDefined: 
            tmpTR = self.getRank('kingdom')
            tmpTR.name = '&K:'+str(self.getRank('phylum').name)
            newRec.add(tmpTR)     
                                         
        if not self.getRank('superkingdom').name == 'N/D' and 'D' in taxF:            
            newRec.add(self.getRank('superkingdom'))
        elif 'D' in taxF and isStarterDefined: 
            tmpTR = self.getRank('superkingdom')
            tmpTR.name = '&D:'+str(self.getRank('kingdom').name)
            newRec.add(tmpTR)                                   
        newRec.add(self.getRank('end'))
        
#        print len(newRec.records)
        doadd = True
        sl = newRec.toStrList()      
        for et in excludedTaxa:
            if et in sl:
                doadd = False
        if doadd:
            return newRec
        else:
            return None 

''' record for one line of taxonomy
basically contains taxID, name and rank '''
class TaxRecordOneLine:
    def __init__(self, _taxID = 0, _name = 'N/D', _rank = ''):
        self.taxID = _taxID
        self.rank = _rank
        self.name = _name
        if _rank == 'species':
            self.name = 'S:'+self.name
        elif _rank == 'genus':
            self.name = 'G:'+self.name
        elif _rank == 'family':
            self.name = 'F:'+self.name
        elif _rank == 'order':
            self.name = 'O:'+self.name
        elif _rank == 'class':
            self.name = 'C:'+self.name
        elif _rank == 'phylum':
            self.name = 'P:'+self.name        
        elif _rank == 'kingdom':
            self.name = 'K:'+self.name
        elif _rank == 'superkingdom':
            self.name = 'D:'+self.name
        elif _rank == 'end':
            self.name = 'END'
            
                    
    def toStr(self):
        return self.name+","+self.rank
    def toStrB(self):
        return self.rank+":"+self.name    
    def toStrS(self):
        return self.name,self.rank
    def toStrAll(self):
        return self.taxID+","+self.name+","+self.rank
    def toStrAllS(self):
        return self.taxID,self.name,self.rank    

''' loads taxonomy nodes from NCBI taxonomy nodes dump
    loads it to hash: taxID -> taxID, rank, name '''
def initTax(pathNames, pathNodes):
    global taxMap
    global taxIDNames
    taxMap = {}
    ''' load nodes '''
    with open(pathNodes) as nodesFile:
        for l in nodesFile:                
            l = l.replace("\t|\n","")
            splitL = l.split('\t|\t')
#            print str(splitL)
            sL = []
            #sL.append(splitL[0])
            sL.append(splitL[1])
            sL.append(splitL[2])
            taxMap[int(splitL[0])] = sL
    ''' load names '''    
    taxIDNames = {}
    with open(pathNames) as namesFile:
#        lines = namesFile.readlines()
        for l in namesFile:
            l = l.replace("\t|\n","")
            splitL = l.split('\t|\t')
#            print str(splitL)
            if splitL[3].strip() == 'scientific name':
                key = int(splitL[0])
                sL = []
                sL.append(splitL[1])
#                sL.append(splitL[1])
                taxIDNames[key] = sL                                                       

''' *****  retrieve taxonomy from TAXid ***** 
    creates taxonomy record for target taxID
                                            '''
def getTaxFromTaxID(taxID):
    global taxIDNames
    global taxMap
    doGo = True
    taxRecord = TaxRecord()
    lvl = 0
    maxLvl = 50
    while doGo:
        lvl += 1
        try:
            l = taxMap[taxID]
            futureTaxID = int(l[0])
            taxRecOne = TaxRecordOneLine(taxID, taxIDNames[taxID][0],l[1].strip())
            taxRecord.add(taxRecOne)        
            taxID = futureTaxID
            if lvl > maxLvl or int(taxID) == 1 or int(taxID) == 131567 or int(taxID) == -1:
                doGo = False
        except:
            print 'warning: taxID:',taxID,'does not exist in taxonomy!'
            doGo = False     
    return taxRecord

''' basically just for specific situations, returns
list of all taxIDs this guy has '''
def getNumericTax (taxID):
    global taxIDNames
    global taxMap
    doGo = True
    lvl = 0
    maxLvl = 50
    ret = []   
    #print taxID 
    while doGo:
        ret.append(taxID)
        lvl +=1
        try:        
            l = taxMap[taxID]
            futureTaxID = int(l[0])
            taxID = futureTaxID    
            if lvl > maxLvl or int(taxID) == 1 or int(taxID) == 131567 or int(taxID) == -1:
                doGo = False            
        except:
            print 'warning: taxID:',taxID,'does not exist in taxonomy!'
            doGo = False
                         
    ret.append(taxID)
    return ret 
   
   
'''
  --- extra functions for id Mapping ---
Created on 7 Oct 2014
connects various sequence IDs
format should be:
 --> gi: GI == GI, so no need to anything here
 --> ref|xyz| for refseq                            (NP_149472.1 or similar)
 --> sp|xyz| for swissprot [annotated uniprotKB] -  usually Q... (ex Q197F2) 
 --> tr|xyz| trembl (not annotated uniprotKB)
 --> pdb|xyz| for protein data bank (PDB)           (number)
 --> gb|xyz| genebank
 --> emb|xyz| embl
grabs headers, grabs all if possible, otherwise just some
then goes with those against DB and returns GI (if possible)
if not, doesn't return anything :)

notes:
 --> if GI is inside, just throws it back, no need to parse anything
'''

''' IMPORTS '''
import string

''' extracts GIs from idMap file 
 -> notes: 
- takes LIST of HEADER hashes produced by parseHeader
so basically something like: 
[{'ref':'12312312','pdb':'12937'},{},{}...]
- input is this stuff to enable extraction in single run
to minimize time required to read huge idMap file
- returns list of gi lists (as header can point to multiple GIs, 
it can return stuff that looks like: 
[[123123],[12312,23423,245345],[],[5435345],[543534,345345]]
'''
def extractFromIdMap(hdrList, inputPath, maxGIsPerHdr=-1, out='gi'):
    # possible out = 'gi' or 'taxid'
    #possibleIDs = ('gi','ref','sp','tr','emb','pdb','gb')
    #idMapperFile, tab split:
    # [0] = UniProtKB-AC {sp, tr or emb} <---; [1] = UniProtKB-ID {useless}
    # [2] = GeneID, entrez (useless?) 
    # [3] = refseq {ref} <--
    # [4] = gi (can be multiple!) {gi}  <--- we grab this one!
    # [5] = pdb (pdb) <--
    # [6] = go (useless?); [7] = ipi (?!); [8-11] = uni stuff? 
    # [12] = pir, [13] = ncbi tax, ...
    # [17] = embl <-- 
    if maxGIsPerHdr == -1: 
        maxGIsPerHdr = 100000
        
    maxGIsPerHdr = 10
    outIDs = []    
    cntFound = 0
    cntHeaders = len(hdrList)
    headerGis = []
    for h in hdrList:
        headerGis.append([])
    c = 0
    print '  --> extracting GIs from IdMapFile <-- '
    with open(inputPath,'r') as idMapper:
        print '    --> opened ',inputPath
        print '    --> searching...'
        potentialRefs = []
        #print 'headers'
        for hdr in hdrList:
            #print hdr
            for refType in hdr.keys():
                for pr in hdr[refType]:
                    potentialRefs.append(pr)
        #d print potentialRefs
        #print ''
        
        for l in idMapper:
            if (c+1) % 5000000 == 0:
                print '      --> ',(c+1),'lines scanned ...'
            l = string.lower(l)
            doSeek = False
            potHit = ''
            for pr in potentialRefs:
                if pr in l:
                    #print 'potential hit: ',pr
                    potHit = pr
                    doSeek = True
                    
            if doSeek:
                idFound = False
                #d print l
                ls = l.split('\t')
                #print hdrList
                for hdr in hdrList:                    
                    for refType in hdr.keys():
                        '''
                         decide on output: 
                         if out='gi', grab GIs
                         if out='taxID', grab TaxIDs
                         '''                                            
                        # strip GIs
                        if out=='gi':
                            outIDs = ls[4].split(';')
                            for i in range(0,len(outIDs)):
                                outIDs[i] = outIDs[i].strip()
                        elif string.lower(out)=='taxid':
                            outIDs = ls[13].split(';')
                            for i in range(0,len(outIDs)):
                                outIDs[i] = outIDs[i].strip()
                        
                        # GI (ls[4]) 
                        if not idFound and refType == 'gi':                         
                            for ids in ls[4].split(';'):
                                ids = ids.strip()
                                for ref in hdr[refType]:
                                    # trim fucker
                                    if not ids.find('.') == -1:
                                        ids = ids[0:ids.find('.')]
                                    #print ' --> trying ref:',ref,' => id',ids                        
                                    if ids == ref:
                                        #print 'found',refType,ref,': returning all',out,str(outIDs)
                                        idFound = True                                    
                                        hc = 0                                         
                                        for h in hdrList:
                                            for k in h.keys():
                                                for hOne in h[k]: 
                                                    if k == refType and hOne == ref:
                                                        if len(outIDs) > maxGIsPerHdr:
                                                            for i in range(0,maxGIsPerHdr):
                                                                headerGis[hc].append(outIDs[i])
                                                        else:
                                                            headerGis[hc] = outIDs
                                                            
                                                        
                                            hc+=1                               
                        # ------ REF (ls[3]) ----------
                        if not idFound and refType == 'ref':
                            for ids in ls[3].split(';'):
                                ids = ids.strip()
                                for ref in hdr[refType]:
                                    # trim fucker
                                    if not ids.find('.') == -1:
                                        ids = ids[0:ids.find('.')]
                                    #print ' --> trying ref:',ref,' => id',ids                                              
                                    #print 'ids:',ids,'ref:',ref                                                               
                                    if ids == ref:                                        
                                        #print 'found',refType,ref,': returning all',out,str(outIDs)
                                        idFound = True                                    
                                        hc = 0                                         
                                        for h in hdrList:
                                            for k in h.keys():
                                                for hOne in h[k]: 
                                                    if k == refType and hOne == ref:
                                                        if len(outIDs) > maxGIsPerHdr:
                                                            for i in range(0,maxGIsPerHdr):
                                                                headerGis[hc].append(outIDs[i])
                                                        else:
                                                            headerGis[hc] = outIDs
                                            hc+=1                               
                        # ------ UNIPROT (sp,tre,uniprot,emb) (ls[0]) ---------
                        if not idFound and (refType == 'sp' or refType == 'uniprot' or refType == 'tre' or refType == 'emb'\
                                         or refType=='tr' or refType=='trembl' or refType =='embl'):
                            for ids in ls[0].split(';'):
                                ids = ids.strip()
                                for ref in hdr[refType]:
                                    # trim fucker
                                    if not ids.find('.') == -1:
                                        ids = ids[0:ids.find('.')]                        
                                    if ids == ref:
                                        #d print 'found',refType,ref,': returning all',out,str(outIDs)
                                        idFound = True                                    
                                        hc = 0                                         
                                        for h in hdrList:
                                            for k in h.keys():
                                                for hOne in h[k]: 
                                                    if k == refType and hOne == ref:
                                                        if len(outIDs) > maxGIsPerHdr:
                                                            for i in range(0,maxGIsPerHdr):
                                                                headerGis[hc].append(outIDs[i])
                                                        else:
                                                            headerGis[hc] = outIDs
                                            hc+=1                               
                        # ------ PDB (pdb) (ls[5]) ---------
                        if not idFound and (refType == 'pdb'):
                            for ids in ls[5].split(';'):
                                ids = ids.strip()
                                for ref in hdr[refType]:
                                    # trim fucker
                                    if not ids.find('.') == -1:
                                        ids = ids[0:ids.find('.')]                        
                                    if ids == ref:
                                        #d print 'found',refType,ref,': returning all',out,str(outIDs)
                                        idFound = True                                    
                                        hc = 0                                         
                                        for h in hdrList:
                                            for k in h.keys():
                                                for hOne in h[k]: 
                                                    if k == refType and hOne == ref:
                                                        if len(outIDs) > maxGIsPerHdr:
                                                            for i in range(0,maxGIsPerHdr):
                                                                headerGis[hc].append(outIDs[i])
                                                        else:
                                                            headerGis[hc] = outIDs
                                            hc+=1                               
                        # ------ GENEBANK (gb) (ls[] ) ---------
                        #debug print ls
                        if not idFound and (refType == 'gb' or refType == 'genebank'):
                            for ids in ls[17].split(';'):
                                ids = ids.strip()
                                for ref in hdr[refType]:    
                                    # trim fucker
                                    if not ids.find('.') == -1:
                                        ids = ids[0:ids.find('.')]                        
                                    if ids == ref:
                                        #d print 'found',refType,ref,': returning all',out,str(outIDs)
                                        idFound = True                                    
                                        hc = 0                                         
                                        for h in hdrList:
                                            for k in h.keys():
                                                for hOne in h[k]: 
                                                    if k == refType and hOne == ref:
                                                        if len(outIDs) > maxGIsPerHdr:
                                                            for i in range(0,maxGIsPerHdr):
                                                                headerGis[hc].append(outIDs[i])
                                                        else:
                                                            headerGis[hc] = outIDs
                                            hc+=1                               

                if not idFound:
                    #print '   --> ',potHit,' not found!'
                    pass                                                        
                if idFound:                                          
#                   print 'header GIs: ',headerGis                            
                    cntFound += 1
                    potentialRefs.remove(potHit)
                    #print '   --> ',potHit,' found! [',cntFound,'/',cntHeaders,']'                    
                    if cntFound == cntHeaders:
                        return headerGis                                  
            c+=1   
#            if c > 10000000: 
#                return headerGis         
        
    #print 'done,cnt=',c        
    return headerGis
   
''' 
loadTaxID returns taxID number from binarized taxonomy
'''
def getSize(fileobject):
    fileobject.seek(0,2) # move the cursor to the end of the file
    size = fileobject.tell()
    return size
def loadTaxID (fileP, gi):
    if gi == -1: 
        print 'warning: something went wrong (got gi = -1), will not retrieve! '
    else: 
        gi = int(gi-1)
        with open(fileP,'rb') as outFB:
            maxP = getSize(outFB)
            if 8*gi < maxP:
                outFB.seek(8*gi)
                out = outFB.read(8)
                #            print out
                return int(out,base=16)    
            else:
                print ' --> ERROR: attempting to reach position > EoF in ',fileP,' gi = ',gi
                return -1
    return -1