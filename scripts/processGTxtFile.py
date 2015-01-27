'''
Created on 25 Oct 2013

@author: ranko

processes prepareNCBItaxFiles
made file (.gtxt)

prepares it for dot drawing

does things such as:
- adds diffrent sizes for different tax units
- adds different colors for different tax units
- makes 'root' and 'dummy taxa' "Invisible"
- generates lables from node names
- ...
'''
import parseNCBItaxonomy as pNT

def grabChild (_inFile, _tax):
    with open (_inFile) as iF:
        taxatarget = '----'
        for l in iF:
            if _tax in l and '->' in l:
                if _tax in l.split('->')[0]:
                    return l.split('->')[1].strip()
    return 'N/F'

'''
    -> inFile = input file
    -> taxN = name of taxa (ex: baceteria, mammalia...)
    -> taxL = level of taxa (ex: 'D:' for domain)
'''
def findChild(_inFile,_taxN,_taxL,trim=False):
    doGo = True
    while doGo:
#d        print _taxN
        c = grabChild (_inFile, _taxN)
#d        print c
        if trim:
            if not c.find('"') == -1:
                c = c[1:]
            if not c.find('"') == -1:
                c = c[0:c.find('"')]
        #d print c
        if c == 'N/F':
            doGo = False
            return 'N/F'
        if _taxL in c:
            return c
            doGo = False
        _taxN = c
    return 'N/F'

def processGTxTFile(_inFile, _outFile):    
    with open (_inFile) as iF:
        with open (_outFile,'w') as oF:
            parsingNodes = True
            parsingEdges = False
            # header
            p = ''
            oF.write(p+'diGraph ncbiTax'+'\n')
            oF.write(p+'{'+'\n')
            oF.write(p+'root="N/D"'+'\n')
            p = '   '
            oF.write(p+'layout=twopi;'+'\n')
            oF.write(p+'node [shape=plaintext];'+'\n')
            oF.write(p+'edge [style=dashed dir=back];'+'\n')
            oF.write(p+'ranksep=6;'+'\n')
            oF.write('# --- NODES FOLLOW --- '+'\n')
            oF.write(p+'"N/D" [label="" height=0.05 width=0.05 shape=circle]'+'\n')
            # now put data here
            for inLn in iF:
                inLn = inLn.strip()
                if not parsingEdges:
                    if '->' in inLn:
                        parsingNodes = False
                        parsingEdges = True
                        oF.write('# --- EDGES FOLLOW ---'+'\n')
                        
                if parsingNodes:                    
                    if len(inLn) > 1:
                        add = '['
                        # grab dummy node
                        #  dummy node shows no text
                        if '&' in inLn or 'N/D' in inLn:
                            add = add + 'label="" height=0.05 width=0.05 shape=circle'
                        else:
                            # open for options                                       
                                                
                            # options depending on tax rank       
                            add = add + ' fontcolor="#505050"'                                                                                               
                            if 'D:' in inLn:                 
                                add = add + ' shape=circle fontsize=16 '                                
                            elif 'K:' in inLn:
                                add = add + ' shape=ellipse fontsize=14 '                                                                                           
                            elif 'P:' in inLn:
                                add = add + ' fontsize=8 '
                            elif 'C:' in inLn:
                                add = add + ' fontsize=6 '
                            elif 'O:' in inLn:
                                add = add + ' fontsize=4 '
                            elif 'F:' in inLn:
                                add = add + ' fontsize=2 '                            
                            # correction for removal of 'X:'
                            #d print inLn
                            add = add+'label="'+inLn[3:]
                            
                        # close it
                        add = add+']'
                        # add tab                                                        
                        inLn = p+inLn+add
                        oF.write(inLn+'\n')
                    
                if parsingEdges:
                    inLn = p+inLn                    
                    add = '['
                    # color (testing)
#                    print findChild(_inFile,inLn.split('->')[0].strip(),'D:')
                    if 'Eukaryota' in findChild(_inFile,inLn.split('->')[0].strip(),'D:'):
                            add = add + ' color = "#800000" '                                                        
                    elif 'Archaea' in findChild(_inFile,inLn.split('->')[0].strip(),'D:'):
                            add = add + ' color = "#008000" '                                                        
                    elif 'Bacteria' in findChild(_inFile,inLn.split('->')[0].strip(),'D:'):
                            add = add + ' color = "#000080" '                                                        
                    add += ']'
                    inLn = inLn + add
                    oF.write(inLn+'\n')                                                            
                    
            # footer
            p = ''
            oF.write(p+'}'+'\n')

'''
   -> appends GTxt File (prepared NCBI taxonomy file)
      to outFileResults (file generated by taxparser)
   -> basically adds NCBI tax stuff unless it already
       exists

'''    
            
            
def grabChildHash (hesh, _tax):
    if _tax in hesh.keys():
        return hesh[_tax]
    else:
        return 'N/F'
'''
    -> inFile = input file
    -> taxN = name of taxa (ex: baceteria, mammalia...)
    -> taxL = level of taxa (ex: 'D:' for domain)
'''
def findChildHash(hesh,_taxN,_taxL):
    doGo = True
    while doGo:
#d      print _taxN
        c = grabChildHash (hesh, _taxN)
#d      print c
        if c == 'N/F':
            doGo = False
            return 'N/F'
        if _taxL in c:
            return c
            doGo = False
        _taxN = c
    return 'N/F' 
               
def appendGTxTFile(_inFile, _inFileR, _outFile, feD='DKPC'):
    ''' PREPARSER '''
    nodesR = set()
    edgesR = {}
    nodesT = set()
    edgesT = {}
    with open (_inFile) as iF:
        addTaxNodes = True
        parsingEdges = False
        for l in iF:
            if not parsingEdges:
                if '->' in l:
                    addTaxNodes = False
                    parsingEdges = True                                 
            if addTaxNodes:
                #d print 'N:',l.strip()
                nodesT.add(l.strip().replace('"',''))
            elif parsingEdges:
                s = l.strip().replace('"','').split('->')[0].strip()
                r = l.strip().replace('"','').split('->')[1].strip()
                edgesT[s] = r
                
    with open (_inFileR) as iF:                
        addTaxNodes = False
        parsingEdges = False
        for l in iF:
            if parsingEdges and '}' in l:
                parsingEdges = False
            if parsingEdges:
                if len(l.split('->')) > 0:
                    l = l.strip()
                    #d print l
                    s = l.strip().split('->')[0].strip()                    
                    s = s[1:s[1:].find('"')+1]
                    r = l.strip().split('->')[1].strip()                    
                    r = r[1:r[1:].find('"')+1]
                    edgesR[s] = r                       
            if not parsingEdges:
                if 'EDGES FOLLOW' in l:
                    parsingEdges = True
                    addTaxNodes = False  
            if addTaxNodes:
                l = l.strip()
                l = l[1:]
                l = l[0:l.find('"')]
#d                print l
                nodesR.add(l)            
            if not addTaxNodes:
                if 'NODES FOLLOW' in l:
                    addTaxNodes = True
                    
# --------------- DONE PREPARSING ---------------                    
#    print nodesT
#    print edgesT
#    print nodesR
#    print edgesR
    
    #grabChildHash(edgesR,'S:F. graminearu')
#    findChildHash(edgesR,'O:Tremellales','D:')
        
    dExists = set()
    kExists = set()
    pExists = set()
    cExists = set()
    oExists = set()
    fExists = set()
    gExists = set()
    
    # generate TARGETS (we check against these 
    # when we decide what to actualy draw from taxonomy
    
    for r in edgesR: 
        if 'D' in feD:
            child = findChildHash(edgesR, r, 'D:')
            if not child == 'N/D' and not child == 'N/F':
                dExists.add(child)              
        if 'K' in feD:
            child = findChildHash(edgesR, r, 'K:')
            if not child == 'N/D' and not child == 'N/F':
                kExists.add(child)              
        if 'P' in feD:
            child = findChildHash(edgesR, r, 'P:')
            if not child == 'N/D' and not child == 'N/F':
                pExists.add(child)              
        if 'C' in feD:
            child = findChildHash(edgesR, r, 'C:')
            if not child == 'N/D' and not child == 'N/F':
                cExists.add(child)              
        if 'O' in feD:
            child = findChildHash(edgesR, r, 'O:')
            if not child == 'N/D' and not child == 'N/F':
                oExists.add(child)              
        if 'F' in feD:
            child = findChildHash(edgesR, r, 'F:')
            if not child == 'N/D' and not child == 'N/F':
                fExists.add(child)              
        if 'G' in feD:
            child = findChildHash(edgesR, r, 'G:')
            if not child == 'N/D' and not child == 'N/F':
                gExists.add(child)              
    '''
    print dExists
    print kExists
    print pExists
    print cExists
    print oExists
    print fExists    
    print gExists
    '''
    print ' --> WRITING MERGED FILE '
    with open (_inFileR) as resLines:
        with open (_outFile,'w') as oF:            
            addTaxNodes = False
            parsingEdges = False
            grabResNodes = True
            grabResEdge = False
            for l in resLines:
                # grab stuff before --- NODES FOLLOW ---
                if grabResNodes:
                    if 'EDGES FOLLOW' in l:
                        grabResNodes = False
                        addTaxNodes = True   
                        grabResEdge = True                        
                        print ' ---> ADDING TAX NODES'
                        oF.write('# ---> TAXONOMY NODES FOLLOW <--- \n')                                                                         
                # now add nodes from taxonomy file
                # (we already have those hashed
                # so just add them, but before adding, 
                # compare them against: 
                # 1) results (hash), if results has 
                #      this node, skip it
                # 2) if we do targeted mapping, 
                #      check if results have child equal to 
                #      appropriate taxa    
                        if addTaxNodes:                    
                            p = '    '                          
                            for taxNode in nodesT:                        
                                doWriteNode = True
                                # check for target(s)
                                # note: go down -> up                            
                                if 'G' in feD:
                                    if not findChildHash(edgesT, taxNode, 'G:') in gExists:
                                        doWriteNode = False
                                elif 'F' in feD:
                                    if not findChildHash(edgesT, taxNode, 'F:') in fExists:
                                        doWriteNode = False                        
                                elif 'O' in feD:
                                    if not findChildHash(edgesT, taxNode, 'O:') in oExists:
                                        doWriteNode = False                        
                                elif 'C' in feD:
                                    if not findChildHash(edgesT, taxNode, 'C:') in cExists:
                                        doWriteNode = False                        
                                elif 'P' in feD:
                                    if not findChildHash(edgesT, taxNode, 'P:') in pExists:
                                        doWriteNode = False                        
                                elif 'K' in feD:
                                    if not findChildHash(edgesT, taxNode, 'K:') in kExists:
                                        doWriteNode = False                        
                                elif 'D' in feD:
                                    if not findChildHash(edgesT, taxNode, 'D:') in dExists:
                                        doWriteNode = False
                                
                                if taxNode in nodesR:                                                                                                 
                                    doWriteNode = False
                                    #d print taxNode,' already in results hash!'  
                                                                                                  
                                if doWriteNode:
                                    add = '['
                                    # grab dummy node
                                    #  dummy node shows no text
                                    if '&' in taxNode or 'N/D' in taxNode:
                                        add = add + 'label="" height=0.05 width=0.05 shape=circle'
                                    else:
                                        # options depending on tax rank       
                                        add = add + ' fontcolor="#505050"'                                                                                               
                                        if 'D:' in taxNode:                 
                                            add = add + ' shape=circle fontsize=16 '                                
                                        elif 'K:' in taxNode:
                                            add = add + ' shape=ellipse fontsize=14 '                                                                                           
                                        elif 'P:' in taxNode:
                                            add = add + ' fontsize=8 '
                                        elif 'C:' in taxNode:
                                            add = add + ' fontsize=6 '
                                        elif 'O:' in taxNode:
                                            add = add + ' fontsize=4 '
                                        elif 'F:' in taxNode:
                                            add = add + ' fontsize=4 '                            
                                        elif 'G:' in taxNode:
                                            add = add + ' fontsize=4 '                            
                                        elif 'S:' in taxNode:
                                            add = add + ' fontsize=4 '                            
                                        # correction for removal of 'X:'
                                        #d print inLn
                                        #add = add+'label="'+inLn[2:]                            
                                        # close it
                                    add = add+']'
                                    # add tab                                                        
                                    inLn = p+'"'+taxNode+'"'+add
                                    oF.write(inLn+'\n')
                                    #print inLn
                    # WRITE LINE (NODES)
                    else:                               
                        oF.write(l)                    
                    # add original (result) edges  
                if grabResEdge:                        
                    # end reached                    
                    if l[0] == '}':
                        print ' ---> ADDING TAX EDGES'
                        oF.write('# ---> TAXONOMY EDGES FOLLOW <--- '+'\n')
                        # go through all the possible edges from tax file
                        # if it exists in results (start and end of edge are the same)
                        # don't write it
                        # if it doesn't exist in targets, don't write it either                        
                        for teS in edgesT.keys():
                            teE = edgesT[teS]
                            writeEdge = True
                            # check for result edges
                            if teS in edgesR.keys():                                
                                reE = edgesR[teS] 
                                if reE == teE:
                                    writeEdge = False
                            # check for target(s)
                            if 'G' in feD:
                                if not findChildHash(edgesT, teS, 'G:') in gExists:
                                    writeEdge = False
                            if 'F' in feD:
                                if not findChildHash(edgesT, teS, 'F:') in fExists:
                                    writeEdge = False
                            if 'O' in feD:
                                if not findChildHash(edgesT, teS, 'O:') in oExists:
                                    writeEdge = False
                            if 'C' in feD:
                                if not findChildHash(edgesT, teS, 'C:') in cExists:
                                    writeEdge = False
                            if 'P' in feD:
                                if not findChildHash(edgesT, teS, 'P:') in pExists:
                                    writeEdge = False
                            if 'K' in feD:
                                if not findChildHash(edgesT, teS, 'K:') in kExists:
                                    writeEdge = False
                            if 'D' in feD:
                                if not findChildHash(edgesT, teS, 'D:') in dExists:
                                    writeEdge = False    
                                                                                                                                                                                        
                            if writeEdge:
                                #d print 'writing: ',inLn                                
                                add = '['
                                # color (testing)
                                if 'Eukaryota' in findChildHash(edgesT, teS, 'D:'):
                                    add = add + ' color = "#800000" '                                                        
                                elif 'Archaea' in findChildHash(edgesT, teS, 'D:'):
                                    add = add + ' color = "#008000" '                                                        
                                elif 'Bacteria' in findChildHash(edgesT, teS, 'D:'):
                                    add = add + ' color = "#000080" '                                                        
                                add += ']'
                                outputTxt = p+'"'+teS+'" -> "'+teE+'" ' + add
                                oF.write(outputTxt+'\n')                                                  
#                    else:
                    oF.write(l)
                    #p = ''
                    #oF.write(p+'}'+'\n')                                                                                               
            # footer
            print ' ---> DONE! <---'


#inFileTax = '/home/ranko/Development/workspace/TaxMapper/src/taxFiles/taxOutD_DKPC.gtxt'
#inFileR = '/home/ranko/Development/workspace/TaxMapper/src/testBlast_out_graph.gv'
#outFileP = '/home/ranko/Development/workspace/TaxMapper/src/TESTMERGED.gv'
#print ' ---> PROCESSING <---- '
#appendGTxTFile(inFileTax, inFileR, outFileP)
'''
inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOutD_DKP.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOutD_DKPC.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOutD_DKPCO.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOutD_DKPCOF.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOut_DKP.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOut_DKPC.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOut_DKPCO.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

inPath = '/home/ranko/Development/workspace/TaxMapper/src/'
inFileTax = 'taxOut_DKPCOF.gtxt'
outFileTax = '__gviz_'+inFileTax[0:-4]+'gv'
processGTxTFile(inPath+inFileTax,inPath+outFileTax)

print ' ----> DONE ! <---- '
'''

