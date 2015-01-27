#!/usr/bin/env python

'''
Created on Aug 29, 2013

@author: ranko
'''

''' FUNCTION(s) THAT CREATES D.gtxt files from NCBI taxonomy
.gtxt files are files used to generate .gv (graphviz) files
and are used with taxmapper to merge results with ncbi taxonomy
'''

import parseNCBItaxonomy as pNT
#taxIDNames = {}
#taxMap = {}

class TaxIdDataDirectedT():
    def __init__(self, _startTax,_targetTax):
        self.targetTax = _targetTax
        self.startTax = _startTax
        self.rank = ''        

def prepareNCBITaxForDrawing(namesFile, nodesFile,drawRanks='SGFOCPKD',doTrimSName = 'Y',doDummies=True,verbose=False):
    names = namesFile
    nodes = nodesFile
    print ' ---> LOADING TAXBASE '
    pNT.initTax(names, nodes)         
    print ' ---> DONE '
    print ' ---> GRABBIN ALL --- '
    allTaxes = {}
    c = 0
    ct = 0
    print 'total taxa: ',len(pNT.taxIDNames.keys())
    for n in pNT.taxIDNames.keys():
        t = pNT.getTaxFromTaxID(n)
        c+=1            
        if not t.getRank('species').name == 'N/D':
            ct+=1
            if doDummies:
                rl = t.getTax9EwDummy(drawRanks,doTrimSName)
                #if 'phylum:P:Bacteroidetes -> superkingdom:D:Bacteria' in rl.toStrRow():
                #    print rl.toStrRow()
                #    print t.toStrRow()
                #    print rl.getTax9EwDummy().toStrList()
                #    exit (-1)
            else:
                rl = t.getTax9E(drawRanks,doTrimSName)
            h = rl.toTupleList()
            taxLvlC = 0  
            for n2 in h:
                taxLvlC +=1            
#                if n2[0] in allTaxes.keys():            
#                    pass
#                else:
                allTaxes[n2[0]] = TaxIdDataDirectedT(n2[0],n2[1])
                allTaxes[n2[0]].rank = rl.records[taxLvlC-1].rank                
            if ct % 10000 == 0 and verbose: 
                print 'grabbed',ct,'taxa', '; now have: ',len(allTaxes.keys()),'taxa'                     
        if c % 10000 == 0 and verbose:         
            print 'parsed',c,'taxa'
        elif c % 50000 == 0:
            print 'parsed',c,'taxa'
#        if c > 5000:
#            break
    print '----> DONE ! <----'
    outName = 'taxOut'
    if doDummies:
        outName = 'taxOutD'
        
    with open (outName+'_'+drawRanks+'.gtxt','w') as out:
        for t in allTaxes.keys():
            out.write('"'+allTaxes[t].startTax+'"'+'\n')
        for t in allTaxes.keys():
            out.write('"'+allTaxes[t].startTax+'"'+' -> '+'"'+allTaxes[t].targetTax+'"'+'\n')
 
#    for l in allTaxes.keys():
#        print allTaxes[l].startTax, '->',allTaxes[l].targetTax


''' ------------------------ MAIN ---------------------------- '''       	

''' run it '''
import sys
if len (sys.argv) < 3 or len(sys.argv) > 3:
    print '------------------- USAGE -----------------------'
    print '> python prepareNCBItaxFiles <namefile> <nodefile>             '
    print '> where: '
    print '>      namefile = names.dmp of NCBI taxonomy database'
    print '>      nodefile = nodes.dmp of NCBI taxonomy database'
    print '-------------------------------------------------'
    exit(-1)

nameF = sys.argv[1]
nodeF = sys.argv[2]
'''
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DK',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKP',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPC',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCO',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCOF',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCOFG',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCOFGS',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCF',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPF',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKF',doTrimSName = 'Y',doDummies=True)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKOF',doTrimSName = 'Y',doDummies=True)
'''
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DK',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKP',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPC',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCO',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCOG',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCOGF',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCOGFS',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPC',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCG',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCGF',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPCGFS',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKP',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPG',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPGF',doTrimSName = 'Y',doDummies=False)
prepareNCBITaxForDrawing(nameF,nodeF,drawRanks='DKPGFS',doTrimSName = 'Y',doDummies=False)

