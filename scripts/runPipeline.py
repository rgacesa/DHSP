#!/usr/bin/env python
# RUNS PIPELINE FOR DISTANT HOMOLOGY SEARCH !

# DEPENDENCIES:
#  LIBS 
#  -> sh subprocess handler (see http://amoffat.github.io/sh/)
#    -> installation: (sudo) pip install sh
#  -> biopython (mainly for BLAST parser)
#  SCRIPTS (LOCAL/INCLUDED)
#  -> following scripts (included): 
#     -> TLM:
#         -> for details how to set up and use, see TaxMapper documentation
#		  -> will not work unless set up properly and binarized and prepared NCBI taxonomy files are generated
#  3RD PARTY          
#  -> following 3rd party programs:
#     -> NCBI BLAST+ (blastp and psiblast are used)
#     -> HMMER (phmmer, jackhmmer and hmmsearch are used)
#     -> HHsuite (hhblits is used); this one can be tricky to install and has lot of dependencies of its own
#        -> check its documentation at http://toolkit.tuebingen.mpg.de/hhpred and ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/hhsuite-userguide.pdf

# HOW TO USE:
#  CONFIG
#   -> make sure paths are correct (see paths below)
#   -> optionally configure parameters (see below)
#
#  input   
#   -> put sequence of interest into ./input/ORIG_SEQ
#   -> put its _KNOWN_ (ideally well characterized by EXPERIMENTAL MEANS)
#   homologs into ./input/H_SEQ (those are optional but will make prediction
#   better by allowing software to cross reference hits by multiple models
#   based on well known, well characterized input
#   NOTE: 
#     -> ideally these sequences are experimentally characterized
#        homologs of target (further refered to as original) as evolutionally 
#        distant as possible, pick between 2 and 10 of such sequences which are
#        not very evolutionally close
#
#  START
#   -> start this script
#   run will take a while (the more homologs in H_SEQ the longer it takes)
#   NOTE: 
#     --> script should be run from scripts folder locally (so go there and start it)
#         otherwise paths will get mixed up and it will die!

# ---------- FOLDERS AND STRUCTURE -------------
# make sure following folders exist:
#./input
#./input/ORIG_SEQ                # original input we want to find homologues of
#./input/H_SEQ                   # known homologues of original used for model generation
#./scripts                       # python scripts and all other scripts (basically code of pipeline is here)
#./results
#./results/COMPARATIVE_results   # this is INITIAL set of results for assessment of BLAST vs vs PSI-BLAST vs HMMER vs JACKHMMER vs HHPRED 
#./results/BLAST                 # initial BLAST results
#./results/F_BLAST               # filtered BLAST results
#./results/M_ALIGNMENT           # multiple alignments
#./results/HMM_PROFILES          # HMM profiles
#./results/HMM_results           # HMMER results
#./results/HHpred_results        # HHpred results
#./results/POTENTIAL_HITS        # filtered and cross referenced HMMER results; also waterman filtered ones
#./results/DH_HMM_results        # final predicted homologues for HMMER prediction path
#./results/DH_HH_results         # final predcited homologues for HHblits prediction path
 
# ------- CONFIG ------------

# --- PATHS ----
# --> PATHS HAVE TO BE EDITED OR CODE WILL NOT WORK!
#
# clustalomega path/command:
CLUSTAL_O='clustalo'
# hmm search (profile vs seq DB)
HMMSEARCH='/home/ranko/Development/tools/hmmer/binaries/hmmsearch'
# hmm search (profile vs seq DB)
HMMBUILD='/home/ranko/Development/tools/hmmer/binaries/hmmbuild'
# phmmer (seq vs seq DB)
PHMMER='/home/ranko/Development/tools/hmmer/binaries/phmmer'
# jackhmmer (seq vs seq DB)
JACKHMMER='/home/ranko/Development/tools/hmmer/binaries/jackhmmer'
# BLAST
BLAST='/home/ranko/Development/tools/blast+/blastp'
# BLAST DB CMD
BLAST_DB_CMD='/home/ranko/Development/tools/blast+/blastdbcmd'
# PSI-BLAST
PSIBLAST='/home/ranko/Development/tools/blast+/psiblast'
# HH-BLITS
HHBLITS='/home/ranko/Development/tools/HHsuite/bin/hhblits'
# NR TAXONOMY BASE FOLDER
NRTAX='/home/ranko/Development/DBs/ncbi_tax'
# NR DB
NRDB='/home/ranko/Development/DBs/nr_fasta/nr.fa'
# WATERMAN (usually installed and put into path)
WATERMAN='water'




# -----------------------------------------------
# ------------- THESE SHOULD BE FINE ------------
# -----------------------------------------------
# ----- WATERMAN CONFIG -----
waterMaxResults = 1000				# pick best X results
# NR TAXONOMY NAMES DUMP; don't edit unless really needed
NRTAXNAMES=NRTAX+'/names.dmp'
# NR TAXONOMY NODES DUMP; don't edit unless really needed
NRTAXNODES=NRTAX+'/nodes.dmp'
# _BINARIZED_ NR Taxonomy; don't edit unless really needed
NRTAXBINARIZED=NRTAX+'/gi_taxid_prot.bin'
# ------ CONFIG PATHS ----------
# this should be fine, dont change unless really needed
TAXMAPPER='./TLM.py'
RES_BLAST='../results/BLAST'
RES_F_BLAST='../results/F_BLAST'
RES_M_ALIGNMENT='../results/M_ALIGNMENT'
RES_HMMS='../results/HMM_PROFILES'
RES_HMMSEARCH='../results/HMM_RESULTS'
RES_POTENTIALHITS = '../results/POTENTIAL_HITS'
RES_COMPARATIVE='../results/COMPARATIVE_RESULTS'
# --- OTHER CONFIG ----

# DETAILS (OR WHAT THE PIPELINE DOES)
# BASICALLY DOES FOLLOWING:
# 0) runs BLAST, PSI-BLAST, HMMER, JACK-HMMER, HHblits with ORIGINAL SEQ
# vs target DB (entered as command line parameter)
# -> puts this one into COMPARATIVE_results
# Note: 
#  -> this is generally low quality sweep mainly here for finding
#  what conventional tools will spit out
#
#  --- analysis itself is as follows ---
# 1) runs python code for generating BLAST results from selected sequences
#    -> original sequence (which is our target) 
#    -> its known homologs (pick 2 - 5 known distant homologs (in data/
#    
#    -> NOTE: generates results for each input sequence in data/SEQS folder
#    -> puts them into BLAST folder
# 1a)    blast results are parsed and filtered and stored into F_BLAST folder
#        -> each file (*.fa) contains set of fastas for filtered blast results 
# 2) generates multiple alignments for each set of filtered BLAST results
# 3) generates HMMER model of each multiple alignment
# 4) does HMMER search against target DB with each model
# 5) filters results as following: 
#  -> sequences with eV better then X are considered (default 0.1)
#  -> sequences hit by at least Y models are considered (default Y = number of models/2)
#      -> these are then processed by smith-waterman against ORIGINAL and only those with match > 20%
#         are considered for 3d modeling (20% is threshold for similar fold)

# ------------ CODE ------------------------------------
# DO NOT MODIFY UNLESS YOU HAVE GOOD REASON TO !!!
# ------------------------------------------------------
# --- IMPORTS ----
import argparse
import sh
import math
import os
import csv
import distutils.core
import copy
from datetime import datetime

from Bio.Blast import NCBIXML
from RBioTools import extractBlastTitleGiTitleMap
from RBioTools import extractHMMERTitleGiTitleMap
from RBioTools import extractAllGIs
from RBioTools import HmmResultRecord
from RBioTools import fetchFastaFromBlastDB
from RBioTools import doWaterman
from RBioTools import getFastaLt
from RBioTools import loadAllFastas 
from RBioTools import getFaSeq 
from RBioTools import getFaHeader 
from RBioTools import addEvToFaHdr 
from parseHHblitz import parseHHResults
from parseHHblitz import filterHHblits
from parseHMMER import parseHMMResults

# ------------------------------------------------------
# CLASSES, FUNCTIONS AND STUFF
# ------------------------------------------------------
#comparative result hit
#basically holds description of which tools hit this sequence and 
#how well; 
#model: 
#   gi -> name (title) -> blast eV -> psiblast eV -> hmmer eV -> jack eV -> HHblits eV
#

def clearFolder(fPath):
	folder = fPath
	for the_file in os.listdir(folder):
		file_path = os.path.join(folder, the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception, e:
			print e

class CompRes():
	def __init__(self):
		self.gi = -1
		self.title = ''
		self.blastEV = -1.0
		self.psiblastEV = -1.0
		self.hmmerEV = -1.0
		self.jackEV = -1.0
		self.HHblitsEV = -1.0
		self.toolsHits = 0
		self.score = 0.0
	def toStr(self):
		ret = str(self.gi).strip().ljust(12)+'\t'+self.title.strip()+'\t'+'{:.1e}'.format(self.blastEV)+'\t'+'{:.1e}'.format(self.psiblastEV)+'\t'+ '{:.1e}'.format(self.hmmerEV)+'\t'+'{:.1e}'.format(self.jackEV)+'\t'+'{:.1e}'.format(self.HHblitsEV)+'\n'
		ret = ret.replace('-1.0e+00','   X   ')
		return ret 
	def toStrS(self):
		ret = str(self.gi).strip().ljust(12)+'\t'+self.title.strip()[:30].ljust(30)+'\t'+'{:.1e}'.format(self.blastEV)+'\t'+'{:.1e}'.format(self.psiblastEV)+'\t'+ '{:.1e}'.format(self.hmmerEV)+'\t'+'{:.1e}'.format(self.jackEV)+'\t'+'{:.1e}'.format(self.HHblitsEV)+'\n'
		ret = ret.replace('-1.0e+00','   X   ')
		return ret 
	def getTitle(self):
		ret = '--- GI ---    \t----------- TITLE ------------ \t--BL --\t--PSI--\t--HMM--\t-JACK -\t--HHb--'
		return ret

# ------------------------------------------------------
# COMMAND LINE PROCESSING FOLLOWS
# ------------------------------------------------------
parser = argparse.ArgumentParser(description='Distant Homology Search Pipeline')
parser.add_argument('-I','--input',nargs='?', required=True, help='input "Original" sequence, path included; should be in input/ORIG_SEQ folder (or code might get angry)')
parser.add_argument('-D','--database',nargs='?', required=True, help='target database, including path to database; keep in mind BLAST db of same name should exist')
parser.add_argument('-H','--useHHblits',nargs='?', default='Y', help='if Y, uses hhblits in comparative results and in distant homology search, otherwise does not!')
parser.add_argument('--hhfilter',nargs='?', default='none',choices=['none','bacteria','archaea','fungi','protista','misc','metazoa','plants'], required=False, help='MANDATORY if -H is Y: hhfilter')
parser.add_argument('--hhdb',nargs='?', default='X',required=False, help='MANDATORY if -H is Y: location of hh database (nr20 db prebuilt by hhsuite)')
parser.add_argument('--doComparativeRun',nargs='?', default='Y', help='if Y, does comparative run for various tools')
parser.add_argument('--threads',nargs='?', type=int, default=4, help='[DEF = 4] number of threads to run per comparative analysis run tool')
# WATERMAN (FINAL FILTER OPTIONS)
parser.add_argument('--wMinLt',nargs='?', type=float, default=0.4, help='[DEF = 0.4] min lt of potential homolog (as ratio of original)')
parser.add_argument('--wMaxLt',nargs='?', type=float, default=2.5, help='[DEF = 2.5] max lt of potential homolog (as ratio of original)')
parser.add_argument('--wMinSim',nargs='?', type=float, default=0.5, help='[DEF = 0.5] minimum waterman coverage/similarity of potential homologs')
parser.add_argument('--wMaxRes',nargs='?', type=float, default=50000, help='[DEF = 50000] maximum results for waterman search')
# PSI BLAST (MODEL BUILDING OPTIONS)
parser.add_argument('--pbEVMax',nargs='?', type=float, default=1.0e-3, help='[DEF = 1.0e-3] max e-Value of local search to consider for model building')
parser.add_argument('--pbCovMin',nargs='?', type=float, default=0.6, help='[DEF = 0.6] min local search coverage to consider for model building')
parser.add_argument('--pbLtMin',nargs='?', type=int, default=0, help='[DEF = 0] minimum local search sequence length (in aa) to consider for model building')
parser.add_argument('--pbLtMax',nargs='?', type=int, default=2000, help='[DEF = 2000] maximum local search sequence length (in aa) to consider for model building')
parser.add_argument('--pbAlLtMin',nargs='?', type=float, default=10, help='[DEF = 10] local search vs original alignment length (in aa) to consider for model building')
parser.add_argument('--pbIterNR',nargs='?', type=int, default=3, help='[DEF = 3] number of local search psi-Blast iterations')
parser.add_argument('--pbMaxSeqs',nargs='?', type=int, default=500, help='[DEF = 500] max numbder of local search sequences to put into model')
# HMMER Search Options
parser.add_argument('--hmrMaxRes',nargs='?', type=float, default=50000, help='[DEF = 50000] number of hmmer search hits')
parser.add_argument('--hmrAvgEV',nargs='?', type=float, default=100, help='[DEF = 100] max avg e-Value of hmmer search')
parser.add_argument('--hmrMinModels',nargs='?', type=float, default=0.35, help='[DEF = 0.35] min number of models seq has to classify for in hmmer search (as ratio of total models)')
parser.add_argument('--hmrMaxEV',nargs='?', type=float, default=100, help='[DEF = 100] hmmer maximum e-Value for seq to be considereds (for each model)')

print " ***** DISTANT HOMOLOGY SEARCH PIPELINE STARTED  ***** "
args = parser.parse_args() 
tarDB = args.database
useHH = False
if args.useHHblits == 'Y':
	useHH = True
	if args.hhdb == 'X': 
		' ERROR: --hhdb is required! '
		exit(-1)
 
# ------ CONFIG: MODEL PREPARATION FILTERS ------
parseEVmax = args.pbEVMax	    # MAX eV (psi-blast)
minCoverage = args.pbCovMin		# 1.0 = 100% (psi-blast)
parseBSmin = 0			# min bitscore
ltMin = args.pbLtMin				# min LT of target
ltMax = args.pbLtMax			# max LT of target
minQL = args.pbAlLtMin 				# min size of alignment
psiIteration = -1		# always grab last iteration if -1
psiIterationToDo = args.pbIterNR
modelResultsMax = args.pbMaxSeqs
# ----- CONFIG: HMMER RESULT vs TARGET DB FILTERS ------
hrMaxEV = args.hmrMaxEV						# MAX EV to consider hit
hrMinHitsForSameTarget = args.hmrMinModels		# min hits per targets to classify as potential homolog (as ratio of models)
hrMaxEVavg = args.hmrAvgEV				    # MAX average EV of hits to consider potential homologs
hrMaxResNR = args.hmrMaxRes					# max number of results to consider for model hit
# WATERMAN OPTIONS (LOADED FROM CL)
waterMinLt = args.wMinLt				    
waterMaxLt = args.wMaxLt				    
waterMinSim = args.wMinSim  
waterMaxResults = args.wMaxRes 
              
# input ...
inputNO = args.input[args.input.rfind("/")+1:]
inputNOChomp = inputNO.replace('.fa','').replace('.faa','').replace('.fasta','')
dbNameChomp = ''
if '/' in tarDB:
	dbNameChomp = tarDB[tarDB.rfind("/")+1:tarDB.rfind(".")]
elif '\\' in tarDB:
	dbNameChomp = tarDB[tarDB.rfind("\\")+1:tarDB.rfind(".")]
else:
	dbNameChomp = tarDB[:tarDB.rfind(".")]		

_inNames = args.hhfilter
inputSeq = args.input
outP = "../results/COMPARATIVE_RESULTS/"
print " -> target database = ", tarDB
print " -> input (Original) sequence = ", inputSeq
doCA = False
if args.doComparativeRun == 'Y':
	doCA = True
	
# *************************************************************************************
# CLEANUP
# *************************************************************************************
clearFolder(RES_BLAST)
clearFolder(RES_F_BLAST)
clearFolder(RES_M_ALIGNMENT)
clearFolder(RES_HMMS)
clearFolder(RES_HMMSEARCH)
clearFolder(RES_POTENTIALHITS)
clearFolder(RES_COMPARATIVE)
# *************************************************************************************
# *************************************************************************************
# COMPARATIVE ANALYSIS CODE FOLLOWS
# NOTE: analysis runs parallel and picks up stuff after all is done
# *************************************************************************************
# *************************************************************************************
if doCA:
	threadsPerJob = max(1,args.threads/8)
	tOutSeq = outP+ inputNO.replace(".fasta","").replace(".faa","").replace(".fa","")
	print " ---------- PERFORMING COMPARATIVE ANALYSIS ----------- "
	# RUN BLAST
	print " ---> RUNNING BLAST vs ", tarDB
	pblast = sh.Command(BLAST)
	procBlast = pblast("-query",inputSeq,"-db",tarDB,"-out",tOutSeq+"_"+dbNameChomp+"_blastR.xml","-parse_deflines","-num_threads",threadsPerJob,"-show_gis","-outfmt","5",_bg=True)
	# RUN PSI BLAST	
	print " ---> RUNNING PSI-BLAST vs ", tarDB
	psiblast = sh.Command(PSIBLAST)
	procPsiBlast = psiblast("-query",inputSeq,"-db",tarDB,"-out",tOutSeq+"_"+dbNameChomp+"_psiblastR.xml","-parse_deflines","-num_threads",threadsPerJob,"-show_gis","-outfmt","5",_bg=True)
	# RUN HMMER (pHMMER)		
	print " ---> RUNNING HMMER (phmmer) vs ", tarDB
	pHmmer = sh.Command(PHMMER)
	procPHmmer = pHmmer("-o",tOutSeq+"_"+dbNameChomp+"_phres.hmr","--notextw","--cpu",threadsPerJob, inputSeq, tarDB, _bg=True)
	# RUN JACK		
	print " ---> RUNNING Jack-HMMER vs ", tarDB
	jHmmer = sh.Command(JACKHMMER)
	procJHmmer = jHmmer("-o",tOutSeq+"_"+dbNameChomp+"_jhres.jhr","--notextw","-N","3","--cpu",threadsPerJob, inputSeq, tarDB, _bg=True)
	# RUN HHblits
	if useHH:
		print " ---> RUNNING HHblits vs ", tarDB
		hhblits = sh.Command(HHBLITS)
		procHH = hhblits("-o",tOutSeq+"_"+dbNameChomp+"_hhres.hhr","-i",inputSeq,"-d",args.hhdb, "-cpu", max(1,args.threads-args.threads/8*4), "-Z","10000","-B","10000","-aliw","10000","-E","10.0","-p","20", _bg=True)					
	# PICK UP PROCESSES AS THE FINISH
	procBlast.wait()
	print "    --> BLAST DONE! "
	procPsiBlast.wait()
	print "    --> PSI BLAST DONE! "
	procPHmmer.wait()
	print "    --> HMMER DONE! "
	procJHmmer.wait()
	print "    --> JACK DONE! "
	if useHH:
		procHH.wait()			
		print "    --> HHblits DONE! "
	# all is done
	print " ---> ALL JOBS DONE!"
	# ---------------------------------------------
	# ------------------- PARSERS -----------------
	# --------------------------------------------- 	
	parseEVmax = 10
	parseBSmin = 0
	ltMax = 10000
	ltMin = 0
	minQL = 0		
	print " ---> PROCESSING RESULTS "
	# comparative results are stored in gi -> compRes format
	# --------------------------------------		
	# BLAST
	# --------------------------------------
	allCompResB = {}	
	bResult = tOutSeq+"_"+dbNameChomp+ "_blastR.xml" 
	#print "     --> BLAST (",bResult,")"
	cnt=0
	with open(bResult) as resultHandle:  
		br = NCBIXML.parse(resultHandle) 
		for b in br:
			for al in b.alignments:
				for hsp in al.hsps:
					cnt+=1					
					gitm = extractBlastTitleGiTitleMap(al.title)
					for k in gitm.keys():
						newHit = CompRes()
						newHit.gi = k						
						newHit.title = gitm[k]
						adder = False
						newTitle = ''
						for a in newHit.title:
							if a == '|' or a == '>':
								adder = False
							if adder == False and a == ' ':
								adder = True
							if adder:
								newTitle += a
						newHit.title = newTitle 											
						newHit.blastEV = hsp.expect  						
						allCompResB[newHit.gi] = newHit
	#print CompRes().getTitle()
	#for k in allCompResB.keys():
	#	print allCompResB[k].toStrS().strip()
		
	# --------------------------------------		
	# PSI BLAST
	# --------------------------------------
	allCompResPB = {}	
	pbResult = tOutSeq+"_"+dbNameChomp+"_psiblastR.xml" 
	#print "     --> PSI-BLAST (",pbResult,")"
	cnt=0
	maxIter = 0
	psiIter = 0
	with open(pbResult) as resultHandle:  
		br = NCBIXML.parse(resultHandle) 
		for b in br:
			maxIter+=1
	
	with open(pbResult) as resultHandle:  
		br = NCBIXML.parse(resultHandle) 
		for b in br:
			psiIter+=1 
			for al in b.alignments:
				for hsp in al.hsps:
					if psiIter == maxIter:
						cnt+=1
						gitm = extractBlastTitleGiTitleMap (al.title)
						for k in gitm.keys():
							newHit = CompRes()
							newHit.gi = k
							newHit.title = gitm[k]
							adder = False
							newTitle = ''							
							for a in newHit.title:
								if a == '|' or a == '>':
									adder = False
								if adder == False and a == ' ':
									adder = True
								if adder:
									newTitle += a
							newHit.title = newTitle 							
							newHit.psiblastEV = hsp.expect  						
							allCompResPB[newHit.gi] = newHit
	#print CompRes().getTitle()
	#for k in allCompResPB.keys():
	#	print allCompResPB[k].toStrS().strip()
		
	# --------------------------------------		
	# HMMER
	# --------------------------------------
	allCompResH = {}		
	phResult = tOutSeq+"_"+dbNameChomp+"_phres.hmr" 
	#print "     --> HMMER (",phResult,")"
	hmmAllResults = parseHMMResults(phResult)
	# go over all results in last iteration (only one in case of non JACK)			 
	hmmQR = hmmAllResults.qry[0]
	hmmQI = hmmQR.getLastIteration()
	cnt = 0	
	for h in hmmQI.hmmTable:
		cnt += 1
		#print '**** HIT # '+str(cnt)+' ****'
		#print 'sequence:','eV:', h.fullSeqEV, h.seqID,
		#print h.description		
		gitm = extractHMMERTitleGiTitleMap(h.description)
		#print gitm
		for k in gitm.keys():
			newHit = CompRes()
			#print k,'->',gitm[k][0:30]
			newHit.gi = k
			newHit.title = gitm[k]
			newHit.hmmerEV = h.fullSeqEV  						
			allCompResH[newHit.gi] = newHit
					
	#print CompRes().getTitle()
	#for k in allCompResH.keys():
	#	print allCompResH[k].toStrS().strip()
	# --------------------------------------   		
	# J-HMMER
	# --------------------------------------
	allCompResJH = {}		
	phResult = tOutSeq+"_"+dbNameChomp+"_jhres.jhr" 
	#print "     --> JACK HMMER (",phResult,")"
	hmmAllResults = parseHMMResults(phResult)
	# go over all results in last iteration (only one in case of non JACK)			 
	hmmQR = hmmAllResults.qry[0]
	hmmQI = hmmQR.getLastIteration()
	cnt = 0	
	for h in hmmQI.hmmTable:
		cnt += 1
		#print '**** HIT # '+str(cnt)+' ****'
		#print 'sequence:','eV:', h.fullSeqEV, h.seqID,
		#print h.description		
		gitm = extractHMMERTitleGiTitleMap(h.description)
		#print gitm
		for k in gitm.keys():
			newHit = CompRes()
			#print k,'->',gitm[k][0:30]
			newHit.gi = k
			newHit.title = gitm[k]
			newHit.jackEV = h.fullSeqEV  						
			allCompResJH[newHit.gi] = newHit
					
	#print CompRes().getTitle()
	#for k in allCompResJH.keys():
	#	print allCompResJH[k].toStrS().strip()
	# --------------------------------------
	# HH (requires additional filtering thingy!)
	# --------------------------------------	
	if useHH:
		allCompResHHb = {}
		#print "     --> HHBLITS (",phResult,")"	
		hhResult = tOutSeq+"_"+dbNameChomp+"_hhres.hhr" 
		#print "     	--> parsing "	
		res = parseHHResults(hhResult)
		#print "         --> filtering ... "
		resF = filterHHblits(res, NRTAXBINARIZED, NRTAXNAMES, NRTAXNODES, targetTaxNames=_inNames)
		for hhRnew in resF.hhResults: 
#			print hhRnew.nr,' ',hhRnew.hHdrShort,'p',hhRnew.prob,'ev:',hhRnew.eV,'score:',hhRnew.rawScore,'ss:',hhRnew.pripredScore, 'm:',hhRnew.nrMatches, 'qry:',hhRnew.qHMMStart,'-',hhRnew.qHMMEnd, 'tar:',hhRnew.hHMMStart,'-',hhRnew.hHMMEnd, 'I:'
#			print 'long header: ', hhRnew.hitHdr
			for k in hhRnew.gis:
				newHit = CompRes()
				newHit.gi = k
				newHit.title = hhRnew.hitHdr
				newHit.HHblitsEV = hhRnew.eV
				# chump on title to remove gis:
				newTitle = ''
				adder = True
				for a in newHit.title:
					if a == '|' or a == '>':
						adder = False
					if adder == False and a == ' ':
						adder = True
					if adder:
						newTitle += a
				newHit.title = newTitle 
				allCompResHHb[k] = newHit
		#print CompRes().getTitle()
		#for k in allCompResHHb.keys():
		#	print allCompResHHb[k].toStrS().strip()	
	# --------------------------------------
	# CONSOLIDATE RESULTS
	# --------------------------------------		
	allRes = {}
	for k in allCompResB.keys():
		if k not in allRes.keys(): 
			allRes[k] = allCompResB[k]
			allRes[k].toolsHits = 1
		else: 
			allRes[k].blastEV = allCompResB[k].blastEV	
			allRes[k].toolsHits += 1	
	for k in allCompResPB.keys():
		if k not in allRes.keys(): 
			allRes[k] = allCompResPB[k]
			allRes[k].toolsHits = 1				
		else: 
			allRes[k].psiblastEV = allCompResPB[k].psiblastEV
			allRes[k].toolsHits += 1								
	for k in allCompResH.keys():
		if k not in allRes.keys(): 
			allRes[k] = allCompResH[k]
			allRes[k].toolsHits = 1				
		else: 
			allRes[k].hmmerEV = allCompResH[k].hmmerEV
			allRes[k].toolsHits += 1						
	for k in allCompResJH.keys():
		if k not in allRes.keys(): 
			allRes[k] = allCompResJH[k]
			allRes[k].toolsHits = 1				
		else: 
			allRes[k].jackEV = allCompResJH[k].jackEV
			allRes[k].toolsHits += 1
	if useHH:							
		for k in allCompResHHb.keys():
			if k not in allRes.keys():
				allRes[k] = allCompResHHb[k]
				allRes[k].toolsHits = 1			
			else: 
				allRes[k].HHblitsEV = allCompResHHb[k].HHblitsEV
				allRes[k].toolsHits += 1				 				
	#print CompRes().getTitle()
	#for k in allRes.keys():
	#		print allRes[k].toStrS().strip()
		
	# calculate scores
	for k in allRes.keys():
		if allRes[k].blastEV > 0:
			allRes[k].score += -1.0 * math.log10(allRes[k].blastEV) 
		if allRes[k].psiblastEV > 0:
			allRes[k].score += -1.0 * math.log10(allRes[k].psiblastEV) 
		if allRes[k].hmmerEV > 0:
			allRes[k].score += -1.0 * math.log10(allRes[k].hmmerEV)
		if allRes[k].jackEV > 0:
			allRes[k].score += -1.0 * math.log10(allRes[k].jackEV)
 		if allRes[k].HHblitsEV > 0:
		 	allRes[k].score += -1.0 * math.log10(allRes[k].HHblitsEV)
		 	
	# write it out and to file
	cROF = tOutSeq+"_"+dbNameChomp+ "_results.csv" 
	with open (cROF,'w') as conROF:
		#d print CompRes().getTitle()
		conROF.write(CompRes().getTitle()+'\n')
		for w in sorted(allRes.items(), key=lambda x: (x[1].toolsHits,x[1].score), reverse=True):
#d			print w[1].toStrS().strip()
			conROF.write(w[1].toStrS().strip()+'\n')
			
	# ----------------- TAXMAPPER -------------
	# PASS comparative results to TaxMapper
    # -----------------------------------------
	print " ---------- GENERATING TAXONOMY MAPS ----------- "
   	inRes = tOutSeq+"_"+dbNameChomp+"_blastR.xml"
	outRes = inRes.replace('.xml','')	            
	pTaxB = sh.Command(TAXMAPPER)
	procTaxB = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes,_bg=True)
	procTaxB2 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKP','--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3','--graphRankSep','4',_bg=True)
	procTaxB3 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKPC','--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3',_bg=True)		
	
   	inRes = tOutSeq+"_"+dbNameChomp+"_psiblastR.xml"
	outRes = inRes.replace('.xml','')	
	pTaxPB = sh.Command(TAXMAPPER)
	procTaxPB = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes,_bg=True)
	procTaxPB2 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKP','--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3','--graphRankSep','4',_bg=True)
	procTaxPB3 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKPC','--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3',_bg=True)			
	
   	inRes = tOutSeq+"_"+dbNameChomp+"_phres.hmr"
	outRes = inRes.replace('.hmr','')	
	pTaxH = sh.Command(TAXMAPPER)
	procTaxH  = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes,_bg=True)
	procTaxH2 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKP','--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3','--graphRankSep','4',_bg=True)
	procTaxH3 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKPC','--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3',_bg=True)					
	
	inRes = tOutSeq+"_"+dbNameChomp+"_jhres.jhr"
	outRes = inRes.replace('.jhr','')	
	pTaxJ = sh.Command(TAXMAPPER)
	procTaxJ  = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes,'--doSpeciesOut','N',_bg=True)
	procTaxJ2 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKP','--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3','--graphRankSep','4',_bg=True)
	procTaxJ3 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKPC','--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3',_bg=True)						
	# HHblits
	if useHH:
		inRes = tOutSeq+"_"+dbNameChomp+"_hhres.hhr"
		outRes = inRes.replace('.hhr','')	
		pTaxHHb  = sh.Command(TAXMAPPER)
		procTaxHHb =  pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes,"--HHblitsFilter",_inNames.strip(),_bg=True)
		procTaxHHb2 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKP',"--HHblitsFilter",_inNames.strip(),'--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3','--graphRankSep','4',_bg=True)
		procTaxHHb3 = pTaxB("-I",inRes,"-X",NRTAX,"-O",outRes+'_DKPC',"--HHblitsFilter",_inNames.strip(),'--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P3,C3',_bg=True)										
	
	procTaxB.wait()
	procTaxB2.wait()
	procTaxB3.wait()
	procTaxPB.wait()
	procTaxPB2.wait()
	procTaxPB3.wait()
	procTaxH.wait()
	procTaxH2.wait()
	procTaxH3.wait()
	procTaxJ.wait()
	procTaxJ2.wait()
	procTaxJ3.wait()
	if useHH:
		procTaxHHb.wait()	
		procTaxHHb2.wait()	
		procTaxHHb3.wait()	
	print ' ---- GRAPHS GENERATED! ---- '	
# *************************************************************************************
# *************************************************************************************
# Distant Homology search CODE FOLLOWS
# *************************************************************************************
# *************************************************************************************		
# 1) do local search against NR for each input seq and original
# --------------------------------------
# "LOCAL SEARCH" - PSI BLAST AGAINST NR
# --------------------------------------
	
pblastI = []
for fn in os.listdir ("../input/H_SEQ"):
	if os.path.isdir("../input/H_SEQ/"+fn):
		pass
	elif '.fa' in fn:
		pblastI.append("../input/H_SEQ/"+fn)
# add original
pblastI.append(args.input)
# not set # cores per job
threadsPerJob = max(1,args.threads/len(pblastI))
# run blast on each of those
print ' -> generating models for:',pblastI
print '    -> psi blasting vs NR (local search) with ',threadsPerJob,' threads / job!'
# run multithreaded psi-blasting for these:
blastJobP = []
psiblast = (sh.Command(PSIBLAST))
c = 0
for fn in pblastI:
	c+=1
	tar = fn.replace('.fa','').replace('.fasta','').replace('.faa','')[fn.rfind('/')+1:]
	print '      -> job',c,':',fn,'vs',NRDB
	blastJobP.append(psiblast("-query",fn,"-db",NRDB,"-out",RES_BLAST+"/"+tar+"_vs_nr_bpsi.xml","-parse_deflines","-num_threads",threadsPerJob,"-num_iterations",str(psiIterationToDo),"-show_gis","-outfmt","5",_bg=True))

for job in blastJobP:
	job.wait()
print '    -> DONE psi-blasting'
# --------------------------------------
# FILTER LOCAL SEARCH RESULTS
# --------------------------------------		
print '    -> filtering local search results'
# open blast result
for filename in os.listdir (RES_BLAST):
    # OPTIONS ----
    maxIter = 0
    inTitle = 'gi'
    ids = ''
    parsedef = True    
    # MISC THINGIES 
    psiBlast = False 
    doRecord = False 
    outputFilePath = RES_F_BLAST+'/'+filename.replace(".xml","").strip()+'_'+inTitle+'_e'+str(parseEVmax)+'_c_'+str(minCoverage)+'.fa'
    blastResultPath = RES_BLAST+'/'+filename    
    queryFile = args.input      
    #print " ***** NEW INPUT ***** "
    if 'bpsi_' in blastResultPath :
        psiBlast = True   
        outputFilePath = '/home/ranko/workspace_e4/PyBlastParsing/output/'+filename.replace(".xml","").strip()+'_'+inTitle+'_e'+str(parseEVmax)+'_c_'+str(minCoverage)+'.fa'
    else:
        doRecord = True      
    #d print "R -> ",blastResultPath, psiBlast
    #d print "Q -> ",queryFile        
    #d print "O -> ",outputFilePath
    #d print " ******************** "
    with open(outputFilePath,'w') as outFile:
        with open(blastResultPath) as resultHandle:
			qryCnt = 0
			blastRecords = NCBIXML.parse(resultHandle)
			cnt = 0
			for blastRecord in blastRecords:
				maxIter += 1    	
		
        if psiIteration == -1: 
			psiIteration = maxIter			    	
        with open(blastResultPath) as resultHandle:
            qryCnt = 0
            #print " --> opened "+blastResultPath
    		# parse it
            blastRecords = NCBIXML.parse(resultHandle)
    		# now iterate over it
            cnt = 0
            for blastRecord in blastRecords:
                #print "QUERY : "
                if psiBlast:
                    qryCnt += 1
                    #print "iteration:",qryCnt
                    if qryCnt == psiIteration:
                        doRecord = True
                    else:
                        doRecord = False
                    if qryCnt > psiIteration:
                        break                
                #print "def:",blastRecord.query
                #print "let:",blastRecord.query_letters
                with open(queryFile) as qryS:
                    qryLngt = 0
                    qryL = qryS.readlines()                    
                    for l in qryL:
                        l = l.strip()
                        outFile.write(l+"\n")
                        if not '>' in l:
                            qryLngt += len(l.strip())
            
                for alignment in blastRecord.alignments:
                    if inTitle in alignment.title:
                        hsps_cnt = 0
                        
                        for hsp in alignment.hsps:
                            hsps_cnt +=1
                            coverage = float(hsp.align_length)/float(qryLngt)                            
                            if hsp.expect <= parseEVmax and hsp.score >= parseBSmin and alignment.length <= ltMax and alignment.length >= ltMin\
                            and len(hsp.sbjct) >= minQL and coverage >= minCoverage and doRecord and cnt <= modelResultsMax:
                                #d print psiIteration,maxIter
                                cnt += 1                                                                                
                                defName = str(alignment.title)
                                if parsedef:
                                    defNameP = defName
                                    r = defName.find('|')                                
                                    r2 = r + defName[r+1:].find('|')
                                
                                    defName = defName[0:r2+2]+defName[defName.find('[')+1:defName.find(']')]
                                    defName = defName.replace(" ","_")
                                    
                                if hsps_cnt == 1:
                                    #print str(cnt)+' -> '+ defName #+" eV: ",hsp.expect," lt: ",len(hsp.sbjct)
                                    outFile.write(">"+defName+"\n")                                
                                    ids = ids + alignment.hit_id+"&"
                                    outFile.write(hsp.sbjct+"\n")                                    
# --------------------------------------
# GENERATE MULTIPLE ALIGNMENTS
# --------------------------------------		                                        
print '    -> generating multiple alignments'
clustalOI = []
for fn in os.listdir (RES_F_BLAST):
	if os.path.isdir(RES_F_BLAST+"/"+fn):
		pass
	elif '.fa' in fn:
		clustalOI.append(RES_F_BLAST+"/"+fn)
# now set # cores per job
threadsPerJob = max(1,args.threads/len(clustalOI))
# run blast on each of those
print '       -> generating multiple alignments for:',clustalOI
print '       -> clustal omega with ',threadsPerJob,' threads / job!'
# run multithreaded clustalO
job = []
clustalO = (sh.Command(CLUSTAL_O))
c = 0
for fn in clustalOI:
	c+=1
	tar = RES_M_ALIGNMENT+"/"+fn.replace('.fa','.fna')[fn.rfind('/')+1:]
	print '      -> job',c,': clustal omega for ',fn
	job.append(clustalO("--in",fn,"--out",tar,"--threads="+str(threadsPerJob),"--force","--auto","--outfmt","fasta",_bg=True))
for j in job:
	j.wait()
print '    -> DONE with MULTIPLE ALIGNMENTS'	
# --------------------------------------
# GENERATE HMM MODELS
# --------------------------------------		                                        
print '    -> generating HMM models'
nrModels = 0
hmmBuildI = []
for fn in os.listdir(RES_M_ALIGNMENT):
	if os.path.isdir(RES_M_ALIGNMENT+"/"+fn):
		pass
	elif '.fna' in fn:
		hmmBuildI.append(RES_M_ALIGNMENT+"/"+fn)
# now set # cores per job
threadsPerJob = max(1,args.threads/len(hmmBuildI))
# run blast on each of those
print '       -> building HMM models for:',hmmBuildI
print '       -> hmmbuild with ',threadsPerJob,' threads / job!'
# run multithreaded clustalO
job = []
hmmBuild = (sh.Command(HMMBUILD))
for fn in hmmBuildI:
	c+=1
	tar = RES_HMMS+"/"+fn.replace('.fna','.hmm')[fn.rfind('/')+1:]
	print '      -> job',c,': hmmbuild for ',fn
	job.append(hmmBuild("--cpu",threadsPerJob,tar,fn,_bg=True))
for j in job:
	j.wait()
	nrModels+=1
print ' -> ',nrModels,'MODELS GENERATED! <-'
# --------------------------------------
# HAMMER TARGET DB 
# --------------------------------------		                                        
print ' -> HMMERING target DB (',tarDB,')'
hmmer = sh.Command(HMMSEARCH)
job = []
hmmerI = []
for fn in os.listdir(RES_HMMS):
	if os.path.isdir(RES_HMMS+"/"+fn):
		pass
	elif '.hmm' in fn:
		hmmerI.append(RES_HMMS+"/"+fn)
# now set # cores per job
threadsPerJob = max(1,args.threads/len(hmmerI))
c = 0
for fn in hmmerI:
	c+=1
	tar = RES_HMMSEARCH+"/"+fn.replace('.hmm','')[fn.rfind('/')+1:]+'_vs_'+dbNameChomp+".hmr"
	print '      -> job',c,': hmmsearch ',fn,'vs',tarDB
	job.append(hmmer("--notextw","-o", tar,"--cpu",threadsPerJob,fn,tarDB,_bg=True))
for j in job:
	j.wait()
print ' -> HMMERING DONE '
# --------------------------------------
# PARSE AND FILTER HMMER RESULTS
# --------------------------------------		
print ' -> PROCESSING RESULTS '
# all potential results
allResultsHash = {} #key = gi !
for f in os.listdir(RES_HMMSEARCH):
	resFile = RES_HMMSEARCH+'/'+f
	if os.path.isfile(resFile):
		print "    -> parsing ",resFile
		hmmAllResults = parseHMMResults(resFile)
		cnt = 0		
		for hmmQR in hmmAllResults.qry:
			for hmmQI in hmmQR.iterationResults:        
				for h in hmmQI.hmmTable:					
					if cnt < hrMaxResNR and h.fullSeqEV <= hrMaxEV:
						#d print h.description
						title = h.seqID+h.description
						giList = extractAllGIs(title)						
						for gi in giList:
							cnt +=1							
							newResult = HmmResultRecord(gi, 1, h.fullSeqEV, '', f, h.description)
							newResult.hmmType = "HMMER"
							if gi in allResultsHash.keys():
								allResultsHash[gi].number += 1
								allResultsHash[gi].eVtotal += float(h.fullSeqEV)
							else:
								allResultsHash[gi] = newResult
		print '       -> added',cnt,'hits'								
# print all out
lr = 0
lh = 0
for r in allResultsHash.values():
	lr+=1
	lh+=r.number
print ' -> processed total of ',lr,'different hits',lh,'total hits'

# ------------------- FILTER RESULTS --------------
# -------------------------------------------------
# -> allResultsHash contains all results from different files
# -> we go through those and drop ones with eV below X and less hits then Y
# (less hits meaning only ones hit by Y or more models pass filter)
# -> then we load sequences for those (from NR)
# and filter out duplicates, filling their header with all the gis of
# that sequence (to minimize junk)

# RESULTS, FILTER I
# -> pick results with eValue <= X and hits# > Y
# fill their sequence in
print '    -> filtering'
allResultsS = []
c = 0
cl = len(allResultsHash.values())
for res in allResultsHash.values():
    res.eVavg = float(res.eVtotal) / float(res.number)
    if float(res.number)/float(nrModels) >= hrMinHitsForSameTarget and res.eVavg <= hrMaxEVavg:	
    	fet = fetchFastaFromBlastDB(res.gi,BLAST_DB_CMD,NRDB)
    	c+=1
    	if c % 10000 == 0: 
    		print 'processed',c,'/',cl,'results'
    	if 'ERROR' not in fet.lower():
			res.fasta = fet
			res.header = getFaHeader(res.fasta)
			res.header = addEvToFaHdr(res.header, res.eVavg)
			res.seq = getFaSeq(res.fasta)
			allResultsS.append(res)
		
# sort by sequence (so we remove dups with one pass instead of O(N^2)  
allResultsSortedD = sorted(allResultsS, key=lambda HmmResultRecord: HmmResultRecord.seq)
allResultsNoDupSeqs = set()

joinedRes = HmmResultRecord()
noDupRes = []
for i in range(0,len(allResultsSortedD)):
	if i == 0:
		joinedRes = copy.copy(allResultsSortedD[i])
	# got duplicate seq:
	elif allResultsSortedD[i].seq == allResultsSortedD[i-1].seq:
		joinedRes.header += copy.copy(allResultsSortedD[i].header)
	else: 
		joinedRes.header = '>'+joinedRes.header[1:].replace('>','_')
		noDupRes.append(joinedRes)
		joinedRes = HmmResultRecord()
		joinedRes = copy.copy(allResultsSortedD[i])
		
joinedRes.header = '>'+joinedRes.header[1:].replace('>','_')		
noDupRes.append(joinedRes)
#dfor p in noDupRes:
#d	print p.header
#d	print p.seq
# write results WITHOUT duplicates		
with open(RES_POTENTIALHITS+'/'+inputNOChomp+'_vs_'+dbNameChomp+'_p_hits_nd.csv' ,'w') as outFile:	
	writer = csv.writer(outFile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_ALL)
#	writer.writerow(res.getCSVheader2())        
	for res in sorted(noDupRes, key=lambda HmmResultRecord: (HmmResultRecord.number, HmmResultRecord.eVavg), reverse=True): 
		writer.writerow(res.toSeq2()) 	
print '    -> DONE, ',len(allResultsS),' potential homologs found!','(',len(noDupRes),'without duplicates)' ,'saved output to',RES_POTENTIALHITS+'/'+inputNO+'_p_hits_nd.csv'
# --------------------------------------
# RUN WATERMAN of potential seqs (hmmer) VS ORIGINAL ! 
# --------------------------------------
print '    -> HMMER results: waterman vs original (',args.input,')'
# -> load stuff
# -> water each one
pRes = RES_POTENTIALHITS+'/'+inputNOChomp+'_vs_'+dbNameChomp+'_p_hits_nd.csv'
# prepare output file
inputNO = RES_POTENTIALHITS+'/'+inputNOChomp+'_vs_'+dbNameChomp+'_p_hits__WF.fa'
with open(inputNO,'w') as outF:
	pass
with open(pRes,'r') as inFile:
	psNR = 0
	tsNR = 0
	ulNR = 0
	passNR = 0
	reader = csv.reader(inFile, delimiter='\t', quotechar='"')
	c = 0
	for row in reader:
		c +=1
		if (c > 1):
			rec = HmmResultRecord()
			rec.loadFromCsv2(row)
			psNR +=1
			if 'Error' not in rec.seq:
				seqLt = len(rec.seq)
				origLt = getFastaLt(args.input) 		
				#d print 'water lt (seqLT/origLT',float(seqLt)/float(origLt)			
				if float(seqLt)/float(origLt) >= float(waterMinLt) and float(seqLt)/float(origLt) <= float(waterMaxLt):
					with open('tmpWaterQry.fa','w') as seqA:
						seqA.write(rec.header+'\n'+rec.seq.strip()+'\n')                
					wtrO = doWaterman(WATERMAN, 'tmpWaterQry.fa', args.input)                              
					#check if similarity is good enough
					#d print ' -> simil:', wtrO.getNonGapCoverage()
					if wtrO.getNonGapCoverage() >= waterMinSim:
						passNR +=1
						#print rec.fasta
						if passNR <= waterMaxResults: 
							with open(inputNO,'a') as outF:
								outF.write(rec.header+'\n'+rec.seq.strip()+'\n')
					else: #nope
						ulNR+=1					
				else:
					tsNR+=1					 
					#d print 'seq too short! [',seqLt,'vs',origLt,']'
print '    -> done!',psNR,'considered; ',tsNR,'too short; ',ulNR,'underlimit; ',passNR,'passed check!'
print '    -> waterman confirmed homologs written to ', inputNO
# ----------------------------------------------
# ----- GENERATE GRAPHS OF RESULTS -------------
# ----------------------------------------------
print '    -> generating taxonomy graphs of results'          
pTaxB = sh.Command(TAXMAPPER)
procTaxB = pTaxB("-I",inputNO,"-X",NRTAX,"-O",inputNO,_bg=True)
procTaxB2 = pTaxB("-I",inputNO,"-X",NRTAX,"-O",inputNO+'_DKP','--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P4,C4','--graphRankSep','4',_bg=True)
procTaxB3 = pTaxB("-I",inputNO,"-X",NRTAX,"-O",inputNO+'_DKPC','--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P4,C4',_bg=True)						
procTaxB.wait()
procTaxB2.wait()
procTaxB3.wait()
print '    --> DONE!'
# ----------------------------------------
# ------------ RUN HHBLITS RUN -----------
# ----------------------------------------
# note: only if option for it is entered!
if useHH:
	print " ---> RUNNING HHblits vs ", tarDB, "(",args.threads,"threads)"
	hhblits = sh.Command(HHBLITS)
	hhRNUF = RES_POTENTIALHITS+"/"+inputNOChomp+"_"+dbNameChomp+"_hhres_unfiltered.hhr"
	hhRNF = RES_POTENTIALHITS+"/"+inputNOChomp+"_"+dbNameChomp+"_hhres_filtered.fa"
	hhRNFND = RES_POTENTIALHITS+"/"+inputNOChomp+"_"+dbNameChomp+"_hhres_filtered_nd.fa"
	procHH = hhblits("-o",hhRNUF,"-i",inputSeq,"-d",args.hhdb, "-cpu",args.threads, "-Z","10000","-B","10000","-aliw","10000","-E","10.0","-p","20", _bg=True)				
	procHH.wait()
	allCompResHHb = {}
	#print "     --> HHBLITS (",phResult,")"	
	#print "     	--> parsing "	
	res = parseHHResults(hhRNUF)
	#print "         --> filtering ... "
	hhseqs = 0
	hhseqsnd = 0
	hhps = []
	print '    -> got ',len(res.hhResults),'results! Filtering! ...'
	resF = filterHHblits(res, NRTAXBINARIZED, NRTAXNAMES, NRTAXNODES, targetTaxNames=_inNames)
	with open (hhRNF,'w') as hhPHOut:
		for hhRnew in resF.hhResults: 
			for gi in hhRnew.gis:
				fet = fetchFastaFromBlastDB(gi,BLAST_DB_CMD,NRDB)
				if 'error' not in fet.lower():
					res_fasta = fet
					res_header = getFaHeader(res_fasta)
					res_header = addEvToFaHdr(res_fasta, hhRnew.eV)				
					res_seq = getFaSeq(res_fasta)
					hhps.append([res_header,res_seq])
					hhPHOut.write(res_header+'\n'+res_seq+'\n')
					hhseqs+=1
	hhpsnd = []
	hhps = sorted(hhps, key=lambda x: x[1],reverse=True)
	joinedRes = ['','']
	for i in range(0,len(hhps)):
		if i == 0:
			joinedRes = copy.copy(hhps[i])
			# got duplicate seq:
		elif hhps[i][1] == hhps[i-1][1]:
			joinedRes[0] += copy.copy(hhps[i][0])
		else: 
			joinedRes[0] = '>'+hhps[i][0][1:].replace('>','_')
			hhpsnd.append(joinedRes)
			joinedRes = ['','']
			joinedRes = copy.copy(hhps[i])
			
	with open (hhRNFND,'w') as hhPHOut:
		for f in hhpsnd:
			hhseqsnd+=1
			hhPHOut.write(f[0]+'\n'+f[1]+'\n')
	print '    -> HHblits:',hhseqs,'totals seqs:',hhseqsnd,'without duplicates:'					
# --------------------------------------
# RUN WATERMAN of potential seqs (hhblits) VS ORIGINAL ! 
# --------------------------------------
	pRes = hhRNFND
	pResO = hhRNFND.replace('.fa','__WF.fa')
	# prepare output file
	with open(pResO,'w') as outF:
		pass
	with open(pRes,'r') as inFile:
		psNR = 0
		tsNR = 0
		tlNR = 0
		ulNR = 0
		passNR = 0
		for f in loadAllFastas(pRes):
			psNR+=1
			seq = getFaSeq(f)
			hdr = getFaHeader(f)
			seqLt = len(seq)
			origLt = getFastaLt(args.input) 	
			#d print 'water lt (seqLT/origLT',float(seqLt)/float(origLt)			
			if float(seqLt)/float(origLt) >= float(waterMinLt) and float(seqLt)/float(origLt) <= float(waterMaxLt):
				with open('tmpWaterQry.fa','w') as seqA:
					seqA.write(hdr+'\n'+seq+'\n')                
				wtrO = doWaterman(WATERMAN, 'tmpWaterQry.fa', args.input)
				#d print ' -> simil:', wtrO.getNonGapCoverage()                             
				if wtrO.getNonGapCoverage() >= waterMinSim:
					passNR +=1
					if passNR <= waterMaxResults: 
						with open(pResO,'a') as outF:
							outF.write(hdr+'\n'+seq+'\n')
					else: #nope
						ulNR+=1					
			else:
				tsNR+=1					 
	print '    -> done!',psNR,'considered; ',tsNR,'too short/too long; ', ulNR,'under sim. limit; ',passNR,'passed check!'
	print '    -> waterman confirmed homologs written to ', pResO											
# ---------------------------------------------------------
# ----- GENERATE GRAPHS OF RESULTS (HHBLITS) -------------
# ---------------------------------------------------------
	print '    -> generating taxonomy graphs of results'            
	pTaxB = sh.Command(TAXMAPPER)
	procTaxB =  pTaxB("-I",pResO,"-X",NRTAX,"-O",pResO,"--HHblitsFilter",_inNames,_bg=True)
	procTaxB2 = pTaxB("-I",inRes,"-X",NRTAX,"-O",pResO+'_DKP',"--HHblitsFilter",_inNames,'--doSpeciesOut','N','--taxMergeDivisions','DKP','--rank','DKP','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P4,C4','--graphRankSep','4',_bg=True)
	procTaxB3 = pTaxB("-I",inRes,"-X",NRTAX,"-O",pResO+'_DKPC',"--HHblitsFilter",_inNames,'--doSpeciesOut','N','--taxMergeDivisions','DKPC','--rank','DKPC','--graphFontSizes',"'D16,K16,P12,C11",'--graphOutStats','D1,K1,P4,C4',_bg=True)						
	
	procTaxB.wait()
	procTaxB2.wait()
	procTaxB3.wait()	
	print '    --> DONE!'

# --------------------------------------
# COPY RESULTS TO OUTPUT
# --------------------------------------
outN=inputNOChomp+'_vs_'+dbNameChomp+'-'+str(datetime.now().strftime('%d-%m-%Y_%H:%M'))
print ' -- COPYING OUTPUT -- '
print '	-> destination folder: ', "../output/"+outN
distutils.dir_util.copy_tree("../results", "../output/"+outN)
distutils.dir_util.copy_tree("../input", "../output/"+outN)
print ' -> DONE! '