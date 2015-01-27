Distant Homology Search Pipeline (DHSP)
======================================================

Pipeline generates and uses multiple HMM models and verification by Smith-Waterman for high 
accuracy distant homology detection. 

## Installation
Download, ensure following dependencies are met: 

### Programs & libraries:
- BioPython is installed (code requires NCBIXML module)
- sh subprocess handler is installed (see http://amoffat.github.io/sh/), 
   - installation: (sudo) pip install sh (if pip is installed)
- NCBI blast+
- HMMER toolkit (v3)
- clustal omega
- (optional) hhsuite (hhblits)
   - use of hhsuite is not necessary and is problematic, see notes
- smith-waterman (recommended: water code, EMBOSS package)

### Following databases:
- NCBI NR database (required for detection of "close homology" of input sequences)
   - note: optionally DHSP can use RefSeq or other general database, but has not been tested with those
- NCBI taxonomy database (required for taxonomy mapping)
   - binarized NCBI taxonomy database files are required (see TLM below)
- database(s) of sequences to be searched; database must be provided in fasta format AND in BLAST format of same name
   - example: proteins.fa (list of proteins in fasta format) AND BLAST DB called proteins.fa made by makeblastDB

## HOW TO USE: 
### Pipeline Configuration:
- edit runPipeline.py code and set paths for appropriate software (lines 69 - 91)
- make sure folder structure is as provided 

### Preparing for run: 
- put sequence of interest into ./input/ORIG_SEQ
- put its _known_ homologs into ./input/H_SEQ; ideally these should be 2 - 10 experimentally confirmed homologs (if possible) which are not closely related (if possible); for example: little is to be gained by using mouse gene as orig_seq, and rat and another rat homologs in H_SEQ

### Running the pipeline
- runPipeline.py runs the code; -h or --help provides help
- example: ./runPipeline.py -I ../input/ORIG_SEQ/Keap1_hs.fa -D ~/Development/DBs/c_elegans/c_elegans_prot.fa --doComparativeRun N -H N --threads 20 will run pipeline using 20 threads and sample sequences included; database of C. elegans proteome and corresponding BLAST database should be placed in path listed after -D

### Results
- results are automatically copied into output/X_vs_Y_time folder (X being input 'original' sequence and Y being target database
- results are organized as follows: 
	- COMPARATIVE_RESULTS: comparison of BLAST, psi-BLAST, HMMER, jack-HMMER and (optionally) HHblits results (if --doComparativeRun is set to Y)
	- BLAST: BLAST results of 'close homology'
	- F_BLAST: sequences of BLAST results
	- M_ALIGNMENT: multiple alignments of BLAST results
	- HMM_PROFILES: HMM profiles generated
	- ORIG_SEQ and H_SEQ: input sequences
	- DH_HH_RESULTS: HHsuite results (if used)
	- HMM_RESULTS: results of HMM model searches
	- POTENTIAL_HITS: final set of results, as follows: 
		- <name>_nd.csv file: results that passed multiple HMM search filter
		- <name>__WF.fa file: results that also passed smith-waterman alignment filter (these are usually results of interest, but if _nd.csv results might also be of interest)
		- .gv files: graphviz files - results of taxonomical mapping against NCBI taxonomy database via TLM (see TLM below)

### Other instructions
- pipeline can run for a while, especially if whole NR is used to generate HMM models and if large number of homologous sequences is used; if multiple searches are required, it can be helpful to queue them in a bash script and let them run (see _RUN_FOR_ALL_DBS.sh and _RUN_for_drosophila.sh scripts for examples)
- code might contain bugs 
- see https://github.com/rgacesa/TLM for detailed instructions on TLM

### Preparing TLM
- in order to function properly, Taxonomy Landscape Mapper (TLM) requires binarized version of NCBI taxonomy and set of prepared files for taxonomy mapping (TaxFiles); TaxFiles are included (scripts/taxFiles), but might be outdated, while binarized NCBI taxonomy files are NOT included and have to be generate
- to create those, download NCBI taxonomy and then: 
	- to generate NCBI taxonomy binarized files, use binarizeNCBItax.py (run it once for gi_taxid_prot.dmp file and once for gi_taxid_nucl.dmp (if needed))
	- to generate taxFiles, use prepareNCBItaxFiles.py; run it with names.dmp and nodes.dmp as arguments

## More notes: 
- last update: 27/01/2015
- if python is not in /usr/bin/python2.7, code might not be able to self-execute; run it as <python> <codename> or edit 1st line
- tested under Ubuntu 12.04.5 LTS
- TLM tries to generate graphviz graph where everything will be clear and visible, but that is not always possible automatically, in which case some manual tweaking with graphviz might be required
- code is early release and might contain bugs
- keep in mind NCBI databases might have problems of their own which result in artefacts in mapping or other strange behavior

## Other info: 
author: Ranko Gacesa

copyright: 2014 King's College London. All rights reserved.

license: Attribution-NonCommercial-ShareAlike 4.0 International (http://creativecommons.org/licenses/by-nc-sa/4.0/)
        Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.        

contact: ranko.gacesa@kcl.ac.uk
last update: 27/01/2015
