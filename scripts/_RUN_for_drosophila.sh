#!/bin/bash
T=72
# --------- 2 B ----------------
# --------- NRF -------------
# copy stuff to input
rm ../input/H_SEQ/*
rm ../input/ORIG_SEQ/*
cp ../data/ORIGS/Nrf2_hs.fa ../input/ORIG_SEQ
cp ~/Dropbox/doktorat/sequences/Nrf2/homologs/* ../input/H_SEQ
IN=Nrf2_hs.fa
# --- RUN FOR C ELEGANS ---
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/c_elegans/c_elegans_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.05 --wMaxLt 20.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 500 --pbCovMin 0.6 --doComparativeRun no
wait
# --- RUN FOR DROSOPHILA ---
echo '--------- RUNNING PIPELINE FOR' $IN
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/drosophila/drosophila_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.05 --wMaxLt 20.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 500 --pbCovMin 0.6 --doComparativeRun no
wait

# -------- KEAP ----------
# copy stuff to input
rm ../input/H_SEQ/*
rm ../input/ORIG_SEQ/*
cp ../data/ORIGS/Keap1_hs.fa ../input/ORIG_SEQ
cp ~/Dropbox/doktorat/sequences/Keap1/homologs/* ../input/H_SEQ
IN=Keap1_hs.fa
# --- RUN FOR DROSOPHILA ---
echo '--------- RUNNING PIPELINE FOR' $IN
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/drosophila/drosophila_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.05 --wMaxLt 20.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 500 --pbCovMin 0.6 --doComparativeRun no
wait
# --- RUN FOR C ELEGANS ---
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/c_elegans/c_elegans_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.05 --wMaxLt 20.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 500 --pbCovMin 0.6 --doComparativeRun no
wait

# -------- RUN 2 (2/A)-----------

# --------- NRF -------------
# copy stuff to input
rm ../input/H_SEQ/*
rm ../input/ORIG_SEQ/*
cp ../data/ORIGS/Nrf2_hs.fa ../input/ORIG_SEQ
cp ~/Dropbox/doktorat/sequences/Nrf2/homologs/* ../input/H_SEQ
IN=Nrf2_hs.fa
# --- RUN FOR C ELEGANS ---
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/c_elegans/c_elegans_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.5 --wMaxLt 2.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 100 --pbCovMin 0.6 --doComparativeRun no
wait
# --- RUN FOR DROSOPHILA ---
echo '--------- RUNNING PIPELINE FOR' $IN
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/drosophila/drosophila_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.5 --wMaxLt 2.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 100 --pbCovMin 0.6 --doComparativeRun no
wait

# -------- KEAP ----------
# copy stuff to input
rm ../input/H_SEQ/*
rm ../input/ORIG_SEQ/*
cp ../data/ORIGS/Keap1_hs.fa ../input/ORIG_SEQ
cp ~/Dropbox/doktorat/sequences/Keap1/homologs/* ../input/H_SEQ
IN=Keap1_hs.fa
# --- RUN FOR DROSOPHILA ---
echo '--------- RUNNING PIPELINE FOR' $IN
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/drosophila/drosophila_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.5 --wMaxLt 2.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 100 --pbCovMin 0.6 --doComparativeRun no
wait
# --- RUN FOR C ELEGANS ---
./runPipeline.py -I ../input/ORIG_SEQ/$IN -D ~/Development/DBs/c_elegans/c_elegans_prot.fa --useHHblits N --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun N --threads $T --wMinLt 0.5 --wMaxLt 2.0 --wMinSim 0.3 --hmrAvgEV 10 --hmrMaxEV 10 --hmrMinModels 0.3 --pbEVMax 1.0e-3 --pbMaxSeqs 100 --pbCovMin 0.6 --doComparativeRun no
wait

