#!/bin/bash
if [ ! "$#" == 4 ]; then
echo "----------- PIPELINE RUNNER! ---------------"
echo " --> runs pipeline against all DBs          "
echo " input parameters = original seq (with path)"
echo "                    min waterman lt"
echo "                    max waterman lt"
echo "                    min waterman coverage"
echo "--------------------------------------------"
fi

T=72
if [ "$#" == 4 ]; then
echo '--------- RUNNING PIPELINE FOR' $1
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/nr.fa --hhfilter none --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/allseqs_nr_archaea_merged.fa --hhfilter archaea --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/allseqs_nr_bacteria_merged.fa --hhfilter bacteria --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/allseqs_nr_protista_merged.fa --hhfilter protista --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/allseqs_nr_plants_merged.fa --hhfilter plants --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/allseqs_nr_fungi_merged.fa --hhfilter fungi --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
./runPipeline.py -I ../input/ORIG_SEQ/$1 -D ~/Development/DBs/nr_fasta/allseqs_nr_misc_merged.fa --hhfilter misc --hhdb ~/Development/DBs/HHsuiteDBs/nr/nr20_12Aug11 --doComparativeRun Y --threads $T --wMinLt $2 --wMaxLt $3 --wMinSim $4
wait
fi
