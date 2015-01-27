#!/bin/bash
#cd output
for D in ./*
do
   echo $D
   cd $D/POTENTIAL_HITS/
   echo ' - removing old graphs!'
   rm *_WG_graph*.gv
   rm *.txt
   echo ' - generating new graphs!'
   for T in ./*__WF.fa
   do
      echo '-->' $T
      ~/Development/workspace/DHPipeline/scripts/taxMapper.py -I $T -X ~/Development/DBs/ncbi_tax/ --inFileTax ~/Development/workspace/DHPipeline/scripts/taxFiles/ -O $T
   done
   cd ..
#   echo $D
#   cd ./COMPARATIVE_RESULTS/
#   echo ' - removing old graphs!'
#   rm *.gv
#   rm *.txt
#   echo ' - generating new graphs!'
#   for T in ./*.xml
#   do
#      echo '-->' $T
#      ~/Development/workspace/DHPipeline/scripts/taxMapper.py -I $T -X ~/Development/DBs/ncbi_tax/ --inFileTax ~/Development/workspace/DHPipeline/scripts/taxFiles/ -O $T
#   done
#   for T in ./*.hhr
#   do
#      echo '-->' $T
#      ~/Development/workspace/DHPipeline/scripts/taxMapper.py -I $T -X ~/Development/DBs/ncbi_tax/ --inFileTax ~/Development/workspace/DHPipeline/scripts/taxFiles/ -O $T
#   done
#   for T in ./*.jhr
#   do
#      echo '-->' $T
#      ~/Development/workspace/DHPipeline/scripts/taxMapper.py -I $T -X ~/Development/DBs/ncbi_tax/ --inFileTax ~/Development/workspace/DHPipeline/scripts/taxFiles/ -O $T
#   done
#   for T in ./*.hmr
#   do
#      echo '-->' $T
#      ~/Development/workspace/DHPipeline/scripts/taxMapper.py -I $T -X ~/Development/DBs/ncbi_tax/ --inFileTax ~/Development/workspace/DHPipeline/scripts/taxFiles/ -O $T
#   done
#   exit
#   cd ..
   cd ..
done
