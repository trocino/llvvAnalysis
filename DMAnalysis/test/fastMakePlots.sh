#!/bin/sh 
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <INPUT DIR> <OUTPUT DIR>" 
    exit 1
fi
input="$1" 
indirfin=${1:(-1)}
if [ $indirfin != "/" ]; then 
    input="$input/" 
fi 
output="$2" 
#output="$input/results" 
main="${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis"
json="$main/data/sample_13TeV_25ns_ZHinv_plot.json"
#json="$main/data/sample_13TeV_25ns_ZHinv_noData_plot.json"
outfile="$output/plotter.root" 
onlyplots="--channel all_ --channel ee --channel mumu --channel emu --channel ll --only raw --only presel --only final --only eventflow "

### Remove data final plots
#onlyplots=" --channel all_ --channel ee --channel mumu --only final --isDataBlind "

Ecm="13"
Lumi="2110.246"
#Lumi="5000.000"
#Lumi="10000.000"
runPlotter --json $json --inDir $input --outDir $output --outFile $outfile $onlyplots --iEcm $Ecm --iLumi $Lumi   
exit
