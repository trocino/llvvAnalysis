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
#onlyplots=" --channel all --channel eeeq0jets  --channel mumueq0jets --channel emueq0jets --channel eeeq1jets  --channel mumueq1jets --channel emueq1jets  --only cut1 --only optim_systs --only mt_shapes "
#onlyplots=" --channel all --channel eelesq1jets  --channel mumulesq1jets --channel emulesq1jets  --only cut1 --only optim_systs --only mt_shapes "
#onlyplots="--channel all_ --channel ee --channel mumu --channel emu --only raw --only final --only presel --only eventflow --only optim"
onlyplots="--channel all_ --channel ee --channel mumu --only pfmet --only mt --only zpt --only eventflow "
#onlyplots="--channel all_ --channel ee --channel mumu --channel emu --only zmass_wwctrl "

### Remove data, only mt, pfmet 
#onlyplots=" --channel all --channel lesq1jets  --only cut1 --only pfmet_ --only mt_ "

Ecm="13"
Lumi="2110.246"
#Lumi="5000.000"
#Lumi="10000.000"
runPlotter --json $json --inDir $input --outDir $output --outFile $outfile $onlyplots --iEcm $Ecm --iLumi $Lumi   
exit
