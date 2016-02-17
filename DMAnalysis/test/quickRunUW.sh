#!/bin/sh 
execute="run2015_WIMPAnalysis"
main="${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis"
templ="$main/test/runWIMPSelAnalysis_cfg.py_templ"
CWD=`pwd`
outfile=$CWD/analysis_$(git describe --always)
runlog=$outfile/LOGFILES
mkdir -p $runlog
json="$main/data/sample_13TeV_25ns_ZHinv_15Feb2016.json"
#json="$main/data/sample_13TeV_25ns_ZHinv_plot.json"
input="/store/user/nsmith/llvvNtuple_15Feb2016/"
wimpweights="$main/data/weights/PileupWeights_Sep182015.root"
queue="2nd"
mkdir -p /tmp/`whoami`
#runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=False @is2011=False @runOptimization=False @wimpweights=${wimpweights}" -s $queue 
#runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=True @is2011=False @runOptimization=False @wimpweights=${wimpweights}" -s $queue 
runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=True @is2011=False @runOptimization=True @wimpweights=${wimpweights}" -s $queue 
exit
