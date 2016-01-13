#!/bin/sh 
execute="run2015_WIMPAnalysis"
main="${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis"
stamp=`date +%Y%m%d_%H%M`
templ="$main/test/runWIMPSelAnalysis_cfg.py_templ"
CWD=`pwd`
#outfile=$CWD/analysis_ZH125_${stamp}
outfile=$CWD/analysis_pfmet100_newTruePUweights_byEyeEWKcorr_MITCuts_noSysts_${stamp}
runlog=$outfile/LOGFILES
mkdir -p $runlog
#json="$main/data/sample_13TeV_25ns_ZHinv.json"
json="$main/data/sample_13TeV_25ns_ZHinv_plot.json"
input="/store/group/phys_exotica/monoZ/llvvNtuple_25ns_29Oct2015"
wimpweights="$main/data/weights/PileupWeights_Sep182015.root"
queue="2nd"
mkdir -p /tmp/`whoami`
#runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=False @is2011=False @runOptimization=False @wimpweights=${wimpweights}" -s $queue 
#runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=True @is2011=False @runOptimization=False @wimpweights=${wimpweights}" -s $queue 
runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=False @is2011=False @runOptimization=False @wimpweights=${wimpweights}" -s $queue 
exit
