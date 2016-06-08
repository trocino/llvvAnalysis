#!/bin/sh 
execute="run2015_WIMPAnalysis"
main="${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis"
stamp=`date +%Y%m%d_%H%M`
templ="$main/test/runWIMPSelAnalysis_cfg.py.templ"
CWD=`pwd`
outfile=$CWD/analysis_attempt_${stamp}
runlog=$outfile/LOGFILES
mkdir -p $runlog
#json="$main/data/sample_13TeV_25ns_testReweighting.json"
json="$main/data/sample_13TeV_25ns_testReweighting_2.json"
input="/store/group/phys_exotica/monoZ/llvvNtuple_ReReco_25Feb2016/"
queue="2nd"
runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=False @is2011=False @runOptimization=False" -s $queue 
#runLocalAnalysisOverSamples.py -g $runlog -e $execute -j $json -o $outfile -d $input -c $templ -p "@runSystematics=True  @is2011=False @runOptimization=True " -s $queue
exit
