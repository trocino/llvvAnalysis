#! /bin/sh 

phase=3 
out="/afs/cern.ch/work/t/trocino/Work/Analysis_2015/ZHtoLLInv/forked/CMSSW_7_4_15_patch1/src/llvvAnalysis/DMAnalysis/test/analysis_pfmet60_newTruePUweights_byEyeEWKcorr_plusOptimization_20151202_1200/plots_for_limits_lessSysts_10invfb/limits_count"
shapeName="mt_shapes" 
inUrl="/afs/cern.ch/work/t/trocino/Work/Analysis_2015/ZHtoLLInv/forked/CMSSW_7_4_15_patch1/src/llvvAnalysis/DMAnalysis/test/analysis_pfmet60_newTruePUweights_byEyeEWKcorr_plusOptimization_20151202_1200/plots_for_limits_lessSysts_10invfb/plotter.root"
jsonUrl="/afs/cern.ch/work/t/trocino/Work/Analysis_2015/ZHtoLLInv/forked/CMSSW_7_4_15_patch1/src/llvvAnalysis/DMAnalysis/data/sample_13TeV_25ns_ZHinv_plot.json" 
queue="2nd" 
otherOpts=""
#otherOpts=" -O --subNRB15 "
landsArg=" _13TeV "
#cl="0.95"
cl="0.90"

### Mode 0 (cut-n-count) 
python optimize_ZH.py -m 0 -p $phase -f CUTS_ZH.txt -o $out -s $shapeName -G phase${phase}_${shapeName}_mode_0.log -i $inUrl -j $jsonUrl -Q $queue -L $landsArg $otherOpts -C $cl

### Mode 1 (shape) -- cut at 120 GeV 
#python optimize_ZH.py -m 1 -p $phase -f CUTS_ZH_shapes.txt -o $out -s $shapeName -G phase${phase}_${shapeName}_mode_1.log -i $inUrl -j $jsonUrl -Q $queue -L $landsArg $otherOpts -C $cl 

### Mode 1 (shape) 
#python optimize_ZH.py -m 1 -p $phase -f CUTS_ZH_shapes.txt -o $out -s $shapeName -G phase${phase}_${shapeName}_mode_1.log -i $inUrl -j $jsonUrl -Q $queue -L $landsArg $otherOpts -C $cl 

###  python optimize_ZH.py -m 1 -p 3 -f CUTS_ZH.txt -o /afs/cern.ch/work/t/trocino/Work/Analysis_2015/ZHtoLLInv/forked/CMSSW_7_4_15_patch1/src/llvvAnalysis/DMAnalysis/test/analysis_allZH_withOptimAndSysts_20151108_1526/limits/ -s mt_shapes -G phase3_mt_shapes_mode_1.log -i /afs/cern.ch/work/t/trocino/Work/Analysis_2015/ZHtoLLInv/forked/CMSSW_7_4_15_patch1/src/llvvAnalysis/DMAnalysis/test/analysis_allZH_withOptimAndSysts_20151108_1526/resultsForLimits/plotter.root -j /afs/cern.ch/work/t/trocino/Work/Analysis_2015/ZHtoLLInv/forked/CMSSW_7_4_15_patch1/src/llvvAnalysis/DMAnalysis/data/sample_13TeV_25ns_ZHinv_plot.json -Q 2nd -L " --bins lesq1jets --subNRB" > & ! log.tmp & 
### combine -M Asymptotic --cl 0.95 --rRelAcc 0.00000001 --rAbsAcc 0.000000001 -m 1 -n  --run expected --expectSignal=1 -t -1 card_combined.dat > COMB_asympt.log; 