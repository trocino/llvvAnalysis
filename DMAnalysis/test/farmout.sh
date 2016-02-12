# source me
inputdir=$1
shift

farmoutAnalysisJobs \
  --infer-cmssw-path \
  --fwklite \
  --input-dir=$(pwd)/$inputdir \
  --match-input-files="*_cfg.py" \
  --job-generates-output-name \
  $@ \
  llvvAnalysis \
  $CMSSW_BASE/bin/$SCRAM_ARCH/run2015_WIMPAnalysis \
  '$inputFileNames'

