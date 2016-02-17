# source me
inputdir=$1
shift

if [[ $@ != *"--resubmit-failed-jobs"* ]]; then
  rm -rf /nfs_scratch/nsmith/llvvAnalysis-run2015_WIMPAnalysis
  gsido rm -rf /hdfs/store/user/nsmith/llvvAnalysis-run2015_WIMPAnalysis
fi

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

