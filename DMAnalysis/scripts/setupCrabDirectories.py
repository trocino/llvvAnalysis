import shutil
import os
import subprocess
import sys



### Test Flag
test = False;
outpath = '/store/user/aalbert/crab/'
workarea = '/disk1/albert/crab/'
filesperjob = 3


### Configuration dependent on the CMSSW version you sourced.
cmssw = os.environ['CMSSW_VERSION']
if( 'CMSSW_7_6' in cmssw ):
    inputpath   = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/skimlist_MC13TeV_76X.txt' )
    configfile  = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_mc_cfg_76X.py' )
    outtag      = 'RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12'
    outpath     = os.path.join( outpath, 'llvv_76' )
    workarea    = os.path.join( workarea, 'llvv_76' )
elif( 'CMSSW_8_0' in cmssw ):
    inputpath   = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/skimlist_MC13TeV_80X.txt' )
    configfile  = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_mc_cfg_80X.py' )
    outtag      = 'RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3'
    outpath     = os.path.join( outpath, 'llvv_80' )
    workarea    = os.path.join( workarea, 'llvv_80' )
else:
    print "Unknown CMSSW version: %s" % cmssw
    print "Exiting."
    sys.exit(1)

templatepath = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/crab_template_cfg.py')

### Read CRAB template config
with open( templatepath, 'r' )as f:
    template = f.read()


### Get datasets from input list
datasets = []
with open(inputpath,'r') as f:
    for l in f.readlines():
        if(':' in l):
            datasets.append( tuple(l.rstrip('\n').split(':')) )

### For each dataset
###     create a separate subfolder
###     put the template config in it
###     fill the template
for ds in datasets:
    tag = ds[0]
    datapath = ds[1]
    try:
        os.mkdir( os.path.join( workarea,'cfg') )
    except OSError: pass
    cfgpath = os.path.join( workarea,'cfg', 'crab_' + tag + '.py' )

    replacements = [ ( '@NAME',         tag              ),
                     ( '@WORKAREA',     workarea         ),
                     ( '@CONFIGFILE',   configfile       ),
                     ( '@INPUTDATASET', datapath         ),
                     ( '@OUTPATH',      outpath          ),
                     ( '@FILESPERJOB',  str(filesperjob) ),
                     ( '@OUTTAG',       outtag           ) ]

    thistemplate = template
    for rep in replacements:
        thistemplate = thistemplate.replace( *rep )

    with( open(cfgpath,'w') ) as f:
        f.write(thistemplate)
