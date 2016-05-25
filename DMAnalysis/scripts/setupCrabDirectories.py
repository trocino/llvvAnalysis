import shutil
import os
import subprocess
inputpath = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/skimlist_MC13TeV.txt' )

### Read CRAB template config
with open( os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/crab_template_cfg.py', 'r') )as f:
    template = f.read()

### Input info for filling the template
outtag = 'RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3'
outpath = '/store/user/aalbert/crab/llvv_80/'
workarea = '/disk1/albert/crab/llvv_80/'
filesperjob = 3

### Config File to run via CRAB
configfile = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_mc_cfg.py' )

### Test Flag
test = False;


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
