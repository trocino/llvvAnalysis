import shutil
import os
import subprocess
import sys
import logging
import optparse

log = logging.getLogger( 'setupCrabDirectories' )



def parseCommandline():
    usage = '%prog [options] SKIMTAG'
    descr = "SKIMTAG may be any string that is suited as a folder name. It will be used to identify your skim."

    parser = optparse.OptionParser( usage=usage, description=descr, version='%prog v0' )
    parser.add_option(       '--debug', metavar='LEVEL', default='INFO',
                        help='Set the debug level. Allowed values: ERROR, WARNING, INFO, DEBUG. [default = %default]' )
    parser.add_option( '-f', '--filesPerJob', type ='int', default=3,
                        help='Set the maximum number of files to be processed by one CRAB job [default = %default].' )
    parser.add_option( '-o', '--outpath', type='str', default='/store/user/aalbert/crab/',
                        help='The path where CRAB will place output files on STORAGE.[default = %default]')
    parser.add_option( '-s', '--storageSite', metavar='STORAGE', type='str', default='T2_DE_RWTH',
                        help='The storage site where CRAB will place output files.[default = %default].')
    parser.add_option( '-t', '--template', type ='str', default='$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/crab_template_cfg.py',
                        help='Path to the CRAB template config. [default = %default].' )
    parser.add_option( '-w', '--workarea', type ='str', default='.',
                        help='Path where CRAB config files will be created. [default = %default].' )

    ( options, args ) = parser.parse_args()

    format = '%(levelname)s (%(name)s) [%(asctime)s]: %(message)s'
    date = '%F %H:%M:%S'
    logging.basicConfig( level=logging._levelNames[ options.debug ], format=format, datefmt=date )

    # Check if number of arguments is OK
    if( len(args) != 1 ):
        log.error( 'Expect 1 argument, got %i instead.' %len(args) )
        log.info ( 'Exiting.' )
        sys.exit( 1 )

    options.skimtag = args[0]

    ### Configuration dependent on the CMSSW version you sourced.
    cmssw = os.environ['CMSSW_VERSION']
    if( 'CMSSW_7_6' in cmssw ):
        options.inputpath       = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/skimlist_MC13TeV_76X.txt' )
        options.configfile_mc   = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_mc_cfg_76X.py' )
        options.configfile_data = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_data_cfg_76X.py' )
        options.outtag_mc       = 'RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12'
        options.splitting_mc    = 'FileBased'
        options.json_mc         = ''
        options.splitting_data  = 'FileBased'
        options.json_data       = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt')
        options.outpath         = os.path.join( options.outpath, 'llvv_76', options.skimtag )
        options.workarea        = os.path.join( options.workarea, 'llvv_76', options.skimtag )
    elif( 'CMSSW_8_0' in cmssw ):
        options.inputpath       = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/skimlist_MC13TeV_80X.txt' )
        options.configfile_mc   = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_mc_cfg_80X.py' )
        options.configfile_data = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/test/run_mainAnalyzer_data_cfg_80X.py' )
        options.outtag_mc       = 'RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1'
        options.splitting_mc    = 'FileBased'
        options.json_mc         = ''
        options.splitting_data  = 'FileBased'
        options.json_data       = os.path.expandvars( '$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/Cert_271036-274421_13TeV_PromptReco_Collisions16_JSON.txt')
        options.outpath         = os.path.join( options.outpath, 'llvv_80', options.skimtag )
        options.workarea        = os.path.join( options.workarea, 'llvv_80', options.skimtag )
    else:
        print "Unknown CMSSW version: %s" % cmssw
        print "Exiting."
        sys.exit(1)

    options.template = os.path.expandvars( options.template )

    # All good to go
    return options


def main():
    ### Read command line arguments
    options = parseCommandline()
    log.info( 'Using work area: %s' % options.workarea )
    log.info( 'Using T2 storage site: %s' % options.storageSite )
    log.info( 'Using T2 storage path: %s' % options.outpath )

    ### Read CRAB template config
    with open( options.template, 'r' )as f:
        template = f.read()


    ### Get datasets from input list
    datasets = []
    with open(options.inputpath,'r') as f:
        for l in f.readlines():
            if(':' in l):
                datasets.append( tuple(l.rstrip('\n').split(':')) )
    log.info( 'Creating CRAB configs for %i data sets.' % len(datasets) )

    ### For each dataset
    ###     create a separate subfolder
    ###     put the template config in it
    ###     fill the template
    for ds in datasets:
        tag = ds[0]
        datapath = ds[1]

        ### Differentiate between MC and data
        if( 'MC13TeV' in tag ):
            configfile = options.configfile_mc
            outtag = options.outtag_mc
            json = options.json_mc
            split = options.splitting_mc
        elif( 'Data13TeV' in tag ):
            configfile = options.configfile_data
            json = options.json_data
            split = options.splitting_data
            outtag = datapath.split('/')[2]
        else:
            log.error( "Cannot decide whether this tag is MC or data: '%s'" % tag )
            log.info ( "Skipping tag: '%s'" % tag )

        ### Create folder
        subdir = os.path.join( options.workarea,'cfg' )
        if not os.path.exists( subdir ):
            os.makedirs( subdir )

        ### Fill the template
        replacements = [ ( '@NAME',        tag              ),
                        ( '@WORKAREA',     options.workarea         ),
                        ( '@CONFIGFILE',   configfile       ),
                        ( '@INPUTDATASET', datapath         ),
                        ( '@OUTPATH',      options.outpath          ),
                        ( '@FILESPERJOB',  str(options.filesPerJob) ),
                        ( '@JSON',  json ),
                        ( '@Method',  split ),
                        ( '@OUTTAG',       outtag           ),
                        ( '@STORAGESITE', options.storageSite ) ]
        thistemplate = template
        for rep in replacements:
            thistemplate = thistemplate.replace( *rep )

        ### Save the filled template
        cfgpath = os.path.join( subdir, 'crab_' + tag + '.py' )
        with( open(cfgpath,'w') ) as f:
            f.write(thistemplate)
    log.info( 'Done.' )

if __name__ == '__main__':
    main()
