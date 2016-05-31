#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import multiprocessing
import subprocess
import json
import sys
import optparse
import logging

log = logging.getLogger( 'mergeSkimOutput.py' )

def parseCommandline():
   usage = '%prog [options] InputDirectory'
   descr = ""

   parser = optparse.OptionParser( usage=usage, description=descr, version='%prog v0' )

   parser.add_option(       '--debug', metavar='LEVEL', default='INFO',
                       help='Set the debug level. Allowed values: ERROR, WARNING, INFO, DEBUG. [default = %default]' )
   parser.add_option( '-d', '--dryrun', action='store_true', default=False,
                       help='Don\'t do anything, just list the files. [default = %default]' )
   parser.add_option( '-j', '--jobs', metavar='NUMJOBS', type ='int', default=1,
                       help='Set the maximum number of parallel jobs to be ' \
                            'started [default = %default].' )
   parser.add_option(       '--json',  default='$CMSSW_BASE/src/llvvAnalysis/DMAnalysis/data/sample_13TeV_25ns.json',
                       help='Set path to the json file to get splitting info from. [default = %default]' )
   parser.add_option( '-o', '--outpath', default='./out',
                       help='Set the output directory. [default = %default]' )
   ( options, args ) = parser.parse_args()

   format = '%(levelname)s (%(name)s) [%(asctime)s]: %(message)s'
   date = '%F %H:%M:%S'
   logging.basicConfig( level=logging._levelNames[ options.debug ], format=format, datefmt=date )

   # Check if number of arguments is OK
   if( len(args) != 1 ):
      log.error( 'Expect 1 argument, got %i instead.' %len(args) )
      log.info ( 'Exiting.' )
      sys.exit( 1 )


   # Expand paths
   inputpath = os.path.expandvars( args[0] )
   options.outpath = os.path.expandvars( options.outpath )
   options.json = os.path.expandvars( options.json )

   # Check if inputpath exists
   if( not os.path.exists( inputpath ) ):
      log.error("Input path does not exist: '%s'" % inputpath )
      log.info('Exiting.')
      sys.exit(1)

   # Check if json file exists
   if( not os.path.exists( inputpath ) ):
      log.error("Input json file does not exist: '%s'" % inputpath )
      log.info('Exiting.')
      sys.exit(1)

   # Check if outpath exists, create if necessary
   if( not os.path.exists( options.outpath ) ):
      os.mkdir( options.outpath )

   # All good to go
   return ( options, inputpath )




"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal

def find_all_matching(substring, path):
    result = []
    for root, dirs, files in os.walk(path):
        for thisfile in files:
            if substring in thisfile:
                result.append(os.path.join(root, thisfile ))
    return result

def get_split_dict( jsonpath  ):
   with open( jsonpath, 'r' ) as f:
      procList=json.load(f,encoding='utf-8')

   split_dict = {};
   for desc in procList["proc"] :
      tag = getByLabel(desc,'tag','')
      for d in desc['data']:
         dtag = getByLabel(d,'dtag','')
         split = getByLabel(d,'split',1)
         split_dict[dtag]=split
   return split_dict

def run( cmd ):
   p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
   (stdout, stderr) = p.communicate()

def main():
   (options, inputpath) = parseCommandline()
   split_dict = get_split_dict( os.path.expandvars( options.json ) )

   for part in inputpath.split('/'):
      if( 'crab_' in part ):
         dtag = part.replace('crab_','')
         break;

   input_files =  find_all_matching('.root',inputpath)

   split_in = len(input_files)                     # Number of input files
   split_out = split_dict[ dtag ]                  # Number of output files

   if( split_in == 0 ):
      log.error( 'No input files found.' )
      log.info( 'Exiting.' )
      sys.exit(1)
   if( split_in < split_out ):
      log.error( "Found less input than output files for dtag '%s':" %dtag )
      log.error( '   %i input files vs %i output files' %(split_in,split_out) )
      log.info( 'Exiting.' )
      sys.exit(1)
   log.info( 'Found %i input files, will create %i output files.' % (split_in,split_out) )

   commands = []
   # Make command stubs containing the hadd call and the output file name
   # So each command corresponds to one output file
   for count_out in range(split_out ):
      cmd = [ 'hadd', os.path.join(options.outpath, '%s_%i.root' % (dtag,count_out) ) ]
      commands.append(cmd)


   # Distribute the input files to the output files
   input_files_iter = iter( input_files )
   commands_iter = iter( commands )
   while(True):
      # Iterate over input files until there are none left
      try:
         thisfile = input_files_iter.next()
      except StopIteration:
         break;

      # Iterate over commands/output files
      # After reaching the last commandoutput file, go back to the first
      try:
         thiscommand = commands_iter.next()
      except StopIteration:
         commands_iter = iter( commands )
         thiscommand = commands_iter.next()

      # Add the input file to the command
      thiscommand.append( thisfile )

   if( options.dryrun ):
      for cmd in commands: print cmd
   else:
      pool = multiprocessing.Pool(options.jobs)
      pool.map_async( run, commands ).get(999999)

if __name__ == '__main__':
   main()
