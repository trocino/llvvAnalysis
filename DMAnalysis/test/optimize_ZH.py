#!/usr/bin/env python
#################################
#
#	ZH->ll+invisible
#
#################################

import os,sys
#import json
import getopt
import commands
import ROOT
from ROOT import TFile, TGraph, TCanvas, TF1, TH1

#default values
shapeBased='1'
shapeName='zpt_rebin_shapes'
#shapeName='redMet_rebin_shapes'
#shapeName='met_rebin_shapes'

inUrl='~/work/Output_invHiggs/Plots/invHiggs_8TeV_Moriond2013_20130403_V3.5/NewPlots_20130404/plotter.root'
#inUrl='~/work/Output_invHiggs/Plots/invHiggs_8TeV_Moriond2013/NewPlots_20130321/plotter.root'
CWD=os.getcwd()
phase=-1
jsonUrl='$CMSSW_BASE/src/CMGTools/HtoZZ2l2nu/data/invHiggs/samples_MCFMxs_invHiggs_8TeV_Moriond2013_RJ.json'
CMSSW_BASE=os.environ.get('CMSSW_BASE')
LandSArg=''
#LandSArg+=' --indexvbf 78 '
#LandSArg+=' --indexvbf 9 '
#LandSArg+=' --bins eq0jets,eq1jets,lesq1jets'
LandSArg+=' --bins eq0jets,eq1jets'
#LandSArg+=' --bins lesq1jets'
#LandSArg+=' --bins eq1jets'
#LandSArg+=' --bins eq0jets'
#LandSArg+=' --systpostfix _8TeV' #put in config file
#LandSArg+=' --systpostfix _7TeV' #put in config file
#LandSArg+=' --subNRB'
#LandSArg+=' --subDY $CMSSW_BASE/src/CMGTools/HtoZZ2l2nu/bin/G/ZH_DataDrivenDY_7TeV_May26_01jetbins_ForLimits_RedMET_Final/gamma_out.root'
#LandSArg+=' --subDY $CMSSW_BASE/src/CMGTools/HtoZZ2l2nu/bin/G/ZH_DataDrivenDY_8TeV_May26_01jetbins_ForLimits_RedMET_Final/gamma_out.root'
#LandSArg+=' --subDY /afs/cern.ch/work/r/rewang/Zhllnunu/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/bin/G/gamma_out_8TeV.root'
#LandSArg+=' --subDY /afs/cern.ch/work/r/rewang/Zhllnunu/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/bin/G/gamma_out_7TeV.root'

DataCardsDir='cards'

#MASS = [ 400, 500, 600 ] 
MASS = [ 110, 125, 150, 200, 300, 400, 500, 600 ] 
#MASS = [ 110, 150, 200, 300, 400, 500, 600 ] 
#MASS = [ 125 ] 
SUBMASS = MASS

cl=0.95

cutList=''
LOGFILE=''
def help() :
   print '\n\033[92m optimize.py \033[0m \n'
   print ' -p phase (no default value is assigned)'
   print '\t 1 --> submit landS jobs for all selection point'
   print '\t 2 --> check the logs to find the optimal selection point'
   print '\t      from the ouptut of the logs you can search for the optimal points yourself ;)'
   print '\t      and edit phase3 of this script with your optimal points (note: change this to be given as input)'
   print '\t 3 --> you are prompted to choose the best cuts for the selection: the limits will be comptued for the list of cuts'
   print '\t       if -f LIST.txt is given the LIST of cuts passed will be used instead'
   print '\t 4 --> once all the final limit jobs have been run, use this phase to build the brazilian flag plot'
   print ' -m mode (default='+shapeBased+')'
   print '\t 0 --> cut and count based analysis'
   print '\t 1 --> shape based analysis'
   print ' -s shapename (default='+shapeName+')'
   print ' -i inputfile (default='+inUrl+')'
   print ' -o output (default='+CWD+')'
   print ' -j jsonfile (default='+jsonUrl+')'
   print ' -G logfile from submitting queue jobs'
   print ' -Q queue for submitting jobs'
   print ' -L LandSArg'
   print ' -D LandSArg subDY'
   print ' -O <other options of the form "--comm">'
   print ' -C confidence level for limits (default: 0.95)'
   print '\nUsage example: \033[92m python optimize.py -m 0 -i ~/work/plotter.root -o ~/work/ -p 1 \033[0m'
   print '\nNote: CMSSW_BASE must be set when launching optimize.py (current values is: ' + CMSSW_BASE + ')\n'

#parse the options
try:
   # retrive command line options
   shortopts  = "p:f:m:i:s:j:o:h:G:Q:L:D:O:C:?"
   opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
   # print help information and exit:
   print "ERROR: unknown options in argument %s" % sys.argv[1:]
   help()
   sys.exit(1)

for o,a in opts:
   if o in("-?", "-h"):
      help()
      sys.exit(1)
   elif o in('-m'): shapeBased = a
   elif o in('-i'): inUrl = a
   elif o in('-p'): phase = int(a)
   elif o in('-o'): CWD=a
   elif o in('-j'): jsonUrl=a
   elif o in('-s'): shapeName=a
   elif o in('-f'): cutList=a
   elif o in('-G'): LOGFILE=a
   elif o in('-Q'): QUEUE=a
   elif o in('-L'): LandSArg = LandSArg+' --systpostfix '+a
   elif o in('-D'): LandSArg = LandSArg+' --subDY '+a
   elif o in('-O'): LandSArg = LandSArg+' '+a
   elif o in('-C'): cl = a

#print "LogFile: "+LOGFILE
#print shapeName
#sys.exit(1)


if(phase<0 or len(CMSSW_BASE)==0):
   help()
   sys.exit(1)


#auxiliary function
def findCutIndex(cut1, hcut1, cut2, hcut2, cut3, hcut3):#, cut4, hcut4):
   for i in range(1, hcut1.GetXaxis().GetNbins()):
      if(hcut1.GetBinContent(i)<cut1-0.001):continue;
      if(hcut2.GetBinContent(i)<cut2-0.001):continue;
      if(hcut3.GetBinContent(i)<cut3-0.001):continue;
      #if(hcut4.GetBinContent(i)<cut4-0.001):continue;
      return i;
   return hcut1.GetXaxis().GetNbins();
#
#def findSideMassPoint(mass):
#   global MASS
#   LMass=0
#   RMass=9999
#   for m in MASS:
#      if(m<=mass and m>=LMass):LMass=m
#      if(m>=mass and m<=RMass):RMass=m
#   return [LMass,RMass]
#
#######################

#prepare the output
OUT = CWD+'/JOBS/'
if(shapeBased=='1'): OUT+='SHAPE/'
else:		     OUT+='COUNT/'
os.system('mkdir -p ' + OUT)

if(shapeBased=='1'): DataCardsDir+='Shape'

#get the cuts
file = ROOT.TFile(inUrl)
cuts1   = file.Get('Data/optim_cut1_MET')
cuts2   = file.Get('Data/optim_cut1_Balance')
cuts3   = file.Get('Data/optim_cut1_DphiZMET')
#cuts4   = file.Get('WW#rightarrow 2l2#nu/optim_cut1_zm')
#cuts4   = file.Get('WW#rightarrow 2l2#nu/optim_cut1_jetthr')



######################################################################

if( phase == 1 ):

   print '########################################'
   print '# PHASE 1                              #'
   print '# RUN LIMITS FOR ALL POSSIBLE CUTS     #'
   print '########################################'

   ALLSCRIPT = open(OUT+"/ALL_SCRIPT.sh","w")
   LOCALSCRIPT = open(OUT+"/LOCAL_SCRIPT.sh","w")
   FILE = open(OUT+"/LIST.txt","w")
   for i in range(1,cuts1.GetNbinsX()):
      #if(shapeBased=='1' and cuts3.GetBinContent(i)<780):continue
      FILE.writelines("index="+str(i) + " --> pfMET>" + str(cuts1.GetBinContent(i)).rjust(5) + " |1-PFMET/qT|<" + str(cuts2.GetBinContent(i)).rjust(5) + " dphi(Z,PFMET)>" + str(cuts3.GetBinContent(i)).rjust(5) + "\n")

      #create wrappper script for each set of cuts ans submit it
      SCRIPT = open(OUT+'script_'+str(i)+'.sh',"w")
      SCRIPT.writelines('echo "TESTING SELECTION : ' + str(i).rjust(5) + ' --> pfMET>' + str(cuts1.GetBinContent(i)).rjust(5) + ' |1-PFMET/qT|<' + str(cuts2.GetBinContent(i)).rjust(5) + ' dphi(Z,PFMET)>'+str(cuts3.GetBinContent(i)).rjust(5)+'";\n')
      SCRIPT.writelines('cd ' + CMSSW_BASE + '/src;\n')
      SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc5_amd64_gcc462")+";\n")
      SCRIPT.writelines("eval `scram r -sh`;\n")
      SCRIPT.writelines('cd ' + CWD + ';\n')
      #SCRIPT.writelines('cd /tmp/;\n')

      pos = 0
      for m in MASS:
	 #if(m == 500): continue ## just run one mass point
         shapeBasedOpt=''
         if(shapeBased=='1') : shapeBasedOpt='--shape'
         cardsdir = 'optim/H'+ str(m);
         if(shapeBased=='0'): cardsdir+='_count_'+str(i)
         if(shapeBased=='1'): cardsdir+='_shape_'+str(i)
         SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
         SCRIPT.writelines("runLimits --m " + str(m) + " --metEXTRAPOL "+ str(cuts1.GetBinContent(i)).rjust(5) + " --respEXTRAPOL " + str(cuts2.GetBinContent(i)).rjust(5) + " --dphiEXTRAPOL " + str(cuts3.GetBinContent(i)).rjust(5) + " --histo " + shapeName + " --in " + inUrl + " --syst " + shapeBasedOpt + " --index " + str(i)     + " --json " + jsonUrl +" " + LandSArg + " ;\n")
         SCRIPT.writelines("sh combineCards.sh;\n")
	 SCRIPT.writelines("combine -M Asymptotic --rRelAcc 0.00000001 --rAbsAcc 0.000000001 -m " + str(m) + " --expectSignal=1 -t -1 --run expected card_combined.dat > COMB_asympt.log;\n") #unblind data
         #SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + " --expectSignal=1 -t -1 --run expected card_combined.dat > COMB_asympt.log;\n")
         #SCRIPT.writelines("combine -M MaxLikelihoodFit -m " +  str(m) + " --saveNormalizations card_combined.dat > COMB_maxlkl.log;\n")
         #SCRIPT.writelines("extractFitNormalization.py mlfit.root hzz2l2v_"+str(m)+"_?TeV.root > fit.txt;\n")
         #SCRIPT.writelines('tail -n 100 COMB.log > ' +OUT+str(m)+'_'+str(i)+'.log;\n')
         SCRIPT.writelines('cd -;\n\n')
         pos += 1
      SCRIPT.close()
      LOGFILE_i = LOGFILE+"_optim"+str(i)+".log"
      ALLSCRIPT.writelines("bsub -sp 99 -G u_zh -q "+QUEUE+" -o "+LOGFILE_i+" -J optim"+str(i)+" 'sh " + OUT+"script_"+str(i)+".sh &> "+OUT+"script_"+str(i)+".log' \n")
      LOCALSCRIPT.writelines("sh "+ OUT+"script_"+str(i)+".sh >& "+OUT+"script_"+str(i)+".log & \n")
   FILE.close()

   ALLSCRIPT.close()
   LOCALSCRIPT.close()


######################################################################
elif(phase == 2):

   print '########################################'
   print '# PHASE 2                              #'
   print '# SCANNING ALL SETS OF CUTS            #'
   print '########################################'

   Listfile = OUT + "/LIST.txt"
   fileName = CWD + "/OPTIM_"
   if(shapeBased=='1'):  fileName+='SHAPE'
   else:                 fileName+='COUNT'

   for m in MASS:
      #if(m == 500): continue ## just run one mass point
      FILE = open(fileName+"_"+str(m)+".txt","w")
      print 'Starting mass ' + str(m)
      FILE.writelines("------------------------------------------------------------------------------------\n")
      BestLimit = []
      fileList = commands.getstatusoutput("ls " + CWD + "/optim/"+"H"+str(m)+"_*/"+"COMB_asympt.log")[1].split();
      for f in fileList:
         exp = commands.getstatusoutput("cat " + f + " | grep \"Expected 50.0%:\"")[1];
         if(len(exp)<=0):continue
	 median = float(exp.split()[4])

         if(float(median)<=0.0):continue
	 index = f.split('/')
	 for item in index:
	 	if "H"+str(m) in item:
			INDEX = int(item.split('_')[2])

	 myindex = "index="+('%d ' % int(INDEX) )
	 myListfile = commands.getstatusoutput("cat " + Listfile + " | grep \""+myindex+"\"")[1];
	 BestLimit.append("mH="+str(m)+ " --> " + ('%f' % float(median)) + "\t" + myListfile)

      #sort the limits for this mass
      BestLimit.sort()
      for s in BestLimit:
         FILE.writelines(s+"\n")

      #all done
      FILE.close()
      print("file: "+fileName+"_"+str(m)+".txt is written: it contains all selection points ordered by exp limit")


######################################################################
elif(phase == 6):

   print '########################################'
   print '# PHASE 6                              #'
   print '# SCANNING BEST CUT FROM 5 MASS POINTS #'
   print '########################################'

   fileName = CWD + "/OPTIM_"
   if(shapeBased=='1'):  fileName+='SHAPE'
   else:                 fileName+='COUNT'
   myFILE = open(CWD+"/Best.txt","w")

   for m in MASS:
      aFILE = fileName+"_"+str(m)+".txt"
      print 'Starting mass ' + str(m)
      myFILE.writelines("------------------------------------------------------------------------------------\n")

      exp = commands.getstatusoutput("cat " + aFILE + " | grep \"mH=\"")[1];
      if(len(exp)<=0):continue
      explimit = exp.split()[2]
      myexp = commands.getstatusoutput("cat " + aFILE + " | grep \""+explimit+"\"")[1];
      myFILE.writelines(myexp+"\n")

   myFILE.close()

######################################################################

elif(phase == 3 ):

   print '########################################'
   print '# PHASE 3                              #'
   print '# FINAL LIMITS                         #'
   print '########################################'

   Gmet  = ROOT.TGraph(len(SUBMASS));
   Gbalance = ROOT.TGraph(len(SUBMASS));
   Gdphi = ROOT.TGraph(len(SUBMASS));
   #Gzm = ROOT.TGraph(len(SUBMASS));

   if(cutList=='') :
      fileName = OUT+"/OPTIM_"
      if(shapeBased=='1'):
         fileName+='SHAPE'
      else:
         fileName+='COUNT'
      fileName+=".txt"

      mi=0
      for m in MASS:

         #if you want to display more than 3 options edit -m3 field
         cut_lines=commands.getstatusoutput("cat " + fileName + " | grep 'mH="+str(m)+"' -m20")[1].split('\n')
         print 'mH='+str(m)+'\tOption \tR \tmin redMET\tBalance\tDPhi\tZMass'
         ictr=1
         for c in cut_lines:
            print '\t #'+ str(ictr) + '\t' + c.split()[2] + '\t' + c.split()[4] + '\t' + c.split()[5] + '\t' + c.split()[6] + '\t' + c.split()[7]
            ictr+=1
         print "Which option you want to keep?"
         opt = int(raw_input(">"))-1

         #save cut chosen
         metCut=float(cut_lines[opt].split()[4])
         balanceCut=float(cut_lines[opt].split()[5])
	 dphiCut=float(cut_lines[opt].split()[6])
         #zmCut=float(cut_lines[opt].split()[7])

         Gmet .SetPoint(mi, m, metCut);
         Gbalance.SetPoint(mi, m, balanceCut);
	 Gdphi.SetPoint(mi, m, dphiCut);
         #Gzm.SetPoint(mi, m, zmCut);
         mi+=1
   else :
      mi=0
      f= open(cutList,'r')
      for line in f :
         vals=line.split(' ')
         Gmet .SetPoint(mi, float(vals[0]), float(vals[1]));
         Gbalance.SetPoint(mi, float(vals[0]), float(vals[2]));
         Gdphi.SetPoint(mi, float(vals[0]), float(vals[3]));
         #Gzm.SetPoint(mi, float(vals[0]), float(vals[4]));
         mi+=1
      f.close()

   Gmet.Set(mi);
   Gbalance.Set(mi);
   Gdphi.Set(mi);
   #Gzm.Set(mi);



   #display cuts chosen
   ROOT.gROOT.SetStyle('Plain')
   ROOT.gStyle.SetOptStat(False);
   ###c1 = ROOT.TCanvas("c1", "c1",600,600);
   ###c1.Divide(2,2);
   ###c1.cd(1);
   Gmet.SetMarkerStyle(20);
   Gmet.SetTitle("MET");
   ###Gmet.Draw("APC");
   Gmet.GetXaxis().SetTitle("mass (GeV/c^{2})");
   Gmet.GetYaxis().SetTitle("met cut");

   ###c1.cd(2);
   Gbalance.SetMarkerStyle(20);
   Gbalance.SetTitle("|1-PFMET/qT|");
   ###Gbalance.Draw("APC");
   Gbalance.GetXaxis().SetTitle("mass (GeV/c^{2})");
   Gbalance.GetYaxis().SetTitle("Balance cut");

   ###c1.cd(3);
   Gdphi.SetMarkerStyle(20);
   Gdphi.SetTitle("dphi(Z,PFMET)");
   ###Gdphi.Draw("APC");
   Gdphi.GetXaxis().SetTitle("mass (GeV/c^{2})");
   Gdphi.GetYaxis().SetTitle("Dphi cut");

   #c1.cd(4);
   #Gzm.SetMarkerStyle(20);
   #Gzm.SetTitle("ZMass");
   #Gzm.Draw("APC");
   #Gzm.GetXaxis().SetTitle("m_{H} (GeV/c^{2})");
   #Gzm.GetYaxis().SetTitle("zmass cut");

   ###c1.cd(0);
   ###c1.Update();
   ###print OUT
   ###c1.SaveAs(OUT+"/OptimizedCuts.png")

   #run limits for the cuts chosen (for intermediate masses use spline interpolation)
   for m in SUBMASS:
        index = findCutIndex(Gmet.Eval(m,0,""), cuts1, Gbalance.Eval(m,0,""), cuts2,  Gdphi.Eval(m,0,""), cuts3 );# Gzm.Eval(m,0,""), cuts4);
        print("mass="+str(m).rjust(3)+ " pfMET>"+str(cuts1.GetBinContent(index)).rjust(5) + " |1-PFMET/qT|<" + str(cuts2.GetBinContent(index)).rjust(5) + "  dphi(Z,PFMET)>"+str(cuts3.GetBinContent(index)).rjust(5) )#+ "  |zll-zmass|<"+str(cuts4.GetBinContent(index)).rjust(5) )

   #while True:
        #ans = raw_input('Use this fit and compute final limits? (y or n)\n')
        #if(ans=='y' or ans == 'Y'): break;
        #else:			    sys.exit(0);
   #print 'YES'


   list = open(OUT+'list.txt',"w")
   listcuts = open(OUT+'cuts.txt',"w")
   pos = 0
   for m in SUBMASS:

      index = findCutIndex(Gmet.Eval(m,0,""), cuts1, Gbalance.Eval(m,0,""), cuts2,  Gdphi.Eval(m,0,""), cuts3);# Gzm.Eval(m,0,""), cuts4);
      SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
      SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
      SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
      SCRIPT.writelines("eval `scram r -sh`;\n")
      SCRIPT.writelines('cd ' + CWD + ';\n')
      shapeBasedOpt=''
      if(shapeBased=='1') : shapeBasedOpt='--shape'
#
#        SideMassesArgs = ' '
#        SideMasses = findSideMassPoint(m)
#        if(not (SideMasses[0]==SideMasses[1])):
#           print "Side Mass for mass " + str(m) + " are " + str(SideMasses[0]) + " and " + str(SideMasses[1])
#           Lindex = findCutIndex(Gmet.Eval(SideMasses[0],0,""), cuts1, Gbalance.Eval(SideMasses[0],0,""), cuts2,  Gdphi.Eval(SideMasses[0],0,""), cuts3);# Gzm.Eval(SideMasses[0],0,""), cuts4);
#           Rindex = findCutIndex(Gmet.Eval(SideMasses[1],0,""), cuts1, Gbalance.Eval(SideMasses[1],0,""), cuts2,  Gdphi.Eval(SideMasses[1],0,""), cuts3);# Gzm.Eval(SideMasses[1],0,""), cuts4);
#           print "cutIndex for sideBand are " + str(Lindex) + " and " + str(Rindex)
#           SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + str(Lindex) +  " --indexR " + str(Rindex) + " "
#

      cardsdir=DataCardsDir+"/"+str(m);
      SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
      ## for SigmaxBr/Sigma_{SM Higgs}
      if(shapeName=="coslZ1Dshape"):
         SCRIPT.writelines("runLimits --shapeMin -1 --m " + str(m) + " --histo " + shapeName + " --in " + inUrl + " " + " --syst " + shapeBasedOpt + " --index " + str(index) + " --json " + jsonUrl + " " + LandSArg +" ;\n")
      else:
         #SCRIPT.writelines("runLimits --m " + str(m) + " --metEXTRAPOL "+ str(cuts1.GetBinContent(index)).rjust(5) + " --respEXTRAPOL " + str(cuts2.GetBinContent(index)).rjust(5) + " --dphiEXTRAPOL " + str(cuts3.GetBinContent(index)).rjust(5) + " --histo " + shapeName + " --in " + inUrl + " " + " --syst " + shapeBasedOpt + " --index " + str(index) + " --json " + jsonUrl + " " + LandSArg +" ;\n")
         SCRIPT.writelines("runLimits --m " + str(m) + " --histo " + shapeName + " --in " + inUrl + " " + " --syst " + shapeBasedOpt + " --index " + str(index) + " --json " + jsonUrl + " " + LandSArg +" ;\n")


      SCRIPT.writelines("sh combineCards.sh;\n")
      #####SCRIPT.writelines("combine -M Asymptotic -m " + str(m) + " card_combined.dat > COMB_asympt.log;\n") #unblind data
      #####SCRIPT.writelines("combine -M Asymptotic -m " + str(m) + " card_combinedLL.dat > COMB_asymptLL.log;\n") #unblind data
      ####SCRIPT.writelines("combine -M Asymptotic -m " + str(m) + " --run expected --expectSignal=1 -t -1 card_combined.dat > COMB_asympt.log;\n")
      #SCRIPT.writelines("combine -M MaxLikelihoodFit -m " + str(m) + " --saveNormalizations --plots card_combined.dat > COMB_maxlkl.log;\n")
      #SCRIPT.writelines("combine -M MaxLikelihoodFit -m " + str(m) + " --plot --rMin=-2 --rMax=4 --robustFit=1 --X-rtd FITTER_DYN_STEP -n HZ card_combined.dat > COMB_maxlkl.log;\n")
      #SCRIPT.writelines("combine -M HybridNew --frequentist  -m 0 --testStat LHC card_combined.dat -H ProfileLikelihood --fork 4 > COMB_hybrid.log;\n")

      ## Blinded 
      #SCRIPT.writelines("combine -M Asymptotic --cl " + str(cl) + " --rRelAcc 0.00000001 --rAbsAcc 0.000000001 -m 1 --run expected --expectSignal=1 -t -1 card_combined.dat > COMB_asympt.log;\n") 
      ## Unblinded 
      SCRIPT.writelines("combine -M Asymptotic --cl " + str(cl) + " --rRelAcc 0.00000001 --rAbsAcc 0.000000001 -m 1 card_combined.dat > COMB_asympt.log;\n") 

      SCRIPT.writelines('cd ..;\n\n')
      SCRIPT.close()
      #os.system("bsub -G u_zh -q 2nd -o "+LOGFILE+" 'sh " + OUT+"script_mass_"+str(m)+".sh'")
      os.system("sh " + OUT+"script_mass_"+str(m)+".sh")
      if(shapeBased=='1'):   list.writelines('H'+str(m)+'_shape_'+str(index)+'\n');
      else:                  list.writelines('H'+str(m)+'_count_'+str(index)+'\n');
      listcuts.writelines(str(m)+' ' + str(Gmet.Eval(m,0,"")) + ' ' + str(Gbalance.Eval(m,0,""))+' '+ str(Gdphi.Eval(m,0,""))+'\n');#' '+ str(Gzm.Eval(m,0,"")) +'\n');
      pos += 1
   list.close();
   listcuts.close();

######################################################################

elif(phase == 4 ):

   print '#            #'
   print '# FINAL PLOT #'
   print '#            #'
   ouputName = 'COUNT'
   if(shapeBased=='1'):ouputName='SHAPE'

   os.system("hadd -f "+ouputName+"_LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.Asymptotic.*.root")
   os.system("root -l -b -q plotLimit.C++'(\""+ouputName+"\",\""+ouputName+"_LimitTree.root\", 7 , 5.035 )'")

######################################################################

elif(phase == 5):

   print '#                                      #'
   print '# SCANNING ALL SETS OF CUTS (ASYMPT)   #'
   print '# (you may want to go for a coffee...) #'
   print '#                                      #'

   fileName = OUT + "/OPTIM_"
   if(shapeBased=='1'):  fileName+='SHAPE'
   else:                 fileName+='COUNT'

   FILE = open(fileName+".txt","w")
   for m in MASS:
      print 'Starting mass ' + str(m)
      FILE.writelines("------------------------------------------------------------------------------------\n")
      BestLimit = []
      fileList = commands.getstatusoutput("ls " + OUT + str(m)+"_*.log")[1].split();
      for f in fileList:
         exp = commands.getstatusoutput("cat " + f + " | grep \"Best fit r: \"")[1];
         if(len(exp)<=0):continue
         median = float(exp.split()[3])
         uncM = float(exp.split()[4].split('/')[0])
         uncP = float(exp.split()[4].split('/')[1])
         unc = (((median+uncP) + (median-uncM))/(2*median))-1
         if(float(median)<=0.0):continue
         index = int(f[f.rfind("_")+1:f.rfind(".log")])
         BestLimit.append("mH="+str(m)+ " --> " + ('%07.3f%%' % float(100.0*unc)) + " " + ('%07.3fpb +- [%07.3f%%/%07.3f%%]' % (float(median*5.9) , 100.0*uncM/median , 100.0*uncP/median) ) + " " + str(index).rjust(5) + " " + str(cuts1.GetBinContent(index)).rjust(5) + " " + str(cuts2.GetBinContent(index)).rjust(5) + " " + str(cuts3.GetBinContent(index)).rjust(5) + " " + str(cuts4.GetBinContent(index)).rjust(5))

      #sort the limits for this mass
      BestLimit.sort()
      for s in BestLimit:
         FILE.writelines(s+"\n")

   #all done
   FILE.close()
   print("file "+fileName+".txt is written: it contains all selection points ordered by exp limit")

######################################################################

else:

   help()



