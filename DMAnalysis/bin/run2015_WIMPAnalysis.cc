#include <iostream>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "llvvAnalysis/DMAnalysis/interface/MacroUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"
#include "llvvAnalysis/DMAnalysis/interface/DMPhysicsEvent.h"
#include "llvvAnalysis/DMAnalysis/interface/SmartSelectionMonitor.h"
#include "llvvAnalysis/DMAnalysis/interface/METUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/BTagUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/EventCategory.h"
#include "llvvAnalysis/DMAnalysis/interface/LeptonEfficiencySF.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"


#include "llvvAnalysis/DMAnalysis/interface/PDFInfo.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
#include "llvvAnalysis/DMAnalysis/interface/BTagCalibrationStandalone.h"


using namespace std;

int main(int argc, char* argv[])
{
    //##################################################################################
    //##########################    GLOBAL INITIALIZATION     ##########################
    //##################################################################################

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    // configure the process
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

    bool isMC       = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");

    TString url=runProcess.getParameter<std::string>("input");
    TString outFileUrl(gSystem->BaseName(url));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);


    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    if ( runProcess.getUntrackedParameter<bool>("saveFinalList", true) ) {
        outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
        printf("TextFile URL = %s\n",outTxtUrl_final.Data());
        fprintf(outTxtFile_final,"run lumi event\n");
    }

    // Embed final list in root file
    int eventList_run, eventList_lumi, eventList_evt, eventList_nJets, eventList_category;
    float eventList_met, eventList_mt, eventList_llpt, eventList_mll;
    bool saveEventList = !isMC || runProcess.getUntrackedParameter<bool>("saveEventList", false);
    TTree * eventList{nullptr};
    if ( saveEventList ) {
        eventList = new TTree("eventList", "List of passed events");
        eventList->Branch("run", &eventList_run);
        eventList->Branch("lumi", &eventList_lumi);
        eventList->Branch("evt", &eventList_evt);
        eventList->Branch("category", &eventList_category);
        eventList->Branch("nJets", &eventList_nJets);
        eventList->Branch("met", &eventList_met);
        eventList->Branch("mt", &eventList_mt);
        eventList->Branch("llpt", &eventList_llpt);
        eventList->Branch("mll", &eventList_mll);
    }

    int fType(0);
    if(url.Contains("DoubleEG")) fType=EE;
    if(url.Contains("DoubleMuon"))  fType=MUMU;
    if(url.Contains("MuonEG"))      fType=EMU;
    if(url.Contains("SingleMuon"))  fType=MUMU;
    if(url.Contains("SingleElectron")) fType=EE;
    bool isSingleMuPD(!isMC && url.Contains("SingleMuon"));
    bool isDoubleMuPD(!isMC && url.Contains("DoubleMuon"));
    bool isSingleElePD(!isMC && url.Contains("SingleElectron"));
    bool isDoubleElePD(!isMC && url.Contains("DoubleEG"));

    bool isMC_ZZ     = isMC && (string(url.Data()).find("TeV_ZZ") != string::npos) && (string(url.Data()).find("TeV_ZZZ") == string::npos);
    bool isMC_WZ     = isMC && (string(url.Data()).find("TeV_WZ") != string::npos) && (string(url.Data()).find("TeV_WZZ") == string::npos);
    //bool isMC_WW     = isMC && (string(url.Data()).find("TeV_WW_")  != string::npos);
    bool isMC_ttbar  = isMC && (string(url.Data()).find("TeV_TT")  != string::npos);
    bool isMC_stop   = isMC && (string(url.Data()).find("TeV_SingleT")  != string::npos);
    //bool isMC_WIMP   = isMC && (string(url.Data()).find("TeV_DM_V_Mx") != string::npos
    //                 || string(url.Data()).find("TeV_DM_A_Mx") != string::npos);
    //bool isMC_ADD    = isMC && (string(url.Data()).find("TeV_ADD_D") != string::npos);
    //bool isMC_Unpart = isMC && (string(url.Data()).find("TeV_Unpart") != string::npos);

    bool isMC_ZH     = isMC && (string(url.Data()).find("TeV_ZH") != string::npos);

    //bool isSignal = (isMC_WIMP || isMC_ADD || isMC_Unpart);
    bool isSignal    = isMC_ZH;


    // https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
    // setup calibration readers
    TString BtagSF = runProcess.getParameter<std::string>("BtagSF");
    gSystem->ExpandPathName(BtagSF);
    cout << "Loading btag sacle factor: " << BtagSF << endl;
    BTagCalibration btagcalib("csvv2", BtagSF.Data());
    BTagCalibrationReader btag_reader(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "central");
    BTagCalibrationReader btag_reader_up(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    BTagCalibrationReader btag_reader_down(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down

    // EWK corrections (from table or from plot) 
    bool useEwkTable = false; 

    // Use MET phi correction?
    bool useMETPhiCorr = runProcess.getUntrackedParameter<bool>("useMETPhiCorr", false);

    // Systematics
    bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
    std::vector<TString> varNames(1,"");
    if(runSystematics) {
        cout << "Systematics will be computed for this analysis" << endl;
        varNames.push_back("_jerup"); 	//1
        varNames.push_back("_jerdown"); //2
        varNames.push_back("_jesup"); 	//3
        varNames.push_back("_jesdown"); //4
        varNames.push_back("_umetup"); 	//5
        varNames.push_back("_umetdown");//6
        varNames.push_back("_lesup"); 	//7
        varNames.push_back("_lesdown"); //8
        varNames.push_back("_puup"); 	//9
        varNames.push_back("_pudown"); 	//10
        varNames.push_back("_btagup"); 	//11
        varNames.push_back("_btagdown");//12
        if(isSignal || isMC_ZZ || isMC_WZ)            {
            varNames.push_back("_pdfup");
            varNames.push_back("_pdfdown");
        }
        if(isSignal || isMC_ZZ || isMC_WZ) {
	  varNames.push_back("_qcdup");
	  varNames.push_back("_qcddown");
        }
	if(isMC_ZZ) { 
            varNames.push_back("_ewkup"); 
            varNames.push_back("_ewkdown"); 
	} 
    }
    size_t nvarsToInclude=varNames.size();

    // Print debug information 
    bool debug(false); 

    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    //jet energy scale uncertainties
    TString uncFile = runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    JetCorrectionUncertainty jecUnc(uncFile.Data());



    //pdf info
    // PDFInfo *mPDFInfo=0;
    // if(isSignal) {
    // 	TString pdfUrl = runProcess.getParameter<std::string>("pdfInput");
    // 	std::string Url = runProcess.getParameter<std::string>("input");
    // 	std::size_t found = Url.find_last_of("/\\");
    //     pdfUrl += '/';
    // 	pdfUrl += Url.substr(found+1);
    //     pdfUrl.ReplaceAll(".root","_pdf.root");

    //     mPDFInfo=new PDFInfo(pdfUrl,"cteq66.LHgrid");
    //     cout << "Readout " << mPDFInfo->numberPDFs() << " pdf variations: " << pdfUrl << endl;
    // }


    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 10,0,10) );
    h->GetXaxis()->SetBinLabel(1,"Trigger && 2 leptons");
    h->GetXaxis()->SetBinLabel(2,"|#it{m}_{ll}-#it{m}_{Z}|<15");
    h->GetXaxis()->SetBinLabel(3,"#it{p}_{T}^{ll}>60");
    h->GetXaxis()->SetBinLabel(4,"3^{rd}-lepton veto");
    h->GetXaxis()->SetBinLabel(5,"b-veto");
    h->GetXaxis()->SetBinLabel(6,"#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss})>2.8");
    h->GetXaxis()->SetBinLabel(7,"#Delta#it{#phi}(#it{l^{+},l^{-}})<#pi/2");
    h->GetXaxis()->SetBinLabel(8,"|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T}<0.4");
    h->GetXaxis()->SetBinLabel(9,"E_{T}^{miss}>100");
    h->GetXaxis()->SetBinLabel(10,"m_{T} > 200");

    mon.addHistogram((TH1F*) h->Clone("eventflow_unweighted"));


    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow", "cutflow", 9, -0.5, 8.5) ) ;

    mon.addHistogram( new TH1F( "nvtx_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "zpt_raw",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_wocorr_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_noHF_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_puppi_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "zmass_raw",    ";#it{m}_{ll} [GeV];Events", 100,40,250) );

    mon.addHistogram( new TH2F( "ptlep1vs2_raw",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 50,-2.6,2.6) );
    mon.addHistogram( new TH1F( "trailep_pt_raw", ";Trailing lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "trailep_eta_raw",";Trailing lepton #eta^{l};Events", 50,-2.6,2.6) );

    mon.addHistogram( new TH1F( "el_reliso_raw",  ";Electron pfRelIsoDbeta;Events",50,0,0.5) );
    mon.addHistogram( new TH1F( "mu_reliso_raw",  ";Muon pfRelIsoDbeta;Events",50,0,2) );

    mon.addHistogram( new TH1F( "leadingjet_pt_raw", ";Leading jet #it{p}_{T}^{j};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "jet_pt_raw", ";all jet #it{p}_{T}^{j};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "jet_eta_raw",";all jet #eta^{j};Events", 50,-2.6,2.6) );


    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nleptons_raw", ";Lepton multiplicity;Events", 3,2,5) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="=";
        label += (ibin+1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }
    mon.addHistogram( (TH1F*) h1->Clone("nGoodleptons_raw") );

    TH1F *h2 = (TH1F *)mon.addHistogram( new TH1F("njets_raw",  ";Jet multiplicity (#it{p}_{T}>30 GeV);Events",5,0,5) );
    for(int ibin=1; ibin<=h2->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h2->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h2->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h3 = (TH1F *)mon.addHistogram( new TH1F("nbjets_raw",    ";b-tag Jet multiplicity;Events",5,0,5) );
    for(int ibin=1; ibin<=h3->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h3->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h3->GetXaxis()->SetBinLabel(ibin,label);
    }

    mon.addHistogram( new TH1F("nSoftMuons_raw", ";Soft muon multiplicity;Events", 4, 0, 4) );

    // Cut flow plots
    mon.addHistogram( new TH1F( "zpt_flow",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( (TH1F*) h1->Clone("nleptons_flow") );
    mon.addHistogram( (TH1F*) h3->Clone("nbjets_flow") );
    mon.addHistogram( new TH1F("nSoftMuons_flow", ";Soft muon multiplicity;Events", 4, 0, 4) );

    // preselection plots
    double METBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
    const int nBinsMET = sizeof(METBins)/sizeof(double) - 1;

    double METBins2[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000};
    const int nBinsMET2 = sizeof(METBins2)/sizeof(double) - 1;

    double MTBins[]= {200,250, 300,400,600,800,1000};
    const int nBinsMT = sizeof(MTBins)/sizeof(double) - 1;

    mon.addHistogram( new TH1F( "zmass_presel",    ";#it{m}_{ll} [GeV];Events", 50,91-15,91+15) );
    mon.addHistogram( new TH1F( "zpt_presel",    ";#it{p}_{T}^{ll} [GeV];Events / 10 GeV", 45,50,500) );
    mon.addHistogram( new TH1F( "pfmet_presel",      ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET, METBins));
    mon.addHistogram( new TH1F( "pfmet2_presel",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, METBins2));
    mon.addHistogram( new TH1F( "pfmetCtrl_presel",     ";E_{T}^{miss} [GeV];Events / 5 GeV", 20, 0, 100));
    mon.addHistogram( new TH1F( "pfmetCtrl_wocorr_presel",     ";E_{T}^{miss} [GeV];Events / 5 GeV", 20, 0, 100));
    mon.addHistogram( new TH1F( "pfmetCtrl_noHF_presel",     ";E_{T}^{miss} [GeV];Events / 5 GeV", 20, 0, 100));
    mon.addHistogram( new TH1F( "pfmetCtrl_puppi_presel",     ";E_{T}^{miss} [GeV];Events / 5 GeV", 20, 0, 100));
    mon.addHistogram( new TH1F( "pfmetParallelZ_presel",     ";E_{T}^{miss} Parallel Z [GeV];Events / 10 GeV", 40, -200, 200));
    mon.addHistogram( new TH1F( "pfmetPerpZ_presel",         ";E_{T}^{miss} Perp. Z [GeV];Events / 5 GeV", 40, -100, 100));
    mon.addHistogram( new TH1F( "dphiZMET_presel",   ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphiLL_presel",     ";#Delta#it{#phi}(#it{l^{+},l^{-}});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "balancedif_presel", ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 5,0,1.0) );
    mon.addHistogram( new TH1F( "mt_presel",         ";#it{m}_{T} [GeV];Events", 12,0,1200) );

    //MET X-Y shift correction
    mon.addHistogram( new TH2F( "pfmetx_vs_nvtx_presel",";Vertices;E_{X}^{miss} [GeV];Events",50,0,50, 200,-75,75) );
    mon.addHistogram( new TH2F( "pfmety_vs_nvtx_presel",";Vertices;E_{Y}^{miss} [GeV];Events",50,0,50, 200,-75,75) );
    mon.addHistogram( new TH2F( "pfmetParallelZ_vs_nvtx_presel",";Vertices;E_{parallel}^{miss} [GeV];Events",50,0,50, 200,-200,200) );
    mon.addHistogram( new TH2F( "pfmetPerpZ_vs_nvtx_presel",";Vertices;E_{perp.}^{miss} [GeV];Events",50,0,50, 200,-100,100) );
    mon.addHistogram( new TH1F( "pfmetphi_wocorr_presel",";#it{#phi}(E_{T}^{miss});Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pfmetphi_wicorr_presel",";#it{#phi}(E_{T}^{miss});Events", 50,-1.*TMath::Pi(),TMath::Pi()) );


    // generator level plots
    mon.addHistogram( new TH1F( "met_Gen", ";#it{p}_{T}(#bar{#chi}#chi) [GeV];Events", 100,0,800) );
    mon.addHistogram( new TH1F( "zpt_Gen", ";#it{p}_{T}(Z) [GeV];Events", 800,0,800) );
    mon.addHistogram( new TH1F( "dphi_Gen", ";#Delta#phi(Z,#bar{#chi}#chi) [rad];Events", 100,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "zmass_Gen", ";#it{m}_{ll} [GeV] [GeV];Events", 250,0,250) );
    mon.addHistogram( new TH2F( "ptlep1vs2_Gen",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );


    // btaging efficiency
    std::vector<TString> CSVkey;
    CSVkey.push_back("CSVL");
    CSVkey.push_back("CSVM");
    CSVkey.push_back("CSVT");
    for(size_t csvtag=0; csvtag<CSVkey.size(); csvtag++) {
        mon.addHistogram( new TH1F( TString("beff_Denom_")+CSVkey[csvtag],      "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("ceff_Denom_")+CSVkey[csvtag],      "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("udsgeff_Denom_")+CSVkey[csvtag],   "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("beff_Num_")+CSVkey[csvtag],        "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("ceff_Num_")+CSVkey[csvtag],        "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("udsgeff_Num_")+CSVkey[csvtag],     "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
    }


    mon.addHistogram( new TH1F( "mt_final",             ";#it{m}_{T} [GeV];Events", nBinsMT, MTBins) );
    mon.addHistogram( new TH1F( "mt_final120",             ";#it{m}_{T} [GeV];Events", nBinsMT, MTBins) );
    mon.addHistogram( new TH1F( "pfmet_final",      ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET, METBins));
    mon.addHistogram( new TH1F( "pfmet2_final",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, METBins2));


    TH1F *hfinalWeights=(TH1F*) mon.addHistogram( new TH1F ("weights_final", ";;Events", 6,0,6) );
    hfinalWeights->GetXaxis()->SetBinLabel(1,"Raw");
    hfinalWeights->GetXaxis()->SetBinLabel(2,"*(gen sign)");
    hfinalWeights->GetXaxis()->SetBinLabel(3,"*(ewk)");
    hfinalWeights->GetXaxis()->SetBinLabel(4,"*(PU)");
    hfinalWeights->GetXaxis()->SetBinLabel(5,"*(btagSF)");
    hfinalWeights->GetXaxis()->SetBinLabel(6,"*(lepSF)");

    mon.addHistogram( (TH1F*) hfinalWeights->Clone("weights_presel") );

    //#################################################
    //############# CONTROL PLOTS #####################
    //#################################################
    // WW/ttbar/Wt/tautau control plots, for k-method (for emu channel)
    mon.addHistogram( new TH1F( "zpt_wwctrl_raw",   ";#it{p}_{T}^{ll} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "zmass_wwctrl_raw", ";#it{m}_{ll} [GeV];Events", 100,20,300) );
    mon.addHistogram( new TH1F( "pfmet_wwctrl_raw", ";E_{T}^{miss} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "mt_wwctrl_raw",";#it{m}_{T}(#it{ll}, E_{T}^{miss}) [GeV];Events", 50,0,300) );






    //##################################################################################
    //########################## STUFF FOR CUTS OPTIMIZATION  ##########################
    //##################################################################################
    //optimization
    std::vector<double> optim_Cuts1_MET;
    std::vector<double> optim_Cuts1_Balance;
    std::vector<double> optim_Cuts1_DphiZMET;

    bool runOptimization = runProcess.getParameter<bool>("runOptimization");
    if(runOptimization) {
        // for optimization
        cout << "Optimization will be performed for this analysis" << endl;
        for(double met=80; met<=120; met+=10) {
            for(double balance=0.2; balance<=0.4; balance+=0.1) {
                for(double dphi=2.7; dphi<2.9; dphi+=0.1) {
                    optim_Cuts1_MET     .push_back(met);
                    optim_Cuts1_Balance .push_back(balance);
                    optim_Cuts1_DphiZMET.push_back(dphi);
                }
            }
        }
    }

    size_t nOptims = optim_Cuts1_MET.size();


    //make it as a TProfile so hadd does not change the value
    TProfile* Hoptim_cuts1_MET      = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_MET",";cut index;met",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_Balance  = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_Balance",";cut index;balance",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_DphiZMET = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_DphiZMET",";cut index;dphi",nOptims,0,nOptims) );

    for(size_t index=0; index<nOptims; index++) {
        Hoptim_cuts1_MET        ->Fill(index, optim_Cuts1_MET[index]);
        Hoptim_cuts1_Balance    ->Fill(index, optim_Cuts1_Balance[index]);
        Hoptim_cuts1_DphiZMET   ->Fill(index, optim_Cuts1_DphiZMET[index]);
    }

    TH1F* Hoptim_systs  =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;



    for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
        Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);

        //1D shapes for limit setting
        mon.addHistogram( new TH2F (TString("mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events",nOptims,0,nOptims,nBinsMT,MTBins) );
        //mon.addHistogram( new TH2F (TString("pfmet_shapes")+varNames[ivar],";cut index; E_{T}^{miss} [GeV];#Events",nOptims,0,nOptims,nBinPFMET,xbinsPFMET) );

        //2D shapes for limit setting
        //
        //
    }

    //##################################################################################
    //#############         GET READY FOR THE EVENT LOOP           #####################
    //##################################################################################

    //open the file and get events tree
    DataEvtSummaryHandler summaryHandler_;
    TFile *file = TFile::Open(url);
    printf("Looping on %s\n",url.Data());
    if(file==0) return -1;
    if(file->IsZombie()) return -1;
    if( !summaryHandler_.attachToTree( (TTree *)file->Get(dirname) ) ) {
        file->Close();
        return -1;
    }


    //check run range to compute scale factor (if not all entries are used)
    const Int_t totalEntries= summaryHandler_.getEntries();
    float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
    if(evEnd<0 || evEnd>summaryHandler_.getEntries() ) evEnd=totalEntries;
    if(evStart > evEnd ) {
        file->Close();
        return -1;
    }

    //MC normalization (to 1/pb)
    float cnorm=1.0;
    if(isMC) {
        //TH1F* cutflowH = (TH1F *) file->Get("mainAnalyzer/llvv/nevents");
        //if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        TH1F* posH = (TH1F *) file->Get("mainAnalyzer/llvv/n_posevents");
        TH1F* negH = (TH1F *) file->Get("mainAnalyzer/llvv/n_negevents");
        if(posH && negH) cnorm = posH->GetBinContent(1) - negH->GetBinContent(1);
        if(rescaleFactor>0) cnorm /= rescaleFactor;
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);

    //pileup weighting
    TString PU_Central = runProcess.getParameter<std::string>("PU_Central");
    gSystem->ExpandPathName(PU_Central);
    cout << "Loading PU weights Central: " << PU_Central << endl;
    TFile *PU_Central_File = TFile::Open(PU_Central);
    TH1F* weight_pileup_Central = (TH1F *) PU_Central_File->Get("pileup");

    TString PU_Up = runProcess.getParameter<std::string>("PU_Up");
    gSystem->ExpandPathName(PU_Up);
    cout << "Loading PU weights Up: " << PU_Up << endl;
    TFile *PU_Up_File = TFile::Open(PU_Up);
    TH1F* weight_pileup_Up = (TH1F *) PU_Up_File->Get("pileup");

    TString PU_Down = runProcess.getParameter<std::string>("PU_Down");
    gSystem->ExpandPathName(PU_Down);
    cout << "Loading PU weights Down: " << PU_Down << endl;
    TFile *PU_Down_File = TFile::Open(PU_Down);
    TH1F* weight_pileup_Down = (TH1F *) PU_Down_File->Get("pileup");

    // event categorizer
    EventCategory eventCategoryInst(1);   //jet(0,1,>=2) binning

    // Lepton scale factors 
    LeptonEfficiencySF lsf(2015);


    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################

    // loop on all the events
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    DuplicatesChecker duplicatesChecker;
    int nDuplicates(0);
    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
            printf("#");
            fflush(stdout);
        }

        //##############################################   EVENT LOOP STARTS   ##############################################
        //load the event content from tree
        summaryHandler_.getEntry(iev);
        DataEvtSummary_t &ev=summaryHandler_.getEvent();
        if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) {
            nDuplicates++;
            cout << "nDuplicates: " << nDuplicates << endl;
            continue;
        }


        //prepare the tag's vectors for histo filling
        std::vector<TString> tags(1,"all");

        //genWeight
        float genWeight = 1.0;
        if(isMC && ev.genWeight<0) genWeight = -1.0;

        // total event weight, accrued over time
        float weight = 1.0;
        if(isMC) weight *= genWeight;

        // pielup reweightiing
        double puWeight = 1.;
        double puWeight_plus = 1.;
        double puWeight_minus = 1.;

        if(isMC) {
            float xval = ev.ngenTruepu;
            int xbins = weight_pileup_Central->GetXaxis()->GetNbins();
            if     (xval > weight_pileup_Central->GetXaxis()->GetBinUpEdge(xbins)    ) xval = weight_pileup_Central->GetXaxis()->GetBinUpEdge(xbins);
            else if(xval < weight_pileup_Central->GetXaxis()->GetBinLowEdge(1)       ) xval = weight_pileup_Central->GetXaxis()->GetBinLowEdge(1);

            int binx = weight_pileup_Central->GetXaxis()->FindBin(xval);

            puWeight       *= weight_pileup_Central->GetBinContent(binx);
            puWeight_plus  *= weight_pileup_Up->GetBinContent(binx);
            puWeight_minus *= weight_pileup_Down->GetBinContent(binx);
        }

        Hcutflow->Fill(1,weight);
        Hcutflow->Fill(2,weight*puWeight);
        Hcutflow->Fill(3,weight*puWeight_minus);
        Hcutflow->Fill(4,weight*puWeight_plus);

        weight *= puWeight;


        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        float ewk_w = 1.; 
        /// Also include an extra 10% to account for gg->ZZ contribution. 
        ///  N.B.:  sigma(gg->ZZ) = sigma(qq->ZZ) BEFORE NLO EWK contribution!!! 
        ///  Therefore it's not  sigma(qq->ZZ) * (1 + ewk_corr) * (1 + gg_contr) 
        ///  but rather  
        ///    sigma(qq->ZZ) + ewk_corr*sigma(qq->ZZ) + gg_contr*sigma(qq->ZZ) 
        ///    = sigma(qq->ZZ) * (1 + ewk_corr + gg_contr) 
        ///  where gg_contr = 0.1 
        float ggZZ_contr = 0.1; 
        
        /////// 
        // EWK correction for ZZ and WZ (ewk_w = 1.0 otherwise) 
        /// 
        // WZ 
        if(false && // turn them off for now (probably negligible) 
           isMC && (url.Contains("MC13TeV_WZ")) && (!url.Contains("MC13TeV_WZZ"))) {
          TLorentzVector wz_z, wz_w;
          TLorentzVector 
            nu(0., 0., 0., 0.), l1(0., 0., 0., 0.), 
            l2(0., 0., 0., 0.), l3(0., 0., 0., 0.), 
            lw(0., 0., 0., 0.), lz1(0., 0., 0., 0.), lz2(0., 0., 0., 0.); 
          int nuid(0), l1id(0), l2id(0), l3id(0); 
          bool foundNeut(false), foundLep1(false), foundLep2(false), foundLep3(false); 
          //std::cout << std::endl; 
          //std::cout << " New Event:" << std::endl; 
          for(Int_t ipart=0; ipart<ev.nmcparticles; ++ipart) { 
            if(ev.mc_status[ipart]!=1) continue; 
            //std::cout << " " << ev.mc_id[ipart] << "/" << ev.mc_status[ipart]; 
            //continue; 
        
            if( !foundNeut && (fabs(ev.mc_id[ipart])==12  || fabs(ev.mc_id[ipart])==14  || fabs(ev.mc_id[ipart])==16 ) ) { 
              nu.SetPxPyPzE(ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              nuid = ev.mc_id[ipart]; 
              foundNeut = true; 
            } 
            else if( !foundLep1 && (fabs(ev.mc_id[ipart])==11  || fabs(ev.mc_id[ipart])==13  || fabs(ev.mc_id[ipart])==15 ) ) { 
              l1.SetPxPyPzE(ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              l1id = ev.mc_id[ipart]; 
              foundLep1 = true; 
            } 
            else if( foundLep1 && !foundLep2 && (fabs(ev.mc_id[ipart])==11  || fabs(ev.mc_id[ipart])==13  || fabs(ev.mc_id[ipart])==15 ) ) { 
              l2.SetPxPyPzE(ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              l2id = ev.mc_id[ipart]; 
              foundLep2 = true; 
            } 
            else if( foundLep1 && foundLep2 && !foundLep3 && (fabs(ev.mc_id[ipart])==11  || fabs(ev.mc_id[ipart])==13  || fabs(ev.mc_id[ipart])==15 ) ) { 
              l3.SetPxPyPzE(ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              l3id = ev.mc_id[ipart]; 
              foundLep3 = true; 
            } 
            if(foundNeut && foundLep3) break; 
          } // end for nmcparticles 
        
          if(nuid!=0 && l1id!=0 && l2id!=0 && l3id!=0) { 
            // Now find the right combinations... good luck! 
            if( (nuid+l1id)==(nuid/abs(nuid)) && (l2id+l3id)==0  ) { 
              lw = l1; lz1 = l2; lz2 = l3; 
            } 
            if( (nuid+l2id)==(nuid/abs(nuid)) && (l1id+l3id)==0  ) { 
              if( fabs((nu+l2).M()-80.) < fabs((nu+lw).M()-80.) ) { 
        	lw = l2; lz1 = l1; lz2 = l3; 
              } 
            } 
            if( (nuid+l3id)==(nuid/abs(nuid)) && (l1id+l2id)==0  ) { 
              if( fabs((nu+l3).M()-80.) < fabs((nu+lw).M()-80.) ) { 
        	lw = l3; lz1 = l1; lz2 = l2; 
              } 
            } 
            wz_z = lz1 + lz2; wz_w = nu + lw; 
        
            if(debug) {
              std::cout << std::endl; 
              std::cout << " *** Z: " << lz1.Pt() << "-" << lz2.Pt() << ";  W: " << nu.Pt() << "-" << lw.Pt() << std::endl; 
              std::cout << "        nu(" << nuid << "): " << nu.Pt() << ", l1(" << l1id << "): " << l1.Pt() << ", l2(" << l2id << "): " 
        		<< l2.Pt() << ", l3(" << l3id << "): " << l3.Pt() << std::endl; 
            } 
        
            float wz_min_pt = wz_z.Pt() < wz_w.Pt() ? wz_z.Pt() : wz_w.Pt(); 
            if( !useEwkTable ) { 
              // Reading by eye from the paper
              if(     wz_min_pt<60.)  ewk_w = 1. - ( 0.9/100.) + ggZZ_contr; 
              else if(wz_min_pt<80.)  ewk_w = 1. - ( 0.9/100.) + ggZZ_contr; 
              else if(wz_min_pt<100.) ewk_w = 1. - ( 1.0/100.) + ggZZ_contr; 
              else if(wz_min_pt<120.) ewk_w = 1. - ( 1.5/100.) + ggZZ_contr; 
              else if(wz_min_pt<140.) ewk_w = 1. - ( 2.0/100.) + ggZZ_contr; 
              else if(wz_min_pt<160.) ewk_w = 1. - ( 2.6/100.) + ggZZ_contr; 
              else if(wz_min_pt<180.) ewk_w = 1. - ( 3.0/100.) + ggZZ_contr; 
              else if(wz_min_pt<200.) ewk_w = 1. - ( 4.9/100.) + ggZZ_contr; 
              else if(wz_min_pt<220.) ewk_w = 1. - ( 5.2/100.) + ggZZ_contr; 
              else if(wz_min_pt<240.) ewk_w = 1. - ( 6.5/100.) + ggZZ_contr; 
              else if(wz_min_pt<260.) ewk_w = 1. - ( 7.5/100.) + ggZZ_contr; 
              else if(wz_min_pt<280.) ewk_w = 1. - ( 8.0/100.) + ggZZ_contr; 
              else if(wz_min_pt<300.) ewk_w = 1. - ( 9.9/100.) + ggZZ_contr; 
              else if(wz_min_pt<320.) ewk_w = 1. - (10.9/100.) + ggZZ_contr; 
              else if(wz_min_pt<340.) ewk_w = 1. - (12.3/100.) + ggZZ_contr; 
              else if(wz_min_pt<360.) ewk_w = 1. - (12.6/100.) + ggZZ_contr; 
              else if(wz_min_pt<380.) ewk_w = 1. - (13.5/100.) + ggZZ_contr; 
              else                    ewk_w = 1. - (14.0/100.) + ggZZ_contr; 
            } // end "if not useEwkTable" 
          } // end "if all gen-particles are identified" 
          else { 
            if(debug) 
              std::cout << " ### WARNING: some gen-particles not identified: nu:" << nuid 
        		<< ", l1:" << l1id << ", l2:" << l2id << ", l3:" << l3id << "###" << std::endl; 
          } 
        } // end "if isMC and is WZ" 
        
        /// 
        // ZZ 
        TLorentzVector l1(0.,0.,0.,0.), l2(0.,0.,0.,0.), v1(0.,0.,0.,0.), v2(0.,0.,0.,0.); 
        if(isMC && (url.Contains("MC13TeV_ZZTo2L2Nu"))) { // No 4l nor 2l2q for now 
          // Neutrinos 
          bool foundNeut1(false), foundNeut2(false); 
          bool foundLept1(false), foundLept2(false); 
          int firstNeut(0), firstLept(0); 
          int v1id(0), v2id(0), l1id(0), l2id(0); 
          for(Int_t ipart=0; ipart<ev.nmcparticles; ++ipart) { 
            if(ev.mc_status[ipart]!=1) continue; 
            if( !foundNeut1 && 
        	(fabs(ev.mc_id[ipart])==12 || fabs(ev.mc_id[ipart])==14 || fabs(ev.mc_id[ipart])==16 ) ) { 
              v1.SetPxPyPzE( ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              foundNeut1 = true; 
              firstNeut = ev.mc_id[ipart]; 
              v1id = ev.mc_id[ipart]; 
            } 
            else if( foundNeut1 && !foundNeut2 && 
        	     //(fabs(ev.mc_id[ipart])==12 || fabs(ev.mc_id[ipart])==14 || fabs(ev.mc_id[ipart])==16) && 
        	     ev.mc_id[ipart] + firstNeut == 0 ) { 
              v2.SetPxPyPzE( ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              foundNeut2 = true; 
              v2id = ev.mc_id[ipart]; 
            } 
            else if( !foundLept1 && 
        	     (fabs(ev.mc_id[ipart])==11 || fabs(ev.mc_id[ipart])==13 || fabs(ev.mc_id[ipart])==15) ) { 
              l1.SetPxPyPzE( ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              foundLept1 = true; 
              firstLept = ev.mc_id[ipart]; 
              l1id = ev.mc_id[ipart]; 
            } 
            else if( foundLept1 && !foundLept2 && 
        	     //(fabs(ev.mc_id[ipart])==11 || fabs(ev.mc_id[ipart])==13 || fabs(ev.mc_id[ipart])==15) && 
        	     ev.mc_id[ipart] + firstLept == 0 ) { 
              l2.SetPxPyPzE( ev.mc_px[ipart], ev.mc_py[ipart], ev.mc_pz[ipart], ev.mc_en[ipart]); 
              foundLept2 = true; 
              l2id = ev.mc_id[ipart]; 
            } 
            if(foundNeut2 && foundLept2) break; 
          } // end for mcparticles 
        
          if(foundNeut2 && foundLept2) { 
            //float z_max_pt=(0.); 
            float z_min_pt(0.);
            float z1_pt = (l1+l2).Pt(); 
            float z2_pt = (v1+v2).Pt();
            //z_max_pt = z1_pt<z2_pt ? z2_pt : z1_pt; 
            z_min_pt = z1_pt>z2_pt ? z2_pt : z1_pt;
        
            if(debug) {
              std::cout << std::endl; 
              std::cout << " *** Zll pt: " << z1_pt << ";  Zvv: " << z2_pt << ";  Zmin: " << z_min_pt << std::endl; 
              std::cout << "        l1(" << l1id << "): " << l1.Pt() << ", l2(" << l2id << "): " << l2.Pt() 
        		<< ", v1(" << v1id << "): " << v1.Pt() << ", v2(" << v2id << "): " << v2.Pt() << std::endl; 
            } 
        
            if( !useEwkTable ) { 
              //Reading by eye from the paper
              if(     z_min_pt<60.)  ewk_w = 1.-(4.0/100.);
              else if(z_min_pt<80.)  ewk_w = 1.-(5.0/100.);
              else if(z_min_pt<100.) ewk_w = 1.-(6.3/100.);
              else if(z_min_pt<120.) ewk_w = 1.-(7.6/100.);
              else if(z_min_pt<140.) ewk_w = 1.-(9.2/100.);
              else if(z_min_pt<160.) ewk_w = 1.-(10.0/100.);
              else if(z_min_pt<180.) ewk_w = 1.-(11.4/100.);
              else if(z_min_pt<200.) ewk_w = 1.-(12.8/100.);
              else if(z_min_pt<220.) ewk_w = 1.-(14.2/100.);
              else if(z_min_pt<240.) ewk_w = 1.-(15.6/100.);
              else if(z_min_pt<260.) ewk_w = 1.-(17.0/100.);
              else if(z_min_pt<280.) ewk_w = 1.-(18.4/100.);
              else if(z_min_pt<300.) ewk_w = 1.-(20.0/100.);
              else if(z_min_pt<320.) ewk_w = 1.-(21.2/100.);
              else if(z_min_pt<340.) ewk_w = 1.-(22.4/100.);
              else if(z_min_pt<360.) ewk_w = 1.-(23.6/100.);
              else if(z_min_pt<380.) ewk_w = 1.-(24.8/100.);
              else                   ewk_w = 1.-(26.0/100.);
            } // end "if not useEwkTable" 
          } // end "foundNeut2 and foundLept2" 
          else { 
            if(debug) 
              std::cout << " ### WARNING: some gen-particles not identified: v1:" << v1id 
        		<< ", v2:" << v2id << ", l1:" << l1id << ", l2:" << l2id << "###" << std::endl; 
          } 
        } // end "if isMC and is ZZ" 
        
        // Apply EWK weight! 
        weight *= ewk_w; 
        
        // End EWK corrections 
        ////// 

        // FIXME need to have a function: loop all leptons, find a Z candidate,
        // can have input, ev.mn, ev.en
        // assign ee,mm,emu channel
        // check if channel name is consistent with trigger
        // store dilepton candidate in lep1 lep2, and other leptons in 3rdleps


        bool hasMMtrigger = ev.triggerType & 0x1;
        bool hasMtrigger  = (ev.triggerType >> 1 ) & 0x1;
        bool hasEEtrigger = (ev.triggerType >> 2 ) & 0x1;
        bool hasEtrigger  = (ev.triggerType >> 3 ) & 0x1;
        bool hasEMtrigger = (ev.triggerType >> 4 ) & 0x1;



        //#########################################################################
        //####################  Generator Level Reweighting  ######################
        //#########################################################################



        //#########################################################################
        //#####################      Objects Selection       ######################
        //#########################################################################

        //
        // MET ANALYSIS
        //
        //apply Jet Energy Resolution corrections to jets (and compute associated variations on the MET variable)
        std::vector<PhysicsObjectJetCollection> variedJets;
        LorentzVectorCollection variedMET;

        //METUtils::computeVariation(phys.jets, phys.leptons, phys.metNoHF, variedJets, variedMET, &jecUnc);
        if ( !isMC ) {
            for(size_t ijet=0; ijet<phys.jets.size(); ijet++) {
                if ( phys.jets[ijet].genPt > 0. )
                    cout << "fuck" << endl;
            }
        }
        METUtils::computeVariation(phys.jets, phys.leptons, phys.met, variedJets, variedMET, &jecUnc);

        LorentzVector metP4=variedMET[0];
        PhysicsObjectJetCollection &corrJets = variedJets[0];

        //
        // LEPTON ANALYSIS
        //

        // looping leptons (electrons + muons)
        int nGoodLeptons(0);
        std::vector<std::pair<int,LorentzVector> > goodLeptons;
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
            if(lep.pt()<20) continue;
            if(abs(lepid)==13 && fabs(lep.eta())> 2.4) continue;
            if(abs(lepid)==11 && fabs(lep.eta())> 2.5) continue;

            bool hasTightIdandIso(true);
            if(abs(lepid)==13) { //muon
                hasTightIdandIso &= phys.leptons[ilep].isMediumMu;
                hasTightIdandIso &= fabs(phys.leptons[ilep].d0)<0.02 && fabs(phys.leptons[ilep].dz)<0.1;
                if(hasTightIdandIso) mon.fillHisto("mu_reliso_raw",   tags, phys.leptons[ilep].m_pfRelIsoDbeta(), weight);

                hasTightIdandIso &= ( phys.leptons[ilep].m_pfRelIsoDbeta() < 0.12 ); // MIT moved to 0.15
            } else if(abs(lepid)==11) { //electron
                hasTightIdandIso &= phys.leptons[ilep].isElpassMedium;
                if(hasTightIdandIso) mon.fillHisto("el_reliso_raw",   tags, phys.leptons[ilep].e_pfRelIsoDbeta(), weight);

                // Electron tight ID from VID already includes isolation
                // hasTightIdandIso &= ( phys.leptons[ilep].e_pfRelIsoDbeta() < 0.1 );
            } else continue;


            if(!hasTightIdandIso) continue;
            nGoodLeptons++;
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(lepid,lep);
            goodLeptons.push_back(goodlep);

        }

        if(nGoodLeptons<2) continue; // 2 tight leptons
        

        float _MASSDIF_(999.);
        int id1(0),id2(0);
        LorentzVector lep1(0,0,0,0),lep2(0,0,0,0);
        for(size_t ilep=0; ilep<goodLeptons.size(); ilep++) {
            int id1_ = goodLeptons[ilep].first;
            LorentzVector lep1_ = goodLeptons[ilep].second;

            for(size_t jlep=ilep+1; jlep<goodLeptons.size(); jlep++) {
                int id2_ = goodLeptons[jlep].first;
                LorentzVector lep2_ = goodLeptons[jlep].second;
                if(id1_*id2_>0) continue; // opposite charge

                LorentzVector dilepton=lep1_+lep2_;
                double massdif = fabs(dilepton.mass()-91.);
                if(massdif < _MASSDIF_) {
                    _MASSDIF_ = massdif;
                    lep1.SetPxPyPzE(lep1_.px(),lep1_.py(),lep1_.pz(),lep1_.energy());
                    lep2.SetPxPyPzE(lep2_.px(),lep2_.py(),lep2_.pz(),lep2_.energy());
                    id1 = id1_;
                    id2 = id2_;
                }
            }
        }


        // ID + ISO scale factors (only muons for the time being) 
        // Need to implement variations for errors (unused for now) 
        float llScaleFactor = 1.0; 
        if ( isMC ) {
            llScaleFactor *= lsf.getLeptonEfficiency( lep1.pt(), lep1.eta(), abs(id1) ).first; 
            llScaleFactor *= lsf.getLeptonEfficiency( lep2.pt(), lep2.eta(), abs(id2) ).first; 
            if(llScaleFactor>0)
                weight *= llScaleFactor; 
            else 
                std::cout << " *** WARNING: ll scale factor = " << llScaleFactor 
                    << ", l1(" << lep1.pt() << ", " << lep1.eta() << ", " << id1 << ")" 
                    << ", l2(" << lep2.pt() << ", " << lep2.eta() << ", " << id2 << ")" 
                    << std::endl; 
        }

        if(id1*id2==0) continue;
        LorentzVector zll(lep1+lep2);
        bool passZmass(fabs(zll.mass()-91)<15);
        bool passZpt(zll.pt()>60);

        double dphiLL{fabs(deltaPhi(lep1.phi(), lep2.phi()))};
        bool passDiLepDphi(dphiLL<M_PI/2);


        TString tag_cat;
        int evcat = getDileptonId(abs(id1),abs(id2));
        switch(evcat) {
        case MUMU :
            tag_cat = "mumu";
            break;
        case EE   :
            tag_cat = "ee";
            break;
        case EMU  :
            tag_cat = "emu";
            break;
        default   :
            continue;
        }


        bool hasTrigger(false);

        if(!isMC) {
            if(evcat!=fType) continue;

            if(evcat==EE   && !(hasEEtrigger||hasEtrigger) ) continue;
            if(evcat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
            if(evcat==EMU  && !hasEMtrigger && !(hasEtrigger && hasMtrigger) ) continue;

            //this is a safety veto for the single mu PD
            if(isSingleMuPD) {
                if(!hasMtrigger) continue;
                if(hasMtrigger && hasMMtrigger) continue;
            }
            if(isDoubleMuPD) {
                if(!hasMMtrigger) continue;
            }

            //this is a safety veto for the single Ele PD
            if(isSingleElePD) {
                if(!hasEtrigger) continue;
                if(hasEtrigger && hasEEtrigger) continue;
            }
            if(isDoubleElePD) {
                if(!hasEEtrigger) continue;
            }

            hasTrigger=true;

        } else {
            if(evcat==EE   && (hasEEtrigger || hasEtrigger) ) hasTrigger=true;
            if(evcat==MUMU && (hasMMtrigger || hasMtrigger) ) hasTrigger=true;
            if(evcat==EMU  && (hasEMtrigger || (hasEtrigger && hasMtrigger)) ) hasTrigger=true;
            if(!hasTrigger) continue;
        }

        tags.push_back(tag_cat); //add ee, mumu, emu category


        // Now we can look at nvtx distribution (not enough QCD MC to allow data comparison before this)
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      genWeight);
        mon.fillHisto("nvtxwgt_raw",   tags, phys.nvtx,      genWeight*puWeight);


        //
        // 3rd LEPTON ANALYSIS
        //

        //loop over all lepton again, check deltaR with dilepton,
        bool pass3dLeptonVeto(true);
        int n3rdLeptons(0);
        vector<LorentzVector> allLeptons;
        allLeptons.push_back(lep1);
        allLeptons.push_back(lep2);
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
            if(lep.pt()<10) continue;
            if(abs(lepid)==13 && fabs(lep.eta())> 2.4) continue;
            if(abs(lepid)==11 && fabs(lep.eta())> 2.5) continue;


            bool isMatched(false);
            isMatched |= (deltaR(lep1,lep) < 0.01);
            isMatched |= (deltaR(lep2,lep) < 0.01);
            if(isMatched) continue;

            bool hasLooseIdandIso(true);
            if(abs(lepid)==13) { //muon
                hasLooseIdandIso &= phys.leptons[ilep].isLooseMu;
                hasLooseIdandIso &= ( phys.leptons[ilep].m_pfRelIsoDbeta() < 0.2 );
            } else if(abs(lepid)==11) { //electron
                hasLooseIdandIso &= phys.leptons[ilep].isElpassVeto;
                // Electron ID from VID already includes isolation
                //hasLooseIdandIso &= ( phys.leptons[ilep].e_pfRelIsoDbeta() < 0.2 );
            } else continue;


            if(!hasLooseIdandIso) continue;

            allLeptons.push_back(lep);
            n3rdLeptons++;
        }

        pass3dLeptonVeto=(n3rdLeptons==0);

        //
        // Also look for soft muons for b veto
        //
        int nSoftMuons(0);
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            if ( lep.pt() < 3 || abs(phys.leptons[ilep].id) != 13 ) continue;
            if ( ! phys.leptons[ilep].isSoftMu ) continue;
            if ( std::find(allLeptons.begin(), allLeptons.end(), lep) != allLeptons.end() ) continue;
            
            nSoftMuons++;
        }

        //
        //JET AND BTAGING ANALYSIS
        //
        PhysicsObjectJetCollection GoodIdJets;
        for(size_t ijet=0; ijet<corrJets.size(); ijet++) {
            if(corrJets[ijet].pt()<30) continue;
            if(fabs(corrJets[ijet].eta())>5.) continue;

            if(!corrJets[ijet].isPFLoose) continue;

            //check overlaps with good leptons
            double minDR(999.);
            for( auto& goodLep : goodLeptons ) {
                double dR = deltaR( corrJets[ijet], goodLep.second );
                if(dR > minDR) continue;
                minDR = dR;
            }
            if(minDR < 0.3) continue;


            GoodIdJets.push_back(corrJets[ijet]);

        }

        bool passBveto(true);
        int nCSVLtags(0),nCSVMtags(0),nCSVTtags(0);
        double BTagScaleFactor(1.0);
        for(size_t ijet=0; ijet<corrJets.size(); ijet++) {
            if(corrJets[ijet].pt()<20) continue;
            if(fabs(corrJets[ijet].eta())>2.4) continue;
            if(!corrJets[ijet].isPFLoose) continue;

            nCSVLtags += (corrJets[ijet].btag0>0.605);
            nCSVMtags += (corrJets[ijet].btag0>0.890);
            nCSVTtags += (corrJets[ijet].btag0>0.970);


            if(!isMC) continue;
            bool isCSVMtagged(corrJets[ijet].btag0>0.890);
            if(abs(corrJets[ijet].flavid)==5) {
                if(isCSVMtagged) BTagScaleFactor *= btag_reader.eval( BTagEntry::FLAV_B, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
            } else if(abs(corrJets[ijet].flavid)==4) {
                if(isCSVMtagged) BTagScaleFactor *= btag_reader.eval( BTagEntry::FLAV_C, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
            } else {
                //no SF available yet for light flavor jet (CSVv2.csv)
                //BTagScaleFactor *= btag_reader.eval( BTagEntry::FLAV_UDSG, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
            }

            // for Btag efficiency
            if(isMC_ttbar||isMC_stop) {
                int flavid = abs(corrJets[ijet].flavid);
                for(size_t csvtag=0; csvtag<CSVkey.size(); csvtag++) {

                    bool isBTag(false);
                    if	   (CSVkey[csvtag]=="CSVL" && (corrJets[ijet].btag0>0.605)) isBTag = true;
                    else if(CSVkey[csvtag]=="CSVM" && (corrJets[ijet].btag0>0.890)) isBTag = true;
                    else if(CSVkey[csvtag]=="CSVT" && (corrJets[ijet].btag0>0.970)) isBTag = true;

                    if(flavid==5) {
                        mon.fillHisto(TString("beff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        if(isBTag) mon.fillHisto(TString("beff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                    } else if(flavid==4) {
                        mon.fillHisto(TString("ceff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        if(isBTag) mon.fillHisto(TString("ceff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                    } else {
                        mon.fillHisto(TString("udsgeff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        if(isBTag) mon.fillHisto(TString("udsgeff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                    }

                }
            }

        }

        passBveto=(nCSVMtags==0)&&(nSoftMuons==0);
        // apply BTag weights
        if(isMC) weight *= BTagScaleFactor;

        LorentzVector metP4_XYCorr = METUtils::applyMETXYCorr(metP4,isMC,phys.nvtx);
        LorentzVector metP4_final = metP4;
        if ( useMETPhiCorr )
            metP4_final = metP4_XYCorr;

        TVector2 met2D(metP4_final.px(), metP4_final.py());
        TVector2 zPt2D(zll.px(), zll.py());
        TVector2 metZBasis = met2D.Rotate(-zPt2D.Phi());
        double metParallelZ = metZBasis.X();
        double metPerpZ = metZBasis.Y();

        double dphiZMET=fabs(deltaPhi(zll.phi(),metP4_final.phi()));
        bool passDphiZMETcut(dphiZMET>2.8);

        //missing ET
        bool passMETcut=(metP4_final.pt()>100);
        bool passMETcut120=(metP4_final.pt()>120);

        //missing ET balance
        double balanceDif = fabs(1-metP4_final.pt()/zll.pt());
        bool passBalanceCut=(balanceDif<0.4);

        //transverse mass
        double MT_massless = METUtils::transverseMass(zll,metP4_final,false);

        //event category
        int eventSubCat  = eventCategoryInst.Get(phys,&GoodIdJets);
        TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);

        tags.push_back(tag_cat+tag_subcat); // add jet binning category
        if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back(tag_cat+"lesq1jets");
        if(tag_cat=="mumu" || tag_cat=="ee") {
            tags.push_back("ll"+tag_subcat);
            if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back("lllesq1jets");
        }


        mon.fillHisto("eventflow",tags,0,weight);
        mon.fillHisto("eventflow_unweighted",tags,0,1.);

        size_t iLeadingJet=0;
        for(size_t ij=0; ij<GoodIdJets.size(); ij++) {
            if ( GoodIdJets[ij].pt() > GoodIdJets[iLeadingJet].pt() )
                iLeadingJet = ij;
            mon.fillHisto("jet_pt_raw",   tags, GoodIdJets[ij].pt(),weight);
            mon.fillHisto("jet_eta_raw",  tags, GoodIdJets[ij].eta(),weight);
        }
        if ( GoodIdJets.size() > 0 ) mon.fillHisto("leadingjet_pt_raw", tags, GoodIdJets[iLeadingJet].pt(), weight);
        mon.fillHisto("njets_raw"                       ,tags, GoodIdJets.size(), weight);

        mon.fillHisto("nGoodleptons_raw",tags, nGoodLeptons, weight);
        mon.fillHisto("nleptons_raw",tags, allLeptons.size(), weight);
        if(lep1.pt()>lep2.pt()) {
            mon.fillHisto("leadlep_pt_raw",   tags, lep1.pt(), weight);
            mon.fillHisto("leadlep_eta_raw",  tags, lep1.eta(), weight);
            mon.fillHisto("trailep_pt_raw",   tags, lep2.pt(), weight);
            mon.fillHisto("trailep_eta_raw",  tags, lep2.eta(), weight);
        } else {
            mon.fillHisto("leadlep_pt_raw",   tags, lep2.pt(), weight);
            mon.fillHisto("leadlep_eta_raw",  tags, lep2.eta(), weight);
            mon.fillHisto("trailep_pt_raw",   tags, lep1.pt(), weight);
            mon.fillHisto("trailep_eta_raw",  tags, lep1.eta(), weight);
        }

        mon.fillHisto("zpt_raw"                         ,tags, zll.pt(),   weight);
        mon.fillHisto("pfmet_raw"                       ,tags, metP4_final.pt(), weight);
        mon.fillHisto("pfmet_wocorr_raw"                ,tags, metP4.pt(), weight);
        mon.fillHisto("pfmet_noHF_raw"                  ,tags, phys.metNoHF.pt(), weight);
        mon.fillHisto("pfmet_puppi_raw"                 ,tags, phys.metPUPPI.pt(), weight);
        mon.fillHisto("zmass_raw"                       ,tags, zll.mass(), weight);
        mon.fillHisto("nbjets_raw"                      ,tags, nCSVMtags, weight);
        mon.fillHisto("nSoftMuons_raw"                      ,tags, nSoftMuons, weight);
        if(lep1.pt()>lep2.pt()) mon.fillHisto("ptlep1vs2_raw"                   ,tags, lep1.pt(), lep2.pt(), weight);
        else                    mon.fillHisto("ptlep1vs2_raw"                   ,tags, lep2.pt(), lep1.pt(), weight);



        //##############################################
        //########  Main Event Selection        ########
        //##############################################

        bool passFinal=false;
        if(passZmass) {
            mon.fillHisto("eventflow",  tags, 1, weight);
            mon.fillHisto("eventflow_unweighted",  tags, 1, 1.);
            mon.fillHisto("zpt_flow", tags, zll.pt(),   weight);

            if(passZpt) {
                mon.fillHisto("eventflow",  tags, 2, weight);
                mon.fillHisto("eventflow_unweighted",  tags, 2, 1.);
                mon.fillHisto("nleptons_flow",tags, allLeptons.size(), weight);

                if(pass3dLeptonVeto) {
                    mon.fillHisto("eventflow",  tags, 3, weight);
                    mon.fillHisto("eventflow_unweighted",  tags, 3, 1.);
                    mon.fillHisto("nbjets_flow",tags, nCSVMtags, weight);
                    mon.fillHisto("nSoftMuons_flow",tags, nSoftMuons, weight);

                    if(passBveto) {
                        mon.fillHisto("eventflow",  tags, 4, weight);
                        mon.fillHisto("eventflow_unweighted",  tags, 4, 1.);

                        //for MET X-Y shift correction
                        mon.fillHisto("pfmetphi_wocorr_presel",tags, metP4.phi(), weight);
                        mon.fillHisto("pfmetphi_wicorr_presel",tags, metP4_XYCorr.phi(), weight);
                        mon.fillHisto("pfmetParallelZ_presel",tags, metParallelZ, weight);
                        mon.fillHisto("pfmetPerpZ_presel",    tags, metPerpZ, weight);
                        mon.fillHisto("pfmetx_vs_nvtx_presel",tags,phys.nvtx,metP4_final.px(), weight);
                        mon.fillHisto("pfmety_vs_nvtx_presel",tags,phys.nvtx,metP4_final.py(), weight);
                        mon.fillHisto("pfmetParallelZ_vs_nvtx_presel",tags,phys.nvtx, metParallelZ, weight);
                        mon.fillHisto("pfmetPerpZ_vs_nvtx_presel",    tags,phys.nvtx, metPerpZ, weight);

                        //preselection plots
                        mon.fillHisto("zmass_presel", tags, zll.mass(), weight);
                        mon.fillHisto("zpt_presel", tags, zll.pt(), weight);
                        mon.fillHisto("pfmet_presel", tags, metP4_final.pt(), weight, true);
                        mon.fillHisto("pfmet2_presel",tags, metP4_final.pt(), weight, true);
                        mon.fillHisto("pfmetCtrl_presel", tags, metP4_final.pt(), weight);
                        mon.fillHisto("pfmetCtrl_wocorr_presel", tags, metP4.pt(), weight);
                        mon.fillHisto("pfmetCtrl_noHF_presel", tags, phys.metNoHF.pt(), weight);
                        mon.fillHisto("pfmetCtrl_puppi_presel", tags, phys.metPUPPI.pt(), weight);
                        mon.fillHisto("mt_presel",   tags, MT_massless, weight);
                        mon.fillHisto("dphiZMET_presel",tags, dphiZMET, weight);
                        mon.fillHisto("dphiLL_presel",tags, dphiLL, weight);
                        mon.fillHisto("balancedif_presel",tags, balanceDif, weight);

                        //check weights
                        mon.fillHisto("weights_presel", tags, 0., 1.);
                        mon.fillHisto("weights_presel", tags, 1., 1.*genWeight);
                        mon.fillHisto("weights_presel", tags, 2., 1.*genWeight*ewk_w);
                        mon.fillHisto("weights_presel", tags, 3., 1.*genWeight*ewk_w*puWeight);
                        mon.fillHisto("weights_presel", tags, 4., 1.*genWeight*ewk_w*puWeight*BTagScaleFactor);
                        mon.fillHisto("weights_presel", tags, 5., 1.*genWeight*ewk_w*puWeight*BTagScaleFactor*llScaleFactor);


                        if(passDphiZMETcut) {
                            mon.fillHisto("eventflow",  tags, 5, weight);
                            mon.fillHisto("eventflow_unweighted",  tags, 5, 1.);

                            if(passDiLepDphi) {
                                mon.fillHisto("eventflow", tags, 6, weight);
                                mon.fillHisto("eventflow_unweighted", tags, 6, 1.);

                                if(passBalanceCut) {
                                    mon.fillHisto("eventflow",  tags, 7, weight);
                                    mon.fillHisto("eventflow_unweighted",  tags, 7, 1.);

                                    if(passMETcut) {
                                        mon.fillHisto("eventflow",  tags, 8, weight);
                                        mon.fillHisto("eventflow_unweighted",  tags, 8, 1.);

                                        if(MT_massless>200) {
                                            passFinal=true;
                                            mon.fillHisto("eventflow",  tags, 9, weight);
                                            mon.fillHisto("eventflow_unweighted",  tags, 9, 1.);

                                            mon.fillHisto("mt_final",   tags, MT_massless, weight);
                                            mon.fillHisto("pfmet_final",tags, metP4_final.pt(), weight);
                                            mon.fillHisto("pfmet2_final",tags, metP4_final.pt(), weight);
                                            if(passMETcut120) mon.fillHisto("mt_final120",   tags, MT_massless, weight);

                                            if(!isMC && outTxtFile_final) fprintf(outTxtFile_final,"%d | %d | %d | pfmet: %f | mt: %f \n",ev.run,ev.lumi,ev.event,metP4_final.pt(), MT_massless);
                                            if(saveEventList) {
                                                eventList_run = ev.run;
                                                eventList_lumi = ev.lumi;
                                                eventList_evt = ev.event;
                                                eventList_nJets = GoodIdJets.size();
                                                eventList_category = evcat; // MUMU=1,EE=2,EMU=3
                                                eventList_met = metP4_final.pt();
                                                eventList_mt = MT_massless;
                                                eventList_llpt = zll.pt();
                                                eventList_mll = zll.mass();
                                                eventList->Fill();
                                            }

                                            mon.fillHisto("weights_final", tags, 0., 1.);
                                            mon.fillHisto("weights_final", tags, 1., 1.*genWeight);
                                            mon.fillHisto("weights_final", tags, 2., 1.*genWeight*ewk_w);
                                            mon.fillHisto("weights_final", tags, 3., 1.*genWeight*ewk_w*puWeight);
                                            mon.fillHisto("weights_final", tags, 4., 1.*genWeight*ewk_w*puWeight*BTagScaleFactor);
                                            mon.fillHisto("weights_final", tags, 5., 1.*genWeight*ewk_w*puWeight*BTagScaleFactor*llScaleFactor);

                                        } //pass MT cut

                                    } //passMETcut

                                } //passBalanceCut

                            } //passDiLepDphi

                        } //passDphiZMETcut

                    } //passBveto

                } //pass3dLeptonVeto

            } //passZpt

        } //passZmass




        //##############################################
        //########  WW/ttbar/etc. Control       ########
        //##############################################

        // WW/ttbar/Wt/tautau control
        mon.fillHisto("zpt_wwctrl_raw"                      ,tags, zll.pt(),   weight);
        mon.fillHisto("zmass_wwctrl_raw"                    ,tags, zll.mass(), weight);
        mon.fillHisto("pfmet_wwctrl_raw"                    ,tags, metP4_final.pt(), weight);
        mon.fillHisto("mt_wwctrl_raw"              	    ,tags, MT_massless, weight);







        //##############################################################################
        //### HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //##############################################################################


        //Fill histogram for posterior optimization, or for control regions
        for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
	  float iweight = weight;                                                //nominal
	  if(varNames[ivar]=="_puup")        iweight *=puWeight_plus/puWeight;   //pu up
	  if(varNames[ivar]=="_pudown")      iweight *=puWeight_minus/puWeight;  //pu down

	  if(isSignal || isMC_ZZ || isMC_WZ) { 
	    ///// 
	    // (varNames[ivar]=="_pdfup" || varNames[ivar]=="_pdfdown") ) {
	    // if(mPDFInfo) {
	    //   float PDFWeight_plus(1.0), PDFWeight_down(1.0);
	    //   std::vector<float> wgts=mPDFInfo->getWeights(iev);
	    //   for(size_t ipw=0; ipw<wgts.size(); ipw++) {
	    // 	PDFWeight_plus = TMath::Max(PDFWeight_plus,wgts[ipw]);
	    // 	PDFWeight_down = TMath::Min(PDFWeight_down,wgts[ipw]);
	    //   }
	    //   if(varNames[ivar]=="_pdfup")    iweight *= PDFWeight_plus;
	    //   if(varNames[ivar]=="_pdfdown")  iweight *= PDFWeight_down;
	    // }
	    // }
	    ///// 

	    //// *** Long if/else if/else *** 
	    /// * 1, 2 * 
	    /// 
	    // // PDF uncertainties' Rosetta stone: 
	    // // NNPDF30_lo_as_0130: 
	    // //   LHAID 263000, pdf uncertainties: 263001-263100 
	    // //                 no alpha_s uncertainties
	    // // NNPDF30_lo_as_0130_nf_4: 
	    // //   LHAID 263400, pdf uncertainties: 263401-263500 
	    // //                 no alpha_s uncertainties
	    // // NNPDF30_nlo_as_0118:  <========= this should be the one for NLO samples (not WZ...) 
	    // //   LHAID 260000, pdf uncertainties: 260001-260100 (index in ntuples: 9-108 [starting from 0]) 
	    // //                 alpha_s uncertainties: 265000,266000 (index in ntuples: 109,110 [starting from 0])
	    // // NPDF30_nlo_as_0118_nf_4: 
	    // //   LHAID 260400, pdf uncertainties: 260401-260500 
	    // //                 alpha_s uncertainties: 265400,266400 
	    // // NNPDF30_nlo_nf_5_pdfas: 
	    // //   LHAID 292200, pdf uncertainties: 292201-292300 
	    // //                 alpha_s uncertainties: 292301,292302 
	    // // NNPDF30_nlo_nf_4_pdfas: 
	    // //   LHAID 292000, pdf uncertainties: 292001-292100 
	    // //                 alpha_s uncertainties: 292101,292102 

	    if(varNames[ivar]=="_pdfup" || varNames[ivar]=="_pdfdown") { 
	      float meanPDFweight(0.); 
	      float rmsPDFweight(0.); 
	      for(size_t i=9; i<=108; ++i) { 
		meanPDFweight += ev.lheWeights[i]/ev.lheOriginalWeight; 
		rmsPDFweight  += ev.lheWeights[i]*ev.lheWeights[i]/pow(ev.lheOriginalWeight,2); 
	      } 
	      meanPDFweight /= 100.; 
	      rmsPDFweight = sqrt( rmsPDFweight/100. - meanPDFweight*meanPDFweight ); // should still be multiplied by 100./99. 

	      float avgASvar = 0.75*(fabs(ev.lheWeights[109]/ev.lheOriginalWeight-1.)+fabs(ev.lheWeights[110]/ev.lheOriginalWeight-1.)); // 0.75 = 1.5/2.
	      float pdfAsVar = sqrt(rmsPDFweight*rmsPDFweight + avgASvar*avgASvar); 
	      pdfAsVar = (pdfAsVar>1. ? 1. : pdfAsVar); 

	      if(varNames[ivar]=="_pdfup") { 
		iweight *= (1. + pdfAsVar); 
		Hcutflow->Fill(5, iweight); 
	      } 
	      else { 
		iweight *= (1. - pdfAsVar); 
		Hcutflow->Fill(6, iweight); 
	      } 
	    } 

	    /// * 3 * 
	    else if(varNames[ivar]=="_qcdup") {  
	      float maxQcdScl(-999.); 
	      if( ev.weight_QCDscale_muR0p5_muF0p5 > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR0p5_muF0p5; 
	      if( ev.weight_QCDscale_muR0p5_muF1   > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR0p5_muF1  ; 
	      //if( ev.weight_QCDscale_muR0p5_muF2   > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR0p5_muF2  ; 
	      if( ev.weight_QCDscale_muR1_muF0p5   > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR1_muF0p5  ; 
	      //if( ev.weight_QCDscale_muR1_muF1     > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR1_muF1    ; // 1 by construction 
	      if( ev.weight_QCDscale_muR1_muF2     > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR1_muF2    ; 
	      //if( ev.weight_QCDscale_muR2_muF0p5   > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR2_muF0p5  ; 
	      if( ev.weight_QCDscale_muR2_muF1     > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR2_muF1    ; 
	      if( ev.weight_QCDscale_muR2_muF2     > maxQcdScl ) maxQcdScl = ev.weight_QCDscale_muR2_muF2    ; 
	      if( maxQcdScl < -998. ) maxQcdScl = 1.000; 
	      iweight *= maxQcdScl; 
	      Hcutflow->Fill(7, iweight); 
	    } 

	    /// * 4 * 
	    else if(varNames[ivar]=="_qcddown") { 
	      float minQcdScl(999.); 
	      if( ev.weight_QCDscale_muR0p5_muF0p5 < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR0p5_muF0p5; 
	      if( ev.weight_QCDscale_muR0p5_muF1   < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR0p5_muF1  ; 
	      //if( ev.weight_QCDscale_muR0p5_muF2   < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR0p5_muF2  ; 
	      if( ev.weight_QCDscale_muR1_muF0p5   < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR1_muF0p5  ; 
	      //if( ev.weight_QCDscale_muR1_muF1     < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR1_muF1    ; // 1 by construction 
	      if( ev.weight_QCDscale_muR1_muF2     < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR1_muF2    ; 
	      //if( ev.weight_QCDscale_muR2_muF0p5   < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR2_muF0p5  ; 
	      if( ev.weight_QCDscale_muR2_muF1     < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR2_muF1    ; 
	      if( ev.weight_QCDscale_muR2_muF2     < minQcdScl ) minQcdScl = ev.weight_QCDscale_muR2_muF2    ; 
	      if( minQcdScl >  998. ) minQcdScl = 1.000; 
	      iweight *= minQcdScl; 
	      Hcutflow->Fill(8, iweight); 
	    }

	    /// * 5 * 
	    else if(isMC_ZZ && varNames[ivar]=="_ewkup"  ) { iweight /= (1. + (1.56 - 1.)*(1. - ewk_w)); } // Ewk up 

	    /// * 6 * 
	    else if(isMC_ZZ && varNames[ivar]=="_ewkdown") { iweight *= (1. + (1.56 - 1.)*(1. - ewk_w)); } // Ewk down

	    // WZ: assign 100% uncertainty (TEMPORARY) 
	    /// * 7 * 
	    //else if(isMC_WZ && varNames[ivar]=="_ewkup"  ) { iweight /= ewk_w; } // Ewk up

	    /// * 8 * 
	    //else if(isMC_WZ && varNames[ivar]=="_ewkdown") { iweight *= ewk_w; } // Ewk down
	    // //if(varNames[ivar]=="_ewkdown") { iweight /= ewk_w; iweight *=(1-(1-EWK_w)*2);} // Ewk down 2012 

	  } // end if(isSignal || isMC_ZZ || isMC_WZ)


            //##############################################
            // recompute MET/MT if JES/JER was varied
            //##############################################
            //LorentzVector vMET = variedMET[ivar>8 ? 0 : ivar];
            //PhysicsObjectJetCollection &vJets = ( ivar<=4 ? variedJets[ivar] : variedJets[0] );


            LorentzVector vMET = variedMET[0];
            if(varNames[ivar]=="_jerup" || varNames[ivar]=="_jerdown" || varNames[ivar]=="_jesup" || varNames[ivar]=="_jesdown" ||
                    varNames[ivar]=="_umetup" || varNames[ivar]=="_umetdown" || varNames[ivar]=="_lesup" || varNames[ivar]=="_lesdown") {
                vMET = variedMET[ivar];
            }
            if ( useMETPhiCorr )
                vMET = METUtils::applyMETXYCorr(vMET,isMC,phys.nvtx);

            PhysicsObjectJetCollection &vJets = variedJets[0];
            if(varNames[ivar]=="_jerup" || varNames[ivar]=="_jerdown" || varNames[ivar]=="_jesup" || varNames[ivar]=="_jesdown") {
                vJets = variedJets[ivar];
            }


            bool passLocalBveto(true);
            for(size_t ijet=0; ijet<vJets.size(); ijet++) {

                if(vJets[ijet].pt()<30) continue;
                if(fabs(vJets[ijet].eta())>5.) continue;

                //jet ID
                if(!vJets[ijet].isPFLoose) continue;
                //if(vJets[ijet].pumva<0.5) continue;

                //check overlaps with selected leptons
                double minDR(999.);
                for(vector<LorentzVector>::iterator lIt = allLeptons.begin(); lIt != allLeptons.end(); lIt++) {
                    double dR = deltaR( vJets[ijet], *lIt );
                    if(dR > minDR) continue;
                    minDR = dR;
                }
                if(minDR < 0.3) continue;
	    }

            for(size_t ijet=0; ijet<vJets.size(); ijet++) {
                if(vJets[ijet].pt()>20 && fabs(vJets[ijet].eta())<2.4) {
                    passLocalBveto &= (vJets[ijet].btag0<0.890);
                    bool isLocalCSVMtagged(vJets[ijet].btag0>0.890);
                    double BTagWeights_Up=1., BTagWeights_Down=1.;
                    if(abs(vJets[ijet].flavid)==5) {
                        if(isLocalCSVMtagged) BTagWeights_Up *= btag_reader_up.eval( BTagEntry::FLAV_B, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        if(isLocalCSVMtagged) BTagWeights_Down *= btag_reader_down.eval( BTagEntry::FLAV_B, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                    } else if(abs(vJets[ijet].flavid)==4) {
                        if(isLocalCSVMtagged) BTagWeights_Up *= btag_reader_up.eval( BTagEntry::FLAV_C, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        if(isLocalCSVMtagged) BTagWeights_Down *= btag_reader_down.eval( BTagEntry::FLAV_C, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                    } else {
                        // No udsg weights available
                    }
                    if(varNames[ivar]=="_btagup") iweight *= BTagWeights_Up;
                    if(varNames[ivar]=="_btagdown")	iweight *= BTagWeights_Down;
                }

            }

            double mt_massless = METUtils::transverseMass(zll,vMET,false); //massless mt
            bool passLocalMt(mt_massless>200.); 
            bool passBaseSelection( passZmass && passZpt && pass3dLeptonVeto && passLocalBveto && (nSoftMuons==0) && passDiLepDphi && passLocalMt );

            double LocalDphiZMET=fabs(deltaPhi(zll.phi(),vMET.phi()));

            //############
            //optimization
            //############
            for(unsigned int index=0; index<nOptims; index++) {


                double minMET = optim_Cuts1_MET[index];
                double minBalance = optim_Cuts1_Balance[index];
                double minDphi = optim_Cuts1_DphiZMET[index];

                bool passLocalMETcut(vMET.pt()>minMET);
                bool passLocalBalanceCut=(vMET.pt()/zll.pt()>(1.-minBalance) && vMET.pt()/zll.pt()<(1.+minBalance));
                bool passLocalDphiZMETcut(LocalDphiZMET>minDphi);

                bool passOptimSelection(passBaseSelection && passLocalMETcut && passLocalBalanceCut && passLocalDphiZMETcut);



                // fill shapes for limit setting
                if( passOptimSelection ) {
                    mon.fillHisto(TString("mt_shapes")+varNames[ivar],tags,index, mt_massless, iweight);
                    if ( index==17 && passFinal==false && ivar==0 ) {
                        cerr << "Found one!" << endl
                            << " MET, Balance, dphi values = " << minMET << ", " << minBalance << ", " << minDphi << endl
                            << "_run = " <<  ev.run << endl
                            << "_lumi = " <<  ev.lumi << endl
                            << "_evt = " <<  ev.event << endl
                            << "_nJets = " <<  GoodIdJets.size() << endl
                            << "_category = " <<  evcat << endl
                            << "_met = " <<  metP4_final.pt() << endl
                            << "_mt = " <<  MT_massless << endl
                            << "_llpt = " <<  zll.pt() << endl
                            << "_mll = " <<  zll.mass() << endl;
                   }
                }


            }//all optimization END






        }//Systematic variation END

















    } // loop on all events END




    printf("\n");
    file->Close();

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################
    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results saved in %s\n", outUrl.Data());

    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();
    if ( eventList ) eventList->Write();

    ofile->Close();

    PU_Central_File->Close();
    PU_Up_File->Close();
    PU_Down_File->Close();

    if(outTxtFile_final)fclose(outTxtFile_final);
}
/* vim: set ts=4 sw=4 tw=0 et :*/
