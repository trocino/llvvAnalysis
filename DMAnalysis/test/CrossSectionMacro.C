#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>

#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TRegexp.h"
#include "TLegend.h"
#include "TLatex.h"

void extendXlimit(TH1D*&, double); 

void changeFirstPoint(TH1D*&, double, double); 

void CrossSectionMacro(bool plot0jets=true, bool isnormalized=false, bool islogy=false) {
  TString indir = "analysis_approval_20160302_1536/limits_shape/cardsShape/125/"; 
  ifstream infile;
  bool saveplots = true; 

  bool stacksignal = false; 
  bool stackatgc = false; 
  bool addcorrelation = true; 
  bool normalizebins = isnormalized; 
  bool dology = islogy; 
  bool only0jets = plot0jets; 
  bool only1jets = !only0jets; 
  if(only0jets && only1jets) {
    std::cout << " *** ERROR: only0jets and only1jets cannot be both true at the same time!" << std::endl; 
    exit(1); 
  } 

  if     (only0jets) 
    infile.open( (indir+"card_0jets.dat").Data() );
  else if(only1jets)
    infile.open( (indir+"card_1jets.dat").Data() );
  else 
    infile.open( (indir+"card_combined.dat").Data() );

  // Exit if file opening failed
  if(!infile.is_open()){
    std::cerr << "Opening datacard failed!" << std::endl;
    exit(1);
  }

  unsigned int nseps = 0;

  unsigned int nbins = 0;
  unsigned int cbins = 0;
  unsigned int nproc = 0;
  //unsigned int nsyst = 0;

  TString rootfile; 
  std::vector<TString> finstates; 
  std::vector<TString> processes; 
  std::vector<int>     processIds; 
  int nprocidx = 0; 
  std::map<TString, bool> systtype; 
  std::vector<TString> systnames; 
  std::map<TString, std::vector<double> > systs; 

  std::map<TString, TString> procnames; 

  double zzXSerrs = 0.0; 

  // procnames["zz2l2nu"] = "ZZ #rightarrow 2l2#nu";
  // procnames["topwwwjetsdata"] = "top, WW, W+jets (data)";
  // procnames["zlldata"] = "Z #rightarrow 2l (data)";
  // procnames["wz3lnu"] = "WZ #rightarrow 3l#nu";
  // procnames["zz2l2nuf4z0002"] = "ZZ #rightarrow 2l2#nu f_{4}^{Z}=0.002"; 
  // procnames["zz2l2nuf4z0005"] = "ZZ #rightarrow 2l2#nu f_{4}^{Z}=0.005"; 
  // procnames["zz2l2nuf4z001"]  = "ZZ #rightarrow 2l2#nu f_{4}^{Z}=0.01"; 
  // procnames["zz2l2nuf4z002"]  = "ZZ #rightarrow 2l2#nu f_{4}^{Z}=0.02"; 

  // procnames["ZZ"] = "ZZ #rightarrow 2l2#nu";
  // procnames["EM"] = "top + WW";
  // procnames["zjets"] = "Z #rightarrow 2l";
  // procnames["WZ"] = "WZ #rightarrow 3l#nu";
  // procnames["vvv"] = "VVV"; 
  // procnames["ZH_hinv"] = "ZH(125)";

  procnames["ZZ"] = "ZZ";
  procnames["EM"] = "top + WW";
  procnames["zjets"] = "Z/#gamma* + jets";
  procnames["WZ"] = "WZ";
  procnames["vvv"] = "VVV"; 
  procnames["ZH_hinv"] = "ZH(125)";


  TString titlex;
  TString titley;

  string line;

  while(infile.good()) {
    getline(infile, line);
    //std::cout << line << std::endl;
    if(line.size()==0) continue; 

    std::vector<TString> lstr; 
    string temp;
    stringstream s(line);
    while(s >> temp) lstr.push_back(temp.c_str()); 
    if(lstr.size()==0) continue; 
    //std::cout << lstr.size() << std::endl;

    // Read file
    if(lstr[0].Contains("--------")) {++nseps; continue;} 
    if(lstr[0].CompareTo("imax")==0) {nbins = cbins = lstr[1].Atoi(); continue;} 
    if(lstr[0].CompareTo("jmax")==0) {nproc = lstr[1].Atoi(); ++nproc; continue;} 
    //if(lstr[0].CompareTo("kmax")==0) {nsyst = lstr[1].Atoi(); continue;} 
    //std::cout << nbins << " " << nproc << " " << nsyst << std::endl;
    if(lstr[0].CompareTo("shapes")==0 && cbins>0) {
      finstates.push_back(lstr[2]); 
      --cbins; 
      if(rootfile.Length()==0 && lstr[3].EndsWith(".root")) rootfile = lstr[3]; 
      continue; 
    }
    //std::cout << finstates.size() << std::endl;
    if(lstr[0].CompareTo("process")==0 && nprocidx==0) {
      ++nprocidx; 
      for(unsigned int i=0; i<nproc; ++i) 
	processes.push_back(lstr[1+i]); 
      continue; 
    }
    if(lstr[0].CompareTo("process")==0 && nprocidx==1) {
      ++nprocidx; 
      for(unsigned int i=0; i<nproc; ++i) 
	processIds.push_back(lstr[1+i].Atoi()); 
      continue; 
    }
    //std::cout << processes.size() << std::endl;
    if(nseps==4) {
      systnames.push_back(lstr[0]); 
      systtype[lstr[0]] = lstr[1].Contains("shapeN2"); 
      //std::cout << lstr[0] << " " << systtype[lstr[0]] << std::endl; 
      systs[lstr[0]] = std::vector<double>(0., nbins*nproc); 
      for(unsigned int i=0; i<nbins; ++i) 
	for(unsigned int j=0; j<nproc; ++j) {
	  systs[lstr[0]].push_back( lstr[2+i*nproc+j].CompareTo("-")==0 ? 0. : ( systtype[lstr[0]] ? 1. : fabs( lstr[2+i*nproc+j].Atof() - 1.0 ) ) ); 
	}
      continue; 
    }
  }

  infile.close();

  // if ATGC, sort ATGC samples...
  map<int, int> remappedAtgcs; 
  if(!stacksignal) {
    vector<double> sortAtgcsl; 
    for(unsigned int isg=0; isg<processIds.size(); ++isg) {
      if(processIds[isg]<=0) {
	TString tmpATGCstr = processes[isg]( processes[isg].Index("0")+1, processes[isg].Length() ); 
	tmpATGCstr.Prepend("0."); 
	sortAtgcsl.push_back(tmpATGCstr.Atof()); 
	//std::cout << tmpATGCstr.Data() << "  " << sortAtgcsl.back() << std::endl;
	remappedAtgcs[isg] = isg; 
      }
    }

    // Re-order
    for(unsigned int isg1=0; isg1<sortAtgcsl.size(); ++isg1) {
      for(unsigned int isg2=isg1+1; isg2<sortAtgcsl.size(); ++isg2) {
	if(sortAtgcsl[isg1]<sortAtgcsl[isg2]) {
	  double tmpsort = sortAtgcsl[isg1]; 
	  sortAtgcsl[isg1] = sortAtgcsl[isg2]; 
	  sortAtgcsl[isg2] = tmpsort; 

	  int tmpidx = remappedAtgcs[isg1]; 
	  remappedAtgcs[isg1] = remappedAtgcs[isg2]; 
	  remappedAtgcs[isg2] = tmpidx; 
	}
      }
      //std::cout << isg1 << "  " << remappedAtgcs[isg1] << std::endl; 
    }
  }

  bool is2011 (false), is2012 (false), is2015(false); 
  if     (rootfile.Contains("7TeV")) is2011 = true; 
  else if(rootfile.Contains("8TeV")) is2012 = true; 
  else                               is2015 = true; 

  TFile *inrootfile = 0; 
  if(rootfile.EndsWith(".root")) inrootfile = TFile::Open( (indir+rootfile).Data()); 
  if(inrootfile==0 || inrootfile->IsZombie()) {
    std::cerr << "Opening ROOT file failed!" << std::endl;
    exit(1);
  }

  std::map<TString, TH1D*> hists; 
  std::map<TString, TH1D*> hists4xs; 
  TH1D *hsum = 0; 
  TH1D *hdata = 0; 
  TGraphAsymmErrors *gaedata = 0; 
  //TH1D *hdata4xs = 0; 
  THStack *hstk = new THStack("hstk", "stacked histos"); 

  double extendedLimit = 1000.; 
//   double oldFirstPoint = 40.; 
//   double newFirstPoint = 45.; 
//   double extendedLimit = 0.; 
  double oldFirstPoint = 0.; 
  double newFirstPoint = 0.; 

  for(unsigned int i=0; i<nproc; ++i) {
    for(unsigned int j=0; j<nbins; ++j) { 
      // if(only0jets && finstates[j].Contains("eq1jets")) continue; 
      // if(only1jets && finstates[j].Contains("eq1jets")) continue; 
      TH1D *htmp = (TH1D*)inrootfile->Get( (finstates[j]+"/"+processes[i]).Data() ); 
      if(htmp==0) {
	std::cerr << "Retrieving histo " << (finstates[j]+"/"+processes[i]).Data() << " failed!" << std::endl;
	exit(1);
      }
      //std::cout << htmp << "  "; 
      if(extendedLimit>0.) extendXlimit(htmp, extendedLimit); 
      if(oldFirstPoint>0. && newFirstPoint>0.) changeFirstPoint(htmp, oldFirstPoint, newFirstPoint); 
      //std::cout << htmp << std::endl; 
      // Separate plots, sum on final states
      if(hists.count(processes[i])==0) {
 	hists[processes[i]] = (TH1D*)htmp->Clone( (processes[i]+"_ll").Data() ); 
 	hists4xs[processes[i]] = (TH1D*)htmp->Clone( (processes[i]+"_ll").Data() ); 
	// if(processes[i].CompareTo("zz2l2nu")==0) { 
	//   hists[processes[i]]->SetFillColor(19); // Fill signal with very pale gray 
	//   //hists[processes[i]]->SetLineWidth(2); // Fill signal with very pale gray 
	//   hists4xs[processes[i]]->SetFillColor(19); // Fill signal with very pale gray 
	//   //hists4xs[processes[i]]->SetLineWidth(2); // Fill signal with very pale gray 
	// } 
	titlex = htmp->GetXaxis()->GetTitle(); 
	titley = htmp->GetYaxis()->GetTitle(); 
      }
      else {
	hists[processes[i]]->Add(htmp); 
	hists4xs[processes[i]]->Add(htmp); 
      }

      if(processes[i].CompareTo("EM")==0) { 
	hists[processes[i]]->SetFillColor(kOrange-2); 
      } 

      // Total background
      if(hsum==0) {
	if(processIds[i]>0 || stacksignal) hsum = (TH1D*)htmp->Clone("sum"); 
      }
      else {
	if(processIds[i]>0 || stacksignal) hsum->Add(htmp); 
      }      
    }

    // Normalize to bin width?
    if(normalizebins) 
      for(int q=1; q<=hists[processes[i]]->GetNbinsX(); ++q) {
	double yval = hists[processes[i]]->GetBinContent(q); 
	double xval = hists[processes[i]]->GetBinWidth(q); 
	hists[processes[i]]->SetBinContent(q, yval/xval); 
      }

    // Feel stuck? Fill stack! 
    if(processIds[i]>0) hstk->Add(hists[processes[i]]); 
  }

  // Stack signal?
  //std::cout << hists[processes[0]] << " " << processes[0] << " " << hists[processes[0]]->GetName() << std::endl;
  if(stacksignal) hstk->Add(hists[processes[0]]); // only if there's one single signal sample, thus "0" index

  // Remove errors from MC 
  for(unsigned int k=0; k<processes.size(); ++k) {
    for(unsigned int h=1; h<=(unsigned int)hists4xs[processes[k]]->GetNbinsX(); ++h) {
      hists4xs[processes[k]]->SetBinError(h, 0.0); 
    }
  }

  // Now data
  for(unsigned int j=0; j<nbins; ++j) {
    // if(only0jets && finstates[j].Contains("eq1jets")) continue; 
    // if(only1jets && finstates[j].Contains("eq1jets")) continue; 
    TH1D *htmp = (TH1D*)inrootfile->Get( (finstates[j]+"/data_obs").Data() ); 
    if(htmp==0) {
      std::cerr << "Retrieving histo " << (finstates[j]+"/data_obs").Data() << " failed!" << std::endl;
      exit(1);
    }
    if(extendedLimit>0.) extendXlimit(htmp, extendedLimit); 
    if(oldFirstPoint>0. && newFirstPoint>0.) changeFirstPoint(htmp, oldFirstPoint, newFirstPoint); 
    if(hdata==0) {
      const unsigned int nbinsdata = htmp->GetNbinsX(); 
      double hdataarr[nbinsdata+1]; 
      for(unsigned int iii=0; iii<=nbinsdata; ++iii) hdataarr[iii] = htmp->GetBinLowEdge(iii+1); 
      hdata = new TH1D("hdata_new", "hdata_new", nbinsdata, hdataarr); 
      //hdata = (TH1D*)htmp->Clone(); 
      for(unsigned int iii=1; iii<=nbinsdata; ++iii) hdata->SetBinContent(iii, htmp->GetBinContent(iii)); 
      hdata->SetMarkerStyle(20); 
      //hdata4xs = (TH1D*)htmp->Clone(); 
    }
    else {
      //hdata->Add(htmp); 
      const unsigned int nbinsdata = htmp->GetNbinsX(); 
      for(unsigned int iii=1; iii<=nbinsdata; ++iii) hdata->SetBinContent(iii, hdata->GetBinContent(iii) +  htmp->GetBinContent(iii)); 
      //hdata4xs->Add(htmp); 
    } 
  }

  //std::cout << " ******** " << hdata->GetBinErrorOption() << std::endl; 
  hdata->SetBinErrorOption(TH1::kPoisson); 
  //std::cout << " ******** " << hdata->GetBinErrorOption() << std::endl; 
  gaedata = new TGraphAsymmErrors(hdata); 
  gaedata->SetMarkerStyle(20); 
  for(int q=0; q<gaedata->GetN(); ++q) {
    gaedata->SetPointEYlow(q, hdata->GetBinErrorLow(q+1)); 
    gaedata->SetPointEYhigh(q, hdata->GetBinErrorUp(q+1)); 
    // if(hdata->GetBinContent(q+1)<1.0E-6) {
    //   gaedata->SetPoint(q, hdata->GetBinCenter(q+1), -1.); 
    //   gaedata->SetPointEYlow(q, 0.); 
    //   gaedata->SetPointEYhigh(q, 0.); 
    // }
  }

  // Normalize to bin width?
  if(normalizebins) 
    for(int q=0; q<gaedata->GetN(); ++q) {
      double xval, yval; 
      gaedata->GetPoint(q, xval, yval); 
      double yerrup = gaedata->GetErrorYhigh(q); 
      double yerrdn = gaedata->GetErrorYlow(q); 
      double xsize  = gaedata->GetErrorXlow(q) + gaedata->GetErrorXhigh(q); 
      // if(yval>0.) {
      gaedata->SetPoint(q, xval, yval/xsize); 
      gaedata->SetPointEYhigh(q, yerrup/xsize); 
      gaedata->SetPointEYlow(q, yerrdn/xsize); 
      // }
    }

  // Errors
  hsum->SetFillStyle(3013); //3001); 
  hsum->SetFillColor(1); //11); 
  hsum->SetLineStyle(0); 
  hsum->SetLineColor(0); 
  hsum->SetLineWidth(0); 
  // Reset bin errors
  for(unsigned int h=1; h<=(unsigned int)hsum->GetNbinsX(); ++h) hsum->SetBinError(h, 0.); 

  for(unsigned int k=0; k<systnames.size(); ++k) {
    TH1D *hsystmp = 0; 
    std::vector<double> & vsystmp = systs[systnames[k]]; 
    if(vsystmp.size()!=nbins*nproc) {
      std::cerr << "Problem with vector size (" << vsystmp.size() << ") vs nbins*nproc (" << nbins*nproc << ")!"<< std::endl;
      exit(1); 
    }

    //std::cout << systnames[k].Data(); 

   // All bins/processes for this systematic
    for(unsigned int i=0; i<nbins*nproc; ++i) {
      double isystval = vsystmp[i]; 
      if(fabs(isystval)<0.0001) continue; 
      unsigned int ibin = i / nproc; 
      unsigned int ipro = i % nproc; 
      //std::cout << systnames[k].Data() << " " << ibin << " " << ipro << " ";
      //std::cout << processIds[ipro] << std::endl;
      if(processIds[ipro]<=0 && stacksignal==false) continue; // Don't sum ATGC signal errors
      //std::cout << std::endl;
      TH1D *hisyst = 0; 
      TString systploti = finstates[ibin]+"/"+processes[ipro]; 
      hisyst = (TH1D*)((TH1D*)inrootfile->Get( systploti.Data() ))->Clone(); 
      if(systtype[systnames[k]]) {
	systploti += "_"+systnames[k]; 
	TH1D *hisystup = (TH1D*)((TH1D*)inrootfile->Get( (systploti+"Up").Data() ))->Clone(); 
	if(hisystup==0) {
	  std::cerr << "Plot " << (systploti+"Up").Data() << " not found!" << std::endl;
	  exit(1); 
	}
	TH1D *hisystdn = (TH1D*)((TH1D*)inrootfile->Get( (systploti+"Down").Data() ))->Clone(); 
	if(hisystdn==0) {
	  std::cerr << "Plot " << (systploti+"Down").Data() << " not found!" << std::endl;
	  exit(1); 
	}
	for(unsigned int h=1; h<=(unsigned int)hisyst->GetNbinsX(); ++h) {
	  double binref = hisyst->GetBinContent(h); 
	  double binup  = fabs(binref-hisystup->GetBinContent(h)); 
	  double bindn  = fabs(binref-hisystdn->GetBinContent(h)); 
	  //if(binup>bindn) hisyst->SetBinContent(h, binup); 
	  //else            hisyst->SetBinContent(h, bindn); 
	  hisyst->SetBinContent(h, (binup+bindn)/2.0); 
	}
	//hisyst->Add( (TH1D*)inrootfile->Get( (systploti+"Down").Data() ), -1); 
	//hisyst->Scale(0.5); 
      }
      else {
	//hisyst = (TH1D*)((TH1D*)inrootfile->Get( systploti.Data() ))->Clone(); 
	if(hisyst==0) {
	  std::cerr << "Plot " << systploti.Data() << " not found!" << std::endl;
	  exit(1); 
	}
	hisyst->Scale( /*fabs(*/ isystval /*-1.0)*/ ); 
      }
      //std::cout << " - " << finstates[ibin] << "/" << processes[ipro] << " (" 
      //	<< fabs(hisyst->Integral()) << ")"; 

      if(hsystmp==0) {
	hsystmp = (TH1D*)hisyst->Clone(); 
	hsystmp->Reset("ICES"); 
      }
      for(unsigned int h=1; h<=(unsigned int)hsystmp->GetNbinsX(); ++h) {
	double val1 = hsystmp->GetBinContent(h); 
	double val2 = hisyst->GetBinContent(h); 
	hsystmp->SetBinContent(h, sqrt(val1*val1 + val2*val2)); 
	double err1 = hists4xs[processes[ipro]]->GetBinError(h); 
	if( processes[ipro].CompareTo("zz2l2nu")==0 && 
	    (systnames[k].Contains("QCDscale") || systnames[k].Contains("pdf") ||
	     systnames[k].Contains("lumi")     || systnames[k].Contains("acc"))  ) {
	  zzXSerrs *= zzXSerrs; 
	  zzXSerrs += val2*val2; 
	  zzXSerrs = sqrt(zzXSerrs);
	  val2 = 0.0; 
	}
	hists4xs[processes[ipro]]->SetBinError(h, sqrt(err1*err1 + val2*val2)); 
      }

      // Add covariance if appropriate (take 100% correlation)
      if(addcorrelation) {
	for(unsigned int j=i+1; j<nbins*nproc; ++j) {
	  double jsystval = vsystmp[j]; 
	  if(fabs(jsystval)<0.0001) continue; 
	  unsigned int jbin = j / nproc; 
	  unsigned int jpro = j % nproc; 
	  //std::cout << systnames[k].Data() << " " << jbin << " " << jpro << " ";
	  //std::cout << processIds[jpro] << std::endl;
	  if(processIds[jpro]<=0 && stacksignal==false) continue; // Don't sum ATGC signal errors
	  //std::cout << std::endl;
	  TH1D *hjsyst = 0; 
	  TString systplotj = finstates[jbin]+"/"+processes[jpro]; 
	  hjsyst = (TH1D*)((TH1D*)inrootfile->Get( systplotj.Data() ))->Clone(); 
	  if(systtype[systnames[k]]) {
	    systplotj += "_"+systnames[k]; 
	    TH1D *hjsystup = (TH1D*)((TH1D*)inrootfile->Get( (systplotj+"Up").Data() ))->Clone(); 
	    if(hjsyst==0) {
	      std::cerr << "Plot " << (systplotj+"Up").Data() << " not found!" << std::endl;
	      exit(1); 
	    }
	    TH1D *hjsystdn = (TH1D*)((TH1D*)inrootfile->Get( (systplotj+"Down").Data() ))->Clone(); 
	    if(hjsystdn==0) {
	      std::cerr << "Plot " << (systplotj+"Down").Data() << " not found!" << std::endl;
	      exit(1); 
	    }
	    for(unsigned int h=1; h<=(unsigned int)hjsyst->GetNbinsX(); ++h) {
	      double binref = hjsyst->GetBinContent(h); 
	      double binup  = fabs(binref-hjsystup->GetBinContent(h)); 
	      double bindn  = fabs(binref-hjsystdn->GetBinContent(h)); 
	      //if(binup>bindn) hjsyst->SetBinContent(h, binup); 
	      //else            hjsyst->SetBinContent(h, bindn); 
	      hjsyst->SetBinContent(h, (binup+bindn)/2.0); 
	    }
	    //hjsyst->Add( (TH1D*)inrootfile->Get( (systplotj+"Down").Data() ), -1); 
	    //hjsyst->Scale(0.5); 
	  }
	  else {
	    //hjsyst = (TH1D*)((TH1D*)inrootfile->Get( systplotj.Data() ))->Clone(); 
	    if(hjsyst==0) {
	      std::cerr << "Plot " << systplotj.Data() << " not found!" << std::endl;
	      exit(1); 
	    }
	    hjsyst->Scale( fabs(jsystval-1.0) ); 
	  }
	  //std::cout << " - " << finstates[ibin] << "/" << processes[ipro] 
	  //    <<  "-"  << finstates[jbin] << "/" << processes[jpro] << " (" 
	  //    << sqrt(2*fabs(hisyst->Integral()*hjsyst->Integral())) << ")"; 

	  for(unsigned int h=1; h<=(unsigned int)hsystmp->GetNbinsX(); ++h) {
	    double val1 = hsystmp->GetBinContent(h); 
	    double val2 = hisyst->GetBinContent(h); 
	    double val3 = hjsyst->GetBinContent(h); 
	    hsystmp->SetBinContent(h, sqrt(val1*val1 + 2*fabs(val2*val3))); 
	    // Keep signal ZZ->2l2nu separate from other processes (don't count correlations between 
	    // ZZ and other processes, but count correlations between other backgrounds)
	    if( (processes[ipro].CompareTo("zz2l2nu")==0 && processes[jpro].CompareTo("zz2l2nu")==0) ||
		(processes[ipro].CompareTo("zz2l2nu")!=0 && processes[jpro].CompareTo("zz2l2nu")!=0) ) {
	      double err1 = hists4xs[processes[ipro]]->GetBinError(h); 
	      if( processes[ipro].CompareTo("zz2l2nu")==0 && 
		  (systnames[k].Contains("QCDscale") || systnames[k].Contains("pdf") ||
		   systnames[k].Contains("lumi")     || systnames[k].Contains("acc"))  ) {
		zzXSerrs *= zzXSerrs; 
		zzXSerrs += 2*fabs(val2*val3); 
		zzXSerrs = sqrt(zzXSerrs);
		val2 = 0.0; 
		val3 = 0.0; 
	      }
	      hists4xs[processes[ipro]]->SetBinError(h, sqrt(err1*err1 + 2*fabs(val2*val3))); 
	    }
	  }
	}
      }
    }
    //std::cout << std::endl;

    // Add contribution from this systematic source to error of sum 
    for(unsigned int h=1; h<=(unsigned int)hsum->GetNbinsX(); ++h) {
      double val1 = hsum->GetBinError(h); 
      double val2 = hsystmp==0 ? 0 : hsystmp->GetBinContent(h); 
      //std::cout << val1 << "  " << val2 << std::endl;
      hsum->SetBinError(h, sqrt(val1*val1 + val2*val2)); 
    }
  }

  double ymax = 0;

  // Normalize to bin width?
  if(normalizebins) 
    for(int q=1; q<=hsum->GetNbinsX(); ++q) {
      double yval = hsum->GetBinContent(q); 
      double yerr = hsum->GetBinError(q); 
      double xval = hsum->GetBinWidth(q); 
      hsum->SetBinContent(q, yval/xval); 
      hsum->SetBinError(q, yerr/xval); 
    }

  TGraphErrors *gsum=new TGraphErrors(hsum);
  for(unsigned int h=0; h<(unsigned int)hsum->GetNbinsX(); ++h) {
    double xerr = hsum->GetBinWidth(h+1)/2; 
    double yerr = hsum->GetBinError(h+1); 
    gsum->SetPointError(h, xerr, yerr); 

    //std::cout << ymax << " - "; 
    ymax = hsum->GetBinContent(h+1)+yerr>ymax ? hsum->GetBinContent(h+1)+yerr : ymax; 
    //std::cout << ymax << " - "; 
    // ymax = hdata->GetBinContent(h+1)+hdata->GetBinError(h+1)>ymax ? 
    //   hsum->GetBinContent(h+1)+hdata->GetBinError(h+1)+hdata->GetBinError(h+1) : ymax; 
    double xtmp, ytmp; 
    gaedata->GetPoint(h, xtmp, ytmp); 
    ymax = ytmp+gaedata->GetErrorYhigh(h)>ymax ? ytmp+gaedata->GetErrorYhigh(h) : ymax; 
    //std::cout << ymax << std::endl; 
  }

  // if(dology) ymax *= 75.; 
  // else       ymax *= 1.05; 

  gsum->SetLineColor(0);
  gsum->SetLineWidth(0);
  gsum->SetFillStyle(3013); //3001);
  gsum->SetFillColor(1); //11);
  gsum->SetMarkerColor(0);
  gsum->SetMarkerStyle(0);
  gsum->SetMarkerSize(0);

  titlex.ReplaceAll("Z", "Dilepton"); 
  titlex.ReplaceAll("red-MET", "Reduced MET"); 
  titlex.ReplaceAll("/c", ""); 

  TRegexp gevs = "[ (]*/[ 1-9]*GeV[ )]*";
  if(normalizebins) {
    if(titley.Contains(gevs)) {
      TString seg = titley(gevs); 
      titley.ReplaceAll(seg.Data(), " / GeV"); 
      titley.Prepend("< "); 
      titley.Append(" >"); 
    }
  }
  else {
    if(titley.Contains(gevs)) {
      TString seg = titley(gevs); 
      titley.ReplaceAll(seg.Data(), ""); 
      titley.ReplaceAll(seg.Data(), " / bin"); 
    }
  }
  titley.ReplaceAll("Entries", "Events"); 

  TCanvas *cc = new TCanvas("cc", "cc", 600, 600); 
  cc->Draw(); 

  double padfrac = 0.25; 

  TPad* t1 = new TPad("t1","t1", 0.0, padfrac, 1.0, 1.0);
  t1->SetTopMargin(0.085); 
  t1->SetBottomMargin(0.015); 
  t1->SetLeftMargin(0.14); 
  t1->SetRightMargin(0.06); 
  //t1->SetGridx(); 
  //t1->SetGridy(); 
  t1->Draw();
  t1->SetLogx(); 
  t1->SetLogy(dology); 

  TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, padfrac);
  t2->SetPad(0.0, 0.0, 1.0, padfrac);
  t2->SetTopMargin(0.05);
  t2->SetBottomMargin(0.42);
  t2->SetLeftMargin(0.14); 
  t2->SetRightMargin(0.06); 
  t2->SetGridx();
  t2->SetGridy();
  t2->Draw();

  t1->cd(); 
  hstk->Draw("HIST"); 
  hstk->GetXaxis()->SetRangeUser(200., 1000.); 
  if(dology) {
    if(normalizebins) {
      if(is2011) {
	hstk->SetMinimum(0.00004); 
	ymax *= 100.; 
      }
      else if(is2012) {
	hstk->SetMinimum(0.001); 
	ymax *= 100.; 
      }
      else {
	hstk->SetMinimum(0.00004); 
	ymax *= 20.; 
      }
    }
    else {
      if(is2011) {
	hstk->SetMinimum(0.02); 
	ymax *= 50.; 
      }
      else if(is2012) {
	hstk->SetMinimum(0.5); 
	ymax *= 100.; 
      }
      else {
	hstk->SetMinimum(0.02); 
	ymax *= 10.; 
      }
    }
  }
  else {
    ymax *= 1.05; 
  }

  //hstk->GetXaxis()->SetMoreLogLabels(); 
  //hstk->GetXaxis()->SetNoExponent(); 
  //hstk->GetXaxis()->SetTitle(titlex.Data()); 
  //hstk->GetXaxis()->SetLabelOffset(1.50); 
  hstk->GetYaxis()->SetTitleSize(0.05 * 1.0/(1.0-padfrac)); 
  hstk->GetYaxis()->SetTitleOffset(1.20 * (1.0-padfrac)); 
  hstk->GetYaxis()->SetLabelSize(0.035 * 1.0/(1.0-padfrac)); 
  hstk->GetYaxis()->SetTitle(titley.Data()); 
  //std::cout << ymax << std::endl; 
  hstk->SetMaximum(ymax); 
  if(!stacksignal) {
    for(unsigned int isg=0; isg<processIds.size(); ++isg) {
      if(processIds[isg]<=0) {
	if(stackatgc) 
	  hists[processes[remappedAtgcs[isg]]]->Add(hsum); 
	hists[processes[remappedAtgcs[isg]]]->SetLineStyle(1+isg); 
	hists[processes[remappedAtgcs[isg]]]->Draw("HISTSAME"); 
      }
      else 
	break;
    }
  }
  gsum->Draw("2");
  //hdata->Draw("E0SAME"); 
  gaedata->Draw("PSAME"); 

  // Ratio
  //TH1D *hratio = (TH1D*)hdata->Clone(); 
  TGraphAsymmErrors *gaeratio = (TGraphAsymmErrors*)gaedata->Clone(); 
  TGraphErrors *gratio = (TGraphErrors*)gsum->Clone(); 

  for(int q=0; q<gratio->GetN(); ++q) {
    double xval1, yval1, 
      xerrhi1, xerrlo1, 
      yerrhi1, yerrlo1; 
    gaeratio->GetPoint(q, xval1, yval1); 
    xerrhi1 = gaeratio->GetErrorXhigh(q); 
    xerrlo1 = gaeratio->GetErrorXlow(q); 
    yerrhi1 = gaeratio->GetErrorYhigh(q); 
    yerrlo1 = gaeratio->GetErrorYlow(q); 

    double xval2, yval2; 
    gratio->GetPoint(q, xval2, yval2); 
    double xerr2 = gratio->GetErrorX(q); 
    double yerr2 = gratio->GetErrorY(q); 

    //std::cout << yval1/yval2 << " +- " << yerr1/yval2 << std::endl;

    if(yval2>0.0) {
      gaeratio->SetPoint(q, xval1, yval1/yval2); 
      gaeratio->SetPointError(q, xerrlo1, xerrhi1, yerrlo1/yval2, yerrhi1/yval2); 
      //if(yval1<10e-4) {hratio->SetBinContent(q+1, -1.); hratio->SetBinError(q+1, 0.);} 
      // if(yval1<1.0E-6) {
      // 	gaeratio->SetPoint(q, xval1, -1.); 
      // 	gaeratio->SetPointError(q, 0., 0., 0., 0.); 
      // } 
      gratio->SetPoint(q, xval2, 1.0); 
      gratio->SetPointError(q, xerr2, yerr2/yval2); 
    }

    //hratio->SetMarkerSize(2.0); 
  }

  t2->cd(); 
  t2->SetLogx(); 
  gratio->Draw("A2"); 
  //gratio->GetXaxis()->SetRangeUser(200., 1000.); 
  gratio->GetXaxis()->SetLimits(199.99, 1000.); 
  gratio->GetXaxis()->SetMoreLogLabels(); 
  gratio->GetXaxis()->SetNoExponent(); 

  gratio->GetXaxis()->SetTitleSize(0.05 * 1.0/padfrac); 
  gratio->GetXaxis()->SetTitleOffset(0.95); 
  gratio->GetXaxis()->SetLabelSize(0.035 * 1.0/padfrac); 
  gratio->GetXaxis()->SetTickLength(0.03 * 1.0/padfrac);

  gratio->GetYaxis()->SetTitleSize(0.05 * 1.0/padfrac); 
  gratio->GetYaxis()->SetTitleOffset(1.20 * padfrac); 
  gratio->GetYaxis()->SetLabelSize(0.035 * 1.0/padfrac); 
  gratio->GetYaxis()->SetRangeUser(-0.40, 2.40); 
  gratio->GetYaxis()->SetNdivisions(3); 
  //gratio->GetYaxis()->SetTitle("Data/MC"); 
  gratio->GetYaxis()->SetTitle("obs/pred"); 

  gratio->GetXaxis()->SetTitle(titlex.Data()); 
  //hratio->Draw("E1SAME"); 
  gaeratio->Draw("PSAME"); 

  t1->cd(); 

  //TLegend *leg = new TLegend(0.55, 0.65, 0.90, 0.90, NULL, "brNDC");
  TLegend *leg;
  // if(dology) {
  if     (is2011) leg = new TLegend(0.15, 0.65, 0.82, 0.90, NULL, "brNDC");
  else if(is2012) leg = new TLegend(0.15, 0.60, 0.82, 0.90, NULL, "brNDC");
  else            leg = new TLegend(0.45, 0.65, 0.95, 0.90, NULL, "brNDC");
  leg->SetNColumns(2); 
  // }
  // else 
  //   leg = new TLegend(0.57, 0.45, 0.91, 0.90, NULL, "brNDC");
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0); //(1001);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);

  if(stacksignal) { // it's SM ZZ
    leg->AddEntry(hists[processes[0]], procnames[processes[0]].Data(), "F"); 
    for(unsigned int h=processes.size()-1; h>0; --h) { 
      TString tmpProcLab(procnames[processes[h]]); 
      tmpProcLab.ReplaceAll("top", "Top"); 
      tmpProcLab.ReplaceAll("data", "Data"); 
      leg->AddEntry(hists[processes[h]], tmpProcLab.Data(), "F"); 
    }
  }
  else {           // it's ATGC
    for(int h=processes.size()-1; h>=0; --h) { 
      TString tmpProcLab(procnames[processes[h]]); 
      tmpProcLab.ReplaceAll("top", "Top"); 
      tmpProcLab.ReplaceAll("data", "Data"); 
      if(processIds[h]>0)
	leg->AddEntry(hists[processes[h]], tmpProcLab.Data(), "F"); 
      else
	leg->AddEntry(hists[processes[remappedAtgcs[h]]], procnames[processes[remappedAtgcs[h]]].Data(), "L"); 
    }
  }
  if((is2012 || is2015) && dology) leg->AddEntry("", "", ""); 
  leg->AddEntry(hdata, "Data", "PEL"); 
  leg->Draw();

  TLatex *tex, *tex2; 
  //if(is2012) tex = new TLatex(0.30, 0.97, "CMS Preliminary 2012,   #scale[0.5]{#int} L=19.6 fb^{-1},   #sqrt{s} = 8 TeV"); 
  //else       tex = new TLatex(0.33, 0.97, "CMS Preliminary 2011,   #scale[0.5]{#int} L=5.1 fb^{-1},   #sqrt{s} = 7 TeV");
  if     (is2011) tex = new TLatex(0.71, 0.98, "5.1 fb^{-1} (7 TeV)"); 
  else if(is2012) tex = new TLatex(0.69, 0.98, "19.6 fb^{-1} (8 TeV)"); 
  else            tex = new TLatex(0.69, 0.98, "2.3 fb^{-1} (13 TeV)"); 
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.052);
  tex->SetLineWidth(2);
  tex->Draw();

  //tex2 = new TLatex(0.20, 0.98, "CMS Preliminary"); // outside frame 
  tex2 = new TLatex(0.17, 0.89, "CMS"); // inside frame 
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetTextFont(62);
  tex2->SetTextSize(0.068);
  tex2->SetLineWidth(2);
  tex2->Draw();
  tex2 = new TLatex(0.17, 0.81, "Preliminary"); // inside frame 
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetTextFont(52);
  tex2->SetTextSize(0.068);
  tex2->SetLineWidth(2);
  tex2->Draw();

  TString nameoutplot = titlex;
  nameoutplot.ReplaceAll(" ", "_");
  nameoutplot.ReplaceAll("GeV/c", "");
  nameoutplot.ReplaceAll("GeV", "");
  nameoutplot.ReplaceAll("[", "");
  nameoutplot.ReplaceAll("]", "");
  nameoutplot.ReplaceAll("{", "");
  nameoutplot.ReplaceAll("}", "");
  nameoutplot.ReplaceAll("#", "");
  if     (only0jets) nameoutplot += "_eq0jets"; 
  else if(only1jets) nameoutplot += "_eq1jets"; 
  else               nameoutplot += "_lesq1jets"; 
  if(normalizebins) nameoutplot += "_norm";
  if     (is2011) nameoutplot += "_7TeV";
  else if(is2012) nameoutplot += "_8TeV";
  else            nameoutplot += "_13TeV";

  if(dology) nameoutplot += "_logy"; 
  else       nameoutplot += "_nology"; 
  nameoutplot += ".png";
  if(saveplots) cc->SaveAs(nameoutplot.Data());
  nameoutplot.ReplaceAll(".png", ".pdf");
  if(saveplots) cc->SaveAs(nameoutplot.Data());
  nameoutplot.ReplaceAll(".pdf", ".root");
  if(saveplots) cc->SaveAs(nameoutplot.Data());
  nameoutplot.ReplaceAll(".root", ".png");

  return; 

  // // Cross section part
  // //TFile *efffile  = TFile::Open("/afs/cern.ch/work/t/trocino/Work/ZZllnn/aTGC/Analysis/2012-11-16/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/test/aTGC/2013-01-14_8TeV/MC8TeV_ZZ_0.root"); 
  // TFile *efffile  = is2012 ? TFile::Open("/afs/cern.ch/work/t/trocino/Work/ZZllnn/aTGC/Analysis/2012-11-16/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/test/aTGC/2013-04-30_onlyZZ_8TeV/MC8TeV_ZZ.root") : TFile::Open("/afs/cern.ch/work/t/trocino/Work/ZZllnn/aTGC/Analysis/2012-11-16/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/test/aTGC/2013-04-30_onlyZZ_7TeV/MC7TeV_ZZ_madgraph.root"); 
  // TFile *mcfmfile = is2012 ? TFile::Open("/afs/cern.ch/user/t/trocino/ZZ2l2n/MCFM/MCFM/Bin/2013-04-17_pdfSyst_withoutAcceptance/set_zz_8000_MSTW2008nlo68cl_90/ZZlept_tota_MSTW200_90__90__test.root") : TFile::Open("/afs/cern.ch/user/t/trocino/ZZ2l2n/MCFM/MCFM/Bin/2013-04-17_pdfSyst_withoutAcceptance/set_zz_7000_MSTW2008nlo68cl_90/ZZlept_tota_MSTW200_90__90__test.root"); 
  // TFile *mcfmfile_MCcuts = is2012 ? TFile::Open("/afs/cern.ch/user/t/trocino/ZZ2l2n/MCFM/MCFM/Bin/2013-04-17_pdfSyst_withoutAcceptance_MCcuts/set_zz_8000_MSTW2008nlo68cl_90/ZZlept_tota_MSTW200_90__90__test.root") : TFile::Open("/afs/cern.ch/user/t/trocino/ZZ2l2n/MCFM/MCFM/Bin/2013-04-17_pdfSyst_withoutAcceptance_MCcuts/set_zz_7000_MSTW2008nlo68cl_90/ZZlept_tota_MSTW200_90__90__test.root"); 
  // TFile *mcfmfile_acc = is2012 ? TFile::Open("/afs/cern.ch/user/t/trocino/ZZ2l2n/MCFM/MCFM/Bin/2013-04-17_pdfSyst_withAcceptance/set_zz_8000_MSTW2008nlo68cl_90/ZZlept_tota_MSTW200_90__90__test.root") : TFile::Open("/afs/cern.ch/user/t/trocino/ZZ2l2n/MCFM/MCFM/Bin/2013-04-17_pdfSyst_withAcceptance/set_zz_7000_MSTW2008nlo68cl_90/ZZlept_tota_MSTW200_90__90__test.root"); 
  // TH1D *hmcfm        = (TH1D*)mcfmfile->Get("id19"); 
  // //TH1D *hmcfm_MCcuts = (TH1D*)mcfmfile_MCcuts->Get("id19"); 
  // //TH1D *hmcfm_acc    = (TH1D*)mcfmfile_acc->Get("id19"); 

  // bool isinclusive = false; 
  // bool correctforacc = true; 
  // bool accFromMCFM = false;

  // double xbins[20]; 
  // double errxbins[20]; 
  // double xstheo[20]; 
  // double errxstheo[20]; 
  // double xsexpe[20]; 
  // double errxsexpe[20]; 

  // std::vector<std::pair<int, int> > binning; 
  // ////binning.push_back(std::pair<int, int>(1, 5)); 
  // binning.push_back(std::pair<int, int>(3, 5)); 
  // binning.push_back(std::pair<int, int>(6, 6)); 
  // binning.push_back(std::pair<int, int>(7, 8)); 

  // //binning.push_back(std::pair<int, int>(3, 8)); 


  // //// 7 TeV
  // // pt(Z) =   0- 20:  7038 /  6752 / 4575 --> 
  // // pt(Z) =  20- 40: 14326 / 13725 / 8962 --> 
  // // pt(Z) =  40- 60: 12824 / 12327 / 7939 --> 
  // // pt(Z) =  60- 80:  8443 /  8083 / 5159 --> 
  // // pt(Z) =  80-100:  5096 /  4836 / 2930 --> 
  // // pt(Z) = 100-200:  7247 /  6848 / 4466 --> 
  // // pt(Z) = 200-400:   989 /   913 /  714 --> 
  // // pt(Z) = 400-800:    64 /    60 /   55 --> 
  // //
  // //// 8 TeV
  // // pt(Z) =   0- 20: 6121 / 3369 / 2173 --> 0.5504   / 0.355007
  // // pt(Z) =  20- 40: 8870 / 6926 / 4340 --> 0.780834 / 0.48929
  // // pt(Z) =  40- 60: 7120 / 6230 / 3887 --> 0.875    / 0.545927
  // // pt(Z) =  60- 80: 4717 / 4258 / 2642 --> 0.902692 / 0.560102
  // // pt(Z) =  80-100: 2738 / 2514 / 1485 --> 0.918188 / 0.542367
  // // pt(Z) = 100-200: 3922 / 3586 / 2284 --> 0.914329 / 0.582356
  // // pt(Z) = 200-400:  536 /  476 /  373 --> 0.88806  / 0.695896
  // // pt(Z) = 400-800:   17 /   14 /   13 --> 0.823529 / 0.764706

  // double alln_7[]     = {7038., 14326., 12824., 8443., 5096., 7247., 989., 64.}; // all generated events (MadGraph)
  // double masscutn_7[] = {6752., 13725., 12327., 8083., 4836., 6848., 913., 60.}; // generated events with mass cut (MadGraph)
  // double massaccn_7[] = {4575.,  8962.,  7939., 5159., 2930., 4466., 714., 55.}; // generated events with mass cut and acceptance (MadGraph)

  // double alln_8[]     = {6121.,  8870.,  7120., 4717., 2738., 3922., 536., 17.}; // all generated events (MadGraph)
  // double masscutn_8[] = {3369.,  6926.,  6230., 4258., 2514., 3586., 476., 14.}; // generated events with mass cut (MadGraph)
  // double massaccn_8[] = {2173.,  4340.,  3887., 2642., 1485., 2284., 373., 13.}; // generated events with mass cut and acceptance (MadGraph)

  // int nalln     = sizeof(alln_7)/sizeof(double); 
  // int nmasscutn = sizeof(masscutn_7)/sizeof(double); 
  // int nmassaccn = sizeof(massaccn_7)/sizeof(double); 

  // double *alln     = is2012 ?     alln_7 : alln_8; 
  // double *masscutn = is2012 ? masscutn_7 : masscutn_8; 
  // double *massaccn = is2012 ? massaccn_7 : massaccn_8; 

  // if( nalln!=hdata->GetNbinsX() || nmasscutn!=hdata->GetNbinsX() || nmassaccn!=hdata->GetNbinsX() ) {
  //   std::cerr << "Size " << hdata->GetNbinsX() << " does not match with " 
  // 	      << nalln << ", " << nmasscutn << ", " << nmassaccn << std::endl;
  //   //exit(1); 
  //   std::cerr << "Only one bin (1-hdata->GetNbinsX()) possible! " << std::endl; 

  //   binning.clear(); 
  //   binning.push_back(std::pair<int, int>(1, hdata->GetNbinsX())); 
  // }

  // const unsigned int nrebins(binning.size()); 

  // std::vector<double> allv(nrebins, 0.); 
  // std::vector<double> masscutv(nrebins, 0.); 
  // std::vector<double> massaccv(nrebins, 0.); 

  // for(unsigned int kk=0; kk<nrebins; ++kk) {
  //   for(int hh=binning[kk].first-1; hh<binning[kk].second; ++hh) {
  //     allv[kk]      += alln[hh]; 
  //     masscutv[kk]      += masscutn[hh]; 
  //     massaccv[kk]      += massaccn[hh]; 
  //   }
  // }

  // // Min-max scale for cross section plot
  // double ymin2 = 9999.; 
  // double ymax2 = 0.; 

  // // Retrieve post fit values
  // ifstream postfit;
  // postfit.open( ("../"+indir+"fit.txt").Data() );

  // // Exit if file opening failed
  // if(!postfit.is_open()){
  //   std::cerr << "Opening fit.txt failed!" << std::endl;
  //   exit(1);
  // }

  // double 
  //   prefit_nrb(0.), postfit_nrb(0.), prepost_nrb(1.), 
  //   prefit_dy(0.),  postfit_dy(0.),  prepost_dy(1.), 
  //   prefit_wz(0.),  postfit_wz(0.),  prepost_wz(1.); 
 
  // while(postfit.good()) {
  //   getline(postfit, line);
  //   //std::cout << line << std::endl;
  //   if(line.size()==0) continue; 

  //   TString tline(line.c_str()); 
  //   tline.ReplaceAll("&",        " "); 
  //   tline.ReplaceAll("$",        " "); 
  //   tline.ReplaceAll("\\mu\\mu", " "); 
  //   tline.ReplaceAll("ee",       " "); 
  //   tline.ReplaceAll("\\pm",     " "); 
  //   tline.ReplaceAll("\\\\",     " "); 
  //   line = tline.Data(); 

  //   std::cout << tline.Data() << std::endl; 

  //   std::vector<TString> lstr; 
  //   string temp;
  //   stringstream s(line);
  //   while(s >> temp) lstr.push_back(temp.c_str()); 
  //   int nlstr = lstr.size(); 
  //   if(nlstr==0) continue; 

  //   if(lstr[0].Contains("NonResonant")) {
  //     TString tmppre = lstr[1].ReplaceAll("\\pm", ""); 
  //     TString tmppos = lstr[nlstr-1].ReplaceAll("$", ""); 
  //     std::cout << "NonResonant " << tmppre.Data() << " " << tmppos.Data() << std::endl; 
  //     prefit_nrb += tmppre.Atof(); 
  //     postfit_nrb += tmppos.Atof(); 
  //   } 

  //   else if(lstr[0].Contains("Z+Jets")) {
  //     TString tmppre = lstr[1].ReplaceAll("\\pm", ""); 
  //     TString tmppos = lstr[nlstr-1].ReplaceAll("$", ""); 
  //     std::cout << "Z+Jets " << tmppre.Data() << " " << tmppos.Data() << std::endl; 
  //     prefit_dy += tmppre.Atof(); 
  //     postfit_dy += tmppos.Atof(); 
  //   } 

  //   else if(lstr[0].Contains("wz3lnu")) {
  //     TString tmppre = lstr[1].ReplaceAll("\\pm", ""); 
  //     TString tmppos = lstr[nlstr-1].ReplaceAll("$", ""); 
  //     std::cout << "wz3lnu " << tmppre.Data() << " " << tmppos.Data() << std::endl; 
  //     prefit_wz += tmppre.Atof(); 
  //     postfit_wz += tmppos.Atof(); 
  //   } 
  // }

  // postfit.close();

  // prepost_nrb = postfit_nrb / prefit_nrb; 
  // prepost_dy  = postfit_dy  / prefit_dy; 
  // prepost_wz  = postfit_wz  / prefit_wz; 

  // for(unsigned int pi=0; pi<nrebins; ++pi) {
  //   int bin1 = binning[pi].first;
  //   int bin2 = binning[pi].second; 

  //   int x1 = hdata4xs->GetXaxis()->GetBinLowEdge(bin1);
  //   int x2 = hdata4xs->GetXaxis()->GetBinUpEdge(bin2);

  //   //if(pi==0) xbins[pi] = x1; 
  //   //xbins[pi+1] = x2; 
  //   xbins[pi] = (x1+x2)/2.; 
  //   errxbins[pi] = (x2-x1)/2.; 

  //   std::cout << std::endl << "------------------------------------------" << std::endl; 
  //   std::cout << " pt = " << x1 << "-" << x2 << std::endl << std::endl; 

  //   int xsbin1 = hmcfm->GetXaxis()->FindBin(x1+0.001); 
  //   int xsbin2 = hmcfm->GetXaxis()->FindBin(x2-0.001); 
  //   double refxs, refxserr; 
  //   refxs = hmcfm->IntegralAndError(xsbin1, xsbin2, refxserr, "width"); // this xs is with NO acceptance cuts!! (pt>0, |eta|<inf)

  //   //allv[pi] = refxs;
  //   //masscutv[pi] = hmcfm_MCcuts->Integral(xsbin1, xsbin2, "width"); ; 
  //   //massaccv[pi] = hmcfm_acc->Integral(xsbin1, xsbin2, "width"); ; 

  //   refxs *= 3; refxserr *= 3; // from single-flavour to three flavours

  //   double N, errN, Nzz, errNzz, Nwz, errNwz, Ndy, errNdy, Nnr, errNnr; 

  //   N   = hdata4xs->IntegralAndError(bin1, bin2, errN); 
  //   Nzz = hists4xs["zz2l2nu"]->IntegralAndError(bin1, bin2, errNzz); 
  //   Nwz = hists4xs["wz3lnu"]->IntegralAndError(bin1, bin2, errNwz); 
  //   Ndy = hists4xs["zlldata"]->IntegralAndError(bin1, bin2, errNdy); 
  //   Nnr = hists4xs["topwwwjetsdata"]->IntegralAndError(bin1, bin2, errNnr); 

  //   Nwz *= prepost_wz;  errNwz *= prepost_wz; 
  //   Ndy *= prepost_nrb; errNdy *= prepost_nrb; 
  //   Nnr *= prepost_dy;  errNnr *= prepost_dy; 

  //   double Nsig = N - Nwz - Ndy - Nnr; 
  //   if(Nsig<0.) Nsig = 10e-8; 
  //   // N.B. These errors include correlations, but they only make sense when summed together
  //   //      (errors on single background sources do not make sense, because they already 
  //   //       include part of the correlations with other backgrounds)
  //   double errNsig = sqrt(errN*errN + errNwz*errNwz + errNdy*errNdy + errNnr*errNnr); 

  //   double Lumi = is2012 ? 19.577 : 5.051; 
  //   if(isinclusive) {
  //     Lumi *= 1000.; // fb-1 --> pb-1
  //     refxs /= 1000.; refxserr /= 1000.; // fb --> pb
  //   }

  //   double br = 1.0; 
  //   if(isinclusive) br = 0.0403896; // BR(ZZ -> 2l2nu)
  //   refxs /= br; refxserr /= br; 

  //   double acc = 0.516366998; // from MCFM:  (pt>20 && |eta|<2.5 && 60<mZ<120) / (60<mZ<120)
  //   if(!correctforacc) {
  //     refxs *= acc; 
  //     refxserr *= acc; 
  //     acc = 1.; // you're not going to use it anymore...
  //   }

  //   //double acc = (8676.+5851.+1687.+970.+13.)/(13695.+9602.+2741.+1321.+14.) // from MadGraph:  (pt>20 && |eta|<2.5 && 60<mZ<120) / (60<mZ<120)  --> 0.628246813  (21.67% greater than acc from MCFM)

  //   // This is to adjust the number of generated events
  //   //double accMG = (8676.+5851.+1687.+970.+13.)/(18916.+10650.+2982.+1476.+17.); // from MadGraph: (pt>20 && |eta|<2.5 && 60<mZ<120) / all
  //   //double accMass = (13695.+9602.+2741.+1321.+14.)/(18916.+10650.+2982.+1476.+17.); // from MadGraph: (60<mZ<120) / all

  //   double accMG   = massaccv[pi]/allv[pi]; // from MadGraph: (pt>20 && |eta|<2.5 && 60<mZ<120) / all
  //   double accMass = masscutv[pi]/allv[pi]; // from MadGraph: (60<mZ<120) / all
  //   if(!correctforacc) accMass = accMG; 

  //   TH1D *hcutflow = (TH1D*)efffile->Get("all_cutflow"); 
  //   //TH1D *hevtflow = (TH1D*)efffile->Get("all_eventflow"); 
  //   TH1D *hzptbeg = (TH1D*)efffile->Get("all_zpt_begin"); 
  //   TH1D *hzptend = (TH1D*)efffile->Get("all_zpt_end"); 
  //   // Add overflow to last bin
  //   int lastbin = hzptbeg->GetNbinsX(); 
  //   double lastbincont = hzptbeg->GetBinContent(lastbin); 
  //   double lastbinerr  = hzptbeg->GetBinError(lastbin); 
  //   lastbincont += hzptbeg->GetBinContent(lastbin+1); 
  //   lastbinerr *= lastbinerr;
  //   lastbinerr =+ hzptbeg->GetBinError(lastbin+1)*hzptbeg->GetBinError(lastbin+1); 
  //   lastbinerr = sqrt(lastbinerr); 
  //   hzptbeg->SetBinContent(lastbin, lastbincont); 
  //   hzptbeg->SetBinError(lastbin, lastbinerr); 
  //   hzptbeg->SetBinContent(lastbin+1, 0); 
  //   hzptbeg->SetBinError(lastbin+1, 0); 
  //   // --- // 
  //   lastbincont = 0.0; lastbinerr = 0.0; 
  //   lastbincont = hzptend->GetBinContent(lastbin); 
  //   lastbinerr  = hzptend->GetBinError(lastbin); 
  //   lastbincont += hzptend->GetBinContent(lastbin+1); 
  //   lastbinerr *= lastbinerr;
  //   lastbinerr =+ hzptend->GetBinError(lastbin+1)*hzptend->GetBinError(lastbin+1); 
  //   lastbinerr = sqrt(lastbinerr); 
  //   hzptend->SetBinContent(lastbin, lastbincont); 
  //   hzptend->SetBinError(lastbin, lastbinerr); 
  //   hzptend->SetBinContent(lastbin+1, 0); 
  //   hzptend->SetBinError(lastbin+1, 0); 

  //   int zptbin1 = hzptbeg->GetXaxis()->FindBin(x1+0.001); 
  //   int zptbin2 = hzptbeg->GetXaxis()->FindBin(x2-0.001); 
  //   double nFinalWghtEvt  = hzptend->Integral(zptbin1, zptbin2); // n. final ZZ events (with weights)
  //   double nPreselWghtEvt = hzptbeg->Integral(zptbin1, zptbin2); // n. preselected ZZ events (with weights)

  //   //double nFinalWghtEvt  = hevtflow->GetBinContent(11); // n. final ZZ events (with weights)
  //   ////double nPreselWghtEvt = hevtflow->GetBinContent(1);  // n. preselected ZZ events (with weights)
  //   //double nPreselWghtEvt = hcutflow->GetBinContent(3);  // n. preselected ZZ events (with weights)

  //   //std::cout << nFinalWghtEvt << " / " << nPreselWghtEvt << " = " << nFinalWghtEvt/nPreselWghtEvt << std::endl; 

  //   double nPreselEvt     = hcutflow->GetBinContent(2);  // n. preselected ZZ events (no weights)
  //   double nGenEvt        = hcutflow->GetBinContent(1);  // n. generated ZZ events (no weights)

  //   double accXeff = 0.;
  //   if(accFromMCFM) {
  //     nGenEvt *= accMG;  // n. generated ZZ events with MadGraph with acceptance cuts (no weights)
  //     double eff = (nFinalWghtEvt/nPreselWghtEvt) * (nPreselEvt/nGenEvt); 
  //     accXeff = eff*acc; 
  //   }
  //   else {
  //     nGenEvt *= accMass;  // n. generated ZZ events with MadGraph within 60<mZ<120 (no weights)
  //     accXeff = (nFinalWghtEvt/nPreselWghtEvt) * (nPreselEvt/nGenEvt); 
  //   }

  //   double errRelTheoXs = 0.;
  //   double errRelLumi = 0.;
  //   double errRelAcc = 0.02;
  //   // Assuming relative error on ZZ as error on efficiency
  //   // -> does not include errors on lumi, acceptance, QCD+PDF
  //   double errRelEff = errNzz/Nzz; 

  //   int zzidx = -1; 
  //   for(unsigned int i=0; i<nproc; ++i) 
  //     if(processes[i].CompareTo("zz2l2nu")==0) {
  // 	zzidx = int(i); 
  // 	break; 
  //     } 

  //   if(zzidx==-1) {
  //     std::cerr << "ERROR: ZZ->2l2nu process not found! Leave theoretical systematics at 0!" << std::endl; 
  //   }
  //   else {
  //     for(unsigned int k=0; k<systnames.size(); ++k) {
  // 	std::vector<double> & vsystmp = systs[systnames[k]]; 

  // 	// Theoretical error
  // 	if(systnames[k].Contains("QCD") || systnames[k].Contains("pdf")) {
  // 	  double tmperr = vsystmp[zzidx]; 
  // 	  errRelTheoXs = sqrt( errRelTheoXs*errRelTheoXs + tmperr*tmperr ); 
  // 	}

  // 	// Luminosity
  // 	else if(systnames[k].Contains("lumi")) {
  // 	  double tmperr = vsystmp[zzidx]; 
  // 	  //errRelLumi = sqrt( errRelLumi*errRelLumi + tmperr*tmperr ); 
  // 	  errRelLumi = tmperr; 
  // 	}

  // 	// Acceptance
  // 	else if(systnames[k].Contains("acc")) {
  // 	  double tmperr = vsystmp[zzidx]; 
  // 	  //errRelAcc = sqrt( errRelAcc*errRelAcc + tmperr*tmperr ); 
  // 	  errRelAcc = tmperr; 
  // 	}

  // 	else {}
  //     }
  //   }

  //   double errRelAccXeff = sqrt(errRelLumi*errRelLumi + errRelAcc*errRelAcc + errRelEff*errRelEff); 

  //   /*
  //   double errRelTheorXs = zzXSerrs/Nzz; // Assuming relative QCD+PDF error on ZZ as error on theoretical cross section (REPLACE!) 
  //   double errTheorXs = refxs*errRelTheorXs; 
  //   refxserr *= refxserr; 
  //   refxserr += errTheorXs*errTheorXs;
  //   refxserr = sqrt(refxserr); 
  //   */
  //   refxserr = sqrt(refxserr*refxserr + errRelTheoXs*refxs*errRelTheoXs*refxs); 

  //   double xs = Nsig / (Lumi*br*accXeff); 
  //   double errRelXs = sqrt( (errNsig/Nsig)*(errNsig/Nsig) + errRelAccXeff*errRelAccXeff ); 
  //   double errXs = xs*errRelXs; 

  //   double statErrXs = errN / (Lumi*br*accXeff); 

  //   double systErrRelNsig = sqrt(errNwz*errNwz + errNdy*errNdy + errNnr*errNnr) / Nsig; 
  //   double systErrRelXs = sqrt( systErrRelNsig*systErrRelNsig + errRelAccXeff*errRelAccXeff ); 
  //   double systErrXs = xs*systErrRelXs; 

  //   //std::cout << Nsig << " " << Lumi << " " << acc << " " << eff << std::endl; 
  //   const char *unit = isinclusive ? "pb" : "fb"; 
  //   std::cout << "    MCFM cross section: (" << refxs << " +- " << /*errRelTheoXs*refxs*/ refxserr << ") " << unit << std::endl; 
  //   std::cout << "Measured cross section: (" << xs << " +- " << errXs << ") " << unit << std::endl; 
  //   std::cout << "                        (" << xs << " +- " << statErrXs << "(stat) +- " << systErrXs << "(syst)) " << unit << std::endl; 
  //   std::cout << "                        (" << xs << " +- " << sqrt(statErrXs*statErrXs + systErrXs*systErrXs) << "(stat+syst)) " << unit << std::endl; 

  //   xstheo[pi] = refxs; 
  //   errxstheo[pi] = refxserr; 
  //   xsexpe[pi] = xs; 
  //   errxsexpe[pi] = errXs; 

  //   ymin2 = refxs-refxserr<ymin2             ? refxs-refxserr : ymin2; 
  //   ymin2 = xs-errXs<ymin2 && xs-errXs>10e-4 ? xs-errXs       : ymin2; 

  //   ymax2 = refxs+refxserr>ymax2 ? refxs+refxserr : ymax2; 
  //   ymax2 = xs+errXs>ymax2       ? xs+errXs       : ymax2; 
  // }

  // ymin2 *= 0.5; 
  // ymax2 *= 1.5; 

  // //double xmin = xbins[0];
  // //double xmax = xbins[nrebins];

  // TGraphErrors *hxs_t = new TGraphErrors(nrebins, xbins, xstheo, errxbins, errxstheo); 
  // TGraphErrors *hxs_e = new TGraphErrors(nrebins, xbins, xsexpe, errxbins, errxsexpe); 

  // TCanvas *cxs = new TCanvas("cxs", "cxs", 600, 600); 
  // cxs->Draw(); 

  // TPad* t0_1 = new TPad("t0_1","t0_1", 0.0, padfrac, 1.0, 1.0);
  // t0_1->SetTopMargin(0.085); 
  // t0_1->SetBottomMargin(0.015); 
  // //t0_1->SetGridx(); 
  // //t0_1->SetGridy(); 
  // t0_1->Draw();
  // t0_1->SetLogx(); 
  // t0_1->SetLogy(); 

  // TPad* t0_2 = new TPad("t0_2","t0_2", 0.0, 0.0, 1.0, padfrac);
  // t0_2->SetPad(0.0, 0.0, 1.0, padfrac);
  // t0_2->SetTopMargin(0.05);
  // t0_2->SetBottomMargin(0.4);
  // t0_2->SetGridx();
  // t0_2->SetGridy();
  // t0_2->Draw();

  // t0_1->cd(); 
  // hxs_t->SetLineWidth(0); 
  // hxs_t->SetLineColor(0); 
  // hxs_t->SetFillColor(kRed-7); 
  // hxs_t->SetFillStyle(3013); //3001); 
  // hxs_t->Draw("A2"); 
  // //double minxaxis = t0_1->GetLogx() ? 40. : 0.; 
  // double minxaxis = 0.; 
  // hxs_t->GetXaxis()->SetLimits(xbins[0]-errxbins[0]+minxaxis, xbins[nrebins-1]+errxbins[nrebins-1]); 
  // hxs_t->GetXaxis()->SetMoreLogLabels(); 
  // hxs_t->GetXaxis()->SetNoExponent(); 
  // //hxs_t->GetXaxis()->SetTitle(titlex.Data()); 
  // hxs_t->GetXaxis()->SetLabelOffset(1.0); 
  // hxs_t->GetYaxis()->SetTitleSize(0.05 * 1.0/(1.0-padfrac)); 
  // hxs_t->GetYaxis()->SetTitleOffset(1.20 * (1.0-padfrac)); 
  // hxs_t->GetYaxis()->SetLabelSize(0.035 * 1.0/(1.0-padfrac)); 
  // hxs_t->GetYaxis()->SetTitle("#sigma(pp #rightarrow ZZ #rightarrow 2l2#nu) [fb]"); 
  // //std::cout << ymax << std::endl; 
  // //hxs_t->SetMaximum(ymax2); 
  // hxs_t->GetYaxis()->SetRangeUser(ymin2, ymax2); 
  // hxs_e->SetMarkerStyle(20); 
  // hxs_e->SetMarkerColor(kBlack); 
  // hxs_e->SetMarkerSize(1.); 
  // hxs_e->Draw("PE1SAME"); 

  // // Ratio
  // TGraphErrors *hxs_tr = (TGraphErrors*)hxs_t->Clone(); 
  // TGraphErrors *hxs_er = (TGraphErrors*)hxs_e->Clone(); 

  // for(unsigned int q=0; q<nrebins; ++q) {
  //   double yval1 = 1.;
  //   double yerr1 = errxstheo[q]/xstheo[q];
  //   double yval2 = xsexpe[q]/xstheo[q];
  //   double yerr2 = errxsexpe[q]/xstheo[q];
  //   if(yval2<10e-4) {yval2 = -1.; yerr2 = 0.;} 

  //   hxs_tr->SetPoint(q, xbins[q], yval1); 
  //   hxs_er->SetPoint(q, xbins[q], yval2); 
  //   hxs_tr->SetPointError(q, errxbins[q], yerr1); 
  //   hxs_er->SetPointError(q, errxbins[q], yerr2); 

  //   //ymin2 = xstheo[q]-errxstheo[q]<ymin2 ? xstheo[q]-errxstheo[q] : ymin2; 
  //   //ymin2 = xsexpe[q]-errxsexpe[q]<ymin2 ? xsexpe[q]-errxsexpe[q] : ymin2; 

  //   //ymax2 = xstheo[q]+errxstheo[q]>ymax2 ? xstheo[q]+errxstheo[q] : ymax2; 
  //   //ymax2 = xsexpe[q]+errxsexpe[q]>ymax2 ? xsexpe[q]+errxsexpe[q] : ymax2; 
  // }

  // //ymin2 *= 0.5; 
  // //ymax2 *= 1.1; 

  // //hxs_t->GetYaxis()->SetRangeUser(0., ymax2); 
  // //hxs_t->GetYaxis()->SetLimits(ymin2, ymax2); 

  // t0_2->cd(); 
  // t0_2->SetLogx(t0_1->GetLogx()); 
  // hxs_tr->Draw("A2"); 
  // hxs_tr->GetXaxis()->SetLimits(xbins[0]-errxbins[0]+minxaxis, xbins[nrebins-1]+errxbins[nrebins-1]); 
  // hxs_tr->GetXaxis()->SetMoreLogLabels(); 
  // hxs_tr->GetXaxis()->SetNoExponent(); 

  // hxs_tr->GetXaxis()->SetTitleSize(0.05 * 1.0/padfrac); 
  // hxs_tr->GetXaxis()->SetTitleOffset(0.95); 
  // hxs_tr->GetXaxis()->SetLabelSize(0.035 * 1.0/padfrac); 
  // hxs_tr->GetXaxis()->SetTickLength(0.03 * 1.0/padfrac);

  // hxs_tr->GetYaxis()->SetTitleSize(0.05 * 1.0/padfrac); 
  // hxs_tr->GetYaxis()->SetTitleOffset(1.20 * padfrac); 
  // hxs_tr->GetYaxis()->SetLabelSize(0.035 * 1.0/padfrac); 
  // hxs_tr->GetYaxis()->SetRangeUser(0., 2.); 
  // hxs_tr->GetYaxis()->SetNdivisions(3); 
  // hxs_tr->GetYaxis()->SetTitle("Exp/Theo"); 

  // hxs_tr->GetXaxis()->SetTitle(titlex.Data()); 
  // hxs_er->Draw("PE1SAME"); 

  // TLegend *leg2 = new TLegend(0.2, 0.04, 0.50, 0.25, NULL, "brNDC");
  // leg2->SetLineColor(0);
  // leg2->SetLineStyle(0);
  // leg2->SetLineWidth(0);
  // leg2->SetFillColor(0);
  // leg2->SetFillStyle(0); //(1001);
  // leg2->SetBorderSize(0);
  // leg2->SetTextFont(42);
  // //leg2->SetNColumns(2); 

  // leg2->AddEntry(hxs_t, "NLO SM (MCFM)", "F"); 
  // leg2->AddEntry(hxs_e, "Measured", "PE"); 

  // t0_1->cd(); 
  // leg2->Draw(); 
  // tex->Draw();
  // tex2->Draw();

  // if(nrebins>1) nameoutplot.ReplaceAll(".png", "_binned.png"); 
  // if(saveplots) cxs->SaveAs(("xs_"+nameoutplot).Data());

  // return;
}

void extendXlimit(TH1D* &h, double extendedLimit) {

  const int nHbins(h->GetNbinsX()+1); 
  double newXbins[nHbins]; 
  for(int u=1; u<nHbins; ++u) {
    newXbins[u-1] = h->GetBinLowEdge(u); 
  }
  newXbins[nHbins-1] = extendedLimit; 

  TH1D *hnew = (TH1D*)h->Clone(); 
  hnew->Reset("ICE"); 
  hnew->SetBins(nHbins-1, newXbins); 
  //hnew->Sumw2(); 

  for(int u=1; u<nHbins; ++u) {
    hnew->SetBinContent(u, h->GetBinContent(u)); 
    hnew->SetBinError(u, h->GetBinError(u)); 
  }

  //std::cout << h << "  ";
  h = hnew; 
  //std::cout << h << "  ";

  return; 
}

void changeFirstPoint(TH1D* &h, double oldPoint, double newPoint) {

  const int nHbins(h->GetNbinsX()+1); 
  double newXbins[nHbins]; 
  for(int u=1; u<nHbins; ++u) {
    newXbins[u-1] = h->GetBinLowEdge(u); 
    if( fabs(newXbins[u-1]-oldPoint)<0.01 ) newXbins[u-1] = newPoint; 
  }
  newXbins[nHbins-1] = newXbins[nHbins-2] + h->GetBinWidth(nHbins); 

  TH1D *hnew = (TH1D*)h->Clone(); 
  hnew->Reset("ICE"); 
  hnew->SetBins(nHbins-1, newXbins); 
  //hnew->Sumw2(); 

  for(int u=1; u<nHbins; ++u) {
    hnew->SetBinContent(u, h->GetBinContent(u)); 
    hnew->SetBinError(u, h->GetBinError(u)); 
  }

  //std::cout << h << "  ";
  h = hnew; 
  //std::cout << h << "  ";

  return; 
}
