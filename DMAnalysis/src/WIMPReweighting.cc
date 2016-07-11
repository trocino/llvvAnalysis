//
//  WIMPReweighting.cpp
//
//
//  Created by RENJIE WANG on 3/8/15.
//
//

#include "llvvAnalysis/DMAnalysis/interface/WIMPReweighting.h"

using namespace std;

//
WIMPReweighting::WIMPReweighting(const edm::ParameterSet &runProcess)
{
    Init(runProcess); 
}

//
WIMPReweighting::WIMPReweighting(const edm::ParameterSet &runProcess, TString &url)
{
    Init(runProcess, url); 
}

//
void WIMPReweighting::Init(const edm::ParameterSet &runProcess)
{

    std::vector<std::string> WIMPFiles, MCFiles; 

    if(runProcess.existsAs<std::vector<std::string> >("genwimpweights")) WIMPFiles = runProcess.getParameter<std::vector<std::string> >("genwimpweights");
    if(runProcess.existsAs<std::vector<std::string> >("MCweights"))      MCFiles = runProcess.getParameter<std::vector<std::string> >("MCweights");

    for(size_t ifile=0; ifile<WIMPFiles.size(); ifile++) {

        TString File(WIMPFiles[ifile].c_str());
	if(!File.Contains(".root")) continue;
        gSystem->ExpandPathName(File);
        TFile *wimpFile=TFile::Open(File);

        if(wimpFile) {
            cout << "[WIMPReweighting] retrieving WIMP weights from: " << File << endl;
            std::vector<TString> WimpKeys1D;
            WimpKeys1D.push_back("genmet");
	    WimpKeys1D.push_back("pt_chichi");
	    //WimpKeys1D.push_back("pileup_weights");
            for(size_t itag=0; itag<WimpKeys1D.size(); itag++) {
                    TString key = WimpKeys1D[itag];
                    cout << "key: " << key << endl;
                    TH1F *h = (TH1F *) wimpFile->Get(key);
		    if(h==0) continue;
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    wimpWeights1DH_[key] = h;
            }//itag

            std::vector<TString> WimpKeys2D;
            WimpKeys2D.push_back("dphi_vs_met");

            for(size_t itag=0; itag<WimpKeys2D.size(); itag++) {
                    TString key = WimpKeys2D[itag];
                    cout << "key: " << key << endl;
                    TH2F *h = (TH2F *) wimpFile->Get(key);
		    if(h==0) continue;
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    wimpWeights2DH_[key] = h;
            }//itag


	    cout << "[WIMPReweighting] close file: " << File << endl;
        }
        wimpFile->Close();
    }


    for(size_t ifile=0; ifile<MCFiles.size(); ifile++) {

        TString File(MCFiles[ifile].c_str());
	if(!File.Contains(".root")) continue;
        gSystem->ExpandPathName(File);
        TFile *wimpFile=TFile::Open(File);

        if(wimpFile) {
            cout << "[WIMPReweighting] retrieving MC weights from: " << File << endl;
            std::vector<TString> WimpKeys1D;
	    WimpKeys1D.push_back("pileup_weights");
            for(size_t itag=0; itag<WimpKeys1D.size(); itag++) {
                    TString key = WimpKeys1D[itag];
                    cout << "key: " << key << endl;
                    TH1F *h = (TH1F *) wimpFile->Get(key);
		    if(h==0) continue;
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    wimpWeights1DH_[key] = h;
            }//itag

	    cout << "[WIMPReweighting] close file: " << File << endl;
        }
        wimpFile->Close();
    }

}



//
bool WIMPReweighting::Init(const edm::ParameterSet &runProcess, TString &url)
{
    // Retrieve names of reference full-sim samples 
    refPoints_ = runProcess.getParameter<std::vector<std::string> >("ReferencePoints"); 

    // Name of target gen-distribution plot 
    TString hnumlabel = exctractDistrNameFromUrl(url); 

    // Switch to URL of reference full-sim sample 
    bool successfullySwitched = selectReferenceWIMPurl(url); 
    if(!successfullySwitched) {
        std::cerr << " *** WARNING: couldn't retrieve switch to reference URL! This should never happen! *** " << std::endl; 
    } 

    // Name of reference gen-distribution plot 
    TString hdenlabel = exctractDistrNameFromUrl(url); 

    std::string genDistrFileName = runProcess.getParameter<std::string>("GenDistributionsFileName"); 
    TFile *genDistrFile = TFile::Open(genDistrFileName.c_str()); 

    if(genDistrFile!=0 && !genDistrFile->IsZombie()) { 
        TH1F *hden = (TH1F *)genDistrFile->Get(hdenlabel.Data()); 
	TH1F *hnum = (TH1F *)genDistrFile->Get(hnumlabel.Data());
        
	if(hden!=0 && hnum!=0) {
	    //hden->Scale(1./hden->Integral(0,-1)); 
	    //hnum->Scale(1./hnum->Integral(0,-1)); 
	    TH1F *h = (TH1F*)hnum->Clone((hnumlabel+"_weight").Data());

	    // Loop on all bins, including underflow and overflow 
	    for(int i=0; i<=h->GetNbinsX()+1; ++i) {
	        if(hden->GetBinContent(i)>0.) {
		    h->SetBinContent(i, hnum->GetBinContent(i)/hden->GetBinContent(i));
		}
		else {
		    std::cerr << " *** WARNING: bin " << i << " (out of " << h->GetNbinsX()
			      << " bins) has denominator 0. Setting the ratio to 1.0"
			      << " (OK only if bin == 0) ***" << std::endl; 
		    h->SetBinContent(i, 1.); 
		} 
	    }
	    h->SetDirectory(0);
	    wimpWeights1DH_["genmet_acc_simplmod"] = h; 
	}
	else { 
	    std::cerr << " *** WARNING: couldn't find histos with gen-distributions to compute weights! *** " << std::endl;
	    return false; 
	} 
    } 
    else {
        std::cerr << " *** WARNING: couldn't open file with gen-distributions to compute weights! *** " << std::endl;
	return false; 
    }

    return true; 
}



bool WIMPReweighting::selectReferenceWIMPurl(TString &url) { 
    std::pair<float, float> mxmv = extractMassesFromUrl(url); 
    if(mxmv.first<1e-6 || mxmv.second<1e-6)
      return false; 

    double mindist(999999.); 
    size_t minref(refPoints_.size()); 
    for(size_t i=0; i<refPoints_.size(); ++i) { 
        std::pair<float, float> imxmv = extractMassesFromUrl( TString(refPoints_[i].c_str()) ); 
	if(imxmv.first<1e-6 || imxmv.second<1e-6) continue; 
	if(imxmv.first<mxmv.first || imxmv.second<mxmv.second || (imxmv.first==mxmv.first && imxmv.second==mxmv.second)) continue; 
	float diff = sqrt( pow(std::log10(imxmv.first/mxmv.first), 2) + pow(std::log10(imxmv.second/mxmv.second), 2) ); 
	if(diff<mindist) { 
	    mindist = diff; 
	    minref = i; 
	} 
    } 
    if(minref>=refPoints_.size()) 
        return false; 

    url = url(0, url.Last('/')+1) + refPoints_[minref].c_str(); 
    return true; 
} 


std::pair<float, float> WIMPReweighting::extractMassesFromUrl(TString aurl) { 
    // Get the values of Mx and Mv 
    size_t mxidx0 = 3 + aurl.Index("_Mx"); 
    size_t mxidx1 =     aurl.Index("Mv", mxidx0); 
    size_t mvidx1 =     aurl.Index(".root", mxidx1); 
    TString mxstr = aurl(mxidx0,   mxidx1-mxidx0  ); 
    TString mvstr = aurl(mxidx1+2, mvidx1-mxidx1-2); 

    return std::make_pair(mxstr.Atof(), mvstr.Atof()); 
} 


TString WIMPReweighting::exctractDistrNameFromUrl(TString aurl) { 
    // Is it vector or axial-vector? (To be updated when new scenaria become available)
    // For now if it's not vector, just assume it's axial-vector 
    int isVect = aurl.Contains("TeV_DM_V_Mx") ? 1 : 0; 

    // Get the value of gQ
    //   For now just assume dtag is
    //   "MC13TeV_DM_V_Mx1000Mv10.root" (for gQ = 1)  or 
    //   "MC13TeV_DM_V_GQ0p25_Mx1000Mv10.root" (for gQ = 0.25) 
    // 
    TString gqstr = aurl.Contains("_GQ0p25_") ? "0.25" : "1";

    // Get the values of Mx and Mv 
    size_t mxidx0 = 3 + aurl.Index("_Mx"); 
    size_t mxidx1 =     aurl.Index("Mv", mxidx0); 
    size_t mvidx1 =     aurl.Index(".root", mxidx1); 
    TString mxstr = aurl(mxidx0,   mxidx1-mxidx0  ); 
    TString mvstr = aurl(mxidx1+2, mvidx1-mxidx1-2); 

    return TString::Format("monoz_weights_cV%d_cA%d_gDM1_gQ%s_Mx%s_Mmed%s",
			   isVect, 1-isVect, gqstr.Data(), mxstr.Data(), mvstr.Data()); 
} 




//
double WIMPReweighting::get1DWeights(double xval, TString key)
{
    double weight_ = 1.;
    TH1F* h_ = wimpWeights1DH_[key];
    if(h_==0) {
	//cout << "cannot find hist: " << key << " weight will be return as 1" << endl;
	return weight_;
    }

    // 
    // For simplified models, cut it short...
    // We have enough events in the overflow,
    // and there shouldn't be events in the underflow
    // in final distributions
    // 
    if(key.EqualTo("genmet_acc_simplmod"))  
      return h_->GetBinContent( h_->GetXaxis()->FindBin(xval) ); 

    // This is only for EWKDM 
    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    weight_ = h_->GetBinContent(binx);

    return weight_;
}


//
double WIMPReweighting::get2DWeights(double xval, double yval, TString key)
{
    double weight_ = 1.;
    TH2F* h_ = wimpWeights2DH_[key];
    if(h_==0) {
	//cout << "cannot find hist: " << key << " weight will be return as 1" << endl;
	return weight_;
    }

    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int ybins = h_->GetYaxis()->GetNbins();
    if(yval > h_->GetYaxis()->GetBinUpEdge(ybins)    ) yval = h_->GetYaxis()->GetBinUpEdge(ybins);
    if(yval < h_->GetYaxis()->GetBinLowEdge(1)       ) yval = h_->GetYaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    int biny = h_->GetYaxis()->FindBin(yval);
    weight_ = h_->GetBinContent(binx,biny);

    return weight_;
}


//
WIMPReweighting::~WIMPReweighting()
{
}

