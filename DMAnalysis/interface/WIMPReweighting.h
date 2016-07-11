//
//  WIMPReweighting.h
//
//
//  Created by RENJIE WANG on 3/8/15.
//
//

#ifndef ____WIMPReweighting__
#define ____WIMPReweighting__

#include <iostream>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"

#include "Math/LorentzVector.h"
#include "TVector2.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"



class WIMPReweighting {


public:

    WIMPReweighting() {};
    WIMPReweighting(const edm::ParameterSet &runProcess);
    WIMPReweighting(const edm::ParameterSet &runProcess, TString &url);
    void Init(const edm::ParameterSet &runProcess);
    bool Init(const edm::ParameterSet &runProcess,  TString &url);
    double get1DWeights(double xval, TString key);
    double get2DWeights(double xval, double yval, TString key);
    bool selectReferenceWIMPurl(TString &url); 
    std::pair<float, float> extractMassesFromUrl(TString url); 
    TString exctractDistrNameFromUrl(TString url); 

    ~WIMPReweighting();

private:
    std::map<TString, TH1F *> wimpWeights1DH_;
    std::map<TString, TH2F *> wimpWeights2DH_;

    std::vector<std::string> refPoints_;
};

#endif /* defined(____WIMPReweighting__) */
