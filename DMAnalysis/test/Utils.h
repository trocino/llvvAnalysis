#include <iostream>

#include <sstream>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TAxis.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>



void addText(double x1, double x2, double y1, double y2, TString TEXT, Color_t color, Float_t angle = 0)
{
    TPaveText* T = new TPaveText(x1,y1,x2,y2, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    T->SetTextColor(color);
    TText *text = T->AddText(TEXT);
    text->SetTextAngle(angle);
    text->SetTextAlign(22);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);
};

std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

std::string Convert (double number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

std::string Convert (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

TString Convert2TString (int number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TString Convert2TString (float number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}


TString Convert2TString (double number){
    std::ostringstream buff;
    buff<<number;
    std::string str = buff.str();
    TString t_str = str;
    return t_str;
}




