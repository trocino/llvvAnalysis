#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TString.h"

#include "Utils.h"

using namespace std;

Double_t findIntersection(UInt_t nxs, Double_t *xs, Double_t *ys)
{
    UInt_t pos = 0;
    ys[0] -= 1;
    for(UInt_t i=1; i<nxs; ++i) {
        ys[i] -= 1;
        if( ys[i-1]*ys[i]<0. ) {
            pos = i;
            break;
        }
    }

    if(pos>0) {
        Double_t tga = (ys[pos] - ys[pos-1]) / (xs[pos] - xs[pos-1]); // can be negative
        Double_t deltax = ys[pos]/tga;                                // must be positive
        return (xs[pos] - deltax);
    }

    return 999.;
}

//main
///***************************/
//void plotInvHiggs(TString Main = "")
//{
//    TString myFile = Main;
//    makeBetterLimitPlot(myFile);
//}
///***************************/


void makeBetterLimitPlot(string folder)
{

    TString tag = "";
    TString myfolder = folder.c_str();
    if(myfolder.Contains("/D5")) tag="D5";
    if(myfolder.Contains("/D8")) tag="D8";
    if(myfolder.Contains("/D9")) tag="D9";
    if(myfolder.Contains("/C3")) tag="C3";
    if(myfolder.Contains("/Unpart")) tag="Unpart";

    string dbpars[] = {"110","125","150","200","300"};
    string parName = "m_{H} [GeV]";
    bool is7TeV  = false;
    bool is8TeV  = false;
    bool is13TeV = true;
    bool savePlots = true;
    bool normalizePlot = true;
    bool showObserved = false;


    if(is7TeV && is8TeV) normalizePlot = true; // doesn't make sense to use xs in combination 7+8 TeV

    UInt_t nPts = sizeof(dbpars)/sizeof(string);
    const UInt_t npts(nPts);


    Double_t xPts[npts];
    Double_t obs[npts];
    Double_t obserr[npts];
    Double_t m2s[npts];
    Double_t m1s[npts];
    Double_t exp[npts];
    Double_t p1s[npts];
    Double_t p2s[npts];

    // Double_t xPts_l[npts];
    // Double_t obs_l[npts];
    // Double_t obserr_l[npts];
    // Double_t obsm_l[npts];
    // Double_t obsp_l[npts];
    // Double_t m2s_l[npts];
    // Double_t m1s_l[npts];
    // Double_t exp_l[npts];
    // Double_t p1s_l[npts];
    // Double_t p2s_l[npts];

    Double_t xPts_r[npts];
    Double_t obs_r[npts];
    Double_t obserr_r[npts];
    Double_t obsm_r[npts];
    Double_t obsp_r[npts];
    Double_t m2s_r[npts];
    Double_t m1s_r[npts];
    Double_t exp_r[npts];
    Double_t p1s_r[npts];
    Double_t p2s_r[npts];

    Double_t miny = 0.9;
    Double_t maxy = 1.1;

    for(UInt_t i=0; i<npts; ++i) {

        if(dbpars[i].compare("0") == 0) {
            xPts[i] = 0.;
            obs[i] = 999999.;
            m2s[i] = 999999.;
            m1s[i] = 999999.;
            exp[i] = 999999.;
            p1s[i] = 999999.;
            p2s[i] = 999999.;
            continue;
        }
        xPts[i] = atof(dbpars[i].c_str());
	//xPts[i] /= 100.;

        //TFile myfile( (folder+"/"+dbpars[i]+"/higgsCombineTest.Asymptotic.mH"+dbpars[i]+".root").c_str() );
        TFile myfile( (folder+"/"+dbpars[i]+"/higgsCombineTest.Asymptotic.mH1.root").c_str() );
        //TFile myfile( (folder+"/"+dbpars[i]+"/higgsCombineZHinv0jets.Asymptotic.mH"+dbpars[i]+".root").c_str() );
        //TFile myfile( (folder+"/"+dbpars[i]+"/higgsCombineZHinv1jets.Asymptotic.mH"+dbpars[i]+".root").c_str() );

        if ( myfile.IsOpen() && !myfile.IsZombie() ) {

            Double_t lim, limerr;
            TTree *t = (TTree*)myfile.Get("limit");
            t->SetBranchAddress("limit",    &lim);
            t->SetBranchAddress("limitErr", &limerr);

            // Expected  2.5%
            t->GetEntry(0);
            m2s[i] = lim;
            if( m2s[i]<miny ) miny = m2s[i];
	    //std::cout << i << ": m2s = " << m2s[i] << std::endl; 

            // Expected  16.0%
            t->GetEntry(1);
            m1s[i] = lim;
	    //std::cout << i << ": m1s = " << m1s[i] << std::endl; 

            // Expected 50.0%
            t->GetEntry(2);
            exp[i] = lim;
	    //std::cout << i << ": exp = " << exp[i] << std::endl; 

            // Expected 84.0%
            t->GetEntry(3);
            p1s[i] = lim;
	    //std::cout << i << ": p1s = " << p1s[i] << std::endl; 

            // Expected 97.5%
            t->GetEntry(4);
            p2s[i] = lim;
            if( p2s[i]>maxy ) maxy = p2s[i];
	    //std::cout << i << ": p2s = " << p2s[i] << std::endl; 

            // Observed limit
            t->GetEntry(5);
            obs[i] = lim;
            obserr[i] = limerr;

            if( obs[i]-obserr[i]<miny ) miny = obs[i]-obserr[i];
            if( obs[i]+obserr[i]>maxy ) maxy = obs[i]+obserr[i];
        }
        myfile.Close();
    }

    UInt_t negpts(0), pospts(0);

    for(UInt_t i=0; i<npts; ++i) {
      // if(xPts[i]<0) {
      //     xPts_l[negpts] = xPts[i];
      //     m2s_l[negpts] = exp[i] - m2s[i];
      //     m1s_l[negpts] = exp[i] - m1s[i];
      //     exp_l[negpts] = exp[i];
      //     p1s_l[negpts] = p1s[i] - exp[i];
      //     p2s_l[negpts] = p2s[i] - exp[i];
      //     obs_l[negpts] = obs[i];
      //     obserr_l[negpts] = obserr[i];
      //     obsm_l[negpts] = obs[i] - obserr[i];
      //     obsp_l[negpts] = obs[i] + obserr[i];
      //     ++negpts;
      // } else if(xPts[i]>0) {
      xPts_r[pospts] =  xPts[i];
      m2s_r[pospts]  =  exp[i] - m2s[i];
      m1s_r[pospts]  =  exp[i] - m1s[i];
      exp_r[pospts]  =  exp[i]; 
      p1s_r[pospts]  = -exp[i] + p1s[i];
      p2s_r[pospts]  = -exp[i] + p2s[i];
      obs_r[pospts]  =  obs[i];
      obserr_r[pospts] = obserr[i];
      obsm_r[pospts] = obs[i] - obserr[i];
      obsp_r[pospts] = obs[i] + obserr[i];

      ++pospts;
      // } else {}
    }

    /*
    cout << "\n =========================================== \n" << endl;
//
    cout << "double unparticle_du[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << xPts_r[j];
            else cout << "," << xPts_r[j];
    }
    cout << "};\n" << endl;
//
//
    cout << "double unparticle_obs[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << obs_r[j];
            else cout << "," << obs_r[j];
    }
    cout << "};\n" << endl;
//
//
    cout << "double unparticle_exp[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << exp_r[j];
            else cout << "," << exp_r[j];
    }
    cout << "};\n" << endl;
//
//
    cout << "double unparticle_m1s[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << m1s_r[j];
            else cout << "," << m1s_r[j];
    }
    cout << "};\n" << endl;
//
//
    cout << "double unparticle_m2s[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << m2s_r[j];
            else cout << "," << m2s_r[j];
    }
    cout << "};\n" << endl;
//
//
    cout << "double unparticle_p1s[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << p1s_r[j];
            else cout << "," << p1s_r[j];
    }
    cout << "};\n" << endl;
//
//
    cout << "double unparticle_p2s[" << npts << "]={";
    for(UInt_t j=0; j<npts; j++){
            if(j==0) cout << p2s_r[j];
            else cout << "," << p2s_r[j];
    }
    cout << "};\n" << endl;
//
    cout << "\n =========================================== \n" << endl;
    */

    // TGraph *grObs_l;
    // TGraphErrors *grObsErr_l;
    // TGraph *grExp_l;
    // TGraphAsymmErrors *gr1s_l;
    // TGraphAsymmErrors *gr2s_l;
    // if(negpts>0) {
    //     grObs_l = new TGraph(negpts, xPts_l, obs_l);
    //     grObs_l->SetLineWidth(2);
    //     grObs_l->SetLineStyle(1);
    //     //grObs_l->SetMarkerStyle(20);
    //     //grObs_l->SetMarkerSize(0);//RENJIE

    //     grObsErr_l = new TGraphErrors(negpts, xPts_l, obs_l, 0, obserr_l);
    //     grObs_l->SetFillColor(kGray);
    //     grObs_l->SetFillStyle(3002);

    //     grExp_l = new TGraph(negpts, xPts_l, exp_l);
    //     grExp_l->SetLineWidth(2);
    //     grExp_l->SetLineStyle(2);
    //     grExp_l->SetMarkerSize(0); //RJ

    //     gr1s_l = new TGraphAsymmErrors(negpts, xPts_l, exp_l, 0, 0, m1s_l, p1s_l);
    //     gr1s_l->SetFillColor(kGreen);

    //     gr2s_l = new TGraphAsymmErrors(negpts, xPts_l, exp_l, 0, 0, m2s_l, p2s_l);
    //     gr2s_l->SetFillColor(kYellow);
    // }


    TGraph *grObs_r = 0;
    TGraphErrors *grObsErr_r = 0;
    TGraph *grExp_r = 0;
    TGraphAsymmErrors *gr1s_r = 0;
    TGraphAsymmErrors *gr2s_r = 0;


    if(pospts>0) {
        grObs_r = new TGraph(pospts, xPts_r, obs_r);
        grObs_r->SetLineWidth(2);
        grObs_r->SetLineStyle(1);
        grObs_r->SetMarkerStyle(20);
        //grObs_r->SetMarkerSize(0); //RENJIE

        grObsErr_r = new TGraphErrors(pospts, xPts_r, obs_r, 0, obserr_r);
        grObs_r->SetFillColor(kGray);
        grObs_r->SetFillStyle(3002);

        grExp_r = new TGraph(pospts, xPts_r, exp_r);
        grExp_r->SetLineWidth(2);
        grExp_r->SetLineStyle(2);
        grExp_r->SetMarkerSize(0); //RJ


        gr1s_r = new TGraphAsymmErrors(pospts, xPts_r, exp_r, 0, 0, m1s_r, p1s_r);
        gr1s_r->SetFillColor(kGreen);

        gr2s_r = new TGraphAsymmErrors(pospts, xPts_r, exp_r, 0, 0, m2s_r, p2s_r);
        gr2s_r->SetFillColor(kYellow);
    }


    TMultiGraph *mg = new TMultiGraph();
    //if(gr2s_l)                     mg->Add(gr2s_l);
    if(gr2s_r)                     mg->Add(gr2s_r);
    //if(gr1s_l)                     mg->Add(gr1s_l);
    if(gr1s_r)                     mg->Add(gr1s_r);
    //if(grObsErr_l && showObserved) mg->Add(grObsErr_l);
    if(grObsErr_r && showObserved) mg->Add(grObsErr_r);
    //if(grExp_l)                    mg->Add(grExp_l, "L");
    if(grExp_r)                    mg->Add(grExp_r, "L");
    //if(grObs_l && showObserved)    mg->Add(grObs_l, "PL");
    if(grObs_r && showObserved)    mg->Add(grObs_r, "PL");


    TCanvas *canv = new TCanvas("canv", "limits canvas", 800., 650.);
    canv->cd();
    canv->SetGridx(0);
    canv->SetGridy(0);
    TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
    t1->Draw();
    t1->cd();
    t1->SetLogy(true);
    t1->SetLogx(false);


    mg->Draw("a3"); //RJ
    Double_t deltaZ = maxy-miny;
    miny /= 1.5;
    maxy *= 2.0;
    //mg->SetMinimum(5e-3);
    //mg->SetMaximum(1000);//45.25);//10.25);
    //if(tag=="D9") mg->SetMaximum(9000);

    mg->GetXaxis()->SetTitle(parName.c_str());
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetTitleOffset(1.1);
    //if(normalizePlot) mg->GetYaxis()->SetTitle("95% CL limit on #Lambda [GeV]");
    //else              mg->GetYaxis()->SetTitle("95% CL limit on #sigma#timesBR [fb]");
    mg->GetYaxis()->SetTitle("95% CL limit on #sigma #times B(inv)/#sigma_{th}");
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleOffset(1.35);

    TLine *ll;
    ll = new TLine(xPts[0], 1., xPts[npts-1], 1.);
    ll->SetLineWidth(2);
    ll->SetLineColor(kRed);
    ll->SetLineStyle(5);
    ll->Draw();

    float posx1 = 0.2;
    float posx2 = 0.48;
    float posy1 = 0.6;
    float posy2 = 0.82;
    TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetBorderSize(0);


    if(showObserved) {
      //if(grObs_l) leg->AddEntry(grObs_l, "Observed", "PL");
      //else 
      leg->AddEntry(grObs_r, "Observed", "PL");
    }

    //if(grExp_l)   leg->AddEntry(grExp_l, "Expected", "PL");
    //else 
    leg->AddEntry(grExp_r, "Expected", "PL");

    //if(gr1s_l)    leg->AddEntry(gr1s_l, "Expected #pm 1#sigma", "F");
    //else 
    leg->AddEntry(gr1s_r, "Expected #pm 1#sigma", "F");
    //if(gr2s_l)    leg->AddEntry(gr2s_l, "Expected #pm 2#sigma", "F");
    //else 
    leg->AddEntry(gr2s_r, "Expected #pm 2#sigma", "F");

    leg->Draw();


    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    char Buffer[1024];
    double iEcm_7   = 7.0;
    double iEcm_8   = 8.0;
    double iEcm_13  = 13.0;
    double iLumi_7  = 5051;
    double iLumi_8  = 19712;
    double iLumi_13 = 2110;


    T = new TPaveText(0.47, 0.85, 0.73, 0.80, "NDC");
    sprintf(Buffer, "pp #rightarrow ZH #rightarrow l^{+}l^{-} + inv");
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->Draw("same");
    T->SetBorderSize(0);


    if(is13TeV) {
        T = new TPaveText(0.45, 0.90, 0.75, 0.85, "NDC");
        sprintf(Buffer, "#sqrt{s} = %.1f TeV, L = %.1f fb^{-1}", iEcm_13, iLumi_13/1000);
        T->AddText(Buffer);
        T->SetTextAlign(12);
        T->SetTextFont(42);
	T->SetFillColor(0);
	T->SetFillStyle(0);
        T->Draw("same");
        T->SetBorderSize(0);
    }

    T = new TPaveText(0.21,0.92,0.36,0.82, "NDC");
    sprintf(Buffer, "#splitline{#bf{CMS}}{#it{Preliminary}}");
    //sprintf(Buffer, "#splitline{#bf{CMS}}{unpublished}");
    T->AddText(Buffer);
    T->SetTextFont(42);
    //T->Draw("same");
    T->SetBorderSize(0);

    addText(0.21,0.36,0.92,0.82,"#splitline{#bf{CMS}}{#it{             }}",kBlack);

    if( savePlots ) {
        string plotName = parName;
        string toReplace[] = {"_", "{", "}", "#", "^"};
        UInt_t nvars = sizeof(toReplace)/sizeof(string);
        for(UInt_t k=0; k<nvars; ++k) {
            int poschar = plotName.find(toReplace[k].c_str());
            while( poschar>-1 ) {
                plotName.replace(poschar, 1, "");
                poschar = plotName.find(toReplace[k].c_str());
            }
        }
        plotName = "ZH_llMET_"+tag;
        if(is8TeV) plotName += "_8TeV";
        if(!showObserved) plotName += "_noObs";
        if(myfolder.Contains("cardsShape")) plotName += "_Shape";
        canv->SaveAs( (folder+"/"+plotName+"_sigmaNorm.png").c_str() );
        //canv->SaveAs( (folder+"/"+plotName+"_sigmaNorm.pdf").c_str() );
        //canv->SaveAs( (folder+"/"+plotName+"_sigmaNorm.eps").c_str() );
        canv->SaveAs( (folder+"/"+plotName+"_sigmaNorm.root").c_str() );
        //canv->SaveAs( (folder+"/"+plotName+".C").c_str() );
    }

    return;
}

