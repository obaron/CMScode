//
// Feb. 11 2015
//
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
using namespace std;

static const int no_radius = 3;
static const int list_radius[no_radius] = {2,3,4};

TCanvas *tcan[no_radius];
TH1F hPbPb_pthat_fine_R[no_radius];

void drawpthat() {
TFile *nocut = new TFile("PbPb_pp_mc_akPuPF_nocut_20150211.root");
TFile *yescut = new TFile("PbPb_pp_mc_akPuPF_005cut_20150211.root");

	for(int k = 3;k<=no_radius;k++)
	{
	std::cout<<"start of loop"<<std::endl;
	tcan[k] = new TCanvas();
	std::cout<<"tcanvas made"<<std::endl;
	nocut->cd();
	std::cout<<"nocut->cd"<<std::endl;
	hPbPb_pthat_fine_R[k].SetMarkerColor(1);
	hPbPb_pthat_fine_R[k].Draw();
	std::cout<<"first hist drawn"<<std::endl;
	yescut->cd();
	hPbPb_pthat_fine_R[k].SetMarkerColor(2);
	hPbPb_pthat_fine_R[k].Draw("same");
	tcan[k]->Draw();
	std::cout<<"end of loop"<<std::endl;
	}
}




