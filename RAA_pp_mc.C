// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN

// 
// read all the MC files for PbPb and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
// 

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TStyle>
#include <TH2F>
#include <TMath.h>
  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

//static const int nAlgos = 9;
//static const int BinLabelN = 11;
//remember to change this to run akPu3PF for pPb and akVs3PF for Pbpb datasets. or just create a separate header file which will be way easier. 
//static const char *algoName[nAlgos] = { "", "icPu5", "akPu2PF", "akPu3PF", "akPu4PF", "akPu5PF" , "akPu2Calo", "akPu3Calo", "akPu4Calo" };
//static const char *algoNamePP[nAlgos] = { "", "icPu5", "ak2PF", "ak3PF", "ak4PF", "ak5PF" , "ak2Calo", "ak3Calo", "ak4Calo" };
//static const char *algoNameGen[nAlgos] = { "", "icPu5", "akPu2PF", "akVs3PF", "akPu4PF", "akPu2PF", "akPu3PF", "akPu4PF" };
//static const char *BinLabel[BinLabelN] = {"100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180", "180-200", "200-240","240-300" };


// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}

class JetData
{
public:
  JetData(char *fileName, char *jetTree, char *genJetTree, bool loadGenJet = 0,bool isPbPb = 0) {
		cout <<"Open "<<fileName<<endl;
		tFile = new TFile(fileName,"read");
		tEvt = (TTree*)tFile->Get("hiEvtAnalyzer/HiTree");
		tSkim = (TTree*)tFile->Get("skimanalysis/HltTree");
		tJet = (TTree*)tFile->Get(jetTree);
		tJet->SetBranchAddress("jtpt" , jtpt );
		tJet->SetBranchAddress("rawpt" , rawpt ); //set branch address for rawpt, like for jtpt
		tJet->SetBranchAddress("trackMax" , trackMax );
		tJet->SetBranchAddress("chargedMax",chargedMax);
		tJet->SetBranchAddress("chargedSum",chargedSum);
		tJet->SetBranchAddress("neutralMax",neutralMax);
		tJet->SetBranchAddress("neutralSum",neutralSum);
		tJet->SetBranchAddress("refpt", refpt);
		tJet->SetBranchAddress("nref" ,&njets);
		tJet->SetBranchAddress("jteta", jteta);
		tJet->SetBranchAddress("jtphi", jtphi);
		tJet->SetBranchAddress("jtm",jtmass);
		tJet->SetBranchAddress("pthat",&pthat);
		if (loadGenJet) tGenJet = (TTree*)tFile->Get(genJetTree);
		if (loadGenJet) tGenJet->SetBranchAddress("ngen" ,&ngen);
		if (loadGenJet) tGenJet->SetBranchAddress("genpt", genpt);
		if (loadGenJet) tGenJet->SetBranchAddress("gensubid", gensubid);
		tEvt->SetBranchAddress("hiBin",&bin);
		tEvt->SetBranchAddress("vz",&vz);
		tSkim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
		if(isPbPb) tSkim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
		else tSkim->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
		tJet->AddFriend(tEvt);
		tJet->AddFriend(tSkim);
	};
	TFile *tFile;
	TTree *tJet;
	TTree *tGenJet;
	TTree *tEvt;
	TTree* tSkim;
	float jtpt[1000];
	float rawpt[1000];
	float refpt[1000];
	float jtphi[1000];
	float jteta[1000];
	float jtmass[1000];
	float trackMax[1000];
	float chargedMax[1000];
    float neutralMax[1000];
    float chargedSum[1000];
    float neutralSum[1000];
	float genpt[1000];
	int gensubid[1000];
	float vz;
	float pthat;
	int njets;
	int ngen;
	int bin;     
	int pHBHENoiseFilter;
	int pPAcollisionEventSelectionPA;
	int pcollisionEventSelection;
};


using namespace std;

void RAA_pp_mc(int radius = 7, char *algo = "Pu"){
   
  TStopwatch timer;
  timer.Start();

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  cout<<"Running for Radius = "<< radius<<" and Algorithm is "<<algo<<endl;
 
  const int nbinsPP_pthat = 11;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
  //Long64_t nentries = 500;
  
//  gStyle->SetOptLogx();
//  gStyle->SetOptLogy();
  
  boundariesPP_pthat[0]=15;
  fileNamePP_pthat[0] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT15_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[0]= 0.2034;
//  entries[0] = 71680;  
  
  boundariesPP_pthat[1]=30;
  fileNamePP_pthat[1] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT30_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[1]= 0.01075;
  //entries[1] = 52160;
  
  boundariesPP_pthat[2]=50;
  fileNamePP_pthat[2] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT50_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[2]= 0.001025;
 // entries[2] = 50240;
  
  boundariesPP_pthat[3]=80;
  fileNamePP_pthat[3] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT80_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[3]= 9.8650e-05;
 // entries[3] = 52160;
  
  boundariesPP_pthat[4]=120;
  fileNamePP_pthat[4] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT120_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[4]= 1.1290e-05;
 // entries[4] = 53760;

  boundariesPP_pthat[5] = 170;
  fileNamePP_pthat[5] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT170_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[5]= 1.4650e-06;
  //entries[5] = 53120;
  
  boundariesPP_pthat[6]=220;
  fileNamePP_pthat[6] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT220_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[6]= 2.8370e-07;
 // entries[6] = 54080;
  
  boundariesPP_pthat[7]=280;
  fileNamePP_pthat[7] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT280_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[7]= 5.3230e-08;
 // entries[7] = 53120;
  
  boundariesPP_pthat[8]=370;
  fileNamePP_pthat[8] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT370_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[8]= 5.9340e-09;
  //entries[8] = 52800;
  
  boundariesPP_pthat[9]=460;
  fileNamePP_pthat[9] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT460_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root";
  xsectionPP[9]= 8.1250e-10;
  //entries[9] = 54080;
  
  boundariesPP_pthat[10]=540;
  fileNamePP_pthat[10] = "/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet29_cut1/QCDpT540_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_1GeVcut_badJEC_forest.root ";
  xsectionPP[10]= 1.4670e-10;
  //entries[10] = 53440;  
  
  xsectionPP[11] = 0;
  boundariesPP_pthat[11]=1000;
  
  
  static const int nbins_cent = 6;
  Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  //now we have to multiply by 5, since centrality goes from 0-200. 
  Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
  
  TH1F* hpp_gen = new TH1F("hpp_gen",Form("gen refpt ak%d%sPF",radius,algo),nbins_pt,boundaries_pt);
  TH1F* hpp_raw = new TH1F("hpp_raw",Form("raw pt ak%d%sPF",radius,algo),nbins_pt,boundaries_pt); //histo for rawpt
  TH1F* hpp_eta = new TH1F("hpp_eta",Form("jteta ak%d%sPF",radius,algo),500,-4,4); //histo for eta
  TH1F* hpp_phi = new TH1F("hpp_phi",Form("jtphi ak%d%sPF",radius,algo),20,-1*TMath::Pi(),TMath::Pi());
  TH1F* hpp_reco = new TH1F("hpp_reco",Form("reco jtpt ak%d%sPF",radius,algo),nbins_pt,boundaries_pt);
  TH2F* hpp_matrix = new TH2F("hpp_matrix",Form("matrix refpt rawpt ak%d%sPF",radius,algo),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  //TH2F* hpp_response = new TH2F("hpp_response","response jtpt refpt",1000,0,1000,1000,0,1000);
  TH1F* hpp_mcclosure_data = new TH1F("hpp_mcclosure_data","data for unfolding mc closure test pp",nbins_pt,boundaries_pt);

  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
  TH1F *hCentData = new TH1F("hCentData","",40,0,40);
  TH1F *hCentMC = new TH1F("hCentMC","",40,0,40);
	
  TH1F *hVzData = new TH1F("hVzData","",60,-15,15);
  TH1F *hVzMC = new TH1F("hVzMC","",60,-15,15);
  TH1F *hVzPPData = new TH1F("hVzPPData","",60,-15,15);
  TH1F *hVzPPMC = new TH1F("hVzPPMC","",60,-15,15);
  hCent->Sumw2();
  hCentData->Sumw2();
  hCentMC->Sumw2();
  hVzData->Sumw2();
  hVzMC->Sumw2();
  hVzPPData->Sumw2();
  hVzPPMC->Sumw2();
  

  // Setup jet data branches
    JetData *dataPP[nbinsPP_pthat];
  
  
  
  for (int i=0;i<nbinsPP_pthat;i++) dataPP[i] = new JetData(fileNamePP_pthat[i],Form("ak%dPFJetAnalyzer/t",radius),Form("ak%dPFJetAnalyzer/t",radius),0,0);
	
	
  TH1F *hPtHatPP = new TH1F("hPtHatPP",Form("p_{T} Hat Distribution ak%d%s",radius,algo),nbinsPP_pthat,boundariesPP_pthat);
  TH1F *hPtHatRawPP = new TH1F("hPtHatRawPP",Form("Raw p_{T} Hat Distribution ak%d%s",radius,algo),nbinsPP_pthat,boundariesPP_pthat);
  
 
  cout<<"reading all the pp mc files"<<endl;
  for (int i=0;i<nbinsPP_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
    dataPP[i]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRawPP->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
  
  hPtHatRawPP->Print("base");
  
  // fill pp MC
  for (int i=0;i<nbinsPP_pthat;i++) {
    if (xsectionPP[i]==0) continue;
    //float scale=(xsectionPP[i]-xsectionPP[i+1])/dataPP[i]->tJet->GetEntries(Form("pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1])); 
    double nevts = dataPP[i]->tJet->GetEntries();
	double nents = (1./nevts);
	cout <<"Loading PP pthat"<<boundariesPP_pthat[i]
	 <<" sample, cross section = "<<xsectionPP[i]
	 << Form(" pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1])<<endl;
	cout <<"Number of events(nevts): "<<nevts<<endl;
	
    //cout<<""<<endl;
    for (Long64_t jentry2=0; jentry2<dataPP[i]->tJet->GetEntries();jentry2++) {
      //for (Long64_t jentry2=0; jentry2<10;jentry2++) {
      dataPP[i]->tEvt->GetEntry(jentry2);
      dataPP[i]->tJet->GetEntry(jentry2);
      //dataPP[i]->tGenJet->GetEntry(jentry2);
      //if(dataPP[i]->pthat<boundariesPP_pthat[i] || dataPP[i]->pthat>boundariesPP_pthat[i+1]) continue;
      //if(dataPP[i]->bin<=28) continue;
	  
      int pthatBin = hPtHatPP->FindBin(dataPP[i]->pthat);
	  //cout<< fileNamePP_pthat[i]<<endl; //DEBUGGING PURPOSES ONLY
      float scalepp = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/hPtHatRawPP->GetBinContent(pthatBin);
      if(fabs(dataPP[i]->vz)>15) continue;
      double weight_cent=1;
      double weight_pt=1;
      double weight_vz=1;
      if(!dataPP[i]->pPAcollisionEventSelectionPA || !dataPP[i]->pHBHENoiseFilter) continue;
      
      //weight_vz = fVzPP->Eval(dataPP[i]->vz); //Commented out for MC - OB (7/01/14)
      //if (weight_vz>5||weight_vz<0.5) cout <<dataPP[i]->vz<<" "<<weight_vz<<endl;
      weight_vz = 1; //uncommented for MC - OB (7/01/14)
      hPtHatPP->Fill(dataPP[i]->pthat,scalepp*weight_vz);
	  
      int hasLeadingJet = 0;
      hVzPPMC->Fill(dataPP[i]->vz,scalepp*weight_vz);
      /*
	for (int k= 0; k < dataPP[i]->njets; k++) { 
	if ( dataPP[i]->jteta[k]  > 2. || dataPP[i]->jteta[k] < -2. ) continue;
	if ( dataPP[i]->jtpt[k]>100) {
	hasLeadingJet = 1;
	}
	break;
	
	}
	if (hasLeadingJet == 0) continue;
      */
      for (int k= 0; k < dataPP[i]->njets; k++) {   //JET LOOP
	int subEvt=-1;
	if ( dataPP[i]->refpt[k]  < 30. ) continue;
	
		hpp_eta->Fill(dataPP[i]->jteta[k]);  //jteta fill, BEFORE eta cut!
		hpp_phi->Fill(dataPP[i]->jtphi[k],nents); //normalize by entries, hardcoded!!
		
	
	if ( dataPP[i]->jteta[k]  > 2. || dataPP[i]->jteta[k] < -2. ) continue;
	//if ( dataPP[i]->chargedMax[k]/dataPP[i]->jtpt[k]<0.01) continue; //COMMENTED OUT FOR BAD JEC
	//if ( dataPP[i]->neutralMax[k]/TMath::Max(dataPP[i]->chargedSum[k],dataPP[i]->neutralSum[k]) < 0.975)continue; //COMMENTED OUT FOR BAD JEC
	//if ( dataPP[i]->neu)

	//if (uhist[nbins_cent]->hMeasMatch!=0) {
	//   int ptBinNumber = uhist[nbins_cent]->hMeasMatch->FindBin(dataPP[i]->jtpt[k]);
	//   int ratio = uhist[nbins_cent]->hMeasMatch->GetBinContent(ptBinNumber);
	//if (ratio!=0) weight_pt = 1./ratio;
	//}
					
	//if (!isMC||jentry2<dataPP[i]->tJet->GetEntries()/2.) {
	//if(!isMC){
	//hpp_response->Fill(dataPP[i]->jtpt[k],dataPP[i]->refpt[k],scalepp*weight_vz);

    hpp_matrix->GetXaxis()->SetTitle("Ref p_{T}");
    hpp_matrix->GetXaxis()->CenterTitle();
    hpp_matrix->GetYaxis()->SetTitle("Raw p_{T}");
    hpp_matrix->GetYaxis()->CenterTitle();

    hpp_reco->GetXaxis()->SetTitle("Reco p_{T}");
    hpp_reco->GetXaxis()->CenterTitle();
    hpp_reco->GetYaxis()->SetTitle("Events");
    hpp_reco->GetYaxis()->CenterTitle();
  
    hpp_raw->GetXaxis()->SetTitle("Raw p_{T}");
    hpp_raw->GetXaxis()->CenterTitle();
    hpp_raw->GetYaxis()->SetTitle("Events");
    hpp_raw->GetYaxis()->CenterTitle();
   
    hpp_eta->GetXaxis()->SetTitle("Eta");
    hpp_eta->GetXaxis()->CenterTitle();
  
    hpp_phi->GetXaxis()->SetTitle("Phi");
    hpp_phi->GetXaxis()->CenterTitle();
    hpp_phi->GetYaxis()->SetTitle("Normalized Units");
    hpp_phi->GetYaxis()->CenterTitle();
  
    hpp_gen->GetXaxis()->SetTitle("Gen Ref p_{T}");
    hpp_gen->GetXaxis()->CenterTitle();
    hpp_gen->GetYaxis()->SetTitle("Events");
    hpp_gen->GetYaxis()->CenterTitle();

    hPtHatPP->GetXaxis()->SetTitle("p_{T} Hat");
    hPtHatPP->GetXaxis()->CenterTitle();
    //hPtHatPP->GetYaxis()->SetTitle("Events");
	
	hpp_matrix->Fill(dataPP[i]->refpt[k],dataPP[i]->rawpt[k],scalepp*weight_vz);
	hpp_raw->Fill(dataPP[i]->rawpt[k],scalepp*weight_vz); //rawpt
	hpp_gen->Fill(dataPP[i]->refpt[k],scalepp*weight_vz);   
	hpp_reco->Fill(dataPP[i]->jtpt[k],scalepp*weight_vz);
	
	
	
	//}	  
	if (jentry2>dataPP[i]->tJet->GetEntries()/2.)
	  hpp_mcclosure_data->Fill(dataPP[i]->jtpt[k],scalepp*weight_vz);
	
	//uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);   
	//uhist[nbins_cent]-> hMeas->Fill(dataPP[i]->jtpt[k],scale*weight_vz); 
	//}
      }//njet loop     
    }//nentry loop
  }//ptbins loop
  
   
  TDatime date;

  //declare the output file 
  TFile f(Form("/net/hisrv0001/home/obaron/CMSSW_5_3_16/drawfiles/output/pp_mc_ak%d_%s_%d_chMax_12003cut.root",radius,algo,date.GetDate()),"RECREATE");
  f.cd();
  
  divideBinWidth(hpp_gen);
  divideBinWidth(hpp_reco);
  divideBinWidth(hpp_mcclosure_data);

  TCanvas *c5=new TCanvas();
  c5->cd();
  hpp_eta->Draw();

  TCanvas *c7=new TCanvas();
  c7->cd();
  hpp_phi->GetYaxis()->SetRangeUser(0,2);
  hpp_phi->Draw();
  
  gStyle->SetOptLogy();

  TCanvas *c8=new TCanvas();
  c8->cd();
  hPtHatPP->Draw();
  
  TCanvas *c6=new TCanvas();
  c6->cd();
  hpp_raw->Draw();   
  
  TCanvas *c1=new TCanvas();
  c1->cd();
  hpp_gen->Draw();
   
  TCanvas *c2=new TCanvas();
  c2->cd();
  hpp_reco->Draw();
  
  /*TCanvas *c3=new TCanvas();
  c3.cd();
  hpp_mcclosure_data.Draw();*/
  
  gStyle->SetOptLogx();
  gStyle->SetOptLogz();
  TCanvas *c4=new TCanvas();
  c4->cd();
  hpp_matrix->Draw("colz");
  cout<<nents<<endl;
  
  hpp_gen->Write();
  hpp_gen->Print("base");
  hpp_reco->Write();
  hpp_reco->Print("base");
  hpp_raw->Write();
  hpp_raw->Print("base");
  hpp_eta->Write();
  hpp_eta->Print("base");
  hpp_phi->Write();
  hpp_phi->Print("base");
  hpp_matrix->Write();
  hpp_matrix->Print("base");
 // hpp_mcclosure_data->Write();
 // hpp_mcclosure_data->Print("base");
  hPtHatPP->Write();
  hPtHatPP->Print("base");
  
  //hpp_response->Write();
  //hpp_response->Print("base");
  
  f.Write();
  f.Close();
  
  c1->SaveAs(Form("output/pp_mc_ak%d_%s_%d_refpt_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  c2->SaveAs(Form("output/pp_mc_ak%d_%s_%d_jtpt_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  // c3->SaveAs(Form("output/pp_mc_ak%d_%s_%d_closure_chMax_10993cut.pdf",radius,algo,date.GetDate()),"RECREATE");
  c4->SaveAs(Form("output/pp_mc_ak%d_%s_%d_matrix_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  c5->SaveAs(Form("output/pp_mc_ak%d_%s_%d_eta_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  c6->SaveAs(Form("output/pp_mc_ak%d_%s_%d_rawpt_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  c7->SaveAs(Form("output/pp_mc_ak%d_%s_%d_phi_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  c8->SaveAs(Form("output/pp_mc_ak%d_%s_%d_pthat_chMax_10993cut.png",radius,algo,date.GetDate()),"RECREATE");
  

  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<timer.CpuTime()<<endl;
  cout<<"Real time (min) = "<<timer.RealTime()<<endl;
  
}