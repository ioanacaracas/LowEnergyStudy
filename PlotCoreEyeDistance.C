#include <TProfile.h>
#include <TProfile2D.h>

#include <Math/Polynomial.h>  
#include <Math/Interpolator.h>
#include <TFile.h>
#include <TObject.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <TStyle.h>
#include <Riostream.h>
#include <TMath.h>

#include <iostream>
#include <string>
#include <cmath>

void PlotCoreEyeDistance(){
	

	gStyle->SetOptStat(0);

	const int nLgEBins = 20; //all Erange
  	const int nH1Bins = 18;
  	const int nRcorebins = 200;


 // =======make EyeCoreDistance histos for each H1 Ecal bin
	TH1D* RHistoForEcalBinForH1Bin[nLgEBins][nH1Bins];
  	string nameEprimHistoForEcalBinForH1Bin[nLgEBins][nH1Bins];

	string nameRhist;
	std::stringstream ssRhist;

  for(int i = 0; i < nLgEBins; ++i){
    for(int ih1 = 0; ih1<nH1Bins; ++ih1){
      ssRhist.clear();
      ssRhist.str(std::string());
      ssRhist <<"Rdist_Ecalbin_"<<i<<"H1bin_"<<ih1;
      nameRhist = ssRhist.str();
      // nameEprimHistoForEcalBinForH1Bin[i][ih1] = nameRhist;
      RHistoForEcalBinForH1Bin[i][ih1] = new TH1D(nameRhist.c_str(), "", 200, 0, 100);
    }
  }

//====make Zenith histos for each H1 Ecal bin-====
  	TH1D* ZenithHistoForEcalBinForH1Bin[nLgEBins][nH1Bins];


	string nameZenithhist;
	std::stringstream ssZenithhist;

  for(int i = 0; i < nLgEBins; ++i){
    for(int ih1 = 0; ih1<nH1Bins; ++ih1){
      ssZenithhist.clear();
      ssZenithhist.str(std::string());
      ssZenithhist <<"Zenith_Ecalbin_"<<i<<"H1bin_"<<ih1;
      nameZenithhist = ssZenithhist.str();
      // nameEprimHistoForEcalBinForH1Bin[i][ih1] = nameRhist;
      ZenithHistoForEcalBinForH1Bin[i][ih1] = new TH1D(nameZenithhist.c_str(), "", 70, 110, 180);
    }
  }

	TH2D* ZenithVsREyeCoreHistoForEcalBinForH1Bin[nLgEBins][nH1Bins];


	string nameZenithVsREyeCorehist;
	std::stringstream ssZenithZenithVsREyeCorehist;

  for(int i = 0; i < nLgEBins; ++i){
    for(int ih1 = 0; ih1<nH1Bins; ++ih1){
      ssZenithZenithVsREyeCorehist.clear();
      ssZenithZenithVsREyeCorehist.str(std::string());
      ssZenithZenithVsREyeCorehist <<"ZenithVsREyeCore_Ecalbin_"<<i<<"H1bin_"<<ih1;
      nameZenithVsREyeCorehist = ssZenithZenithVsREyeCorehist.str();
      // nameEprimHistoForEcalBinForH1Bin[i][ih1] = nameRhist;
      ZenithVsREyeCoreHistoForEcalBinForH1Bin[i][ih1] = new TH2D(nameZenithVsREyeCorehist.c_str(), "", 7, 110, 180, 50, 0, 100);
    }
  }

  //=====make Zenith Vs R_EyeCpre for each H1 Ecal bin===

	TFile* InDataPCGFupFile = new TFile("AnalysisPCGFupSelected.root", "READ");

	InDataPCGFupFile->cd();
	TTree* treeDataPCGF;

	 InDataPCGFupFile->GetObject("TPCGFup", treeDataPCGF);
	 //= (TTree*)InDataPCGFupFile->Get("TPCGFup");
 	
 	float LgEMC ;
 	float LgEcalMC;
	float D1;
	float H1;
	double Zenith;
	double GenerateddistanceEyeCore;

	treeDataPCGF->Print();
	


	treeDataPCGF->SetBranchAddress("LgEMC", &LgEMC);
	treeDataPCGF->SetBranchAddress("LgEcalMC", &LgEcalMC);
	treeDataPCGF->SetBranchAddress("D1", &D1);
	treeDataPCGF->SetBranchAddress("H1", &H1);
	treeDataPCGF->SetBranchAddress("Zenith", &Zenith);
	treeDataPCGF->SetBranchAddress("GenerateddistanceEyeCore", &GenerateddistanceEyeCore);
	
	int nEvents = treeDataPCGF->GetEntries();


	TH2D* EcalVsH1EventsHist = new TH2D("EcalVsH1","EcalVsH1", nLgEBins, 16.5, 19, nH1Bins, 0, 9);
	TProfile2D *EcalVsH1VsRcoreEye = new TProfile2D("EcalVsH1VsRcoreEye","EcalVsH1VsRcoreEye",  16.5, 19, 0, 9, 0,100);

	double H1ForHistoTitle;
	double lgEcalForHistoTitle;


	cout<<" event: "<<nEvents<<endl;
	for(int iEv = 0; iEv < nEvents; iEv++){

		treeDataPCGF->GetEntry(iEv);

		//cout<<"H1 = "<<H1<<" D1: "<<D1<<" Zenith: "<< Zenith<<" R" << GenerateddistanceEyeCore<<" Ecal: "<<LgEcalMC<<" Egen: "<< LgEMC<<endl;
		if(LgEcalMC<=18.5){
			EcalVsH1EventsHist->Fill(LgEcalMC, H1/1000);
			EcalVsH1VsRcoreEye->Fill(LgEcalMC, H1/1000, GenerateddistanceEyeCore/1000);
			double EnergyDiff;
      		int EcalBin;

			//=========This works for 20 Ecal bins, with Ecal in 16.5 - 19(for everything > 18.5 content = 0) as in Massimo's DoubleDiff)
                                  //=====EcalBinWidth = 0.125====
			if( LgEcalMC != 18.5 ){
				EnergyDiff = (LgEcalMC - 16.5) * 8;
				EcalBin = (int) EnergyDiff;
			}
			else
				{EcalBin = nLgEBins - 1; }//put E=18.5 in last bin

			//cout<<"Ecal bin: "<<EcalBin<<" e mc: "<<LgEcalMC<<endl;

		//=========This works for 18 H1 bins, with H1 in 0 - 9 (as in Massimo's DoubleDiff)
                                  //=====H1BinWidth = 0.5====
			int H1Bin;

			if(H1 != 9000){
				H1Bin = (int) (H1/1000 * 2);
			}
			else 
				H1Bin = nH1Bins - 1;


			RHistoForEcalBinForH1Bin[EcalBin][H1Bin]->Fill(GenerateddistanceEyeCore/1000);
			ZenithHistoForEcalBinForH1Bin[EcalBin][H1Bin]->Fill(Zenith);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[EcalBin][H1Bin]->Fill(Zenith, GenerateddistanceEyeCore/1000);

						

		}

	}



	TCanvas* EcalVsH1Canvas = new TCanvas("EcalVsH1Canvas", "EcalVsH1Canvas", 200, 10, 700,530);
	EcalVsH1Canvas->cd();
	EcalVsH1Canvas->SetBottomMargin(0.18);
	EcalVsH1Canvas->SetTopMargin(0.07);
	EcalVsH1Canvas->SetFrameLineWidth(2);
	EcalVsH1Canvas->SetLineWidth(2);
	
	EcalVsH1Canvas->SetLeftMargin(0.14);
	EcalVsH1Canvas->SetRightMargin(0.18);

	EcalVsH1EventsHist->GetXaxis()->SetRangeUser(16.5,18.5);
	EcalVsH1EventsHist->SetTitle("EcalVsH1 selected events");
	EcalVsH1EventsHist->GetXaxis()->SetRangeUser(16.5,18.5);
	EcalVsH1EventsHist->GetXaxis()->SetTitle("lg E_{cal} / eV");
	EcalVsH1EventsHist->GetXaxis()->SetNdivisions(9,10,2);
	EcalVsH1EventsHist->GetXaxis()->SetRangeUser(16.5, 18.5);
	EcalVsH1EventsHist->GetYaxis()->SetTitleOffset(0.75);
	EcalVsH1EventsHist->GetYaxis()->SetTitle("H_{1} [km]");
	EcalVsH1EventsHist->GetXaxis()->SetTitleOffset(1.25);
	EcalVsH1EventsHist->GetYaxis()->SetLabelSize(.06);
	EcalVsH1EventsHist->GetYaxis()->SetTitleSize(.065);
	EcalVsH1EventsHist->GetXaxis()->SetLabelSize(.06);
	EcalVsH1EventsHist->GetXaxis()->SetTitleSize(.065);
	EcalVsH1EventsHist->GetXaxis()->SetLabelFont(132);
	EcalVsH1EventsHist->GetYaxis()->SetLabelFont(132);
	EcalVsH1EventsHist->GetXaxis()->SetTitleFont(132);
	EcalVsH1EventsHist->GetYaxis()->SetTitleFont(132);
	EcalVsH1EventsHist->GetXaxis()->CenterTitle();
	EcalVsH1EventsHist->GetYaxis()->CenterTitle();
	
	EcalVsH1EventsHist->GetZaxis()->SetTitle("selected events");
	EcalVsH1EventsHist->GetZaxis()->SetTitleOffset(0.85);
	EcalVsH1EventsHist->GetZaxis()->SetTitleSize(0.065);
	EcalVsH1EventsHist->GetZaxis()->SetLabelSize(0.06);
	EcalVsH1EventsHist->GetZaxis()->SetTitleFont(132);
	EcalVsH1EventsHist->GetZaxis()->SetLabelFont(132);
	EcalVsH1EventsHist->GetZaxis()->CenterTitle();


	EcalVsH1EventsHist->Draw("COLZ");
	cout<<"E bin width: "<<EcalVsH1EventsHist->GetXaxis()->GetBinWidth(1)<<endl;

	TCanvas* EcalVsH1VsRCanvas = new TCanvas("EcalVsH1VsRCanvas", "EcalVsH1VsRCanvas", 800,600);
	//EcalVsH1VsRcoreEye->GetXaxis()->SetLimits(16.5,18.5);
	//EcalVsH1VsRcoreEye->SetMinimum(16.5);

	EcalVsH1VsRcoreEye->Draw("COLZ");

	double nSelectedEvents;
	double Rmean;
	double meanZenith;

	TH2D* EcalVsH1VsRZaxisHist = new TH2D("EcalVsH1VsRZaxisHist","EcalVsH1VsRZaxisHist", nLgEBins, 16.5, 19, nH1Bins, 0, 9);
	TH2D* EcalVsH1VsZenithHist = new TH2D("EcalVsH1VsZenithHist", "EcalVsH1VsZenithHist", nLgEBins, 16.5, 19, nH1Bins, 0, 9);

	for(int iLgeBin = 0; iLgeBin<nLgEBins; iLgeBin++){
    	for(int iH1bin = 0; iH1bin < nH1Bins; iH1bin++){


    		nSelectedEvents = EcalVsH1EventsHist->GetBinContent(iLgeBin+1, iH1bin+1);
    		Rmean = RHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetMean();
    		meanZenith = ZenithHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetMean();

    		cout<<"Ecalbin = "<<iLgeBin<< "  H1 bin = "<<iH1bin
    			<< " H1 = "<< EcalVsH1EventsHist->GetYaxis()->GetBinCenter(iH1bin + 1)
    			<< "Ecal = "<< EcalVsH1EventsHist->GetXaxis()->GetBinCenter(iLgeBin + 1)
    			<<" events = "<<nSelectedEvents<<" <R> "<< Rmean<<endl;

    		//cout<<"Ecalbin = "iLgeBin<< "  H1 bin = "<<iH1bin<< " H1 = "<< EcalVsH1EventsHist->GetYaxis()->GetBinCenter(iH1bin + 1)<<" Ecal = "<<EcalVsH1EventsHist->GetXaxis()->GetBinCenter(iLgeBin + 1)<<endl;
    		if(nSelectedEvents!=0){
    			EcalVsH1VsRZaxisHist->Fill(EcalVsH1EventsHist->GetXaxis()->GetBinCenter(iLgeBin + 1), EcalVsH1EventsHist->GetYaxis()->GetBinCenter(iH1bin + 1),1*Rmean);
    			EcalVsH1VsZenithHist->Fill(EcalVsH1EventsHist->GetXaxis()->GetBinCenter(iLgeBin + 1), EcalVsH1EventsHist->GetYaxis()->GetBinCenter(iH1bin + 1),1*meanZenith);
    		}

    		//set title with values of Ecal and H1

    		H1ForHistoTitle = EcalVsH1EventsHist->GetYaxis()->GetBinCenter(iH1bin + 1);
    		lgEcalForHistoTitle = EcalVsH1EventsHist->GetXaxis()->GetBinCenter(iLgeBin + 1);

    		ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->SetTitle(Form("lg Ecal = %1.2f / eV, H1 =%1.2f km ", lgEcalForHistoTitle, H1ForHistoTitle));

    	}
    }


	
	TCanvas* EcalVsH1VsRonZCanvas = new TCanvas("EcalVsH1VsRonZCanvas", "EcalVsH1VsRonZCanvas", 200, 10, 700,530);
	EcalVsH1VsRonZCanvas->cd();
	EcalVsH1VsRonZCanvas->SetBottomMargin(0.18);
	EcalVsH1VsRonZCanvas->SetTopMargin(0.07);
	EcalVsH1VsRonZCanvas->SetFrameLineWidth(2);
	EcalVsH1VsRonZCanvas->SetLineWidth(2);
	
	EcalVsH1VsRonZCanvas->SetLeftMargin(0.14);
	EcalVsH1VsRonZCanvas->SetRightMargin(0.18);

	EcalVsH1VsRZaxisHist->SetTitle("EcalVsH1VsR_EyeCore");
	EcalVsH1VsRZaxisHist->GetXaxis()->SetRangeUser(16.5,18.5);
	EcalVsH1VsRZaxisHist->GetXaxis()->SetTitle("lg E_{cal} / eV");
	EcalVsH1VsRZaxisHist->GetXaxis()->SetNdivisions(9,10,2);
	EcalVsH1VsRZaxisHist->GetXaxis()->SetRangeUser(16.5, 18.5);
	EcalVsH1VsRZaxisHist->GetYaxis()->SetTitleOffset(0.75);
	EcalVsH1VsRZaxisHist->GetYaxis()->SetTitle("H_{1} [km]");
	EcalVsH1VsRZaxisHist->GetXaxis()->SetTitleOffset(1.25);
	EcalVsH1VsRZaxisHist->GetYaxis()->SetLabelSize(.06);
	EcalVsH1VsRZaxisHist->GetYaxis()->SetTitleSize(.065);
	EcalVsH1VsRZaxisHist->GetXaxis()->SetLabelSize(.06);
	EcalVsH1VsRZaxisHist->GetXaxis()->SetTitleSize(.065);
	EcalVsH1VsRZaxisHist->GetXaxis()->SetLabelFont(132);
	EcalVsH1VsRZaxisHist->GetYaxis()->SetLabelFont(132);
	EcalVsH1VsRZaxisHist->GetXaxis()->SetTitleFont(132);
	EcalVsH1VsRZaxisHist->GetYaxis()->SetTitleFont(132);
	EcalVsH1VsRZaxisHist->GetXaxis()->CenterTitle();
	EcalVsH1VsRZaxisHist->GetYaxis()->CenterTitle();
	
	EcalVsH1VsRZaxisHist->GetZaxis()->SetTitle("R_EyeCore [km]");
	EcalVsH1VsRZaxisHist->GetZaxis()->SetTitleOffset(0.85);
	EcalVsH1VsRZaxisHist->GetZaxis()->SetTitleSize(0.065);
	EcalVsH1VsRZaxisHist->GetZaxis()->SetLabelSize(0.06);
	EcalVsH1VsRZaxisHist->GetZaxis()->SetTitleFont(132);
	EcalVsH1VsRZaxisHist->GetZaxis()->SetLabelFont(132);
	EcalVsH1VsRZaxisHist->GetZaxis()->CenterTitle();
	// EcalVsH1VsRZaxisHist->GetZaxis()->SetRangeUser(110, 180);
	EcalVsH1VsRZaxisHist->Draw("COLZ");
	
	TCanvas* EcalVsH1VsZenithCanvas = new TCanvas("EcalVsH1VsZenithCanvas", "EcalVsH1VsZenithCanvas",200, 10, 700,530);
	EcalVsH1VsZenithCanvas->cd();
	EcalVsH1VsZenithCanvas->SetBottomMargin(0.18);
	EcalVsH1VsZenithCanvas->SetTopMargin(0.07);
	EcalVsH1VsZenithCanvas->SetFrameLineWidth(2);
	EcalVsH1VsZenithCanvas->SetLineWidth(2);
	
	EcalVsH1VsZenithCanvas->SetLeftMargin(0.14);
	EcalVsH1VsZenithCanvas->SetRightMargin(0.18);


	EcalVsH1VsZenithHist->GetXaxis()->SetTitle("lg E_{cal} / eV");
	EcalVsH1VsZenithHist->GetXaxis()->SetNdivisions(9,10,2);
	EcalVsH1VsZenithHist->GetXaxis()->SetRangeUser(16.5, 18.5);
	EcalVsH1VsZenithHist->GetYaxis()->SetTitleOffset(0.75);
	EcalVsH1VsZenithHist->GetYaxis()->SetTitle("H_{1} [km]");
	EcalVsH1VsZenithHist->GetXaxis()->SetTitleOffset(1.25);
	EcalVsH1VsZenithHist->GetYaxis()->SetLabelSize(.06);
	EcalVsH1VsZenithHist->GetYaxis()->SetTitleSize(.065);
	EcalVsH1VsZenithHist->GetXaxis()->SetLabelSize(.06);
	EcalVsH1VsZenithHist->GetXaxis()->SetTitleSize(.065);
	EcalVsH1VsZenithHist->GetXaxis()->SetLabelFont(132);
	EcalVsH1VsZenithHist->GetYaxis()->SetLabelFont(132);
	EcalVsH1VsZenithHist->GetXaxis()->SetTitleFont(132);
	EcalVsH1VsZenithHist->GetYaxis()->SetTitleFont(132);
	EcalVsH1VsZenithHist->GetXaxis()->CenterTitle();
	EcalVsH1VsZenithHist->GetYaxis()->CenterTitle();
	cout<<"it should work"<<endl;
	EcalVsH1VsZenithHist->GetZaxis()->SetTitle("#theta [#circ]");
	EcalVsH1VsZenithHist->GetZaxis()->SetTitleOffset(0.85);
	EcalVsH1VsZenithHist->GetZaxis()->SetTitleSize(0.065);
	EcalVsH1VsZenithHist->GetZaxis()->SetLabelSize(0.06);
	EcalVsH1VsZenithHist->GetZaxis()->SetTitleFont(132);
	EcalVsH1VsZenithHist->GetZaxis()->SetLabelFont(132);
	EcalVsH1VsZenithHist->GetZaxis()->CenterTitle();
	EcalVsH1VsZenithHist->GetZaxis()->SetRangeUser(110, 180);



	// gPad->Update();
	// TPaletteAxis *palette=(TPaletteAxis*)EcalVsH1VsZenithHist->FindObject("palette");
	
	// palette->SetX1NDC (0.859);
	// palette->SetX2NDC (0.889);


	// EcalVsH1VsZenithHist->GetYaxis()->SetTitleOffset(1);
	// EcalVsH1VsZenithHist->GetXaxis()->SetTitleOffset(1.2);
	// EcalVsH1VsZenithHist->GetYaxis()->SetLabelSize(.06);
	// EcalVsH1VsZenithHist->GetYaxis()->SetTitleSize(.065);
	// EcalVsH1VsZenithHist->GetXaxis()->SetLabelSize(.06);
	// EcalVsH1VsZenithHist->GetXaxis()->SetTitleSize(.065);
	// EcalVsH1VsZenithHist->GetXaxis()->SetLabelFont(132);
	// EcalVsH1VsZenithHist->GetYaxis()->SetLabelFont(132);
	// EcalVsH1VsZenithHist->GetXaxis()->SetTitleFont(132);
	// EcalVsH1VsZenithHist->GetYaxis()->SetTitleFont(132);
	// EcalVsH1VsZenithHist->GetYaxis()->SetTitle("H_{1} [km]");
	// EcalVsH1VsZenithHist->GetXaxis()->SetTitle("E_{sh}");
	// EcalVsH1VsZenithHist->GetXaxis()->CenterTitle();
	// EcalVsH1VsZenithHist->GetYaxis()->CenterTitle();
	
	// EcalVsH1VsZenithHist->GetXaxis()->SetRangeUser(16.5,18.5);
	// EcalVsH1VsZenithHist->GetZaxis()->SetTitle("#theta [#circ]");
	// EcalVsH1VsZenithHist->GetZaxis()->SetTitleOffset();
	// EcalVsH1VsZenithHist->GetZaxis()->SetTitleSize(0.065);
	// EcalVsH1VsZenithHist->GetZaxis()->SetLabelSize(0.06);
	// EcalVsH1VsZenithHist->GetZaxis()->SetTitleFont(132);
	// EcalVsH1VsZenithHist->GetZaxis()->SetLabelFont(132);
	// EcalVsH1VsZenithHist->GetZaxis()->CenterTitle();
	
	EcalVsH1VsZenithHist->Draw("COLZtext");
	

	TFile* FileRhistos = new TFile("FileRCoreEyeHistos.root", "RECREATE");
	for(int iLgeBin = 0; iLgeBin<nLgEBins; iLgeBin++){
    	for(int iH1bin = 0; iH1bin < nH1Bins; iH1bin++){
    		RHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->Write();
			ZenithHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->Write();
			
    	}
	}
	TFile* FileRhistos = new TFile("FileZenithVsRCoreEyeHistosFixedEcalH1Bin.root", "RECREATE");
	
	for(int iLgeBin = 0; iLgeBin<nLgEBins; iLgeBin++){
    	for(int iH1bin = 0; iH1bin < nH1Bins; iH1bin++){
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->SetTitle("#theta [#circ]");
			
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetYaxis()->SetTitle("R_EyeCore [km]");
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->SetTitleOffset(1.25);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetYaxis()->SetLabelSize(.06);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetYaxis()->SetTitleSize(.065);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->SetLabelSize(.06);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->SetTitleSize(.065);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->SetLabelFont(132);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetYaxis()->SetLabelFont(132);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->SetTitleFont(132);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetYaxis()->SetTitleFont(132);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetXaxis()->CenterTitle();
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetYaxis()->CenterTitle();

			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetZaxis()->SetTitle("events");
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetZaxis()->SetTitleOffset(0.85);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetZaxis()->SetTitleSize(0.065);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetZaxis()->SetLabelSize(0.06);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetZaxis()->SetTitleFont(132);
			ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetZaxis()->SetLabelFont(132);

			if(ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->GetEntries()!=0)
				ZenithVsREyeCoreHistoForEcalBinForH1Bin[iLgeBin][iH1bin]->Write();
		}
	}


	treeDataPCGF->ResetBranchAddresses();

}