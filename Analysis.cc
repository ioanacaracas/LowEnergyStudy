#include <algorithm>

#include <RecEvent.h>
#include <RecEventFile.h>
#include <DetectorGeometry.h>
#include <FileInfo.h>
#include <EventInfo.h>
#include <FDSelection.h>
#include <UtilityFunctions.h>
#include <Shower.h>


#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <THistPainter.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TF1.h>
#include <TFormula.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TVector3.h>
#include <TVector.h>
#include "TROOT.h"

#include <vector>
#include <math.h>
#include <stdio.h> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include "Analysis.h"

#include "CompADST.h"

using namespace std;


int main(int argc, char** argv) {

  int nOptions = getOptions(argc, argv);

  if (nOptions < 0 || argc - nOptions != 1) {
    Usage(argv[0]);
    return 1;
  }
  static const string inFileName = argv[nOptions];

  TFile* outFile = new TFile("TestAnalysisPCGFupSelected.root", "RECREATE");

  const int UniqueEvents = ReadData(inFileName);
  if (UniqueEvents == 0) {
    cout << " 0 events read in. This is wrong..." << endl;
    return 0;
  }
  cout << " " << UniqueEvents
       << " unique events input file. Starting Analysis..." << endl;

  //const bool Analysis();

  outFile->Write();
  outFile->Close();

  // TFile* outFileWithTree = new TFile("AnalysisPCGFupSelected.root", "RECREATE");
  
  // outFileWithTree->Write();
  // outFileWithTree->Close();


  return 1;
}



int ReadData(string inFileName) {

	const double pi = 3.14159265358979;

	unsigned int nUnique = 0;

	EventInfo eventInfo;
	RecEventFile* inDataFile;
	RecEvent* recEvent;

	inDataFile = new RecEventFile(inFileName.c_str());
	recEvent = new RecEvent();
	DetectorGeometry* const theGeometry = new DetectorGeometry();
	  
	cout << " opening file " << inFileName;


	inDataFile->SetBranchStatus("event.fSDEvent*", 0);
	inDataFile->SetBranchStatus("event.fMDEvent*", 0);
	inDataFile->SetBranchStatus("event.fRdEvent*", 0);



	inDataFile->SetBuffers(&(recEvent));
	inDataFile->ReadDetectorGeometry(*theGeometry);

	unsigned int nEv = inDataFile->GetNEvents();
	int combinedCount = 0;

	cout << " reading " << nEv << " events... " << endl;




	//FileParamXmaxVsChi0.open ("XmaxVsChi0.txt");
	TH1F* GeneratedCorePositionHisto = new TH1F("GeneratedCorePositionHisto", "GeneratedCorePositionHisto", 1000, 0, 100);
	TH1F* TriggeredCorePositionHisto = new TH1F("TriggeredCorePositionHisto", "TriggeredCorePositionHisto", 1000, 0, 100);
	TH1F* TriggeredCorePositionCoihuecoHisto = new TH1F("TriggeredCorePositionCoihuecoHisto", "TriggeredCorePositionCoihuecoHisto", 1000, 0, 100);
	TH1F* TriggeredCorePositionHEATHisto = new TH1F("TriggeredCorePositionHEATHisto", "TriggeredCorePositionHEATHisto", 1000, 0, 100);

	TVector3 EyeCoord;
	TVector3 DistanceEyeCoreVector;

	TVector3 Eye4Coord;
	Eye4Coord.SetX(-31895.8);
	Eye4Coord.SetY(15026.1);
	Eye4Coord.SetZ(214.902);

	TVector3 Eye5Coord;
	Eye5Coord.SetX(-31741.1);
	Eye5Coord.SetY(15095.6);
	Eye5Coord.SetZ(210.548);

	TVector3 GeneratedDistanceEye4CoreVector;





	vector<double> VectorGenerateddistanceEyeCore; //vector to hold the value of distance eye to core for all generated events
	vector<double> VectorTriggereddistanceEyeCore; //vector to hold the value of distance eye to core for all triggered events
	vector<double> VectorRp;
	vector<double> VectorlgEMC;
	vector<double> VectorlgERec;

	TTree *TPCGFup = new TTree("TPCGFup", "myTree");
 

	float LgEMC;
	float LgEcalMC;
	float D1;
	float H1;
	double Zenith;
	double GenerateddistanceEyeCore;
	float XMaxMC;
	float XMaxGH;

	TPCGFup->Branch("LgEMC", &LgEMC, "LgEMC/F");
	TPCGFup->Branch("LgEcalMC", &LgEcalMC, "LgEcalMC/F");
	TPCGFup->Branch("D1", &D1, "D1/F");
	TPCGFup->Branch("H1", &H1, "H1/F");
	TPCGFup->Branch("Zenith", &Zenith, "Zenith/D");
	TPCGFup->Branch("GenerateddistanceEyeCore", &GenerateddistanceEyeCore, "GenerateddistanceEyeCore/D");
	TPCGFup->Branch("XMaxMC", &XMaxMC, "XMaxMC/F");
	TPCGFup->Branch("XMaxGH", &XMaxGH, "XMaxGH/F");


	for (unsigned int i = 0; i < nEv; i++) {
	    inDataFile->GetEventInfo(i, &eventInfo);
	    inDataFile->ReadEvent(i);

	    const std::vector<FDEvent>& fdEvents = recEvent->GetFDEvents();
	    const GenShower& genEvent = recEvent->GetGenShower(); //MC shower

		GeneratedDistanceEye4CoreVector.SetX(Eye4Coord.X()-genEvent.GetCoreAtAltitudeSiteCS(Eye4Coord.Z()).X());
		GeneratedDistanceEye4CoreVector.SetY(Eye4Coord.Y()-genEvent.GetCoreAtAltitudeSiteCS(Eye4Coord.Z()).Y());
		GeneratedDistanceEye4CoreVector.SetZ(Eye4Coord.Z()-genEvent.GetCoreAtAltitudeSiteCS(Eye4Coord.Z()).Z());

		// cout<< "check core altitude?"<<genEvent.GetCoreAtAltitudeSiteCS(Eye4Coord.Z()).Y()<<endl;
		//GenerateddistanceEyeCore = GeneratedDistanceEye4CoreVector.Mag();


		VectorGenerateddistanceEyeCore.push_back(GenerateddistanceEyeCore/1000);

		    if (recEvent->GetFDEvents().size() > 1) {
		    	    combinedCount++;
		      //combineFDs(recEvent, 0);  //i need to use each eye separately, that's why it;s commented
      		}

//cout<<"i ev= "<<i<<endl;
    nUnique++;    


    	int nEyes = recEvent->GetFDEvents().size();

    
    //Read data
    	LgEMC = log10(genEvent.GetEnergy());
    	LgEcalMC = log10(genEvent.GetElecEnergy());
    	D1 = genEvent.GetD1();
    	Zenith = genEvent.GetZenith() * 180 / M_PI;
    	H1 = -D1*cos(Zenith*M_PI/180.);

    	// cout<<"d1 = "<<D1<<" H1 = "<<H1<<endl;

    	

    	VectorlgEMC.push_back(LgEMC);

	   	XMaxMC = genEvent.GetXmaxInterpolated();
	    float dEdXmaxMC = genEvent.GetdEdXmaxInterpolated();
	    XMaxGH = genEvent.GetXmaxGaisserHillas();
	    float NMaxGH = genEvent.GetNmaxGaisserHillas();
	    double X1gen = genEvent.GetX1(); //the x1 from the tree
	    
	    double Azimuth = genEvent.GetAzimuth();
	    float Dec = genEvent.GetDeclination() * 180 / M_PI;


		int eventID = 0; //it will be 0 for the generated events which have no signal in FD
		int eventIDHEAT = 0;
		float LgE = -1; //lgErec = -1 for the generated events which were not triggered
		float LgEHEAT = -1;
		double Rp = -1;
		double RpHEAT = -1;
		double DistanceFdForVector = -1; //generated distance eye core that will be saved in the comp of a tvectorD for further use. When -1 -> Event not triggered
		double DistanceFdForVectorHEAT = -1;
		double distanceRec = -1; //reconstructed distance eye core, when = -1, event iss nto triggered
		double distanceRecHEAT = -1; //reconstructed distance eye core, when = -1, event iss nto triggered

	    const std::vector<double>& slantdepthMC = genEvent.GetDepth();
	    const std::vector<double>& energydepositMC = genEvent.GetEnergyDeposit();
    


	    for (vector<FDEvent>::const_iterator iEye = fdEvents.begin();
	        iEye != fdEvents.end(); ++iEye) {
			const FDEvent& eye = *iEye;

	    	int eyeID = eye.GetEyeId();
	    	if(genEvent.GetDistanceOfShowerMaximum() >= 0)
	    	cout<<"dist of shower max from core: "<<genEvent.GetDistanceOfShowerMaximum()<<endl;
	    // if(eyeID==5)
	    //cout<<"telescope: "<<theGeometry->GetEye(eyeID).GetTelescope(1).GetFADCBinning()<<endl;
	 //for now I only care about Coihueco!!
		    	if(eyeID==7){ 
		    		eventID = eye.GetEventId();
		    		if(eventID==79 )
		   			cout<<"ev ID: "<<eventID<<endl;
		    		//Eye Coordinates
				    EyeCoord.SetX(theGeometry->GetEye(eyeID).GetEyePos().X());
				    EyeCoord.SetY(theGeometry->GetEye(eyeID).GetEyePos().Y());
				    EyeCoord.SetZ(theGeometry->GetEye(eyeID).GetEyePos().Z());




				    // cout<<" eye x = "<<theGeometry->GetEye(eyeID).GetEyePos().X()
				    // 	<<" eye Y = "<<theGeometry->GetEye(eyeID).GetEyePos().Y()
				    // 	<<" eye Z = "<<theGeometry->GetEye(eyeID).GetEyePos().Z()<<endl;



				    double EyeAltitude = theGeometry->GetEye(eyeID).GetEyePos().Z();

				   

				    // cout<<" Core X = "<<genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).X()
				    // 	<<" Core Y = "<<genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).Y()
				    // 	<<" Core Z = "<<genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).Z()<<endl;


				    DistanceEyeCoreVector.SetX(EyeCoord.X()-genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).X());
				    DistanceEyeCoreVector.SetY(EyeCoord.Y()-genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).Y());
				    DistanceEyeCoreVector.SetZ(EyeCoord.Z()-genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).Z());
				    
				    GenerateddistanceEyeCore = DistanceEyeCoreVector.Mag();
				    double distanceEyeCore = GenerateddistanceEyeCore;
				    TPCGFup->Fill();

				    VectorTriggereddistanceEyeCore.push_back(distanceEyeCore/1000);
				    DistanceFdForVector = distanceEyeCore/1000; // generated distance eye core in km (for triggered events)
				  
					LgE = log10(eye.GetFdRecShower().GetEnergy()); //reconstructed log10E


				    VectorlgERec.push_back(LgE);
					double GpsSec = eye.GetGPSSecond();

					distanceRec = eye.GetFdRecGeometry().GetCoreEyeDistance(); //reconstructed distance eye core in km
				    double distance = eye.GetGenGeometry().GetCoreEyeDistance();



				    Rp = eye.GetFdRecGeometry().GetRp();
				    VectorRp.push_back(Rp/1000);

					double Chi0 = eye.GetFdRecGeometry().GetChi0();

					double NPix = eye.GetFdRecGeometry().GetSDPNdF(); //nr of triggered pixels
					double NPixGeo = eye.GetFdRecGeometry().GetNTimeFitPixels();
					float XMaxFD = eye.GetFdRecShower().GetXmax();  

				    double nSDPfitPixel = eye.GetFdRecPixel().GetNumberOfSDPFitPixels(); //total nr of triggered pixels for each event
				    double nTriggeredPixel = eye.GetFdRecPixel().GetNumberOfTriggeredPixels();
				    double nTotalPixel = eye.GetFdRecPixel().GetNumberOfPixels();
				    //cout<<"ntrigpixels = "<<nrtrigpixel<<endl;

				    vector<double>pixelT;
				    vector<double>pixelidvar;
				    vector<double>pixelChivar;
				    
				    //;

				    const std::vector<UShort_t>& pixelID = eye.GetFdRecPixel().GetID(); 
				    const std::vector<double>& pixelTime = eye.GetFdRecPixel().GetTime();
				    const std::vector<EPixelStatus>& pixStatus = eye.GetFdRecPixel().GetStatus();
				    const std::vector<double>& ChiPixel = eye.GetFdRecPixel().GetChi();
				    double TriggtelID;
				    double TriggEyeId;
				    double countnSDPFitPixels = 0;
				    const int nTelescopes = 27;
				    double nPixelsAtTelescop = 0;
				    double tempTrigTel;
				    


				    vector<double> TelIdTrigEventVector;


				    if(eventID==79 && nSDPfitPixel>=40){
				    	cout<<" eye: "<<theGeometry->GetEye(eyeID).GetEyeName()<<" Total triggered pixels: " << nSDPfitPixel<<" Triggered telescopes: "<<endl;

					    for(int ipixit = 0 ; ipixit<pixStatus.size(); ipixit++){
					    	if(pixStatus[ipixit] >= 3){ //get pixels in nSDPFit 
					    	
						    		// cout<<"i = "<< i<<" tel id from pixel i: "<<eye.GetFdRecPixel().GetTelescopeId(i)<<" status pix: "<< pixStatus[i] <<endl;
						    		TelIdTrigEventVector.push_back(eye.GetFdRecPixel().GetTelescopeId(ipixit));

						    	countnSDPFitPixels+=1;
					    	}
					    }
					    cout<<"sum sdp fit pixesl: "<<countnSDPFitPixels<<endl;
				    
					}
					countnSDPFitPixels = 0;
					int itOverTelescopesTrig = 0;
					int TotTelescopesTrigPerEvent = 0;
					int TotTelescopesTrigPerEventFinal;
					


					for(int itVector = 0; itVector <TelIdTrigEventVector.size(); itVector++ ){
						
						
						if(TelIdTrigEventVector[itVector]!= TelIdTrigEventVector[itVector+1]){
							
								TotTelescopesTrigPerEvent+=1;
						}
							if(itVector == TelIdTrigEventVector.size()-1){
								//cout<<"tot telescopes: "<<TotTelescopesTrigPerEvent<<endl;
								TotTelescopesTrigPerEventFinal = TotTelescopesTrigPerEvent;
							}
								
						
					}

					if(TelIdTrigEventVector.size()!=0){
						cout<<"=============tot telescopes: "<<TotTelescopesTrigPerEventFinal<<endl;
						//TotalNrTelescopesPerEvent[i] = TotTelescopesTrigPerEvent;
						double nPixelsAtTelescopArray[TotTelescopesTrigPerEvent];
						double nEyeTrigArray[TotTelescopesTrigPerEvent];

					

						for(int itVector = 0; itVector <TelIdTrigEventVector.size(); itVector++ ){
							TriggtelID = TelIdTrigEventVector[itVector];
							

							nPixelsAtTelescop += 1;

							if(TelIdTrigEventVector[itVector]!= TelIdTrigEventVector[itVector+1]){

								
								
								itOverTelescopesTrig+=1;
													
							

								if(TriggtelID <=6 ){//LL eye Id = 1
			    					TriggEyeId = 1;
			    					cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
			    						<<" telescope nr: "<<TriggtelID
			    						<<" pixels trigg in eye: "<<nPixelsAtTelescop<<" array val pix: "
			    						<<" it: "<<itOverTelescopesTrig	<<endl;
			    				}
			    				else if(TriggtelID<=13){// LM eye ID = 2
			    					TriggEyeId = 2;
			    					cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
			    						<<" telescope nr: "<<TriggtelID
			    						<<" pixels trigg in eye: "<<nPixelsAtTelescop 
			    						<<" it: "<<itOverTelescopesTrig	<<endl;
			    				}
			    				else if(TriggtelID<=18){// LA eye ID = 3
			    					TriggEyeId = 3;
			    					cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
			    						<<" telescope nr: "<<TriggtelID
			    						<<" pixels trigg in eye: "<<nPixelsAtTelescop
			    						<<" it: "<<itOverTelescopesTrig	<<endl;
			    				}
			    				else if(TriggtelID<=24){//CO eye ID =4
			    					TriggEyeId = 4;
			    					cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
			    						<<" telescope nr: "<<TriggtelID
			    						<<" pixels trigg in eye: "<<nPixelsAtTelescop <<" array val pix: "
			    						<<" it: "<<itOverTelescopesTrig	<<endl;
			    						
			    				}
			    				else{//HEAT eye ID =5
			    					TriggEyeId = 5;
			    					cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
			    						<<" telescope nr: "<<TriggtelID
			    						<<" pixels trigg in eye: "<<nPixelsAtTelescop <<" array val pix: "
			    						<<" it: "<<itOverTelescopesTrig	<<endl;
			    				}
									
			    					
			    						nPixelsAtTelescopArray[itOverTelescopesTrig-1] = nPixelsAtTelescop;
			    						nEyeTrigArray[itOverTelescopesTrig-1] = TriggEyeId;
			    					
			    				
								nPixelsAtTelescop = 0;
									
							}



						}

						for( int itrigtelit = 0; itrigtelit <TotTelescopesTrigPerEventFinal; itrigtelit++ ){
							cout<<"it: "<<itrigtelit<<" TrigPixels: "<<nPixelsAtTelescopArray[itrigtelit]<<" at Eye = " << nEyeTrigArray[itrigtelit]<<endl;
							if(itrigtelit != (TotTelescopesTrigPerEventFinal-1) && nEyeTrigArray[itrigtelit] != nEyeTrigArray[itrigtelit+1])
								cout<<"it: "<<itrigtelit<<"===different eyes trigger=== TrigPixels: "<<nPixelsAtTelescopArray[itrigtelit]<<" at Eye = " << nEyeTrigArray[itrigtelit]<<endl;

						}

						cout<<"max SPF pixels: "<<*max_element(nPixelsAtTelescopArray, nPixelsAtTelescopArray+TotTelescopesTrigPerEventFinal)<<endl;

						//nPixelsAtTelescopArray.clear();
						// delete nPixelsAtTelescopArray;

					}
					
					TelIdTrigEventVector.clear();
				    //cout<<"n sdp fit pixels "<<nSDPfitPixel<<" pix stauts size: "<<pixStatus.size()<<endl;

				    //cout<<"n sdp fit pixels "<<nrtrigpixel<<" n pixels: "<<nTotalPixel<<" trig pixels: "<< eye.GetFdRecPixel().GetNumberOfTriggeredPixels()<<endl;
				    // " telescope for last pixel : "<< eye.GetFdRecPixel().GetTelescopeId(nrtrigpixel)<<endl;
				    

				    //for(int iPixelSDPfit = 0; iPixelSDPfit < nSDPfitPixel; iPixelSDPfit++){
				    	//telID = eye.GetFdRecPixel().GetTelescopeId(iPixelSDPfit);
				    	//if(eventID==79 && nSDPfitPixel>=40){

					    	// cout<<" pixel nr: "<< iPixelSDPfit<<" telescope id: "<<eye.GetFdRecPixel().GetTelescopeId(iPixelSDPfit)<<" telIdcheck: "<< telID
					    	// 	// <<" nr of mirrors: "<<theGeometry->GetEye(eyeID).GetNumberOfMirrors()
					    	// 	<<" pixel id: "<<pixelID[iPixelSDPfit] <<endl;

					    	//cout<<" eye: "<<theGeometry->GetEye(eyeID).GetEyeName()<<" Total triggered pixels: " << nSDPfitPixel<<" Triggered telescopes: "<<endl;
					    		// for(int itel = 0; itel <=theGeometry->GetEye(eyeID).GetNumberOfMirrors(); ++itel){
					    		// 	if(eye.MirrorIsInEvent(itel+1)){
					    				
					    		// 		TriggtelID = itel + 1;
					    		// 		cout<<" triggtelID: "<< TriggtelID<<endl;

					    		// 		if(TriggtelID <=6 ){//LL eye Id = 1
					    		// 			TriggEyeId = 1;
					    		// 			cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
					    		// 				<<" pixels trigg in eye: "<<eye.GetFdRecPixel().GetNumberOfSDPFitPixels() <<endl;
					    		// 		}
					    		// 		else if(TriggtelID<=13){// LM eye ID = 2
					    		// 			TriggEyeId = 2;
					    		// 			cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
					    		// 				<<" pixels trigg in eye: "<<eye.GetFdRecPixel().GetNumberOfSDPFitPixels() <<endl;
					    		// 		}
					    		// 		else if(TriggtelID<=18){// LA eye ID = 3
					    		// 			TriggEyeId = 3;
					    		// 			cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
					    		// 				<<" pixels trigg in eye: "<<eye.GetFdRecPixel().GetNumberOfSDPFitPixels() <<endl;
					    		// 		}
					    		// 		else if(TriggtelID<=24){//CO eye ID =4
					    		// 			TriggEyeId = 4;
					    		// 			cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
					    		// 				<<" pixels trigg in eye: "<<eye.GetFdRecPixel().GetNumberOfSDPFitPixels() <<endl;
					    		// 		}
					    		// 		else{//HEAT eye ID =5
					    		// 			TriggEyeId = 5;
					    		// 			cout<<"eye triggered: "<<theGeometry->GetEye(TriggEyeId).GetEyeName()
					    		// 				<<" pixels trigg in eye: "<<eye.GetFdRecPixel().GetNumberOfSDPFitPixels() <<endl;
					    		// 		}



					    		// 	}
					    		// }

					    		//cout<<"last tel Id of the event: "<<TriggtelID<< " last Eye id: "<<TriggEyeId<<endl;

					    		
				    	//}
				    //}



					// for(int k=0; k<pixelTime.size(); k++){ 
					// 	if(pixelTime[k]!=0 && pixStatus[k]==4){ //we only keep the pixel part of time fit:pixStatus[i]==4
						   
					// 		double testtimevar=pixelTime[k];
					// 	    pixelT.push_back(testtimevar); //vector with time of pixels
					// 	    double testpixelvar = pixelID[k];
					// 	    // cout<< "tel id: "<< eye.GetFdRecPixel().GetTelescopeId(k)<<endl;
					// 	    // cout<<"pixel: "<<pixelID[k]<<endl;
					// 	    pixelidvar.push_back(testpixelvar); //vector with id of pixel
					// 	    double pixelchivar = ChiPixel[k];
					// 	   pixelChivar.push_back(180-pixelchivar*180/3.14); //vector with Chi of pixels		

					// 	   //cout<<" telescope: "<<telID[k]<<endl;
						   		   
						  
					// 	//cout<<"k.."<<k<<" id: "<<pixelID[k]<<" time: "<<pixelTime[k] <<" chi:" <<180-ChiPixel[k]*180/3.14<<" eye: "
					// 	//<<eyeID<<"size: "<<pixelTime.size()<<endl;
					// 	}
					// }
				          



				//cout<<NPix<<endl;
				    // data.xmaxFD.push_back(XMaxFD);
				    // //data.xmaxMC.push_back(XMaxMC);
				    // data.dedxmaxMC.push_back(dEdXmaxMC);
				    // data.xmaxGH.push_back(XMaxGH);
				    // data.nmaxGH.push_back(NMaxGH);
				    // //data.x1gen.push_back(X1gen);
				      
				    // data.SlantDepthMC.push_back(slantdepthMC);
				    // data.EnergyDepositMC.push_back(energydepositMC);
				      
				    // data.dec.push_back(Dec);
				      
				    // data.Chi0RecGeom.push_back(Chi0);
				    // data.distFD.push_back(distance);
				    // data.EyeID.push_back(eyeID);
				    // data.TimePixels.push_back(pixelT);
				    // data.IdPixel.push_back(pixelidvar);
				    // data.CHIPixel.push_back(pixelChivar);

				    	
		    
				}
			  //   else if(eyeID==5){ //HEAT
			  //   	eventIDHEAT = eye.GetEventId();
			  //   	// cout << "eye: "<<eyeID<<" x pos: "<<theGeometry->GetEye(eyeID).GetEyePos().X()<<"  "<<theGeometry->GetEye(eyeID).GetEyePos().Y()<<"  "<<theGeometry->GetEye(eyeID).GetEyePos().Z()<<endl;

			  //       //Eye Coordinates
			  //       EyeCoord.SetX(theGeometry->GetEye(eyeID).GetEyePos().X());
			  //       EyeCoord.SetY(theGeometry->GetEye(eyeID).GetEyePos().Y());
			  //       EyeCoord.SetZ(theGeometry->GetEye(eyeID).GetEyePos().Z());

			  //       double EyeAltitude = theGeometry->GetEye(eyeID).GetEyePos().Z();


			  //       DistanceEyeCoreVector.SetX(EyeCoord.X()-genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).X());
			  //       DistanceEyeCoreVector.SetY(EyeCoord.Y()-genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).Y());
			  //       DistanceEyeCoreVector.SetZ(EyeCoord.Z()-genEvent.GetCoreAtAltitudeSiteCS(EyeAltitude).Z());

			  //       double distanceEyeCore = DistanceEyeCoreVector.Mag();
			  //       VectorTriggereddistanceEyeCore.push_back(distanceEyeCore/1000);
			  //       DistanceFdForVectorHEAT = distanceEyeCore/1000; // in km

			  //   	LgEHEAT = log10(eye.GetFdRecShower().GetEnergy()); //reconstructed log10E


			  //       VectorlgERec.push_back(LgEHEAT);
			  //   	double GpsSec = eye.GetGPSSecond();

			  //   	distanceRecHEAT = eye.GetFdRecGeometry().GetCoreEyeDistance();
			  //       double distance = eye.GetGenGeometry().GetCoreEyeDistance();


			  //       RpHEAT = eye.GetFdRecGeometry().GetRp();
			  //       VectorRp.push_back(RpHEAT/1000);

		   //  		double Chi0 = eye.GetFdRecGeometry().GetChi0();

		   //  		double NPix = eye.GetFdRecGeometry().GetSDPNdF(); //nr of triggered pixels
		   //  		double NPixGeo = eye.GetFdRecGeometry().GetNTimeFitPixels();
		   //  		float XMaxFD = eye.GetFdRecShower().GetXmax();

			  //       double nrtrigpixel = eye.GetFdRecPixel().GetNumberOfSDPFitPixels(); //total nr of triggered pixels for each event
			  //       //cout<<"ntrigpixels = "<<nrtrigpixel<<endl;

			  //       vector<double>pixelT;
			  //       vector<double>pixelidvar;
			  //       vector<double>pixelChivar;

			  //       //;

			  //       const std::vector<UShort_t>& pixelID = eye.GetFdRecPixel().GetID();
			  //       const std::vector<double>& pixelTime = eye.GetFdRecPixel().GetTime();
			  //       const std::vector<EPixelStatus>& pixStatus = eye.GetFdRecPixel().GetStatus();
			  //       const std::vector<double>& ChiPixel = eye.GetFdRecPixel().GetChi();



			  //   	for(int k=0; k<pixelTime.size(); k++){
				 //    	if(pixelTime[k]!=0 && pixStatus[k]==4){ //we only keep the pixel part of time fit:pixStatus[i]==4
					//       	double testtimevar=pixelTime[k];
					//         pixelT.push_back(testtimevar); //vector with time of pixels
					//         double testpixelvar = pixelID[k];
					//         pixelidvar.push_back(testpixelvar); //vector with id of pixel
					//         double pixelchivar = ChiPixel[k];
					//     	pixelChivar.push_back(180-pixelchivar*180/3.14); //vector with Chi of pixels

					// 	    //cout<<"k.."<<k<<" id: "<<pixelID[k]<<" time: "<<pixelTime[k] <<" chi:" <<180-ChiPixel[k]*180/3.14<<" eye: "
					// 	    //<<eyeID<<"size: "<<pixelTime.size()<<endl;
				 //        }
				 //    }




			  //   //cout<<NPix<<endl;
			  //       data.xmaxFD.push_back(XMaxFD);
			  //       //data.xmaxMC.push_back(XMaxMC);
			  //       data.dedxmaxMC.push_back(dEdXmaxMC);
			  //       data.xmaxGH.push_back(XMaxGH);
			  //       data.nmaxGH.push_back(NMaxGH);
			  //       //data.x1gen.push_back(X1gen);

			  //       data.SlantDepthMC.push_back(slantdepthMC);
			  //       data.EnergyDepositMC.push_back(energydepositMC);

					// data.dec.push_back(Dec);

			  //       data.Chi0RecGeom.push_back(Chi0);
			  //       data.distFD.push_back(distance);
			  //       data.EyeID.push_back(eyeID);
			  //       data.TimePixels.push_back(pixelT);
			  //       data.IdPixel.push_back(pixelidvar);
			  //       data.CHIPixel.push_back(pixelChivar);

		   //  }



	    }

    	// data.EVENTid.push_back(eventID);
    	// data.lgEMC.push_back(LgEMC);
    	// data.lgErec.push_back(LgE); //rec energy only for Coihueco
    	// data.RpFD.push_back(Rp); //Rp only for Coihueco
    	// data.fddistance.push_back(DistanceFdForVector);
    	// data.zenith.push_back(Zenith);
    	// data.azimuth.push_back(Azimuth);
    	// data.generateddistance.push_back(GenerateddistanceEyeCore/1000);
    	// data.distFDRec.push_back(distanceRec);

    	// data.EVENTidHEAT.push_back(eventIDHEAT);
    	// data.lgEMCHEAT.push_back(LgEMC);
    	// data.lgErecHEAT.push_back(LgEHEAT);
    	// data.RpFDHEAT.push_back(RpHEAT);
    	// data.fddistanceHEAT.push_back(DistanceFdForVectorHEAT);
    	// data.zenithHEAT.push_back(Zenith);
    	// data.generateddistanceHEAT.push_back(GenerateddistanceEyeCore/1000);
    	// data.distFDRecHEAT.push_back(distanceRecHEAT);

    	// data.x1gen.push_back(X1gen);

    	// data.xmaxMC.push_back(XMaxMC);

    	UtilityFunctions::ShowProgress(i, nEv - 1);
	
	}


	// TVectorD RootVectorEVENTidCO(data.EVENTid.size());
	// TVectorD RootVectorlgEMCCO(data.EVENTid.size());
	// TVectorD RootVectorlgErecCO(data.EVENTid.size());
	// TVectorD RootVectorRpCO(data.EVENTid.size());
	// TVectorD RootVectorDistanceFdCO(data.EVENTid.size());
	// TVectorD RootVectorRecDistanceFdCO(data.EVENTid.size());
	// TVectorD RootVectorGenDistanceFdCO(data.EVENTid.size());
	// TVectorD RootVectorZenith(data.EVENTid.size());
	// TVectorD RootVectorAzimuth(data.EVENTid.size());

	// TVectorD RootVectorEVENTidHEAT(data.EVENTid.size());
	// TVectorD RootVectorlgEMCHEAT(data.EVENTid.size());
	// TVectorD RootVectorlgErecHEAT(data.EVENTid.size());
	// TVectorD RootVectorRpHEAT(data.EVENTid.size());
	// TVectorD RootVectorDistanceFdHEAT(data.EVENTid.size());
	// TVectorD RootVectorGenDistanceFdHEAT(data.EVENTid.size());
	// TVectorD RootVectorRecDistanceFdHEAT(data.EVENTid.size());
	// TVectorD RootVectorZenithHEAT(data.EVENTid.size());

	// //for checking x1 which will be further used to get the height of the first interaction
	// TVectorD RootVectorX1(data.EVENTid.size());

	// //xmax MC (generated from corsika)
	// TVectorD RootVectorXmaxMC(data.EVENTid.size());


	// for(int icheck = 0 ; icheck<  data.EVENTid.size(); icheck++){
	// 	RootVectorEVENTidCO[icheck] = data.EVENTid.at(icheck);
	//     RootVectorlgEMCCO[icheck] = data.lgEMC.at(icheck);
	//     RootVectorlgErecCO[icheck] = data.lgErec.at(icheck);
	//     RootVectorRpCO[icheck] = data.RpFD.at(icheck);
	//     RootVectorDistanceFdCO[icheck] = data.fddistance.at(icheck);
	//     RootVectorGenDistanceFdCO[icheck] = data.generateddistance.at(icheck);
	//     RootVectorRecDistanceFdCO[icheck] = data.distFDRec.at(icheck);
	//     RootVectorZenith[icheck] = data.zenith.at(icheck);
	//     RootVectorAzimuth[icheck] = data.azimuth.at(icheck);

	//     RootVectorEVENTidHEAT[icheck] = data.EVENTidHEAT.at(icheck);
	//     RootVectorlgEMCHEAT[icheck] = data.lgEMCHEAT.at(icheck);
	//     RootVectorlgErecHEAT[icheck] = data.lgErecHEAT.at(icheck);
	//     RootVectorRpHEAT[icheck] = data.RpFDHEAT.at(icheck);
	//     RootVectorDistanceFdHEAT[icheck] = data.fddistanceHEAT.at(icheck);
	//     RootVectorGenDistanceFdHEAT[icheck] = data.generateddistanceHEAT.at(icheck);
	//     RootVectorRecDistanceFdHEAT[icheck] = data.distFDRecHEAT.at(icheck);
	//     RootVectorZenithHEAT[icheck] = data.zenithHEAT.at(icheck);


	//     RootVectorX1[icheck] = data.x1gen.at(icheck);

	//     RootVectorXmaxMC[icheck] = data.xmaxMC.at(icheck);
	// //cout<<"x1 nr elements: "<<RootVectorX1.GetNoElements()<<endl;
	   

	//     //cout<<"check id? "<<icheck<<" EvID: "<<data.EVENTid.at(icheck)<<" EMC: "<<data.lgEMC.at(icheck)<<" Erec :"<<data.lgErec.at(icheck)<<" Rp: "<<
	//     //data.RpFD.at(icheck)<<endl;
	// }

	//for(int icheck = 0 ; icheck<  data.EVENTid.size(); icheck++)
	  //cout<<"ev id: "<<RootVectorEVENTidCO[icheck]<<"  "<<" Emc: "<< RootVectorlgEMCCO[icheck]<<" e rec: "<<RootVectorlgErecCO[icheck]<<" rp "<<
	    //  RootVectorRpCO[icheck]<<endl;


	//==============TEMPORARY COMMENTED==========

	// RootVectorEVENTidCO.Write("RootVectorEVENTidCO");
	// RootVectorlgEMCCO.Write("RootVectorlgEMCCO");
	// RootVectorlgErecCO.Write("RootVectorlgErecCO");
	// RootVectorRpCO.Write("RootVectorRpCO");
	// RootVectorDistanceFdCO.Write("RootVectorDistanceFdCO");
	// RootVectorGenDistanceFdCO.Write("RootVectorGenDistanceFdCO");
	// RootVectorRecDistanceFdCO.Write("RootVectorRecDistanceFdCO");
	// RootVectorZenith.Write("RootVectorZenith");
	// RootVectorAzimuth.Write("RootVectorAzimuth");


	// RootVectorEVENTidHEAT.Write("RootVectorEVENTidHEAT");
	// RootVectorlgEMCHEAT.Write("RootVectorlgEMCHEAT");
	// RootVectorlgErecHEAT.Write("RootVectorlgErecHEAT");
	// RootVectorRpHEAT.Write("RootVectorRpHEAT");
	// RootVectorDistanceFdHEAT.Write("RootVectorDistanceFdHEAT");
	// RootVectorGenDistanceFdHEAT.Write("RootVectorGenDistanceFdHEAT");
	// RootVectorRecDistanceFdHEAT.Write("RootVectorRecDistanceFdHEAT");
	// RootVectorZenithHEAT.Write("RootVectorZenithHEAT");

	// RootVectorX1.Write("RootVectorX1");

	// RootVectorXmaxMC.Write("RootVectorXmaxMC");


 
	// double min;

	// for(int i=0; i<data.TimePixels.size(); ++i){
 //    //if(data.TimePixels[i].size()==0)
 //      //    continue;
 //    for(int j=0; j<data.TimePixels[i].size(); ++j){
 //        min=data.TimePixels[i][0];
 //    	//cout<<data.TimePixels[i][j]<<"  ";
 //        if(data.TimePixels[i][j]<min)
 //            min = data.TimePixels[i][j];
	// 		//cout<<min<<" ";
 //    }
 //    data.minTime.push_back(min);
 //  //  cout<<endl;
 //    } 


	// double minchivar;

	// for(int i=0; i<data.CHIPixel.size(); ++i){
	//     //if(data.TimePixels[i].size()==0)
	//       //    continue;
	//     for(int j=0; j<data.CHIPixel[i].size(); ++j){
	//         minchivar=data.CHIPixel[i][0];
	          
	//         if(data.CHIPixel[i][j]<minchivar)
	//             minchivar = data.CHIPixel[i][j];
	// 	}
	//     data.minChi.push_back(minchivar);	    
	// } 



	// for(int ivar=0; ivar<VectorGenerateddistanceEyeCore.size();ivar++){ 
	// 	const double max_value_GenerateddistanceEyeCore = *max_element(VectorGenerateddistanceEyeCore.begin(), VectorGenerateddistanceEyeCore.end());
	// 	int nbins = floor(max_value_GenerateddistanceEyeCore/0.1);
	// 	GeneratedCorePositionHisto->Fill(VectorGenerateddistanceEyeCore[ivar]);
	// }

	// GeneratedCorePositionHisto->Write("MyGeneratedCorePositionHisto");

	// for(int ivar2=0; ivar2<VectorTriggereddistanceEyeCore.size(); ivar2++){ 
	// 	const double max_value_TriggereddistanceEyeCore = *max_element(VectorTriggereddistanceEyeCore.begin(), VectorTriggereddistanceEyeCore.end());
	// 	int nbins = floor(max_value_TriggereddistanceEyeCore/0.1);
	// 	TriggeredCorePositionHisto->Fill(VectorTriggereddistanceEyeCore[ivar2]);
	// }

	// TriggeredCorePositionHisto->Write("MyTriggeredCorePostionHisto");


	cout<<"vector size trig: "<<VectorTriggereddistanceEyeCore.size()<<endl;



	cout << " " << nUnique << " of " << nEv << " read in, " << combinedCount
	    << " were stereo." << endl;

   TPCGFup->Write();
	
	return nUnique;
}
