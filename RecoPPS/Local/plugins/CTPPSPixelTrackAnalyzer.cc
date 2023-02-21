#include "RecoPPS/Local/interface/CTPPSPixelTrackAnalyzer.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"

#include <iostream>
#include <string>

//#define _CMS_coord_

using namespace std;

CTPPSPixelTrackAnalyzer:: CTPPSPixelTrackAnalyzer(const edm::ParameterSet& pset)
// : theRPixDetTopology_(pset)
{

  _outFile = new TFile("myFile.root","RECREATE");
  _verbosity = pset.getUntrackedParameter<unsigned int> ("Verbosity");
  if(_outFile->IsOpen()) cout<<"file open!"<<endl;
  else cout<<"*** Error in opening file ***"<<endl;
  
  auto tagPixelDigi = pset.getParameter<edm::InputTag>("tagPixelDigi"); 
  pixelDigiToken_ = consumes<edm::DetSetVector<CTPPSPixelDigi> >( tagPixelDigi);

  auto tagPixelRecHit = pset.getParameter<edm::InputTag>("tagPixelRecHit");
  tokenCTPPSPixelRecHit_ = consumes<edm::DetSetVector<CTPPSPixelRecHit> >(tagPixelRecHit);

  auto tagPixelCluster = pset.getParameter<edm::InputTag>("tagPixelCluster");
  tokenCTPPSPixelCluster_ = consumes<edm::DetSetVector<CTPPSPixelCluster> >(tagPixelCluster);

  auto tagPixelTrack = pset.getParameter<edm::InputTag>("tagPixelTrack"); 
  if (not tagPixelTrack.label().empty()){
    pixelTrackToken_   = consumes< edm::DetSetVector<CTPPSPixelLocalTrack> >  (tagPixelTrack);
  }
 		    
#ifdef _SIMU_
  psim_token = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "CTPPSPixelHits"));
#endif
  tracks45posYtot = 0;
  tracks45negYtot = 0;
  tracks45recovered = 0;
}

CTPPSPixelTrackAnalyzer::~CTPPSPixelTrackAnalyzer(){
}

void CTPPSPixelTrackAnalyzer::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>( "tagPixelTrack"  , edm::InputTag( "ctppsPixelLocalTracks"   ) );
  desc.add<edm::InputTag>( "tagPixelDigi"  , edm::InputTag( "ctppsPixelDigis"   ) );
  desc.add<edm::InputTag>( "tagPixelRecHit"  , edm::InputTag( "ctppsPixelRecHits"   ) );
  desc.add<edm::InputTag>( "tagPixelCluster"  , edm::InputTag( "ctppsPixelClusters"   ) );
#ifdef _SIMU_
  desc.add<std::string>("mixLabel", "mix");
  desc.add<std::string>("InputCollection", "g4SimHitsCTPPSPixelHits");
#endif


  desc.addUntracked<unsigned int>("Verbosity",0);
  descriptions.add( "ctppsPixelTrackAnalyzer", desc );
}

void CTPPSPixelTrackAnalyzer::beginJob(){
  for(unsigned int i=0; i<10000; i++){
    tracks45_per_ls[i]=0;
    tracks56_per_ls[i]=0;
    events_per_ls[i]=0;
  }
  ls_max=10000;
  _h_45 = new TH2F("Tracks45","Tracks45",300,-0.5,299.5,200,-0.5,199.5);
  _h_56 = new TH2F("Tracks56","Tracks56",300,-0.5,299.5,200,-0.5,199.5);
  _h_D_45 = new TH2F("Digis45","Digis45",300,-0.5,299.5,200,-0.5,199.5);
  _h_D_56 = new TH2F("Digis56","Digis56",300,-0.5,299.5,200,-0.5,199.5);
  _h_D_ADC = new TH1F("DigisADC","DigisADC",257,-0.5,256.5);
  _h_RH_45 = new TH2F("RH45","RH45",100,-15,15,100,-15,15);
  _h_RH_45_220_0 = new TH2F("RH45f0","RH45f0",100,-15,15,100,-15,15);
  _h_RH_45_220_1 = new TH2F("RH45f1","RH45f1",100,-15,15,100,-15,15);
  _h_RH_45_220_2 = new TH2F("RH45f2","RH45f2",100,-15,15,100,-15,15);
  _h_RH_45_220_3 = new TH2F("RH45f3","RH45f3",100,-15,15,100,-15,15);
  _h_RH_45_220_4 = new TH2F("RH45f4","RH45f4",100,-15,15,100,-15,15);
  _h_RH_45_220_5 = new TH2F("RH45f5","RH45f5",100,-15,15,100,-15,15);
  _h_RH_45_210_0 = new TH2F("RH45n0","RH45n0",100,-15,15,100,-15,15);
  _h_RH_56_220_0 = new TH2F("RH56f0","RH56f0",100,-15,15,100,-15,15);
  _h_RH_56_210_0 = new TH2F("RH56n0","RH56n0",100,-15,15,100,-15,15);

  _h_RH_45_210_1 = new TH2F("RH45n1","RH45n1",100,-15,15,100,-15,15);
  _h_RH_56_220_1 = new TH2F("RH56f1","RH56f1",100,-15,15,100,-15,15);
  _h_RH_56_210_1 = new TH2F("RH56n1","RH56n1",100,-15,15,100,-15,15);

  _h_RH_45_210_2 = new TH2F("RH45n2","RH45n2",100,-15,15,100,-15,15);
  _h_RH_56_220_2 = new TH2F("RH56f2","RH56f2",100,-15,15,100,-15,15);
  _h_RH_56_210_2 = new TH2F("RH56n2","RH56n2",100,-15,15,100,-15,15);

  _h_RH_45_210_3 = new TH2F("RH45n3","RH45n3",100,-15,15,100,-15,15);
  _h_RH_56_220_3 = new TH2F("RH56f3","RH56f3",100,-15,15,100,-15,15);
  _h_RH_56_210_3 = new TH2F("RH56n3","RH56n3",100,-15,15,100,-15,15);

  _h_RH_45_210_4 = new TH2F("RH45n4","RH45n4",100,-15,15,100,-15,15);
  _h_RH_56_220_4 = new TH2F("RH56f4","RH56f4",100,-15,15,100,-15,15);
  _h_RH_56_210_4 = new TH2F("RH56n4","RH56n4",100,-15,15,100,-15,15);

  _h_RH_45_210_5 = new TH2F("RH45n5","RH45n5",100,-15,15,100,-15,15);
  _h_RH_56_220_5 = new TH2F("RH56f5","RH56f5",100,-15,15,100,-15,15);
  _h_RH_56_210_5 = new TH2F("RH56n5","RH56n5",100,-15,15,100,-15,15);

  _h_hits_in_track_45n = new TH1F("hits_in_track_45n","hits_in_track_45n",10,-0.5,9.5);
  _h_hits_in_track_45f = new TH1F("hits_in_track_45f","hits_in_track_45f",10,-0.5,9.5);
  _h_hits_in_track_56n = new TH1F("hits_in_track_56n","hits_in_track_56n",10,-0.5,9.5);
  _h_hits_in_track_56f = new TH1F("hits_in_track_56f","hits_in_track_56f",10,-0.5,9.5);

  _h_detid_track  = new TH1F("detidTrack","detidTrack",300,1900,2100);

  _h_RH_56 = new TH2F("RH56","RH56",100,-15,15,100,-15,15);
#ifdef _SIMU_
  _h_SH_45 = new TH2F("SH45","SH45",100,-15,15,100,-15,15);
  _h_SH_56 = new TH2F("SH56","SH56",100,-15,15,100,-15,15);
#endif
#ifdef _CMS_coord_
//CMS coordinates 
  _h_fRH_45xz = new TH2F("fRH45xz","fRH45xz",50000,210000,221000,3000,-40,0);
  _h_fRH_56xz = new TH2F("fRH56xz","fRH56xz",50000,-221000,-210000,3000,-40,0);
  _h_fRH_45yz = new TH2F("fRH45yz","fRH45yz",50000,210000,221000,3000,-30,30);
  _h_fRH_56yz = new TH2F("fRH56yz","fRH56yz",50000,-221000,-210000,3000,-30,30);
#else
// LHC coordinates
  _h_fRH_45xz = new TH2F("fRH45xz","fRH45xz",50000,-221000,-210000,3000,0,60);
  _h_fRH_56xz = new TH2F("fRH56xz","fRH56xz",50000,210000,221000,3000,0,60);
  _h_fRH_45yz = new TH2F("fRH45yz","fRH45yz",50000,-221000,-210000,3000,-30,30);
  _h_fRH_56yz = new TH2F("fRH56yz","fRH56yz",50000,210000,221000,3000,-30,30);
#endif

  _h_trXY0_45_210 = new TH2F("trXY0_45_210","trXY0_45_210",200,0,30,200,-15,15);
  _h_trXY0_45_220 = new TH2F("trXY0_45_220","trXY0_45_220",200,0,30,200,-15,15);
  _h_trXY0_56_210 = new TH2F("trXY0_56_210","trXY0_56_210",200,0,30,200,-15,15);
  _h_trXY0_56_220 = new TH2F("trXY0_56_220","trXY0_56_220",200,0,30,200,-15,15);


  _h_CL_size = new TH1F("CLsize","CLsize",10,-0.5,9.5);
  _h_ELE_size1 = new TH1F("ELEsize1","ELEsize1",200,-0.5,100000);
  _h_ELE_size2 = new TH1F("ELEsize2","ELEsize2",200,-0.5,100000);
  _h_ELEsum_size2 = new TH1F("ELEsumsize2","ELEsumsize2",200,-0.5,100000);
  _h_ELEmin_size2 = new TH1F("ELEminsize2","ELEminsize2",200,-0.5,100000);
  _h_ELEmax_size2 = new TH1F("ELEmaxsize2","ELEmaxsize2",200,-0.5,100000);


}

void CTPPSPixelTrackAnalyzer::endJob(){

  float avg_tks_45[10000];
  float avg_tks_56[10000];
  float ls[10000];
  
  for(unsigned int i=0; i<10000; i++){
    avg_tks_45[i]=0;
    avg_tks_56[i]=0;
    ls[i]=float(i);
    if(events_per_ls[i]>0){
      avg_tks_45[i]=float(tracks45_per_ls[i])/float(events_per_ls[i]);
      avg_tks_56[i]=float(tracks56_per_ls[i])/float(events_per_ls[i]);
    }

    

  }
  _gr_avg_45 = new TGraph(10000,ls,avg_tks_45);
  _gr_avg_56 = new TGraph(10000,ls,avg_tks_56);




  _outFile->cd();



   _gr_avg_45->SetName("AvgPxlTrk45");
   _gr_avg_45->SetTitle("Avg pxl trk number 45");
   _gr_avg_45->SetMarkerColor(2); 
   _gr_avg_45->SetMarkerStyle(3);
   _gr_avg_45->GetXaxis()->SetLimits(0, ls_max+10);
   _gr_avg_45->GetXaxis()->SetTitle("LS");
//   _gr_avg_45->GetYaxis()->SetLimits(ylo,yhi);
//   _gr_avg_45->GetYaxis()->SetTitle(axis);
//   _gr_avg_45->Draw("A*");

   _gr_avg_56->SetName("AvgPxlTrk56");
   _gr_avg_56->SetTitle("Avg pxl trk number 56");
   _gr_avg_56->SetMarkerColor(4); 
   _gr_avg_56->SetMarkerStyle(3);
   _gr_avg_56->GetXaxis()->SetLimits(0, ls_max+10);
   _gr_avg_56->GetXaxis()->SetTitle("LS");

   _h_45->SetName("PxlTrk45");
   _h_45->SetTitle("pxl trk number 45");
   _h_45->GetXaxis()->SetTitle("LS");
   _h_45->GetYaxis()->SetTitle("N tracks");

   _h_56->SetName("PxlTrk56");
   _h_56->SetTitle("pxl trk number 56");
   _h_56->GetXaxis()->SetTitle("LS");
   _h_56->GetYaxis()->SetTitle("N tracks");

   _h_CL_size->SetName("CLsize");

  _gr_avg_45->Write();
  _gr_avg_56->Write();
  _h_45->Write();
  _h_56->Write();

  _h_D_45->Write();
  _h_D_56->Write();
  _h_D_ADC->Write();
  _h_detid_track->Write();
  _h_RH_45->Write();
  _h_RH_45_220_0->Write();
  _h_RH_45_220_1->Write();
  _h_RH_45_220_2->Write();
  _h_RH_45_220_3->Write();
  _h_RH_45_220_4->Write();
  _h_RH_45_220_5->Write();
  _h_RH_56->Write();
  _h_RH_45_210_0->Write();
  _h_RH_56_210_0->Write();
  _h_RH_56_220_0->Write();

  _h_RH_45_210_1->Write();
  _h_RH_56_210_1->Write();
  _h_RH_56_220_1->Write();

  _h_RH_45_210_2->Write();
  _h_RH_56_210_2->Write();
  _h_RH_56_220_2->Write();

  _h_RH_45_210_3->Write();
  _h_RH_56_210_3->Write();
  _h_RH_56_220_3->Write();

  _h_RH_45_210_4->Write();
  _h_RH_56_210_4->Write();
  _h_RH_56_220_4->Write();

  _h_RH_45_210_5->Write();
  _h_RH_56_210_5->Write();
  _h_RH_56_220_5->Write();

#ifdef _SIMU_ 
  _h_SH_45->Write();
  _h_SH_56->Write();
#endif
  _h_fRH_45xz->Write();
  _h_fRH_56xz->Write();
  _h_fRH_45yz->Write();
  _h_fRH_56yz->Write();

  _h_trXY0_45_210->Write();
  _h_trXY0_45_220->Write();
  _h_trXY0_56_210->Write();
  _h_trXY0_56_220->Write();

  _h_hits_in_track_45n->Write();
  _h_hits_in_track_45f->Write();
  _h_hits_in_track_56n->Write();
  _h_hits_in_track_56f->Write();

  _h_CL_size->Write();
  _h_ELE_size1->Write();
  _h_ELE_size2->Write();
  _h_ELEsum_size2->Write();
  _h_ELEmin_size2->Write();
  _h_ELEmax_size2->Write();

 _outFile->Close();
  delete _outFile;
}

void  CTPPSPixelTrackAnalyzer::analyze(const edm::Event & event, const edm::EventSetup& eventSetup){

#ifdef _SIMU_
  cout << " SIMU defined " << endl;
#endif
#ifndef _SIMU_
  cout << " SIMU NOT defined " << endl;
#endif

  if(_verbosity)
    cout << "--- Run: " << event.id().run()
       << " Event: " << event.id().event()
       << " Bunch crossing: " << event.bunchCrossing()
       << " Lumi block: " << event.luminosityBlock() << endl;
  

// if(event.luminosityBlock()>10000)throw cms::Exception("Track analyzer") << "lumi block > 10000";

  _bunchCrossing = event.bunchCrossing();
  _ls = event.luminosityBlock();
  edm::Handle<edm::DetSetVector<CTPPSPixelDigi> > digis;
  event.getByToken(pixelDigiToken_, digis);

  unsigned int digis45=0;
  unsigned int digis56=0;

 // Loop on digis
  edm::DetSetVector<CTPPSPixelDigi>::const_iterator digiDSViter = digis->begin();
  for (; digiDSViter != digis->end(); digiDSViter++) {
    CTPPSPixelDetId detIdObject(digiDSViter->detId());
    edm::DetSet<CTPPSPixelDigi>::const_iterator begin = (*digiDSViter).begin();
    edm::DetSet<CTPPSPixelDigi>::const_iterator end = (*digiDSViter).end();
    for (edm::DetSet<CTPPSPixelDigi>::const_iterator di = begin; di != end; di++) {
    // Detector ID
      _h_D_ADC->Fill((*di).adc());
      if(detIdObject.arm()==0){
	digis45++;
	_h_D_45->Fill((*di).row(),(*di).column());
      }else{
	digis56++;
	_h_D_56->Fill((*di).row(),(*di).column());
      }
    }
  }

//------------------------------------------------------
  edm::Handle<edm::DetSetVector<CTPPSPixelCluster> > clusters;
  event.getByToken(tokenCTPPSPixelCluster_, clusters);
  edm::DetSetVector<CTPPSPixelCluster>::const_iterator clDSViter = clusters->begin();

  for (; clDSViter != clusters->end(); clDSViter++) {
    CTPPSPixelDetId detIdObject(clDSViter->detId());
    edm::DetSet<CTPPSPixelCluster>::const_iterator begin = (*clDSViter).begin();
    edm::DetSet<CTPPSPixelCluster>::const_iterator end = (*clDSViter).end();
    for (edm::DetSet<CTPPSPixelCluster>::const_iterator cl = begin; cl != end; cl++) {
      _h_CL_size->Fill((*cl).size());

      if( (*cl).size()<=2){
//	cout << (*cl).pixelADC(0) << endl;
	if( (*cl).size()==1) _h_ELE_size1->Fill((*cl).pixelADC(0));
	else{
	  _h_ELE_size2->Fill((*cl).pixelADC(0));
	  _h_ELE_size2->Fill((*cl).pixelADC(1));
	  _h_ELEsum_size2->Fill((*cl).pixelADC(1)+(*cl).pixelADC(0));
	  _h_ELEmin_size2->Fill(fmin((*cl).pixelADC(1),(*cl).pixelADC(0)));
	  _h_ELEmax_size2->Fill(fmax((*cl).pixelADC(1),(*cl).pixelADC(0)));
	}
      }

    }
  }


//-------------------------------------------------------
  edm::Handle<edm::DetSetVector<CTPPSPixelRecHit> > recHits;
  event.getByToken(tokenCTPPSPixelRecHit_, recHits);
  edm::DetSetVector<CTPPSPixelRecHit>::const_iterator rhDSViter = recHits->begin();
  for (; rhDSViter != recHits->end(); rhDSViter++) {
    CTPPSPixelDetId detIdObject(rhDSViter->detId());
    edm::DetSet<CTPPSPixelRecHit>::const_iterator begin = (*rhDSViter).begin();
    edm::DetSet<CTPPSPixelRecHit>::const_iterator end = (*rhDSViter).end();
    for (edm::DetSet<CTPPSPixelRecHit>::const_iterator rh = begin; rh != end; rh++) {
    // Detector ID

      if(detIdObject.arm()==0){
	//	std::cout <<detIdObject.rawId() << " " << detIdObject.station() << " " << detIdObject.rp() << " " <<detIdObject.plane() << std::endl; 
	_h_RH_45->Fill((*rh).point().x(),(*rh).point().y());
	if(detIdObject.station() == 2 && detIdObject.plane() == 0){
	  	_h_RH_45_220_0->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 0){
	  	_h_RH_45_210_0->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 1){
	  	_h_RH_45_220_1->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 1){
	  	_h_RH_45_210_1->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 2){
	  	_h_RH_45_220_2->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 2){
	  	_h_RH_45_210_2->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 3){
	  	_h_RH_45_220_3->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 3){
	  	_h_RH_45_210_3->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 4){
	  	_h_RH_45_220_4->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 4){
	  	_h_RH_45_210_4->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 5){
	  	_h_RH_45_220_5->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 5){
	  	_h_RH_45_210_5->Fill((*rh).point().x(),(*rh).point().y());
	}
      }else{
	_h_RH_56->Fill((*rh).point().x(),(*rh).point().y());
	if(detIdObject.station() == 2 && detIdObject.plane() == 0){
	  	_h_RH_56_220_0->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 0){
	  	_h_RH_56_210_0->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 1){
	  	_h_RH_56_220_1->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 1){
	  	_h_RH_56_210_1->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 2){
	  	_h_RH_56_220_2->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 2){
	  	_h_RH_56_210_2->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 3){
	  	_h_RH_56_220_3->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 3){
	  	_h_RH_56_210_3->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 4){
	  	_h_RH_56_220_4->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 4){
	  	_h_RH_56_210_4->Fill((*rh).point().x(),(*rh).point().y());
	}

	if(detIdObject.station() == 2 && detIdObject.plane() == 5){
	  	_h_RH_56_220_5->Fill((*rh).point().x(),(*rh).point().y());
	}
	if(detIdObject.station() == 0 && detIdObject.plane() == 5){
	  	_h_RH_56_210_5->Fill((*rh).point().x(),(*rh).point().y());
	}






      }
    }
  }

#ifdef _SIMU_

  edm::Handle<edm::PSimHitContainer> simHits;
  event.getByToken(psim_token, simHits);

  for (vector<PSimHit>::const_iterator hit = simHits->begin(); hit != simHits->end(); hit++) {
    CTPPSPixelDetId detIdObject((*hit).detUnitId());
    if(detIdObject.arm()==0){
      _h_SH_45->Fill((*hit).entryPoint().x(),(*hit).entryPoint().y());
    }else{
      _h_SH_56->Fill((*hit).entryPoint().x(),(*hit).entryPoint().y());
    }
  }
#endif




  edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > inputPixelTracks;

  unsigned int tracks45=0;
  unsigned int tracks45posY=0;
  unsigned int tracks45negY=0;
  unsigned int tracks56=0;
  if (not pixelTrackToken_.isUninitialized()){
    event.getByToken( pixelTrackToken_, inputPixelTracks );
    
// process tracks from pixels
    unsigned int Ntra=0;
    for ( const auto& rpv : *inputPixelTracks ) {
      CTPPSPixelDetId detIdObject(rpv.detId());
      const uint32_t rpId = rpv.detId();
            std::cout << rpId << " " ;
	    _h_detid_track->Fill(double(rpId)/1e6);
      for ( const auto& trk : rpv ) {
        if ( !trk.isValid() ) continue;
	
	if(detIdObject.arm() == 0 && detIdObject.station() == 0){
	  _h_trXY0_45_210->Fill(trk.x0(),trk.y0());
	}
	else if(detIdObject.arm() == 0 && detIdObject.station() == 2){
	  _h_trXY0_45_220->Fill(trk.x0(),trk.y0());
	}
	else if(detIdObject.arm() == 1 && detIdObject.station() == 0){
	  _h_trXY0_56_210->Fill(trk.x0(),trk.y0());
	}
	else if(detIdObject.arm() == 1 && detIdObject.station() == 2){
	  _h_trXY0_56_220->Fill(trk.x0(),trk.y0());
	}
	else{}


	Ntra++;
	if(trk.z0() == 0 )tracks45recovered++;
	if(trk.z0() < 0){ // 45
	  if(trk.y0()<0){tracks45negY++;tracks45negYtot++;}
	  if(trk.y0()>0){tracks45posY++;tracks45posYtot++;}
	}
//	std::cout << Ntra << std::endl;
	int ii45=0;
	int ii56=0;
	int ii45n=0;
	int ii56n=0;
	int ii45f=0;
	int ii56f=0;
	if(detIdObject.arm() == 1){ 
	  for ( const auto& fRH : trk.hits() ){
	    for (const auto& ffRH : fRH ){
	      if(ffRH.isRealHit() ){
		if( detIdObject.station() == 0){
		  ii56n++;
		}else{
		  ii56f++;
		}
		ii56++;

		_h_fRH_56xz->Fill(ffRH.globalCoordinates().z(),ffRH.globalCoordinates().x());
		_h_fRH_56yz->Fill(ffRH.globalCoordinates().z(),ffRH.globalCoordinates().y());
	      }
	    }
	  }
	}else{
	  std::cout << "      Track in sector 45 " << std::endl;
	  for ( const auto& fRH : trk.hits() ){
	    for (const auto& ffRH : fRH ){
	      if(ffRH.isRealHit()  ){
		if( detIdObject.station() == 0){
		  ii45n++;
		}else{
		  ii45f++;
		}
		ii45++;
		std::cout  << ffRH.globalCoordinates().x() << " "
			   << ffRH.globalCoordinates().y() << " "
			   << ffRH.globalCoordinates().z() 
			   << std::endl;
	      
		_h_fRH_45xz->Fill(ffRH.globalCoordinates().z(),ffRH.globalCoordinates().x());
		_h_fRH_45yz->Fill(ffRH.globalCoordinates().z(),ffRH.globalCoordinates().y());
	      }
	    }
	  }
	}

	if(rpId & 0x1000000) tracks56++;
	else tracks45++;
	 std::cout << "ii45 " << ii45 << "  ii56 " << ii56 <<std::endl;


	 if(ii45n>0)_h_hits_in_track_45n->Fill(float(ii45n));
	 if(ii45f>0)_h_hits_in_track_45f->Fill(float(ii45f));
	 if(ii56n>0)_h_hits_in_track_56n->Fill(float(ii56n));
	 if(ii56f>0)_h_hits_in_track_56f->Fill(float(ii56f));


      }
     
    }
  }
  
  if(_verbosity){
    cout << "           Tracks 45 - 56 : "<< tracks45 << " - " << tracks56 << endl;
    cout << "           Tracks 45 pos - neg : "<< tracks45posY << " - " << tracks45negY << endl;
    cout << "           Tracks 45 pos - neg - recoevered TOT : "<< tracks45posYtot << " - " << tracks45negYtot << " - " << tracks45recovered << endl;
    cout << "           Digis  45 - 56 : "<< digis45 << " - " << digis56 << endl;
  }

  if(tracks45==0 && digis45>20)  _h_45->Fill(_ls,20.);
  else if(tracks45>19)_h_45->Fill(_ls,19.);
  else _h_45->Fill(_ls,tracks45);
 
  if(tracks56==0 && digis56>20)  _h_56->Fill(_ls,20.);
  else if(tracks56>19)_h_56->Fill(_ls,19.);
  else _h_56->Fill(_ls,tracks56);
 


//incrementing total number of tracks per lumi
  tracks45_per_ls[_ls] += tracks45;
  tracks56_per_ls[_ls] += tracks56;
  
  events_per_ls[_ls]++;

/*

 // Loop on digis
  edm::DetSetVector<CTPPSPixelDigi>::const_iterator digiDSViter = digis->begin();
  for (; digiDSViter != digis->end(); digiDSViter++) {
    CTPPSPixelDetId detIdObject(digiDSViter->detId());
    edm::DetSet<CTPPSPixelDigi>::const_iterator begin = (*digiDSViter).begin();
    edm::DetSet<CTPPSPixelDigi>::const_iterator end = (*digiDSViter).end();
    for (edm::DetSet<CTPPSPixelDigi>::const_iterator di = begin; di != end; di++) {
      // Detector ID
      _arm_digi.push_back(detIdObject.arm());
      _station_digi.push_back(detIdObject.station());
      _rp_digi.push_back(detIdObject.rp());
      _plane_digi.push_back(detIdObject.plane());
      // Pixel data
      _row.push_back(di->row());
      _column.push_back(di->column());
      _adc.push_back(di->adc());
    }
  }

 */
}


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(CTPPSPixelTrackAnalyzer);
