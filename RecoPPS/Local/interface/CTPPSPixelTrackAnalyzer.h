#ifndef CTPPS_RPIX_Track_Analyzer_h
#define CTPPS_RPIX_Track_Analyzer_h

//#define _SIMU_


#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigi.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSPixelDigiCollection.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelCluster.h"
//#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
//#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "FWCore/Utilities/interface/InputTag.h"


#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"

//#include "Geometry/VeryForwardGeometry/interface/CTPPSPixelSimTopology.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TFile.h"

class TGraph;

//namespace edm {
//  class ParameterSet; class Event; class EventSetup;}

class CTPPSPixelTrackAnalyzer : public edm::one::EDAnalyzer<>{
  
 public:
  explicit CTPPSPixelTrackAnalyzer(const edm::ParameterSet& pset);
  virtual ~CTPPSPixelTrackAnalyzer();
  void endJob();
  void beginJob();
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

  TGraph *_gr_avg_45;
  TGraph *_gr_avg_56;
  TH2F *_h_45;
  TH2F *_h_56;

  TH2F *_h_D_45;
  TH2F *_h_D_56;
  TH1F *_h_D_ADC;

  TH2F *_h_RH_45;
  TH2F *_h_RH_45_220_0;
  TH2F *_h_RH_45_220_1;
  TH2F *_h_RH_45_220_2;
  TH2F *_h_RH_45_220_3;
  TH2F *_h_RH_45_220_4;
  TH2F *_h_RH_45_220_5;
  TH2F *_h_RH_56;
  TH2F *_h_RH_45_210_0;
  TH2F *_h_RH_56_210_0;
  TH2F *_h_RH_56_220_0;
  TH2F *_h_RH_45_210_1;
  TH2F *_h_RH_56_210_1;
  TH2F *_h_RH_56_220_1;
  TH2F *_h_RH_45_210_2;
  TH2F *_h_RH_56_210_2;
  TH2F *_h_RH_56_220_2;
  TH2F *_h_RH_45_210_3;
  TH2F *_h_RH_56_210_3;
  TH2F *_h_RH_56_220_3;
  TH2F *_h_RH_45_210_4;
  TH2F *_h_RH_56_210_4;
  TH2F *_h_RH_56_220_4;
  TH2F *_h_RH_45_210_5;
  TH2F *_h_RH_56_210_5;
  TH2F *_h_RH_56_220_5;

#ifdef _SIMU_
  TH2F *_h_SH_45;
  TH2F *_h_SH_56;
#endif

  TH2F *_h_fRH_45xz;
  TH2F *_h_fRH_56xz;
  TH2F *_h_fRH_45yz;
  TH2F *_h_fRH_56yz;

  TH2F *_h_trXY0_45_210;
  TH2F *_h_trXY0_45_220;
  TH2F *_h_trXY0_56_210;
  TH2F *_h_trXY0_56_220;

  TH1F *_h_hits_in_track_45n;
  TH1F *_h_hits_in_track_45f;
  TH1F *_h_hits_in_track_56n;
  TH1F *_h_hits_in_track_56f;

  TH1F *_h_CL_size;
  TH1F *_h_ELE_size1;
  TH1F *_h_ELE_size2;
  TH1F *_h_ELEsum_size2;
  TH1F *_h_ELEmin_size2;
  TH1F *_h_ELEmax_size2;

  TH1F *_h_detid_track;

 private:
  std::string _outFileName;

  edm::EDGetTokenT< edm::DetSetVector<CTPPSPixelLocalTrack> > pixelTrackToken_;
  edm::EDGetTokenT< edm::DetSetVector<CTPPSPixelDigi> > pixelDigiToken_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelRecHit>> tokenCTPPSPixelRecHit_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelCluster>> tokenCTPPSPixelCluster_;

#ifdef _SIMU_
  edm::EDGetTokenT<edm::PSimHitContainer> psim_token;
#endif

  TFile* _outFile;
  unsigned int _verbosity;
  int _bunchCrossing;
  int _ls;

  unsigned int tracks45_per_ls[10000];
  unsigned int tracks56_per_ls[10000];
  unsigned int events_per_ls[10000];
  
  unsigned int ls_max;

  unsigned int tracks45posYtot;
  unsigned int tracks45negYtot;
  unsigned int tracks45recovered;

};

#endif    
