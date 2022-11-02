// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Laura's code
/// \author
/// \since

/*
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "EventFiltering/filterTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/HistogramRegistry.h"

#include <TMath.h>
#include <cmath>
#include <string>
*/
#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//This is an example of a conveient declaration of "using"
using MyCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

const int Nharm = 7;

//STEP 5: Use to write two-particle correlation but with combination
//This is a simple two-particle correlation function filler
//that makes use of both filters and partitions.
//The core part of the 2pc filling now utilises a combination declaration
//that is in principle more efficient.
struct LaurasCode {

  // all defined filters are applied
  Filter trackFilter = nabs(aod::track::eta) < 0.8f && aod::track::pt > 1.0f;
  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f;

  using MyFilteredTracks = soa::Filtered<MyCompleteTracks>;

  Partition<MyFilteredTracks> triggerTracks = aod::track::pt > 4.0f;
  Partition<MyFilteredTracks> assocTracks = aod::track::pt < 4.0f;
  
  //Configurables for the jets
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0,  
                                    "Accepted z-vertex range"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8, 
                                    "eta range"}; //Eta range 0.8
  Configurable<float> cfgJetR{"cfgJetR", 0.4,
                              "jet resolution parameter"}; // jet cone radius
  Configurable<float> cfgJetPtMin{"cfgJetPtMin", 0.05,
    "minimum jet pT constituent cut"}; // minimum jet constituent pT

  double fiducialVolume;                       // 0.8 - jetR

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histograms are defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {
      {"hZvertex", "hZvertex", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hChargedEta", "hChargedEta", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"hChargedPt", "hChargedPt", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"hChargedJetPt", "hChargedJetPt", {HistType::kTH1F, {{nBins, 0., 10.0}}}},//
      {"hChargedJetEta","charged jet eta",{HistType::kTH1F, {{nBins, -1.0, 1.0}}}}, //

      {"hDeltaEta", "hDeltaEta", {HistType::kTH1F, {{nBins, -2., +2}}}},
      {"hDeltaPhi", "hDeltaPhi", {HistType::kTH1F, {{nBins, 0., 2.0*TMath::Pi()}}}},
      {"hDeltaPhiDeltaEta", "hDeltaPhiDeltaEta", {HistType::kTH2F, {{nBins, 0., 2.0*TMath::Pi()},{nBins, -2., +2}}}},
      {"hTPcosDeltaPhi", "hTPcosDeltaPhi", {HistType::kTH1F, {{nBins, -1.0, 1.0}}}},
      {"hTrigg", "hTrigg", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"etaHistogramTrigger", "etaHistogramTrigger", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramTrigger", "ptHistogramTrigger", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"etaHistogramAssoc", "etaHistogramAssoc", {HistType::kTH1F, {{nBins, -1., +1}}}},
      {"ptHistogramAssoc", "ptHistogramAssoc", {HistType::kTH1F, {{nBins, 0., 10.0}}}},
      {"correlationFunction", "correlationFunction", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI}}}}
    }
  };


 /* void init(o2::framework::InitContext&)
  {
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = cfgJetR;

    fiducialVolume = cfgEtaCut - cfgJetR;

  }  
*/

  Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) {
      //To be completely sure, use inner products
      Double_t x1, y1, x2, y2;
      x1 = TMath::Cos( phi1 );
      y1 = TMath::Sin( phi1 );
      x2 = TMath::Cos( phi2 );
      y2 = TMath::Sin( phi2 );
      Double_t lInnerProd  = x1*x2 + y1*y2;
      Double_t lVectorProd = x1*y2 - x2*y1;
      
      Double_t lReturnVal = 0;
      if( lVectorProd > 1e-8 ) {
          lReturnVal = TMath::ACos(lInnerProd);
      }
      if( lVectorProd < -1e-8 ) {
          lReturnVal = -TMath::ACos(lInnerProd);
      }
      
      if( lReturnVal < -TMath::Pi()/2. ) {
          lReturnVal += 2.*TMath::Pi();
      }
      
      return lReturnVal;
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVertexCut); 

  void process(aod::Collision const& collision, MyFilteredTracks const& tracks) 
  //void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::JetFilters>>::iterator const& collision, TrackCandidates const& tracks)
  {

    //jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = cfgJetR;

    jetReclusterer.etaMin = -0.8;
    jetReclusterer.etaMax = 0.8;

    fiducialVolume = cfgEtaCut - cfgJetR;

    jetConstituents.clear();
    jetReclustered.clear();

    double leadingJetPt = -1.0; 
    double leadingTrackPt = -1.0;


    //Fill the event counter
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hZvertex"))->Fill(collision.posZ());

    //Fill the histograms that are based on track info
    for (auto& track : tracks){
      registry.get<TH1>(HIST("hChargedEta"))->Fill(track.eta());
      registry.get<TH1>(HIST("hChargedPt"))->Fill(track.pt());

      if (fabs(track.eta())<cfgEtaCut) // eta cut of constituents
      {
        if (track.pt() > cfgJetPtMin) //jet constituents
        {
          fillConstituents(track, jetConstituents); // ./PWGJE/Core/JetFinder.h, recomb scheme assumed as Escheme with pion mass
        }
        if (track.pt() > leadingTrackPt)
        {
          leadingTrackPt = track.pt();
        }
      }
    }

    //Reconstruct jet from tracks
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    
    //Now we find the leading jet pT in the eta range 
    for (auto& jet : jetReclustered)
    {
      //if (fabs(jet.eta()) < fiducialVolume)
      //{        
        registry.fill(HIST("hChargedJetEta"), jet.eta());
        registry.fill(HIST("hChargedJetPt"), jet.pt());
        if (jet.pt() > leadingJetPt)
        {
          leadingJetPt = jet.pt();
        }
      //}
    }



    
    //The correlation part
    //partitions are not grouped by default
    auto triggerTracksGrouped = triggerTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto assocTracksGrouped = assocTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());

    //Inspect the trigger and associated populations
    for (auto& track : triggerTracksGrouped) { //<- only for a subset
      registry.get<TH1>(HIST("etaHistogramTrigger"))->Fill(track.eta()); //<- this should show the selection
      registry.get<TH1>(HIST("ptHistogramTrigger"))->Fill(track.pt());
      registry.get<TH1>(HIST("hTrigg"))->Fill(track.pt());
    }
    for (auto& track : assocTracksGrouped) { //<- only for a subset
      registry.get<TH1>(HIST("etaHistogramAssoc"))->Fill(track.eta()); //<- this should show the selection
      registry.get<TH1>(HIST("ptHistogramAssoc"))->Fill(track.pt());
    }
    
    //Now we do two-particle correlations, using "combinations"
    for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(triggerTracksGrouped, assocTracksGrouped))) {
      if(trackTrigger.tpcNClsCrossedRows() < 70 ) continue; //can't filter on dynamic
      if(trackAssoc.tpcNClsCrossedRows() < 70 ) continue; //can't filter on dynamic
      registry.get<TH1>(HIST("correlationFunction"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
      registry.get<TH1>(HIST("hDeltaEta"))->Fill( trackTrigger.eta()-trackAssoc.eta() );
      registry.get<TH1>(HIST("hDeltaPhi"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
      registry.get<TH2>(HIST("hDeltaPhiDeltaEta"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ), trackTrigger.eta()-trackAssoc.eta());
      
      registry.get<TH1>(HIST("hTPcosDeltaPhi"))->Fill(TMath::Cos(2.0*ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi())));
      
    } 

  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LaurasCode>(cfgc)
  };
}
