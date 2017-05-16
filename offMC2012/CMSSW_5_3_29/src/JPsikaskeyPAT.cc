// -*- C++ -*-
//
// Package:    JPsiomegaPAT
// Class:      JPsiomegaPAT

//
//  Author:  Jhovanny Andres Mejia, Eduard de la Cruz Burelo
//        
//
//


// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/interface/JPsikaskeyPAT.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <utility>
#include <string>
//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

//
// static data member definitions
//

//
// constructors and destructor
//
JPsikaskeyPAT::JPsikaskeyPAT(const edm::ParameterSet& iConfig)
  :
  hltTriggerResults_ (iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT")) ),
  vtxSample ( iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices")) ),
  genParticles_ ( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  v0Collection_ ( iConfig.getUntrackedParameter<std::string>("V0Collection",std::string("generalV0Candidates")) ),
  muonType ( iConfig.getUntrackedParameter<std::string>("MuonType",std::string("cleanPatMuons")) ),
  muonTypeForPAT ( iConfig.getUntrackedParameter<std::string>("MuonTypeForPAT",std::string("muons")) ),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  // onlyCount_ ( iConfig.getUntrackedParameter<bool>("onlyCount",false) ),
  tree_(0), 

  priRfVtxX(0), priRfVtxY(0), priRfVtxZ(0), priRfVtxXE(0), priRfVtxYE(0), priRfVtxZE(0), priRfVtxCL(0),
  priRfVtxXYE(0), priRfVtxXZE(0), priRfVtxYZE(0),
  priRfNTrkDif(0),
  
  /*
  bctau(0), bctau2D(0), bctauBS(0), bctauBS2D(0), bctauRf(0), bctau2DRf(0),
  bctauE(0), bctau2DE(0), bctauBSE(0), bctauBS2DE(0), bctauRfE(0), bctau2DRfE(0),
  bctau_kaskey(0), bctau2D_kaskey(0),
  bctauE_kaskey(0), bctau2DE_kaskey(0),
  bctau_lam(0), bctau2D_lam(0),
  bctauE_lam(0), bctau2DE_lam(0),
  */

  mumC2(0), mumCat(0), mumAngT(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupCat(0), mupAngT(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  priVtxXBS(0), priVtxYBS(0), priVtxZBS(0), priVtxXBSE(0), priVtxYBSE(0), priVtxZBSE(0), priVtxCLBS(0),
  priVtxXYBSE(0), priVtxXZBSE(0), priVtxYZBSE(0),

  // todos los vertices primarios 
  //pVtxX(0),   pVtxY(0),   pVtxZ(0),   pVtxXE(0),   pVtxYE(0),   pVtxZE(0),   pVtxCL(0),
  //pVtxXYE(0), pVtxXZE(0), pVtxYZE(0),

  // vertice primario CON mejor pointin-angle
  pVtxIPX(0),   pVtxIPY(0),   pVtxIPZ(0), pVtxIPXE(0),   pVtxIPYE(0),   pVtxIPZE(0), pVtxIPCL(0),
  pVtxIPXYE(0),   pVtxIPXZE(0),   pVtxIPYZE(0),

  // todos los vertices primarios CON constrain de Beamspot
  //pVtxBSX(0),   pVtxBSY(0),   pVtxBSZ(0),   pVtxBSXE(0),   pVtxBSYE(0),   pVtxBSZE(0),   pVtxBSCL(0),
  //pVtxBSXYE(0), pVtxBSXZE(0), pVtxBSYZE(0),

  // vertice primario CON constrain de Beamspot y mejor pointin-angle
  pVtxBSIPX(0),   pVtxBSIPY(0),   pVtxBSIPZ(0), pVtxBSIPXE(0),   pVtxBSIPYE(0),   pVtxBSIPZE(0), pVtxBSIPCL(0),
  pVtxBSIPXYE(0),   pVtxBSIPXZE(0),   pVtxBSIPYZE(0),

 // ************ esta es la informacion concerniente al beamspot *************
  PVXBS(0), PVYBS(0), PVZBS(0), PVXBSE(0), PVYBSE(0), PVZBSE(0),
  PVXYBSE(0), PVXZBSE(0), PVYZBSE(0),

  // ************************ ****************************************************

  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),

  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0),  VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 

  V1DecayVtxX(0), V1DecayVtxY(0), V1DecayVtxZ(0), V1DecayVtxXE(0), V1DecayVtxYE(0), V1DecayVtxZE(0), 
  V1DecayVtxXYE(0), V1DecayVtxXZE(0), V1DecayVtxYZE(0),

  JDecayVtxX(0), JDecayVtxY(0), JDecayVtxZ(0), JDecayVtxXE(0), JDecayVtxYE(0), JDecayVtxZE(0), 
  JDecayVtxXYE(0), JDecayVtxXZE(0), JDecayVtxYZE(0),

  // *******************************************************
 
  nB(0), nMu(0),

  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_kaskeym_mass(0), B_kaskeym_px(0), B_kaskeym_py(0), B_kaskeym_pz(0), B_kaskeym_charge(0),
  
  B_kaskeym_pt1(0), B_kaskeym_px1(0), B_kaskeym_py1(0), B_kaskeym_pz1(0), 
  B_kaskeym_pt2(0), B_kaskeym_px2(0), B_kaskeym_py2(0), B_kaskeym_pz2(0), 
  B_kaskeym_charge1(0), B_kaskeym_charge2(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0), 
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0),B_J_charge2(0),

  B_kas_lambda_mass(0), B_kas_lambda_px(0), B_kas_lambda_py(0), B_kas_lambda_pz(0),
  B_kas_lambda_pt1(0), B_kas_lambda_px1(0), B_kas_lambda_py1(0), B_kas_lambda_pz1(0), 
  B_kas_lambda_pt2(0), B_kas_lambda_px2(0), B_kas_lambda_py2(0), B_kas_lambda_pz2(0), 
  B_kas_lambda_charge1(0), B_kas_lambda_charge2(0),

  B_J_parentId1(0), B_J_parentId2(0),  B_J_muId1(0), B_J_muId2(0),
  B_lam_parentId1(0), B_lam_parentId2(0), B_lam_PId1(0), B_lam_piId2(0),
  B_kaskeym_kId3(0),
  B_kaskey_parentId1(0),B_kaskey_parentId2(0),B_kaskey_parentId3(0),
  B_parentId1(0), B_parentId2(0),B_parentId3(0), B_parentId4(0), B_parentId5(0),
  //B_kaskeym_parentId1(0), B_kaskeym_parentId2(0),
  //B_kaskeym_pId1(0), B_kaskeym_pId2(0),


  //B_kas_lambda_trkP_d0(0), B_kas_lambda_trkP_d0E(0),B_kas_lambda_trkM_d0(0),B_kas_lambda_trkM_d0E(0), 
  //B_kas_pion_d0(0), B_kas_pion_d0E(0),

  B_kaskeym_chi2(0), B_J_chi2(0), B_chi2(0), B_kas_lambda_chi2(0),
  B_Prob(0), B_J_Prob(0), B_kaskey_Prob(0), B_kas_lambda_Prob(0),


  triggersMuPL(0), triggersMuML(0), 
  triggersL1L2_MuPL(0), triggersL1L2_MuML(0),
  

  run(0), event(0)

{
   //now do what ever initialization is needed
}


JPsikaskeyPAT::~JPsikaskeyPAT()
{

}


//
// member functions
//

void JPsikaskeyPAT::CheckHLTTriggers(const std::vector<std::string>& TrigList){

    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace helper;
    
    
    string AllTrg="";
    string tmptrig;

    int ntrigs=TrigList.size();
    if (ntrigs==0)
        std::cout << "No trigger name given in TriggerResults of the input " << endl;
    
    for (int itrig=0; itrig< ntrigs; itrig++) {
        //TString trigName = triggerNames_.triggerName(itrig);
        string trigName = TrigList.at(itrig);
    
	     tmptrig = (string) trigName; tmptrig +=" ";
	     AllTrg += tmptrig;
    }

   int m = sprintf(triggersL,"%s","");
   m = sprintf(triggersL,"%s",AllTrg.c_str());
   //cout<<" INFO: Triggers :  "<<m<<"  "<<n<<"  "<<triggersL<<endl;

   //nTrgL = AllTrg.size();
   nTrgL = m;
   //cout<<" INFO: Triggers :  "<<m<<"  "<<nTrgL<<"  "<<triggersL<<endl; 

   return;
}




void JPsikaskeyPAT::CheckL1Triggers(const edm::Event& iEvent, const edm::EventSetup& iSetup,  std::string &TrigListNameL1Tmp)
{
  using namespace std;
   // get L1 trigger info
    
     edm::ESHandle<L1GtTriggerMenu> menuRcd;
     iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
     const L1GtTriggerMenu* menu = menuRcd.product();
   
     edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
     iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
     const DecisionWord dWord = gtRecord->decisionWord();
    
    string AllTrgL1="";
    string tmptrigL1;

    for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
        
        string trigNameL1 = (algo->second).algoName();
        
        tmptrigL1 = (string) trigNameL1; tmptrigL1 +=" ";
        AllTrgL1 += tmptrigL1;
        
        //cout << "Name: " << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias() << std::endl;
    }
    
    TrigListNameL1Tmp   = AllTrgL1.c_str();

    //  if ( menu->gtAlgorithmResult( "L1_SingleMu3", dWord) )  l1_mu3 = 1;

     return;   
}



// aca se hace el machint si se corre en MonteCarlo
/*
int JPsikaskeyPAT::PdgIDatTruthLevel(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles, int &ParentID, int &PParentID)
 {
   double pXTrack = Track->momentum().x();
   double pYTrack = Track->momentum().y();
   double pZTrack = Track->momentum().z();

   //std::cout<<" *** "<<pXTrack<<" -- "<<pYTrack<<" --  "<<pZTrack<<std::endl;

   std::vector<double> chi2v;
   std::vector<int> idPDG;
   std::vector<int> idPDGP;
   std::vector<int> idPDGPP;

   for( size_t k = 0; k < genParticles->size(); k++ ) 
     {

       reco::GenParticle PCand =  (*genParticles)[k];
       double truePx = PCand.px();
       double truePy = PCand.py();
       double truePz = PCand.pz();
       //if(PCand.pdgId()==fabs(211)) std::cout<<truePx<<"  "<<truePy<<"  "<<truePz<<std::endl;
             double dpx = pXTrack - truePx;
       double dpy = pYTrack - truePy;
       double dpz = pZTrack - truePz;
       double chi2 = dpx*dpx + dpy*dpy + dpz*dpz;
       //if(chi2<1.) std::cout<<chi2<<"  id: "<<PCand.pdgId()<<std::endl;
      //if(chi2<1.) std::cout<<truePx<<"  "<<truePy<<"  "<<truePz<<" id: "<<PCand.pdgId()<<std::endl;

        if(chi2>10.0)continue;

       if(PCand.numberOfMothers()!=1) continue;
       const reco::Candidate * PCandM = PCand.mother();

       if(PCandM->numberOfMothers()!=1) continue;
       const reco::Candidate * PCandMM = PCandM->mother();
       //std::cout<<PCandM->pdgId()<<std::endl;

       //B_kaskeym_chi2->push_back(chi2);
       chi2v.push_back(chi2);
       idPDG.push_back(PCand.pdgId());
       idPDGP.push_back(PCandM->pdgId());
       idPDGPP.push_back(PCandMM->pdgId());

     }

   double chi2min_temp=100000;
   //double chi2min=-1;
   int idtemp = -1;
   int idptemp = -1;
   int idpptemp = -1;
   for(int i=0;i<(int)chi2v.size();i++)
     {
       if(chi2min_temp>chi2v.at(i))
 	{
 	  chi2min_temp = chi2v.at(i);
 	  idtemp = idPDG.at(i);
 	  idptemp = idPDGP.at(i);
	  idpptemp = idPDGPP.at(i);
 	}
     }
 
   //std::cout<<chi2min_temp<<"  id: "<<idtemp<<std::endl;
   ParentID = idptemp;
   PParentID = idpptemp;
   return idtemp;
 }

*/

// aca es funcion de ID para el proton y pion del lambda (esto porque son tres padres)

/*
int JPsikaskeyPAT::PdgIDatTruthLevel4(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles, int &ParentID, int &PParentID, int &GPParentID)
 {
   double pXTrack = Track->momentum().x();
   double pYTrack = Track->momentum().y();
   double pZTrack = Track->momentum().z();

   //std::cout<<" *** "<<pXTrack<<" -- "<<pYTrack<<" --  "<<pZTrack<<std::endl;

   std::vector<double> chi2v;
   std::vector<int> idPDG;
   std::vector<int> idPDGP;
   std::vector<int> idPDGPP;
   std::vector<int> idPDGPPP;
   
for( size_t k = 0; k < genParticles->size(); k++ ) 
     {

       reco::GenParticle PCand =  (*genParticles)[k];
       double truePx = PCand.px();
       double truePy = PCand.py();
       double truePz = PCand.pz();
       //if(PCand.pdgId()==fabs(211)) std::cout<<truePx<<"  "<<truePy<<"  "<<truePz<<std::endl;
             double dpx = pXTrack - truePx;
       double dpy = pYTrack - truePy;
       double dpz = pZTrack - truePz;
       double chi2 = dpx*dpx + dpy*dpy + dpz*dpz;
       //if(chi2<1.) std::cout<<chi2<<"  id: "<<PCand.pdgId()<<std::endl;
      //if(chi2<1.) std::cout<<truePx<<"  "<<truePy<<"  "<<truePz<<" id: "<<PCand.pdgId()<<std::endl;

        if(chi2>10.0)continue;

       if(PCand.numberOfMothers()!=1) continue;
       const reco::Candidate * PCandM = PCand.mother();

       if(PCandM->numberOfMothers()!=1) continue;
       const reco::Candidate * PCandMM = PCandM->mother();
       //std::cout<<PCandM->pdgId()<<std::endl;


       if(PCandMM->numberOfMothers()!=1) continue;
       const reco::Candidate * PCandMMM = PCandMM->mother();

       //B_kaskeym_chi2->push_back(chi2);
       chi2v.push_back(chi2);
       idPDG.push_back(PCand.pdgId());
       idPDGP.push_back(PCandM->pdgId());
       idPDGPP.push_back(PCandMM->pdgId());
       idPDGPPP.push_back(PCandMMM->pdgId());

     }

   double chi2min_temp=100000;
   //double chi2min=-1;
   int idtemp = -1;
   int idptemp = -1;
   int idpptemp = -1;
   int idppptemp = -1;

   for(int i=0;i<(int)chi2v.size();i++)
     {
       if(chi2min_temp>chi2v.at(i))
 	{
 	  chi2min_temp = chi2v.at(i);
 	  idtemp = idPDG.at(i);
 	  idptemp = idPDGP.at(i);
	  idpptemp = idPDGPP.at(i);
	  idppptemp = idPDGPPP.at(i);

 	}
     }
 
   //std::cout<<chi2min_temp<<"  id: "<<idtemp<<std::endl;
   ParentID = idptemp;
   PParentID = idpptemp;
   GPParentID = idppptemp;
   return idtemp;
 }
*/


// sta la funcion que hace el machint de los muones

void  JPsikaskeyPAT::MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp){
    
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace helper;
 
    
  string AllTrg="";
  string tmptrig;

  int ntrigs=TrigList.size();
  if (ntrigs==0)
    std::cout << "No trigger name given in TriggerResults of the input " << endl;
      
//   cout<<" Trigger in Muon List :";
	
//   for (int itrig1=0; itrig1< ntrigs; itrig1++) {
//     string trigName = TrigList.at(itrig1);
//     cout<<" "<< trigName; 
//   }
//   cout<<endl;

//   const pat::TriggerObjectStandAloneCollection    muMatch = iMuon.triggerObjectMatches();
//    cout<<"*********************************** "<< endl;
//    cout<<"****** "<<muMatch.size()<< " " <<endl; 
//     for (unsigned int kt = 0; kt < muMatch.size(); kt++){
//       std::vector<std::string>  pathLabels   = muMatch[kt].pathNames(1,0);

//       for(unsigned int j=0;j<pathLabels.size();j++){
// 	//cout<<"*********************************** "<<pathLabels[j].c_str()<<endl;
// 	cout<<" "<<pathLabels.at(j);
//       }
//     }
//     cout<< endl;

  for (int itrig=0; itrig< ntrigs; itrig++) {

      string trigName = TrigList.at(itrig);

      if (iMuon.triggerObjectMatchesByPath(trigName.c_str(),(unsigned int)1,(unsigned int)0).empty()==false){
	   // cout<<"empty" << endl;
	    //continue;
	// cout<<"In muon: "<<trigName << " " << endl;
      tmptrig = (string) trigName; tmptrig +=" ";
      AllTrg += tmptrig;
      }

  }
		     
    //int nMuonTrgLtmp  = AllTrg.size();
    TrigListNameTmp   = AllTrg.c_str();
    

  //////////////////////////
}


void  JPsikaskeyPAT::MatchMuonWithL1L2(const pat::Muon &iMuon, const std::vector<std::string>& TrigListL1L2, std::string &TrigListNameL1L2Tmp){
    //void JPsiOmPATv2::MatchMuonWithTriggers(const pat::Muon &iMuon, const edm::Event& iEvent, const edm::EventSetup& iSetup){
    
    
    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace helper;
    
    
    string AllTrgL1L2="";
    string tmptrigL1L2;
    
    int ntrigsL1L2 = TrigListL1L2.size();

    
    for (int itrigL1L2=0; itrigL1L2< ntrigsL1L2; itrigL1L2++) {
        
        string trigNameL1L2 = TrigListL1L2.at(itrigL1L2);
        
        //cout<<"i = "<<itrig<<"  "<<trigName.c_str() <<" has " << iMuon.triggerObjectMatchesByPath(trigName.c_str(),(unsigned int)0,(unsigned int)0).empty() << endl;
        
        if (iMuon.triggerObjectMatchesByFilter(trigNameL1L2.c_str()).empty()==false){
            //continue;
            tmptrigL1L2 = (string) trigNameL1L2; tmptrigL1L2 +=" ";
            AllTrgL1L2 += tmptrigL1L2;
        }
        
    }
    
    //int nMuonTrgLtmp  = AllTrg.size();
    TrigListNameL1L2Tmp   = AllTrgL1L2.c_str();
    //cout<< "L1/L2: " <<AllTrgL1L2.c_str() << endl;
    
    //return nMuonTrgLtmp;
    //////////////////////////
}





////////////////
/*
float JPsikaskeyPAT::Myctau(const RefCountedKinematicParticle &CandMC, const RefCountedKinematicVertex &DecayVertexMC, 
			const GlobalPoint &PVPtmp, const GlobalError &PVEtmp, float mass_tmp, 
			float & ctau2Dtmp, float &ctauEtemp, float &ctauE2Dtemp){
  using std::vector;
  using namespace edm;
  using namespace reco;

  // measure 3D ctau
  float betagamma = ( CandMC->currentState().globalMomentum().mag()/ mass_tmp );
  float bPt = sqrt(CandMC->currentState().globalMomentum().x()*CandMC->currentState().globalMomentum().x()+
		   CandMC->currentState().globalMomentum().y()*CandMC->currentState().globalMomentum().y() );
  float betagammaT = bPt/mass_tmp  ;
  
  GlobalPoint BVP = GlobalPoint( DecayVertexMC->position() );
  
  GlobalVector sep3D = BVP-PVPtmp;
  
  GlobalVector pBV = CandMC->currentState().globalMomentum();
  float ctau_temp = (mass_tmp  *(sep3D.dot(pBV)))/(pBV.dot(pBV));

  // measure 2D ctau
  GlobalVector pTBV = GlobalVector( CandMC->currentState().globalMomentum().x(), CandMC->currentState().globalMomentum().y(), 0);	       
  float ctau2D_temp = (mass_tmp*(sep3D.dot(pTBV)))/(pTBV.dot(pTBV));

  //ctau2Dtmp = ctau2D_temp;


  // measure 3D Error

  // calculate ctau error.
  // Momentum error is negligible compared to the vertex errors
  GlobalError BVE = DecayVertexMC->error();
  VertexDistance3D theVertexDistance3D; 

  Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVPtmp, PVEtmp) );
  double myError = TheMeasurement.error();	       
  
  //  ctau is defined by the portion of the flight distance along the compoenent of the B momementum, so only
  // consider the error of that component, too, which is accomplished by scaling by ((VB-VP)(dot)PB)/|VB-VP|*|PB|
	       
  float scale = fabs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );    	       
  float ctauE_temp =  (myError*scale)/betagamma;


  // measure 2D Error
  VertexDistanceXY theVertexDistanceXY; 
  Measurement1D TheMeasurementXY = theVertexDistanceXY.distance( VertexState(BVP, BVE), VertexState(PVPtmp, PVEtmp) );
  double myErrorXY = TheMeasurementXY.error();	       
  GlobalVector sep2D = GlobalVector( DecayVertexMC->position().x()-PVPtmp.x(), DecayVertexMC->position().y()-PVPtmp.y(), 0);

  float scaleXY = fabs( (sep2D.dot(pTBV))/(sep2D.mag()*pTBV.mag()) );

  float ctauE2D_temp = (myErrorXY*scaleXY)/betagammaT;

  ctau2Dtmp   = ctau2D_temp;
  ctauEtemp   = ctauE_temp;
  ctauE2Dtemp = ctauE2D_temp;


  return ctau_temp;
  
}
*/


// ------------ method called to for each event  ------------
void JPsikaskeyPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;


  // bool debug = true;

  
  //*********************************
  // Get event content information
  //*********************************  
  Handle<vector<VertexCompositeCandidate> > theV0Handle;
  iEvent.getByLabel(v0Collection_, "Lambda", theV0Handle);

  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  Handle< vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);

  Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel(muonType, thePATMuonHandle);
  //iEvent.getByLabel("cleanPatMuonsTriggerMatch", thePATMuonHandle);

 

 //get genParticles  () esto tambin es para el Monte Carlo
  /*  
 Handle<GenParticleCollection> genParticles;
  if (doMC_) 
    {
   iEvent.getByLabel(genParticles_, genParticles);
         }
  */




 // Get HLT results
  edm::Handle<edm::TriggerResults> hltresults1;
  try {
    std::string const &trig = std::string("TriggerResults::")+hltTriggerResults_;
    iEvent.getByLabel(edm::InputTag(trig),hltresults1);
  }
  catch ( ... ) 
    {
      std::cout << "Couldn't get handle on HLT Trigger!" << endl;
    }
    
    //HLTConfigProvider hltConfig_;
    
    std::vector<std::string> TrigTable; TrigTable.clear();
    // Get hold of trigger names - based on TriggerResults object
    const edm::TriggerNames& triggerNames1_ = iEvent.triggerNames(*hltresults1);
    
    for (unsigned int itrig = 0; itrig < hltresults1->size(); ++itrig){
        if ((*hltresults1)[itrig].accept() == 1){
            std::string trigName1 = triggerNames1_.triggerName(itrig);
            //int trigPrescale = hltConfig_.prescaleValue(itrig, trigName1);
            TrigTable.push_back(trigName1);
        }
    }

    /*
    cout<< "Trigger table: ";
    for( unsigned int i = 0; i<TrigTable.size(); ++i)
        cout<<TrigTable.at(i) << " ";
    
    cout<< endl;
*/
    

  
    
    
    std::vector<string> ListTriggerL1L2;
    
    ListTriggerL1L2.push_back("hltL1MuOpenL1Filtered0");
    ListTriggerL1L2.push_back("hltL2Mu0L2Filtered0");
    ListTriggerL1L2.push_back("hltSingleMu3L2Filtered3");
    ListTriggerL1L2.push_back("hltSingleMu3L3Filtered3");
    ListTriggerL1L2.push_back("hltSingleMu5L3Filtered5");
    
    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    ListTriggerL1L2.push_back("hltDiMuonL2PreFiltered0");
    ListTriggerL1L2.push_back("hltDiMuonL3PreFiltered0");
    ListTriggerL1L2.push_back("hltDiMuonL3PreFiltered");
    ListTriggerL1L2.push_back("hltMu0L1MuOpenL3Filtered0");
    
    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    ListTriggerL1L2.push_back("hltMu3L1MuOpenL3Filtered3");
    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    ListTriggerL1L2.push_back("hltMu5L1MuOpenL3Filtered5");
    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    
    ListTriggerL1L2.push_back("hltL2Mu0L2Filtered0");
    ListTriggerL1L2.push_back("hltMu0TrackJpsiTrackMassFiltered");
    ListTriggerL1L2.push_back("hltMu3TrackJpsiTrackMassFiltered");
    ListTriggerL1L2.push_back("hltMu5TrackJpsiTrackMassFiltered");

    /////////////////////////////////////
    




  //*********************************
  //Now we get the primary vertex 
  //*********************************


  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;

  // get primary vertex
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(vtxSample, recVtxs);
  //unsigned int nVtxTrks = 0;
  bestVtx = *(recVtxs->begin());
  
  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  //priVtxXE = bestVtx.xError();
  //priVtxYE = bestVtx.yError();
  //priVtxZE = bestVtx.zError();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 

  nVtx = recVtxs->size();  
  
  //get primary with beamspot constraint
  Handle<reco::VertexCollection> recVtxsBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", recVtxsBS);
  
  // nVtxTrks = 0;
  bestVtxBS = *(recVtxsBS->begin());

  priVtxXBS = bestVtxBS.x();
  priVtxYBS = bestVtxBS.y();
  priVtxZBS = bestVtxBS.z();
  priVtxXBSE = bestVtxBS.covariance(0, 0);
  priVtxYBSE = bestVtxBS.covariance(1, 1);
  priVtxZBSE = bestVtxBS.covariance(2, 2);
  priVtxXYBSE = bestVtxBS.covariance(0, 1);
  priVtxXZBSE = bestVtxBS.covariance(0, 2);
  priVtxYZBSE = bestVtxBS.covariance(1, 2);  
  priVtxCLBS = ChiSquaredProbability((double)(bestVtxBS.chi2()),(double)(bestVtxBS.ndof())); 


 
  //*********************************
  //Now we get the Beam Spot  
  //*********************************

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
  else std::cout << "No beam spot available from EventSetup" << endl;
 /*
  else
  {
    edm::LogInfo("MyAnalyzer")
      << "No beam spot available from EventSetup \n";
}
  */

  PVXBS = beamSpot.x0();
  PVYBS = beamSpot.y0();
  PVZBS = beamSpot.z0();
  PVXBSE = beamSpot.covariance()(0, 0);
  PVYBSE = beamSpot.covariance()(1, 1);
  PVZBSE = beamSpot.covariance()(2, 2);
  PVXYBSE = beamSpot.covariance()(0,1);
  PVXZBSE = beamSpot.covariance()(0,2);
  PVYZBSE = beamSpot.covariance()(1,2);

 
  //*****************************
  // Let's check triggers
  //*****************************
  //CheckHLTTriggers(iEvent,iSetup);
  //CheckL1Triggers(iEvent,iSetup);

  // Do we want some cuts on triggers?

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-


 int run1   =  iEvent.id().run();
  int event1 =  iEvent.id().event();


 

  // aca empiezan los for para cada uno de lso track que vamos a recosntruir

  unsigned int nMu_tmp = thePATMuonHandle->size();
 
 for( std::vector<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for( std::vector<pat::Muon>::const_iterator iMuon2 = iMuon1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------

	  const pat::Muon *patMuonP  = 0;
	  const pat::Muon *patMuonM  = 0;
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ patMuonP = &(*iMuon1); glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){patMuonM = &(*iMuon1); glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {patMuonP = &(*iMuon2); glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){patMuonM = &(*iMuon2); glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  //Let's check the vertex and mass

	  if(iMuon1->track()->pt()<3.0) continue;
	  if(iMuon2->track()->pt()<3.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

	  TransientTrack muon1TT(glbTrackP, &(*bFieldHandle) );
	  TransientTrack muon2TT(glbTrackM, &(*bFieldHandle) );

	  // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************
	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;   

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }

	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }

	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	   if(psi_vFit_noMC->currentState().mass()<2.95 || psi_vFit_noMC->currentState().mass()>3.25) continue;
	   
	   //More quality cuts?


//ahora veamos los  for para formar el kaskey menos y el lambda usando el contenedor V0


 if ( theV0Handle->size()>0 && thePATMuonHandle->size()>=2 )
    {

     for ( vector<VertexCompositeCandidate>::const_iterator iVee = theV0Handle->begin();   iVee != theV0Handle->end(); ++iVee )
        {
	  //get Lam tracks from V0 candidate
	  vector<RecoChargedCandidate> v0daughters;
	  vector<TrackRef> theDaughterTracks;
	  v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
       			(iVee->daughter(0))) );
	  v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
       			(iVee->daughter(1))) );


	  for(unsigned int j = 0; j < v0daughters.size(); ++j)
	    {
	      theDaughterTracks.push_back(v0daughters[j].track());
	    }

	  // *********************
	  // aca empieza el for para "machar" los traks del contenedor v0 con el contenedor de los traks
	  //para el caso seran el pion y el proton que forman el lambda
	  // ******************
       pat::GenericParticle patTrack1;
       pat::GenericParticle patTrack2;
       pat::GenericParticle patTrack3;

       //lets loop through the pat tracks to find the match for this reco::track
      for (vector<pat::GenericParticle>::const_iterator iTrack = thePATTrackHandle->begin();iTrack != thePATTrackHandle->end(); ++iTrack )
         {

     // how to access the reco::Track object
	 TrackRef hadTrack = iTrack->track();
	 if ( hadTrack.isNull() ) {
	   cout << "continue due to no track ref" << endl;
	   continue;
	 }
	 if ( hadTrack==theDaughterTracks[0] )
	   patTrack1 = *iTrack;
	 if ( hadTrack==theDaughterTracks[1] )
	   patTrack2 = *iTrack;
	 }

          // aca esta el tercer for que es para el kaon que pegaremos al lambda para formar el kaskey menos       





     for(std::vector<pat::GenericParticle>::const_iterator iTrack1 = thePATTrackHandle->begin(); 
            iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) 
       {
		  


     
	 if (iTrack1->track()==theDaughterTracks[0] || iTrack1->track()==theDaughterTracks[1] )
	   {
	     continue;
	   }

	 //if(iTrack1->track()->charge() != -1) continue;
	 

	 patTrack3 = *iTrack1;

	 if(!(patTrack3.track()->quality(reco::TrackBase::highPurity))) continue;
	 //if(!(patTrack2.track()->quality(reco::TrackBase::highPurity))) continue;
	 //if(!(patTrack3.track()->quality(reco::TrackBase::highPurity))) continue;


 //Now let's checks if our muons do not use the same tracks as we are using now

		   bool matchflag = false;

		   const reco::CandidatePtrVector & mu1P_overlaps = patTrack1.overlaps(muonTypeForPAT);
		   if ( mu1P_overlaps.size() > 0 ) //std::cout << "patTrack1 overlaps with a muon." << endl;
		   for (size_t i = 0; i < mu1P_overlaps.size(); ++i) {
		     const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu1P_overlaps[i]);
		     if (mu) {
		       // check here that muon match isn't the same as a muon used in the reco...
		       if (mu==patMuonP || mu==patMuonM) 
			 {
			   //std::cout << "match between patTrack1 and patMuonP/patMuonM" << endl;
			   matchflag=true;
			 }
		     }
		   }
	       
		   const reco::CandidatePtrVector & mu2P_overlaps = patTrack2.overlaps(muonTypeForPAT);
		   if ( mu2P_overlaps.size() > 0 ) //std::cout << "patTrack2 overlaps with a muon." << endl;
		   for (size_t i = 0; i < mu2P_overlaps.size(); ++i) {
		     const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu2P_overlaps[i]);
		     if (mu) {
		       // check here that muon match isn't the same as a muon used in the reco...
		       if (mu==patMuonP || mu==patMuonM)
			 { 
			   //std::cout << "match between patTrack2 and patMuonP/patMuonM" << endl;
			   matchflag = true;
			 }
		     }
		   }
	
	    const reco::CandidatePtrVector & mu3P_overlaps = patTrack3.overlaps(muonTypeForPAT);
		   if ( mu3P_overlaps.size() > 0 ) //std::cout << "patTrack3 overlaps with a muon." << endl;
		   for (size_t i = 0; i < mu3P_overlaps.size(); ++i) {
		     const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu3P_overlaps[i]);
		     if (mu) {
		       // check here that muon match isn't the same as a muon used in the reco...
		       if (mu==patMuonP || mu==patMuonM)
			 { 
			   //std::cout << "match between patTrack2 and patMuonP/patMuonM" << endl;
			   matchflag = true;
			 }
		     }
		   }

		   if(matchflag) continue;



    
		   // **********
		   //  ahora tenemos sinco traks incluyendo dos muones de cargas opuestas y un kaon negativo
		   //  los otros dos son los que forman el lambda

	      // check which daughter track has high momentum and call it the proton
	       int indexP = 0; int indexPi = 1;
	       if( theDaughterTracks[1]->p() > theDaughterTracks[0]->p() ) {
	         indexP = 1;
		 indexPi = 0;
	       } 

	       TransientTrack pTT(theDaughterTracks[indexP], &(*bFieldHandle) );
	       TransientTrack piTT(theDaughterTracks[indexPi], &(*bFieldHandle) );
	       TransientTrack kaMTT(patTrack3.track(), &(*bFieldHandle) );
	      
	       
	       
	       //The mass of a muon and the insignificant mass sigma 
	       //to avoid singularities in the covariance matrix.
	       ParticleMass pion_mass = 0.13957018;
	       //ParticleMass kaon_mass = 0.493677;
	       ParticleMass p_mass = 0.938272;	       
	       ParticleMass lam_mass = 1.115683;
	       ParticleMass kaskeym_mass = 1.32171;
	       float pion_sigma = pion_mass*1.e-6;
	       //float kaon_sigma = kaon_mass*1.e-6;
	       float p_sigma = p_mass*1.e-6;
	       float lam_sigma = 0.000006;
	       float kaskeym_sigma = 0.00007;
	       
	       
	      
	     
	     
	       
	       // primero el vertice para el lambda


	       


	       vector<RefCountedKinematicParticle> pionParticles;
	       try {
		   pionParticles.push_back(pFactory.particle(pTT,p_mass,chi,ndf,p_sigma));
		   pionParticles.push_back(pFactory.particle(piTT,pion_mass,chi,ndf,pion_sigma));
	       }

	       catch(...) { 
		 std::cout<<" Exception caught ... continuing 3 "<<std::endl; 
		 continue;
	       }


	      
	       KinematicParticleVertexFitter fitter;  
 
	       RefCountedKinematicTree lamVertexFitTree;
	       try {
		     lamVertexFitTree = fitter.fit(pionParticles); 
	       }
	       catch (...) { 
		 std::cout<<" Exception caught ... continuing 4 "<<std::endl; 
		 continue;
	       }

	       if (!lamVertexFitTree->isValid()) {
		 //std::cout << "invalid vertex from the lam vertex fit" << std::endl;
		 continue; 
	       }
	       lamVertexFitTree->movePointerToTheTop();


      
	       RefCountedKinematicParticle lam_vFit_noMC = lamVertexFitTree->currentParticle();
	       RefCountedKinematicVertex lam_vFit_vertex_noMC = lamVertexFitTree->currentDecayVertex();

	       if( lam_vFit_vertex_noMC->chiSquared() < 0 )
		     { 
		       //std::cout << "negative chisq from ks fit" << endl;
		       continue;
		     }

	       if(lam_vFit_noMC->currentState().mass()<1.085683 || lam_vFit_noMC->currentState().mass()>1.145683) continue;

	       //if ( lam_vFit_vertex_noMC->chiSquared() < 0 ) cout << "negative chisq from lam fit" << endl;	       
	       lamVertexFitTree->movePointerToTheFirstChild();
	       RefCountedKinematicParticle lamTrk1 = lamVertexFitTree->currentParticle();
	       lamVertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle lamTrk2 = lamVertexFitTree->currentParticle();
	       lamVertexFitTree->movePointerToTheTop();

	       

	       // Do MC for Lam cand and do mass constrained vertex fit
	       // creating the constraint with a small sigma to put in the resulting covariance 
	       // matrix in order to avoid singularities

	       // aca hago la constricion de la masa para el lambda

	       KinematicParticleFitter csFitterLam;
	       KinematicConstraint * lam_c = new MassKinematicConstraint(lam_mass,lam_sigma);
	       // add mass constraint to the lam fit to do a constrained fit:  
 
	       lamVertexFitTree = csFitterLam.fit(lam_c,lamVertexFitTree);
	       if (!lamVertexFitTree->isValid()){
		 //std::cout << "caught an exception in the lam mass constraint fit" << std::endl;
		 continue; 
	       }
	       
	       lamVertexFitTree->movePointerToTheTop();
	       RefCountedKinematicParticle lam_vFit_withMC = lamVertexFitTree->currentParticle();


	       

// aca estara el vertice del kaskey menos formado por el lambda y el otro pion
	       vector<RefCountedKinematicParticle> LambdaParticles;	       
	       try {
	           LambdaParticles.push_back(pFactory.particle(kaMTT,pion_mass ,chi,ndf,pion_sigma));
		   LambdaParticles.push_back(lam_vFit_withMC);
	       }
	       catch(...) { 
		          std::cout<<" Exception caught ... continuing 5 "<<std::endl; 
	                  continue;
	                }
	       RefCountedKinematicTree kaskeymenosFitTree;
		try {
		      kaskeymenosFitTree = fitter.fit(LambdaParticles); 
		    }

		catch (...) { 
		 std::cout<<" Exception caught ... continuing 6 "<<std::endl; 
		 continue;
	       }

		if (!kaskeymenosFitTree->isValid()) {
		  //std::cout << "invalid vertex from the lam vertex fit" << std::endl;
		 continue; 
		  }
		kaskeymenosFitTree->movePointerToTheTop();



		RefCountedKinematicParticle kaskeymenos_vFit_noMC = kaskeymenosFitTree->currentParticle();
		RefCountedKinematicVertex kaskeymenos_vFit_vertex_noMC = kaskeymenosFitTree->currentDecayVertex();

		   if( kaskeymenos_vFit_vertex_noMC->chiSquared() < 0 )
		     { 
		       //std::cout << "negative chisq from ks fit" << endl;
		       continue;
		     }

		   //some loose cuts go here
		   
		   //cout<<"Mass 1 : "<<kaskey_vFit_noMC->currentState().mass()<<endl;

		   if(kaskeymenos_vFit_vertex_noMC->chiSquared()>50) continue;
	     //if(kaskeymenos_vFit_noMC->currentState().mass()<1.62245 || kaskeymenos_vFit_noMC->currentState().mass()>1.72245) continue;
		   if(kaskeymenos_vFit_noMC->currentState().mass()>1.5) continue;

		  kaskeymenosFitTree ->movePointerToTheFirstChild();
		   RefCountedKinematicParticle kaskeyTrk1 = kaskeymenosFitTree->currentParticle();
		   kaskeymenosFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle kaskeyTrk2 = kaskeymenosFitTree->currentParticle();
		   kaskeymenosFitTree->movePointerToTheTop();	







 // Kaskey menos mass constraint 

	       KinematicParticleFitter csFitterkaskey;
	       KinematicConstraint * kaskey_c = new MassKinematicConstraint(kaskeym_mass,kaskeym_sigma);
	       // add mass constraint to the lam fit to do a constrained fit:  
 
	       kaskeymenosFitTree = csFitterkaskey.fit(kaskey_c,kaskeymenosFitTree);
	       if (!kaskeymenosFitTree->isValid()){
		 //std::cout << "caught an exception in the lam mass constraint fit" << std::endl;
		 continue; 
	       }
	       
	       kaskeymenosFitTree->movePointerToTheTop();
	       RefCountedKinematicParticle kaskeymenos_vFit_withMC = kaskeymenosFitTree->currentParticle();


 // JPsi mass constraint is applied in the final Kaskeyb fit,

	       vector<RefCountedKinematicParticle> vFitMCParticles;
	       vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	       vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	       vFitMCParticles.push_back(kaskeymenos_vFit_withMC);


	       MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
	       KinematicConstrainedVertexFitter kcvFitter;
	       RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
	       if (!vertexFitTree->isValid()) {
		 //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		 continue;
	       }
	       vertexFitTree->movePointerToTheTop();

	       RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
	       RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
	       if (!bDecayVertexMC->vertexIsValid()){
		 // cout << "B MC fit vertex is not valid" << endl;
		 continue;
	       }
	       
	       if ( bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50) {
		 //if ( bDecayVertexMC->chiSquared()<0 ) cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		 //if (debug) cout << " continue from bad chi2 = " << bDecayVertexMC->chiSquared() << endl;
		 continue;
	       }
	       
	       
	       if ( (bCandMC->currentState().mass() > 6.0) || (bCandMC->currentState().mass() < 5.6) ) {
	         // (debug) cout << "continue from bmass > 6.5 or < 4.5 = " << bCandMC->currentState().mass() << endl;
	         continue;
	       }

 // get children from final B fit
	       vertexFitTree->movePointerToTheFirstChild();
	       RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
	      
	       vertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
	       
	       vertexFitTree->movePointerToTheNextChild();
	       RefCountedKinematicParticle kaskeymCandMC = vertexFitTree->currentParticle();
	      

	      
// en esta parte llamaremos los Id de las particulas


	       int mId1=0; int mId2=0; int PId1=0; int piId2=0;
	       int kId3;   int lamId1=0; int lamId2=0;
	       int JId1=0; int JId2=0;  int JId3=0; int JId4=0; int JId5=0; 
	       


	       
	       int ngenT1 = 0;	      
	       int ngenT2 = 0;
	       int ngenT3 = 0;// esta sera para el pion que se une al lambda para formar el kaskey menos
	       int ngenT4 = 0;
	       int ngenT5 = 0;
	       


	       /*
	       int ngenT1 = PdgIDatTruthLevel(iMuonP->track(), genParticles, mId1, JId1);
	       int ngenT2 = PdgIDatTruthLevel(iMuonM->track(), genParticles, mId2, JId2);
	       int ngenT3 = PdgIDatTruthLevel(iTrack1->track(), genParticles, kId3, JId3);
	       int ngenT4 = PdgIDatTruthLevel4(theDaughterTracks[indexP], genParticles, PId1, lamId1, JId4 );
	       int ngenT5 = PdgIDatTruthLevel4(theDaughterTracks[indexPi],genParticles, piId2,lamId2, JId5 );
	       */

	       /* 
	     cout<< "Bd Id: "<< JId1 << "\t J/psi Id:" << mId1 << "\t J Daughter Id: " << ngenT1 <<  endl;
	     cout<< "Bd Id: " << JId2 << "\t J/psi Id:" << mId2 << "\t J Daughter Id: " << ngenT2 <<  endl;
	       
	     cout<< "Bd Id: "<< JId4 << "\t kmenos Id:" << lamId1 << "\t lambda Id: " << PId1 << "\t lamda Daughter Id: " << ngenT4 << endl;
	     cout<< "Bd Id: "<< JId5 << "\t kmenos Id:" << lamId2 << "\t lambda Id: " << piId2 << "\t lamda Daughter Id: " << ngenT5 << endl;
	       */
	       






 // get mu+ and mu- momenta from final B fit

	       KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
	       KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
	       KinematicParameters psiMupKP;
	       KinematicParameters psiMumKP;
	       
	       if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
	       if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
	       if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
	       if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 
	       
	       //mupCategory = getMuCat( *iMuonP );
	       //mumCategory = getMuCat( *iMuonM );

	       if(iMuon1->charge() == 1) mupCategory = getMuCat( *iMuon1 );
	       if(iMuon1->charge() == -1) mumCategory = getMuCat( *iMuon1 );
	       if(iMuon2->charge() == 1) mupCategory = getMuCat( *iMuon2 );
	       if(iMuon2->charge() == -1) mumCategory = getMuCat( *iMuon2 );

               const reco::Muon *recoMuonM = patMuonM;
	       const reco::Muon *recoMuonP = patMuonP;

               mumME1Clean = HasGoodME11(*recoMuonM,2.);
               mupME1Clean = HasGoodME11(*recoMuonP,2.);

	       //KinematicParameters VCandKP = kaskeymCandMC->currentState().kinematicParameters();
	       



 // get momenta from final kaskeym fit
	       KinematicParameters kasTrk1KP = kaskeyTrk1->currentState().kinematicParameters();
	       KinematicParameters kasTrk2KP = kaskeyTrk2->currentState().kinematicParameters();
	       KinematicParameters kasTrknKP;  // este es el trak neutro, es decir el del lambda
	       KinematicParameters kasTrkmKP;  // este es el trak cargado negativo, es decir el del Kaon

	       //this lamTrk1KP momentum is defined at the lambda fit vertex.

	       if ( kaskeyTrk1->currentState().particleCharge() == 0 ) kasTrknKP = kasTrk1KP;
	       if ( kaskeyTrk1->currentState().particleCharge() < 0 )  kasTrkmKP = kasTrk1KP;
	       if ( kaskeyTrk2->currentState().particleCharge() == 0 ) kasTrknKP = kasTrk2KP;
	       if ( kaskeyTrk2->currentState().particleCharge() < 0 )  kasTrkmKP = kasTrk2KP;



  // get momenta from final lambda fit
	       KinematicParameters lamTrk1KP = lamTrk1->currentState().kinematicParameters();
	       KinematicParameters lamTrk2KP = lamTrk2->currentState().kinematicParameters();
	       KinematicParameters lamTrkpKP;
	       KinematicParameters lamTrkmKP;

	       //this lamTrk1KP momentum is defined at the lambda fit vertex.

	       if ( lamTrk1->currentState().particleCharge() > 0 ) lamTrkpKP = lamTrk1KP;
	       if ( lamTrk1->currentState().particleCharge() < 0 ) lamTrkmKP = lamTrk1KP;
	       if ( lamTrk2->currentState().particleCharge() > 0 ) lamTrkpKP = lamTrk2KP;
	       if ( lamTrk2->currentState().particleCharge() < 0 ) lamTrkmKP = lamTrk2KP;



	           GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				       mu1CandMC->currentState().globalMomentum().y(),
 				       mu1CandMC->currentState().globalMomentum().z());


 		   GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				       mu2CandMC->currentState().globalMomentum().y(),
 				       mu2CandMC->currentState().globalMomentum().z());


 		   GlobalVector kaskeymp1vec(kaskeyTrk1->currentState().globalMomentum().x(),
					    kaskeyTrk1->currentState().globalMomentum().y(),
					    kaskeyTrk1->currentState().globalMomentum().z());


 		   GlobalVector kaskeymp2vec(kaskeyTrk2->currentState().globalMomentum().x(),
					    kaskeyTrk2->currentState().globalMomentum().y(),
					    kaskeyTrk2->currentState().globalMomentum().z());

		   GlobalVector lamp1vec(lamTrk1->currentState().globalMomentum().x(),
					 lamTrk1->currentState().globalMomentum().y(),
					 lamTrk1 ->currentState().globalMomentum().z());


 		   GlobalVector lamp2vec(lamTrk2->currentState().globalMomentum().x(),
					 lamTrk2 ->currentState().globalMomentum().y(),
					 lamTrk2->currentState().globalMomentum().z());




// ************
		   
		   // Only save the first time
		   if(nB==0){


		     CheckHLTTriggers(TrigTable);
		     // cout<<"*Trigger List "<< triggersL << endl;
		     // Save number of Muons



		     // Get L1 trigger to level Event
		     string ListTriggL1_tmp="";
		     CheckL1Triggers(iEvent, iSetup, ListTriggL1_tmp);
		     //cout<<"L1: " << ListTriggL1_tmp<<endl;
		     //int ntriggersL1 = sprintf(triggersL1L,"%s",ListTriggL1_tmp.c_str());

		     nTrgL1L = ListTriggL1_tmp.size();
		     //triggersL1L->push_back(triggersL1);


		     nMu  = nMu_tmp;
		     // cout<< "*Number of Muons : " << nMu_tmp << endl;
		     
		   // Get number of run and event
		     run   =  run1;
		     event =  event1;
		     
		     // cout<< "* Run: "<< run << "  event: " << event << endl;

		   } // end nB==0

  // ********************* todos los vertices primarios y escogemos el de mejor pointing angle **************** 
		   //reco::Vertex bestVtxIP;


		   Double_t pVtxIPX_temp = -10000.0;
		   Double_t pVtxIPY_temp = -10000.0;
		   Double_t pVtxIPZ_temp = -10000.0;
		   Double_t pVtxIPXE_temp = -10000.0;
		   Double_t pVtxIPYE_temp = -10000.0;
		   Double_t pVtxIPZE_temp = -10000.0;
		   Double_t pVtxIPXYE_temp = -10000.0;
		   Double_t pVtxIPXZE_temp = -10000.0;
		   Double_t pVtxIPYZE_temp = -10000.0;
		   Double_t pVtxIPCL_temp = -10000.0;	
		   Double_t lip1 = -1000000.0;
		     for(size_t i = 0; i < recVtxs->size(); ++i) {
		       const Vertex &vtx = (*recVtxs)[i];
		       
		       Double_t dx1 = (*bDecayVertexMC).position().x() - vtx.x(); 
		       Double_t dy1 = (*bDecayVertexMC).position().y() - vtx.y();
		       Double_t dz1 = (*bDecayVertexMC).position().z() - vtx.z();
		       float cosAlphaXYb1 = ( bCandMC->currentState().globalMomentum().x() * dx1 + bCandMC->currentState().globalMomentum().y()*dy1 + bCandMC->currentState().globalMomentum().z()*dz1  )/( sqrt(dx1*dx1+dy1*dy1+dz1*dz1)* bCandMC->currentState().globalMomentum().mag() );

		       if(cosAlphaXYb1>lip1)
			 {
			   lip1 = cosAlphaXYb1 ;
			   pVtxIPX_temp = vtx.x();
			   pVtxIPY_temp = vtx.y();
			   pVtxIPZ_temp = vtx.z();
			   pVtxIPXE_temp = vtx.covariance(0, 0);
			   pVtxIPYE_temp = vtx.covariance(1, 1);
			   pVtxIPZE_temp = vtx.covariance(2, 2);
			   pVtxIPXYE_temp = vtx.covariance(0, 1);
			   pVtxIPXZE_temp = vtx.covariance(0, 2);
			   pVtxIPYZE_temp = vtx.covariance(1, 2);
			   pVtxIPCL_temp = (TMath::Prob(vtx.chi2(),(int)vtx.ndof()) );
			 
			 }
                 
		     }
		   pVtxIPX->push_back( pVtxIPX_temp);
		   pVtxIPY->push_back(  pVtxIPY_temp);	    
		   pVtxIPZ->push_back(  pVtxIPZ_temp);
		   pVtxIPXE->push_back( pVtxIPXE_temp);
		   pVtxIPYE->push_back( pVtxIPYE_temp);	    
		   pVtxIPZE->push_back( pVtxIPZE_temp);
		   pVtxIPXYE->push_back( pVtxIPXYE_temp);
		   pVtxIPXZE->push_back( pVtxIPXZE_temp);	    
		   pVtxIPYZE->push_back( pVtxIPYZE_temp);
		   pVtxIPCL->push_back(  pVtxIPCL_temp);
		     


     // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle **************** 

		   reco::Vertex bestVtxBSIP;


		   Double_t pVtxBSIPX_temp = -10000.0;
		   Double_t pVtxBSIPY_temp = -10000.0;
		   Double_t pVtxBSIPZ_temp = -10000.0;
		   Double_t pVtxBSIPXE_temp = -10000.0;
		   Double_t pVtxBSIPYE_temp = -10000.0;
		   Double_t pVtxBSIPZE_temp = -10000.0;
		   Double_t pVtxBSIPXYE_temp = -10000.0;
		   Double_t pVtxBSIPXZE_temp = -10000.0;
		   Double_t pVtxBSIPYZE_temp = -10000.0;
		   Double_t pVtxBSIPCL_temp = -10000.0;	 
		   
		   Double_t lip = -1000000.0;

		   for(size_t i = 0; i < recVtxsBS->size(); ++i) {
		     const Vertex &vtxBS = (*recVtxsBS)[i];
		     

		     Double_t dx = (*bDecayVertexMC).position().x() - vtxBS.x();
		     Double_t dy = (*bDecayVertexMC).position().y() - vtxBS.y();
		     Double_t dz = (*bDecayVertexMC).position().z() - vtxBS.z();
		     Double_t cosAlphaXYb = ( bCandMC->currentState().globalMomentum().x() * dx + bCandMC->currentState().globalMomentum().y()*dy + bCandMC->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* bCandMC->currentState().globalMomentum().mag() );

		     //float cosAlphaXYb = ( bCandMC->currentState().globalMomentum().x() * dx + bCandMC->currentState().globalMomentum().y()*dy + bCandMC->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)*sqrt( bCandMC->currentState().globalMomentum().x()*bCandMC->currentState().globalMomentum().x() + bCandMC->currentState().globalMomentum().y()*bCandMC->currentState().globalMomentum().y() + bCandMC->currentState().globalMomentum().z()*bCandMC->currentState().globalMomentum().z() ) );
		     if(cosAlphaXYb>lip)
		       {
			 lip = cosAlphaXYb ;
			 
			 pVtxBSIPX_temp = vtxBS.x();
			 pVtxBSIPY_temp = vtxBS.y();
			 pVtxBSIPZ_temp = vtxBS.z();
			 pVtxBSIPXE_temp = vtxBS.covariance(0, 0);
			 pVtxBSIPYE_temp = vtxBS.covariance(1, 1);
			 pVtxBSIPZE_temp = vtxBS.covariance(2, 2);
			 pVtxBSIPXYE_temp = vtxBS.covariance(0, 1);
			 pVtxBSIPXZE_temp = vtxBS.covariance(0, 2);
			 pVtxBSIPYZE_temp = vtxBS.covariance(1, 2);
			 pVtxBSIPCL_temp = (TMath::Prob(vtxBS.chi2(),(int)vtxBS.ndof()) );
			   
			 bestVtxBSIP = vtxBS;

		       }
		     
		   }
 
		   pVtxBSIPX->push_back( pVtxBSIPX_temp);
		   pVtxBSIPY->push_back(  pVtxBSIPY_temp);	    
		   pVtxBSIPZ->push_back(  pVtxBSIPZ_temp);
		   pVtxBSIPXE->push_back( pVtxBSIPXE_temp);
		   pVtxBSIPYE->push_back( pVtxBSIPYE_temp);	    
		   pVtxBSIPZE->push_back( pVtxBSIPZE_temp);
		   pVtxBSIPXYE->push_back( pVtxBSIPXYE_temp);
		   pVtxBSIPXZE->push_back( pVtxBSIPXZE_temp);	    
		   pVtxBSIPYZE->push_back( pVtxBSIPYZE_temp);
		   pVtxBSIPCL->push_back(  pVtxBSIPCL_temp);
		     

		   //cout << "PV bestVtxBSIP: " <<bestVtxBSIP.x()<< " "<<bestVtxBSIP.y()<<" "<<bestVtxBSIP.z()<< endl;
		   //cout << "como se guarda: " <<pVtxBSIPX_temp<< " "<<pVtxBSIPY_temp<<" "<<pVtxBSIPZ_temp<< endl;


		 // Get Matching Muon to HLT 
		 string ListTriggMuP_tmp="";
                 string ListTriggMuM_tmp="";

                 if(iMuon1->charge()== 1) MatchMuonWithTriggers(*iMuon1, TrigTable, ListTriggMuP_tmp);
                 if(iMuon1->charge()==-1) MatchMuonWithTriggers(*iMuon1, TrigTable, ListTriggMuM_tmp);
 
                 if(iMuon2->charge()== 1) MatchMuonWithTriggers(*iMuon2, TrigTable, ListTriggMuP_tmp);
                 if(iMuon2->charge()==-1) MatchMuonWithTriggers(*iMuon2, TrigTable, ListTriggMuM_tmp);
                 
                 
                 //nMuonPTrgL = ListTriggMuP_tmp.size();
                 //nMuonMTrgL = ListTriggMuM_tmp.size();
                 

                 int nMuonP = sprintf(triggersMuP,"%s",ListTriggMuP_tmp.c_str());
                 int nMuonM = sprintf(triggersMuM,"%s",ListTriggMuM_tmp.c_str());

		 nMuonPTrgL = nMuonP;
                 nMuonMTrgL = nMuonM;

		 //cout<<"MuP = "<<nMuonPTrgL<<" "<<nMuonP<<" "<<triggersMuP<<endl;
		 //cout<<"MuM = "<<nMuonMTrgL<<" "<<nMuonM<<" "<<triggersMuM<<endl;

		 triggersMuPL ->push_back(triggersMuP);
		 triggersMuML ->push_back(triggersMuM);
                 
                 //cout<<"nMuonPTrgL = "<<nMuonPTrgL<< " nMuonMTrgL = "<<nMuonMTrgL<<endl;
                 
		 // Get List for L1/L2 triggers matching to Muon
                 string ListTriggL1L2_MuP="";
                 string ListTriggL1L2_MuM="";
                 if(iMuon1->charge()== 1)  MatchMuonWithL1L2(*iMuon1, ListTriggerL1L2, ListTriggL1L2_MuP);
                 if(iMuon1->charge()== -1) MatchMuonWithL1L2(*iMuon1, ListTriggerL1L2, ListTriggL1L2_MuM);
                 
                 if(iMuon2->charge()== 1)  MatchMuonWithL1L2(*iMuon2, ListTriggerL1L2, ListTriggL1L2_MuP);
                 if(iMuon2->charge()== -1) MatchMuonWithL1L2(*iMuon2, ListTriggerL1L2, ListTriggL1L2_MuM);

		 // cout << "MuonP: " << ListTriggL1L2_MuP << endl;
                 //cout << "MuonM: " << ListTriggL1L2_MuM << endl;
		 int nL1L2MuP = sprintf(triggersL1L2_MuP,"%s",ListTriggL1L2_MuP.c_str());
		 int nL1L2MuM = sprintf(triggersL1L2_MuM,"%s",ListTriggL1L2_MuM.c_str());

		 //ntriggersL1L2_MuP = ListTriggL1L2_MuP.size();
		 //ntriggersL1L2_MuM = ListTriggL1L2_MuM.size();

		 triggersL1L2_MuPL->push_back(triggersL1L2_MuP);
		 triggersL1L2_MuML->push_back(triggersL1L2_MuM);

		 ntriggersL1L2_MuP = nL1L2MuP;
		 ntriggersL1L2_MuM = nL1L2MuM;

                 //cout << "nmuonlil2P: " <<nL1L2MuP<<" "<<triggersL1L2_MuP<< endl;
                 //cout << "nmuonlil2M: " <<nL1L2MuM<<" "<<triggersL1L2_MuM<< endl;





  
// fill candidate variables now

		   // ************
	           B_mass->push_back(bCandMC->currentState().mass());
		   B_px->push_back(bCandMC->currentState().globalMomentum().x());
		   B_py->push_back(bCandMC->currentState().globalMomentum().y());
		   B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		   B_parentId1->push_back(JId1);
		   B_parentId2->push_back(JId2);
		   B_parentId3->push_back(JId3);
		   B_parentId4->push_back(JId4);
		   B_parentId5->push_back(JId5);

		   B_kaskeym_mass->push_back(kaskeymenos_vFit_noMC->currentState().mass() );
		   B_kaskeym_px->push_back(kaskeymenos_vFit_noMC ->currentState().globalMomentum().x() );
		   B_kaskeym_py->push_back(kaskeymenos_vFit_noMC ->currentState().globalMomentum().y() );
		   B_kaskeym_pz->push_back(kaskeymenos_vFit_noMC ->currentState().globalMomentum().z() );
		   B_kaskeym_charge-> push_back(kaskeymCandMC ->currentState().particleCharge() );
		   B_kaskey_parentId1->push_back(kId3);
		   B_kaskey_parentId2->push_back(lamId1);
		   B_kaskey_parentId3->push_back(lamId2);

		   B_J_mass->push_back(psi_vFit_noMC->currentState().mass() );
		   B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		   B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		   B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
		   B_J_parentId1->push_back(mId1);
		   B_J_parentId2->push_back(mId2);

		   B_kaskeym_pt1->push_back(kaskeymp1vec.perp());
		   B_kaskeym_px1->push_back(kasTrk1KP.momentum().x());
		   B_kaskeym_py1->push_back(kasTrk1KP.momentum().y());
		   B_kaskeym_pz1->push_back(kasTrk1KP.momentum().z());
		   B_kaskeym_charge1->push_back(kaskeyTrk1->currentState().particleCharge());

		   B_kaskeym_kId3->push_back(ngenT3);
		   B_kaskeym_pt2->push_back(kaskeymp2vec.perp());
		   B_kaskeym_px2->push_back(kasTrk2KP.momentum().x());
		   B_kaskeym_py2->push_back(kasTrk2KP.momentum().y());
		   B_kaskeym_pz2->push_back(kasTrk2KP.momentum().z());
		   B_kaskeym_charge2->push_back(kaskeyTrk2->currentState().particleCharge());

		   B_J_muId1->push_back(ngenT1);
		   B_J_pt1->push_back(Jp1vec.perp());
		   B_J_px1->push_back(psiMu1KP.momentum().x());
		   B_J_py1->push_back(psiMu1KP.momentum().y());
		   B_J_pz1->push_back(psiMu1KP.momentum().z());
		   B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

		   B_J_muId2->push_back(ngenT2);
		   B_J_pt2->push_back(Jp2vec.perp());
		   B_J_px2->push_back(psiMu2KP.momentum().x());
		   B_J_py2->push_back(psiMu2KP.momentum().y());
		   B_J_pz2->push_back(psiMu2KP.momentum().z());
		   B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

		   B_kaskeym_chi2->push_back(kaskeymenos_vFit_vertex_noMC->chiSquared());
		   B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		   B_chi2->push_back(bDecayVertexMC->chiSquared());





		   double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   double Omb_J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   double Omb_kaskey_Prob_tmp  = TMath::Prob(kaskeymenos_vFit_vertex_noMC->chiSquared(),(int)kaskeymenos_vFit_vertex_noMC->degreesOfFreedom());
		   double Omb_lamb_Prob_tmp  = TMath::Prob(lam_vFit_vertex_noMC->chiSquared(),(int)lam_vFit_vertex_noMC->degreesOfFreedom());		  

		   B_Prob    ->push_back(B_Prob_tmp);
		   B_J_Prob  ->push_back(Omb_J_Prob_tmp);
		   B_kaskey_Prob ->push_back(Omb_kaskey_Prob_tmp);
		   B_kas_lambda_Prob ->push_back(Omb_lamb_Prob_tmp);



		   B_kas_lambda_mass->push_back(lam_vFit_noMC->currentState().mass() );
		   B_kas_lambda_px->push_back(lam_vFit_noMC->currentState().globalMomentum().x() );
		   B_kas_lambda_py->push_back(lam_vFit_noMC->currentState().globalMomentum().y() );
		   B_kas_lambda_pz->push_back(lam_vFit_noMC->currentState().globalMomentum().z() );
		   B_lam_parentId1->push_back(PId1);
		   B_lam_parentId2->push_back(piId2);

		   B_lam_PId1->push_back(ngenT4);
		   B_kas_lambda_pt1->push_back(lamp1vec.perp());
		   B_kas_lambda_px1->push_back(lamTrk1KP.momentum().x());
		   B_kas_lambda_py1->push_back(lamTrk1KP.momentum().y());
		   B_kas_lambda_pz1->push_back(lamTrk1KP.momentum().z());
		   B_kas_lambda_charge1->push_back(lamTrk1->currentState().particleCharge());

		   B_lam_piId2->push_back(ngenT5);
		   B_kas_lambda_pt2->push_back(lamp2vec.perp());
		   B_kas_lambda_px2->push_back(lamTrk2KP.momentum().x());
		   B_kas_lambda_py2->push_back(lamTrk2KP.momentum().y());
		   B_kas_lambda_pz2->push_back(lamTrk2KP.momentum().z());
		   B_kas_lambda_charge2->push_back(lamTrk2->currentState().particleCharge());

		   B_kas_lambda_chi2->push_back(lam_vFit_vertex_noMC->chiSquared());

	   // ************

//  bDecayVtx se refiere a la posicion del decaimiento del kaskey b

		   bDecayVtxX->push_back((*bDecayVertexMC).position().x());
		   bDecayVtxY->push_back((*bDecayVertexMC).position().y());
		   bDecayVtxZ->push_back((*bDecayVertexMC).position().z());
		   bDecayVtxXE->push_back(bDecayVertexMC->error().cxx());
		   bDecayVtxYE->push_back(bDecayVertexMC->error().cyy());
		   bDecayVtxZE->push_back(bDecayVertexMC->error().czz());
		   bDecayVtxXYE->push_back(bDecayVertexMC->error().cyx());
		   bDecayVtxXZE->push_back(bDecayVertexMC->error().czx());
		   bDecayVtxYZE->push_back(bDecayVertexMC->error().czy());

//  VDecayVtx se refiere a la posicion del decaimiento del kaskey menos

		   VDecayVtxX->push_back(kaskeymenos_vFit_vertex_noMC->position().x() );
		   VDecayVtxY->push_back(kaskeymenos_vFit_vertex_noMC->position().y() );
		   VDecayVtxZ->push_back(kaskeymenos_vFit_vertex_noMC->position().z() );
		   VDecayVtxXE->push_back( kaskeymenos_vFit_vertex_noMC->error().cxx() );
		   VDecayVtxYE->push_back( kaskeymenos_vFit_vertex_noMC->error().cyy() );
		   VDecayVtxZE->push_back( kaskeymenos_vFit_vertex_noMC->error().czz() );
		   VDecayVtxXYE->push_back( kaskeymenos_vFit_vertex_noMC->error().cyx() );
		   VDecayVtxXZE->push_back( kaskeymenos_vFit_vertex_noMC->error().czx() );
		   VDecayVtxYZE->push_back( kaskeymenos_vFit_vertex_noMC->error().czy() );


   //  V1DecayVtx se refiere a la posicion del decaimiento del lambda

		   V1DecayVtxX->push_back( lam_vFit_vertex_noMC->position().x() );
		   V1DecayVtxY->push_back( lam_vFit_vertex_noMC->position().y() );
		   V1DecayVtxZ->push_back( lam_vFit_vertex_noMC->position().z() );
		   V1DecayVtxXE->push_back( lam_vFit_vertex_noMC->error().cxx() );
		   V1DecayVtxYE->push_back( lam_vFit_vertex_noMC->error().cyy() );
		   V1DecayVtxZE->push_back( lam_vFit_vertex_noMC->error().czz() );
		   V1DecayVtxXYE->push_back( lam_vFit_vertex_noMC->error().cyx() );
		   V1DecayVtxXZE->push_back( lam_vFit_vertex_noMC->error().czx() );
		   V1DecayVtxYZE->push_back( lam_vFit_vertex_noMC->error().czy() );

//  JDecayVtx se refiere a la posicion del decaimiento del JPsi

		   //JMass->push_back( psi_vFit_noMC->currentState().mass() ); 
		   JDecayVtxX->push_back( psi_vFit_vertex_noMC->position().x() );
		   JDecayVtxY->push_back( psi_vFit_vertex_noMC->position().y() );
		   JDecayVtxZ->push_back( psi_vFit_vertex_noMC->position().z() );
		   JDecayVtxXE->push_back( psi_vFit_vertex_noMC->error().cxx() );
		   JDecayVtxYE->push_back( psi_vFit_vertex_noMC->error().cyy() );
		   JDecayVtxZE->push_back( psi_vFit_vertex_noMC->error().czz() );
		   JDecayVtxXYE->push_back( psi_vFit_vertex_noMC->error().cyx() );
		   JDecayVtxXZE->push_back( psi_vFit_vertex_noMC->error().czx() );
		   JDecayVtxYZE->push_back( psi_vFit_vertex_noMC->error().czy() );

	       // JVtxCL->push_back( ChiSquaredProbability((double)(psi_vFit_vertex_noMC->chiSquared()),(double)(psi_vFit_vertex_noMC->degreesOfFreedom())) );

		  
		   mu1soft->push_back(iMuon1->isSoftMuon(bestVtxBSIP) );
		   mu2soft->push_back(iMuon2->isSoftMuon(bestVtxBSIP) );
		   mu1tight->push_back(iMuon1->isTightMuon(bestVtxBSIP) );
		   mu2tight->push_back(iMuon2->isTightMuon(bestVtxBSIP) );
		   mu1PF->push_back(iMuon1->isPFMuon());
		   mu2PF->push_back(iMuon2->isPFMuon());
		   mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		   mu2loose->push_back(muon::isLooseMuon(*iMuon2));

		   mumC2->push_back( glbTrackM->normalizedChi2() );
		   mumCat->push_back( mumCategory );
		   //mumAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMLastStationAngTight) ); // este no se  que es
		   mumAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMOneStationTight) ); // este es para poner la condicion si es o no softmuon
		   mumNHits->push_back( glbTrackM->numberOfValidHits() );
		   mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		   mupC2->push_back( glbTrackP->normalizedChi2() );
		   mupCat->push_back( mupCategory );
		   //mupAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMLastStationAngTight) );  // este no se  que es
		   mupAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMOneStationTight) ); // este es para poner la condicion si es o no softmuon
		   mupNHits->push_back( glbTrackP->numberOfValidHits() );
		   mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
                   mumdxy->push_back(glbTrackM->dxy(bestVtxBSIP.position()) );// el dxy del Muon negatico respetcto del PV con BSc (el de mayor pt)
		   mupdxy->push_back(glbTrackP->dxy(bestVtxBSIP.position()) );// el dxy del Muon positivo respetcto del PV con BSc (el de mayor pt)
		   mumdz->push_back(glbTrackM->dz(bestVtxBSIP.position()) );
		   mupdz->push_back(glbTrackP->dz(bestVtxBSIP.position()) );
		   muon_dca->push_back(dca);
		   //cout<<" closest approach  "<<dca<<endl;
		  
		 

		  

// aca esta el tiempo de vida del  Kaskey b
		   //////////////////////////////////////////////////////////////////////////////////////////
		   //calculate ctau and ctau2D with ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)
		   /*
		   float mb = 5.795;
		   GlobalPoint PVP = GlobalPoint( bestVtx.position().x(), bestVtx.position().y(), bestVtx.position().z() );
		   GlobalError PVE = GlobalError( bestVtx.error() );
		   float bctau2D_temp  = 0.;
		   float bctauE_temp   = 0.;
		   float bctau2DE_temp = 0.;
		   float bctau_temp = Myctau(bCandMC, bDecayVertexMC, PVP, PVE, mb, bctau2D_temp, bctauE_temp, bctau2DE_temp);

		   bctau    -> push_back( bctau_temp    );
		   bctau2D  -> push_back( bctau2D_temp  );
		   bctauE   -> push_back( bctauE_temp   );
		   bctau2DE -> push_back( bctau2DE_temp );
		   */
		   /////////////////////////////////////////////////////////////////////////////////////////////////////////


 // aca esta el tiempo de vida del  Kaskey menos
		   //////////////////////////////////////////////////////////////////////////////////////////
		   //calculate ctau and ctau2D with ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)
		   /*
		   float mkaskey = 1.32171;		 
		   GlobalPoint BVP = GlobalPoint( bDecayVertexMC->position() );
		   GlobalError BVE = bDecayVertexMC->error();
		   float bctau2D_temp_kaskey  = 0.;
		   float bctauE_temp_kaskey   = 0.;
		   float bctau2DE_temp_kaskey = 0.;
		   float bctau_temp_kaskey = Myctau(kaskeymenos_vFit_noMC, kaskeymenos_vFit_vertex_noMC , BVP, BVE, mkaskey, bctau2D_temp_kaskey, bctauE_temp_kaskey, bctau2DE_temp_kaskey);


		   bctau_kaskey    ->push_back( bctau_temp_kaskey    );	
		   bctau2D_kaskey  ->push_back(bctau2D_temp_kaskey   ); 
		   bctauE_kaskey   -> push_back( bctauE_temp_kaskey  );
		   bctau2DE_kaskey -> push_back( bctau2DE_temp_kaskey);
		   */
		   /////////////////////////////////////////////////////////////////////////////////////////////////////////


// aca esta el tiempo de vida del  Lambda
		   //////////////////////////////////////////////////////////////////////////////////////////
		   //calculate ctau and ctau2D with ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)
		   /*
		   float mlam = 1.115683 ;
		   GlobalPoint OVP = GlobalPoint( kaskeymenos_vFit_vertex_noMC->position() );
		   GlobalError OVE = kaskeymenos_vFit_vertex_noMC->error();
		   float bctau2D_temp_lam  = 0.;
		   float bctauE_temp_lam   = 0.;
		   float bctau2DE_temp_lam = 0.;
		   float bctau_temp_lam = Myctau(lam_vFit_noMC, lam_vFit_vertex_noMC, OVP, OVE, mlam, bctau2D_temp_lam, bctauE_temp_lam, bctau2DE_temp_lam);


		   bctau_lam    ->push_back( bctau_temp_lam   );
		   bctau2D_lam  ->push_back( bctau2D_temp_lam );
		   bctauE_lam   ->push_back( bctauE_temp_lam  );
		   bctau2DE_lam ->push_back( bctau2DE_temp_lam);
		   */
		   /////////////////////////////////////////////////////////////////////////////////////////////////////////


		   //calculate ctau 3D with BS constraint
		   /*
		   GlobalPoint PVBSP = GlobalPoint( bestVtxBS.position().x(), bestVtxBS.position().y(), bestVtxBS.position().z() );
		   GlobalError PVBSE = GlobalError( bestVtxBS.error() );
		   float bctau2DBS_temp  = 0.;
		   float bctauEBS_temp   = 0.;
		   float bctau2DEBS_temp = 0.;
		   float bctauBS_temp = Myctau(bCandMC, bDecayVertexMC, PVBSP, PVE, mb, bctau2DBS_temp, bctauEBS_temp, bctau2DEBS_temp);

		   bctauBS    -> push_back(bctauBS_temp);
		   bctauBSE   ->push_back(bctauEBS_temp);
		   bctauBS2D  -> push_back(bctau2DBS_temp);
		   bctauBS2DE -> push_back(bctau2DEBS_temp);
		   */
		  
		   // try refitting the primary without the tracks in the B reco candidate
		   
		   // first get tracks from the original primary
		   vector<reco::TransientTrack> vertexTracks;
		   
		   for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestVtxBSIP.tracks_begin();
			 iTrack != bestVtxBSIP.tracks_end(); ++iTrack) {
		     // compare primary tracks to check for matches with B cand
		     TrackRef trackRef = iTrack->castTo<TrackRef>();
		     
		     // the 5 tracks in the B cand are theDaughterTracks[0] theDaughterTracks[1] glbTrackP glbTrackM patTrack3
		     if (  !( (theDaughterTracks[0]==trackRef) || (theDaughterTracks[1]==trackRef) ||
			      (patTrack3.track()==trackRef) || (glbTrackP==trackRef) || (glbTrackM==trackRef) ) ) {
		       TransientTrack tt(trackRef, &(*bFieldHandle) );
		       vertexTracks.push_back(tt);
		     } //else { std::cout << "found track match with primary" << endl;}
		   }
		   
		   priRfNTrkDif->push_back( bestVtxBSIP.tracksSize() - vertexTracks.size() );
		   
		   // if no tracks in primary or no reco track included in primary then don't do anything
		   // if so, then update bctau_temp and bctauMPV_temp
		   
		   reco::Vertex bestVtxRf = bestVtxBSIP;
		   
		   //float bctau2DRf_temp  = -1000.0;
		   //float bctauRfE_temp   = -1000.0;
		   //float bctau2DRfE_temp =-1000.0;
		   //float bctauRf_temp = -1000.0; 

		   if (  vertexTracks.size()>0 && (bestVtxBSIP.tracksSize()!=vertexTracks.size()) ) {
		     
		     AdaptiveVertexFitter theFitter;
		     TransientVertex v = theFitter.vertex(vertexTracks);
		     if ( v.isValid() ) {
		       //calculate ctau with the new vertex to compare to the old one.
		       //GlobalPoint PVRfP = GlobalPoint( v.position().x(), v.position().y(), v.position().z() );
			 reco::Vertex recoV = (reco::Vertex)v;
			 //GlobalError PVRfE = GlobalError( recoV.error() );		   
			 //bctauRf_temp = Myctau(bCandMC, bDecayVertexMC, PVRfP, PVRfE, mb, bctau2DRf_temp, bctauRfE_temp, bctau2DRfE_temp);
		       
		       //set bestVtxRf as new best vertex to fill variables for ntuple
		       bestVtxRf = reco::Vertex(v);
		       
		     } 
		   }
		   /*
		     else {
		     //bctauRf_temp = Myctau(bCandMC, bDecayVertexMC, PVP, PVE, mb, bctau2DRf_temp, bctauRfE_temp, bctau2DRfE_temp);
		     bctauRf_temp = bctau_temp;
		     bctau2DRf_temp =  bctau2D_temp;
		     bctauRfE_temp = bctauE_temp;
		     bctau2DRfE_temp = bctau2DE_temp;
		     } 
		   } 
		   else {
		     //bctauRf_temp = Myctau(bCandMC, bDecayVertexMC, PVP, PVE, mb, bctau2DRf_temp, bctauRfE_temp, bctau2DRfE_temp);
		     bctauRf_temp = bctau_temp;
		     bctau2DRf_temp =  bctau2D_temp;
		     bctauRfE_temp = bctauE_temp;
		     bctau2DRfE_temp = bctau2DE_temp;
		   } 
		   
		  bctauRf->push_back( bctauRf_temp );
		  bctauRfE->push_back( bctauRfE_temp );
		  bctau2DRf->push_back( bctau2DRf_temp );
		  bctau2DRfE->push_back( bctau2DRfE_temp );
		   */
		   
		   priRfVtxX->push_back( bestVtxRf.x() );
		   priRfVtxY->push_back( bestVtxRf.y() );
		   priRfVtxZ->push_back( bestVtxRf.z() );
		   priRfVtxXE->push_back( bestVtxRf.covariance(0, 0) );
		   priRfVtxYE->push_back( bestVtxRf.covariance(1, 1) );
		   priRfVtxZE->push_back( bestVtxRf.covariance(2, 2) );
		   priRfVtxXYE->push_back( bestVtxRf.covariance(0, 1) );
		   priRfVtxXZE->push_back( bestVtxRf.covariance(0, 2) );
		   priRfVtxYZE->push_back( bestVtxRf.covariance(1, 2) );		   
		   priRfVtxCL->push_back( ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );
		   


		   nB++;	       


		   
		   /////////////////////////////////////////////////



  
	       
	       LambdaParticles.clear();
               pionParticles.clear();
	       muonParticles.clear();
	       vFitMCParticles.clear();
				 
	       //////////////////////////////
	       // Would Check PAT truth match here, but PAT truth match doesn't work because the V0 tracks need to have the momentum defined at the V0 vertex
	       //////////////////////////////


	
	       // }// este es el for para el segundo muon condicion de carga negativa

	       // }// este es el for para el segundo muon

	     }// este es el for para el primer  muon condicion de carga positiva
	     
           }// este es el for para el primer muon

         }// este cierra el for del kaon que pegaremos al lambda para formar el kaskey menos

	   // }// este cierra el primer for del contenedor de los traks 
	        //(en este caso el que usamos para el pion y el proton que forman el lambda). Realmente este no esta.

    }// este cierra el primer for, es de los traks en el contenedor v0

}// este cierra el if para las condicines del contenedor V0



//fill the tree and clear the vectors
   if (nB > 0 ) 
     {
       //std::cout << "filling tree" << endl;
       tree_->Fill();
     }
 
  nB = 0; nMu = 0;

   //triggersL1L->clear();
   triggersMuPL->clear(); triggersMuML->clear();
   triggersL1L2_MuPL->clear(); triggersL1L2_MuML->clear(); 

 // *******************************************************
  B_mass->clear(); B_px->clear(); B_py->clear(); B_pz->clear();
  B_kaskeym_mass->clear(); B_kaskeym_px->clear(); B_kaskeym_py->clear(); B_kaskeym_pz->clear(); B_kaskeym_charge->clear();
  B_J_mass->clear(); B_J_px->clear(); B_J_py->clear(); B_J_pz->clear();

  B_kaskeym_pt1->clear(); B_kaskeym_px1->clear(); B_kaskeym_py1->clear(); B_kaskeym_pz1->clear(); B_kaskeym_charge1->clear();
  B_kaskeym_pt2->clear(); B_kaskeym_px2->clear(); B_kaskeym_py2->clear(); B_kaskeym_pz2->clear(); B_kaskeym_charge2->clear();

  B_J_pt1->clear(); B_J_px1->clear(); B_J_py1->clear(); B_J_pz1->clear(); B_J_charge1->clear();
  B_J_pt2->clear(); B_J_px2->clear(); B_J_py2->clear(); B_J_pz2->clear(); B_J_charge2->clear();

  B_kas_lambda_mass->clear(); B_kas_lambda_px->clear(); B_kas_lambda_py->clear(); B_kas_lambda_pz->clear();
  B_kas_lambda_pt1->clear(); B_kas_lambda_px1->clear(); B_kas_lambda_py1->clear(); B_kas_lambda_pz1->clear(); B_kas_lambda_charge1->clear();
  B_kas_lambda_pt2->clear(); B_kas_lambda_px2->clear(); B_kas_lambda_py2->clear(); B_kas_lambda_pz2->clear(); B_kas_lambda_charge2->clear();
  //B_kaskeym_parentId1->clear(); B_kaskeym_parentId2->clear();
  //B_kaskeym_pId1->clear(); B_kaskeym_pId2->clear();


  B_J_parentId1->clear(); B_J_parentId2->clear(); B_J_muId1->clear(); B_J_muId2->clear();
  B_lam_parentId1->clear(); B_lam_parentId2->clear(); B_lam_PId1->clear(); B_lam_piId2->clear();
  B_kaskeym_kId3->clear();
  B_kaskey_parentId1->clear(); B_kaskey_parentId2->clear(); B_kaskey_parentId3->clear();
  B_parentId1->clear(); B_parentId2->clear(); B_parentId3->clear(); B_parentId4->clear(); B_parentId5->clear();

  B_kaskeym_chi2->clear(); B_J_chi2->clear(); B_chi2->clear(); B_kas_lambda_chi2->clear();
  B_Prob->clear(); B_J_Prob->clear(); B_kaskey_Prob->clear(); B_kas_lambda_Prob->clear();

// ********* 

   nVtx = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   priVtxXBS = 0;     priVtxYBS = 0;     priVtxZBS = 0; 
   priVtxXBSE = 0;    priVtxYBSE = 0;    priVtxZBSE = 0; priVtxCLBS = 0;
   priVtxXYBSE = 0;   priVtxXZBSE = 0;   priVtxYZBSE = 0;

  // todos los vertices primarios 
   /*
   pVtxX->clear();  pVtxY->clear();  pVtxZ->clear();
   pVtxXE->clear(); pVtxYE->clear(); pVtxZE->clear();
   pVtxCL->clear();
   pVtxXYE->clear(); pVtxXZE->clear(); pVtxYZE->clear();
   */

   // vertice primario CON mejor pointin-angle
   pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
   pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
   pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear(); 

   // todos los vertices primarios CON constrain de Beamspot  
 /*
   pVtxBSX->clear();  pVtxBSY->clear();  pVtxBSZ->clear();
   pVtxBSXE->clear(); pVtxBSYE->clear(); pVtxBSZE->clear();
   pVtxBSCL->clear();
   pVtxBSXYE->clear(); pVtxBSXZE->clear(); pVtxBSYZE->clear();
   */

   // vertice primario CON constrain de Beamspot y mejor pointin-angle
   pVtxBSIPX->clear();  pVtxBSIPY->clear();  pVtxBSIPZ->clear();
   pVtxBSIPXE->clear();  pVtxBSIPYE->clear();  pVtxBSIPZE->clear();  pVtxBSIPCL->clear();
   pVtxBSIPXYE->clear();  pVtxBSIPXZE->clear();  pVtxBSIPYZE->clear(); 


   priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
   priRfVtxZE->clear(); priRfVtxXYE->clear(); priRfVtxXZE->clear(); priRfVtxYZE->clear(); priRfVtxCL->clear(); 
   priRfNTrkDif->clear();

   /*
   bctau->clear(); bctau2D->clear(); bctauBS->clear(); bctauBS2D->clear(); bctauRf->clear(); bctau2DRf->clear();
   bctauE->clear(); bctau2DE->clear(); bctauBSE->clear(); bctauBS2DE->clear(); bctauRfE->clear(); bctau2DRfE->clear();
   bctau_kaskey->clear(); bctau2D_kaskey->clear();
   bctauE_kaskey->clear(); bctau2DE_kaskey->clear();
   bctau_lam->clear(); bctau2D_lam->clear();
   bctauE_lam->clear(); bctau2DE_lam->clear();
   */

   PVXBS = 0;     PVYBS = 0;     PVZBS = 0; 
   PVXBSE = 0;    PVYBSE = 0;    PVZBSE = 0; 
   PVXYBSE = 0;   PVXZBSE = 0;   PVYZBSE = 0; 
   
   
   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear(); 
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear(); 
   bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();  

   VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   VDecayVtxXYE->clear(); VDecayVtxXZE->clear(); VDecayVtxYZE->clear();

   V1DecayVtxX->clear(); V1DecayVtxY->clear(); V1DecayVtxZ->clear();
   V1DecayVtxXE->clear(); V1DecayVtxYE->clear(); V1DecayVtxZE->clear();
   V1DecayVtxXYE->clear(); V1DecayVtxXZE->clear(); V1DecayVtxYZE->clear();
   
   JDecayVtxX->clear(); JDecayVtxY->clear(); JDecayVtxZ->clear();
   JDecayVtxXE->clear(); JDecayVtxYE->clear(); JDecayVtxZE->clear(); 
   JDecayVtxXYE->clear(); JDecayVtxXZE->clear(); JDecayVtxYZE->clear(); 

   
   mumC2->clear();
   mumCat->clear(); mumAngT->clear(); mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupCat->clear(); mupAngT->clear(); mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 

}// este cierra la funcion analize


// estas dos primeras funciones son para los vertices y los pt a nivel generacion, es decir, para cuando estamos corriendo el MC.
// en particular la primera es para el J/Psi y la segunda es para el Lambda_0
/*
void JPsikaskeyPAT::fillPsi(const reco::Candidate& genpsi) {
  
    for (uint i=0; i<genpsi.numberOfDaughters(); i++) {
    if (genpsi.daughter(i)->pdgId()==13) { //13 is a mu-
      trueMumPx = genpsi.daughter(i)->px();
      trueMumPy = genpsi.daughter(i)->py();
      trueMumPz = genpsi.daughter(i)->pz();
    }
    if (genpsi.daughter(i)->pdgId()==-13) { //-13 is a mu+
      trueMupPx = genpsi.daughter(i)->px();
      trueMupPy = genpsi.daughter(i)->py();
      trueMupPz = genpsi.daughter(i)->pz();
    }
  }
}

void JPsikaskeyPAT::fillV0(const reco::Candidate& genv0) {
  
  for (uint i=0; i<genv0.numberOfDaughters(); i++) {
    if (genv0.daughter(i)->charge()>0 && genv0.numberOfDaughters()==2) {
      trueVTrkpPx = genv0.daughter(i)->px();
      trueVTrkpPy = genv0.daughter(i)->py();
      trueVTrkpPz = genv0.daughter(i)->pz();
      trueVTrkpMass = genv0.daughter(i)->mass();
      trueVDecayVtxX = genv0.daughter(i)->vx();
      trueVDecayVtxY = genv0.daughter(i)->vy();
      trueVDecayVtxZ = genv0.daughter(i)->vz();
    }
    if (genv0.daughter(i)->charge()<0 && genv0.numberOfDaughters()==2) {
      trueVTrkmPx = genv0.daughter(i)->px();
      trueVTrkmPy = genv0.daughter(i)->py();
      trueVTrkmPz = genv0.daughter(i)->pz();
      trueVTrkmMass = genv0.daughter(i)->mass();
    }
  }
}
*/


int const JPsikaskeyPAT::getMuCat(reco::Muon const& muon) const{
  int muCat = 0;
  if (muon.isGlobalMuon()) {
    if (muon.isTrackerMuon()) muCat = 1;
    else muCat = 10;
  }
  else if (muon.isTrackerMuon()) {
    if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)) {
      if (muon::isGoodMuon(muon, muon::TMLastStationAngTight)) muCat = 6;
      else if (muon::isGoodMuon(muon, muon::TMOneStationTight)) muCat = 5;
      else if (muon::isGoodMuon(muon, muon::TMOneStationLoose)) muCat = 4;
      else muCat = 3;
    } else muCat = 2;
  }
  else if (muon.isStandAloneMuon()) muCat = 7;
  else if (muon.isCaloMuon()) muCat = 8;
  else muCat = 9;
  
  if ( !(muon::isGoodMuon(muon, muon::TMOneStationLoose)) && muon::isGoodMuon(muon, muon::TMOneStationTight) )
    std::cout << "inconsistent muon cat 1" << std::endl;
  if ( !(muon::isGoodMuon(muon, muon::TMOneStationTight)) && muon::isGoodMuon(muon, muon::TMLastStationAngTight) )
    std::cout << "inconsistent muon cat 2" << std::endl;

  return muCat;
}




bool const JPsikaskeyPAT::HasGoodME11(reco::Muon const& muon, double const dxdzCut) const{
  bool retVal = false;
  for(std::vector<reco::MuonChamberMatch>::const_iterator mcm = muon.matches().begin();
    mcm != muon.matches().end(); ++mcm) {
    DetId const& chamberId = mcm->id;
    if (chamberId.det() != DetId::Muon) continue;
    if (chamberId.subdetId() != MuonSubdetId::CSC) continue;
    CSCDetId id(chamberId.rawId());
    if (id.station() != 1) continue;
    if (fabs(mcm->dXdZ) > dxdzCut) continue;
    retVal = true;
  }
  return retVal;
}



// ------------ method called once each job just before starting event loop  ------------

void 
JPsikaskeyPAT::beginJob()
{


  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->J/psi kaskey menos ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_kaskeym_mass", &B_kaskeym_mass);
  tree_->Branch("B_kaskeym_px", &B_kaskeym_px);
  tree_->Branch("B_kaskeym_py", &B_kaskeym_py);
  tree_->Branch("B_kaskeym_pz", &B_kaskeym_pz);
  tree_->Branch("B_kaskeym_charge", &B_kaskeym_charge);

  //tree_->Branch("B_kaskeym_parentId1", &B_kaskeym_parentId1);
  //tree_->Branch("B_kaskeym_parentId2", &B_kaskeym_parentId2);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_kas_lambda_mass", &B_kas_lambda_mass);
  tree_->Branch("B_kas_lambda_px", &B_kas_lambda_px);
  tree_->Branch("B_kas_lambda_py", &B_kas_lambda_py);
  tree_->Branch("B_kas_lambda_pz", &B_kas_lambda_pz);


  tree_->Branch("B_kaskeym_pt1", &B_kaskeym_pt1);
  tree_->Branch("B_kaskeym_px1", &B_kaskeym_px1);
  tree_->Branch("B_kaskeym_py1", &B_kaskeym_py1);
  tree_->Branch("B_kaskeym_pz1", &B_kaskeym_pz1);
  tree_->Branch("B_kaskeym_charge1", &B_kaskeym_charge1); 
  //tree_->Branch("B_kaskeym_pId1", &B_kaskeym_pId1);
 
  tree_->Branch("B_kaskeym_pt2", &B_kaskeym_pt2);
  tree_->Branch("B_kaskeym_px2", &B_kaskeym_px2);
  tree_->Branch("B_kaskeym_py2", &B_kaskeym_py2);
  tree_->Branch("B_kaskeym_pz2", &B_kaskeym_pz2);
  tree_->Branch("B_kaskeym_charge2", &B_kaskeym_charge2);
  // tree_->Branch("B_kaskeym_pId2", &B_kaskeym_pId2);

  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);


  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_kas_lambda_pt1", &B_kas_lambda_pt1);
  tree_->Branch("B_kas_lambda_px1", &B_kas_lambda_px1);
  tree_->Branch("B_kas_lambda_py1", &B_kas_lambda_py1);
  tree_->Branch("B_kas_lambda_pz1", &B_kas_lambda_pz1);
  tree_->Branch("B_kas_lambda_charge1", &B_kas_lambda_charge1);

  tree_->Branch("B_kas_lambda_pt2", &B_kas_lambda_pt2);
  tree_->Branch("B_kas_lambda_px2", &B_kas_lambda_px2);
  tree_->Branch("B_kas_lambda_py2", &B_kas_lambda_py2);
  tree_->Branch("B_kas_lambda_pz2", &B_kas_lambda_pz2);
  tree_->Branch("B_kas_lambda_charge2", &B_kas_lambda_charge2);

  tree_->Branch("B_J_parentId1", &B_J_parentId1);
  tree_->Branch("B_J_parentId2", &B_J_parentId2);
  tree_->Branch("B_J_muId1", &B_J_muId1);
  tree_->Branch("B_J_muId2", &B_J_muId2);
  tree_->Branch("B_lam_parentId1", &B_lam_parentId1);
  tree_->Branch("B_lam_parentId2", &B_lam_parentId2);
  tree_->Branch("B_lam_PId1", &B_lam_PId1);
  tree_->Branch("B_lam_piId2", &B_lam_piId2);
  tree_->Branch("B_kaskeym_kId3", &B_kaskeym_kId3);
  tree_->Branch("B_kaskey_parentId1", &B_kaskey_parentId1);
  tree_->Branch("B_kaskey_parentId2", &B_kaskey_parentId2);
  tree_->Branch("B_kaskey_parentId3", &B_kaskey_parentId3);
  tree_->Branch("B_parentId1", &B_parentId1);
  tree_->Branch("B_parentId2", &B_parentId2);
  tree_->Branch("B_parentId3", &B_parentId3);
  tree_->Branch("B_parentId4", &B_parentId4);
  tree_->Branch("B_parentId5", &B_parentId5);
  



  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("B_kaskeym_chi2", &B_kaskeym_chi2);
  tree_->Branch("B_J_chi2", &B_J_chi2);
  tree_->Branch("B_kas_lambda_chi2", &B_kas_lambda_chi2);
  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_kaskey_Prob", &B_kaskey_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
  tree_->Branch("B_kas_lambda_Prob", &B_kas_lambda_Prob);
  // *************************

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");
 
  tree_->Branch("priVtxXBS",&priVtxXBS, "priVtxXBS/f");
  tree_->Branch("priVtxYBS",&priVtxYBS, "priVtxYBS/f");
  tree_->Branch("priVtxZBS",&priVtxZBS, "priVtxZBS/f");
  tree_->Branch("priVtxXBSE",&priVtxXBSE, "priVtxXBSE/f");
  tree_->Branch("priVtxYBSE",&priVtxYBSE, "priVtxYBSE/f");
  tree_->Branch("priVtxZBSE",&priVtxZBSE, "priVtxZBSE/f");
  tree_->Branch("priVtxXYBSE",&priVtxXYBSE, "priVtxXYBSE/f");
  tree_->Branch("priVtxXZBSE",&priVtxXZBSE, "priVtxXZBSE/f");
  tree_->Branch("priVtxYZBSE",&priVtxYZBSE, "priVtxYZBSE/f");
  tree_->Branch("priVtxCLBS",&priVtxCLBS, "priVtxCLBS/f");

  /*
  tree_->Branch("pVtxX",     &pVtxX);
  tree_->Branch("pVtxY",     &pVtxY);
  tree_->Branch("pVtxZ",     &pVtxZ);
  tree_->Branch("pVtxXE",     &pVtxXE);
  tree_->Branch("pVtxYE",     &pVtxYE);
  tree_->Branch("pVtxZE",     &pVtxZE);
  tree_->Branch("pVtxXYE",     &pVtxXYE);
  tree_->Branch("pVtxXZE",     &pVtxXZE);
  tree_->Branch("pVtxYZE",     &pVtxYZE);
  tree_->Branch("pVtxCL",     &pVtxCL);
  */

  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);

  /*
  tree_->Branch("pVtxBSX",     &pVtxBSX);
  tree_->Branch("pVtxBSY",     &pVtxBSY);
  tree_->Branch("pVtxBSZ",     &pVtxBSZ);
  tree_->Branch("pVtxBSXE",     &pVtxBSXE);
  tree_->Branch("pVtxBSYE",     &pVtxBSYE);
  tree_->Branch("pVtxBSZE",     &pVtxBSZE);
  tree_->Branch("pVtxBSXYE",     &pVtxBSXYE);
  tree_->Branch("pVtxBSXZE",     &pVtxBSXZE);
  tree_->Branch("pVtxBSYZE",     &pVtxBSYZE);
  tree_->Branch("pVtxBSCL",     &pVtxBSCL);
  */

  tree_->Branch("pVtxBSIPX",     &pVtxBSIPX);
  tree_->Branch("pVtxBSIPY",     &pVtxBSIPY);
  tree_->Branch("pVtxBSIPZ",     &pVtxBSIPZ);
  tree_->Branch("pVtxBSIPXE",     &pVtxBSIPXE);
  tree_->Branch("pVtxBSIPYE",     &pVtxBSIPYE);
  tree_->Branch("pVtxBSIPZE",     &pVtxBSIPZE);
  tree_->Branch("pVtxBSIPXYE",     &pVtxBSIPXYE);
  tree_->Branch("pVtxBSIPXZE",     &pVtxBSIPXZE);
  tree_->Branch("pVtxBSIPYZE",     &pVtxBSIPYZE);
  tree_->Branch("pVtxBSIPCL",     &pVtxBSIPCL);


  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");

  /*
  tree_->Branch("nTrgL",      &nTrgL,     "nTrgL/I");
  tree_->Branch("triggers",   &triggers,  "triggers[nTrgL]/C");
  tree_->Branch("nMuonPTrgL",  &nMuonPTrgL,  "nMuonPTrgL/I");
  tree_->Branch("triggersMuP", &triggersMuP, "triggersMuP[nMuonPTrgL]/C");
  
  tree_->Branch("nMuonMTrgL",  &nMuonMTrgL,  "nMuonMTrgL/I");
  tree_->Branch("triggersMuM", &triggersMuM, "triggersMuM[nMuonMTrgL]/C");
  tree_->Branch("triggersL1", &triggersL1, "triggersL1/C");
  tree_->Branch("ntriggersL1L2_MuP", &ntriggersL1L2_MuP, "ntriggersL1L2_MuP/I");
  tree_->Branch("triggersL1L2_MuP",  &triggersL1L2_MuP,  "triggersL1L2_MuP/C");
  
  tree_->Branch("ntriggersL1L2_MuM", &ntriggersL1L2_MuM, "ntriggersL1L2_MuM/I");
  tree_->Branch("triggersL1L2_MuM",  &triggersL1L2_MuM,  "triggersL1L2_MuM/C");
  */
  tree_->Branch("nTrgL",      &nTrgL,    "nTrgL/I");  
  tree_->Branch("nTrgL1L",    &nTrgL1L,  "nTrgL1L/I");  

  tree_->Branch("triggersL",         &triggersL,  "triggersL[nTrgL]/C");
  tree_->Branch("triggersL1L",       &triggersL1L,"triggersL1L[nTrgL1L]/C");
  tree_->Branch("triggersMuPL",      &triggersMuPL);
  tree_->Branch("triggersMuML",      &triggersMuML);
  tree_->Branch("triggersL1L2_MuPL", &triggersL1L2_MuPL);
  tree_->Branch("triggersL1L2_MuML", &triggersL1L2_MuML);


 
  tree_->Branch("priRfVtxX",&priRfVtxX);
  tree_->Branch("priRfVtxY",&priRfVtxY);
  tree_->Branch("priRfVtxZ",&priRfVtxZ);
  tree_->Branch("priRfVtxXE",&priRfVtxXE);
  tree_->Branch("priRfVtxYE",&priRfVtxYE);
  tree_->Branch("priRfVtxZE",&priRfVtxZE);
  tree_->Branch("priRfVtxXYE",&priRfVtxXYE);
  tree_->Branch("priRfVtxXZE",&priRfVtxXZE);
  tree_->Branch("priRfVtxYZE",&priRfVtxYZE);
  tree_->Branch("priRfVtxCL",&priRfVtxCL);
  tree_->Branch("priRfNTrkDif",&priRfNTrkDif);
 
  /*
  tree_->Branch("bctau",&bctau);
  tree_->Branch("bctau2D",&bctau2D);
  tree_->Branch("bctauBS",&bctauBS);
  tree_->Branch("bctauBS2D",&bctauBS2D);
  tree_->Branch("bctauE",&bctauE);
  tree_->Branch("bctau2DE",&bctau2DE);
  tree_->Branch("bctauBSE",&bctauBSE);
  tree_->Branch("bctauBS2DE",&bctauBS2DE);
  tree_->Branch("bctauRf",&bctauRf);
  tree_->Branch("bctauRfE",&bctauRfE);
  tree_->Branch("bctau2DRf",&bctau2DRf);
  tree_->Branch("bctau2DRfE",&bctau2DRfE);
  tree_->Branch("bctau_kaskey",&bctau_kaskey);
  tree_->Branch("bctau2D_kaskey",&bctau2D_kaskey);
  tree_->Branch("bctauE_kaskey",&bctauE_kaskey);
  tree_->Branch("bctau2DE_kaskey",&bctau2DE_kaskey);
  tree_->Branch("bctau_lam",&bctau_lam);
  tree_->Branch("bctau2D_lam",&bctau2D_lam);
  tree_->Branch("bctauE_lam",&bctauE_lam);
  tree_->Branch("bctau2DE_lam",&bctau2DE_lam);
  */

  tree_->Branch("PVXBS",&PVXBS, "PVXBS/D");
  tree_->Branch("PVYBS",&PVYBS, "PVYBS/D");
  tree_->Branch("PVZBS",&PVZBS, "PVZBS/D");
  tree_->Branch("PVXBSE",&PVXBSE, "PVXBSE/D");
  tree_->Branch("PVYBSE",&PVYBSE, "PVYBSE/D");
  tree_->Branch("PVZBSE",&PVZBSE, "PVZBSE/D");
  tree_->Branch("PVXYBSE",&PVXYBSE, "PVXYBSE/D");
  tree_->Branch("PVXZBSE",&PVXZBSE, "PVXZBSE/D");
  tree_->Branch("PVYZBSE",&PVYZBSE, "PVYZBSE/D");


  tree_->Branch("bDecayVtxX",&bDecayVtxX);
  tree_->Branch("bDecayVtxY",&bDecayVtxY);
  tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
  tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
  tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
  tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
  tree_->Branch("bDecayVtxXYE",&bDecayVtxXYE);
  tree_->Branch("bDecayVtxXZE",&bDecayVtxXZE);
  tree_->Branch("bDecayVtxYZE",&bDecayVtxYZE);
 

  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("VDecayVtxXYE",&VDecayVtxXYE);
  tree_->Branch("VDecayVtxXZE",&VDecayVtxXZE);
  tree_->Branch("VDecayVtxYZE",&VDecayVtxYZE);
 
  tree_->Branch("V1DecayVtxX",&V1DecayVtxX);
  tree_->Branch("V1DecayVtxY",&V1DecayVtxY);
  tree_->Branch("V1DecayVtxZ",&V1DecayVtxZ);
  tree_->Branch("V1DecayVtxXE",&V1DecayVtxXE);
  tree_->Branch("V1DecayVtxYE",&V1DecayVtxYE);
  tree_->Branch("V1DecayVtxZE",&V1DecayVtxZE);
  tree_->Branch("V1DecayVtxXYE",&V1DecayVtxXYE);
  tree_->Branch("V1DecayVtxXZE",&V1DecayVtxXZE);
  tree_->Branch("V1DecayVtxYZE",&V1DecayVtxYZE);


  tree_->Branch("JDecayVtxX",&JDecayVtxX);
  tree_->Branch("JDecayVtxY",&JDecayVtxY);
  tree_->Branch("JDecayVtxZ",&JDecayVtxZ);
  tree_->Branch("JDecayVtxXE",&JDecayVtxXE);
  tree_->Branch("JDecayVtxYE",&JDecayVtxYE);
  tree_->Branch("JDecayVtxZE",&JDecayVtxZE);
  tree_->Branch("JDecayVtxXYE",&JDecayVtxXYE);
  tree_->Branch("JDecayVtxXZE",&JDecayVtxXZE);
  tree_->Branch("JDecayVtxYZE",&JDecayVtxYZE); 
 
  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumCat",&mumCat);
  tree_->Branch("mumAngT",&mumAngT);
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupCat",&mupCat);
  tree_->Branch("mupAngT",&mupAngT);
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

}


// ------------ method called once each job just after ending the event loop  ------------
void JPsikaskeyPAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsikaskeyPAT);