#ifndef _JPsikaskeyPAT_h
#define _JPsikaskeyPAT_h

// user include files
//#include "myAnalyzers/JPsiKsPAT/interface/JPsif0PAT.h"

#include <memory>

// user include files
//#include "myAnalyzers/JPsiLambdaPAT/interface/JPsiLambdaPAT.h"
#include "myAnalyzers/JPsiKsPAT/interface/JPsikaskeyPAT.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"


#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"


#include "TFile.h"
#include "TTree.h"




//
// class decleration
//

class 	JPsikaskeyPAT	 : public edm::EDAnalyzer {
public:
  explicit 	JPsikaskeyPAT	(const edm::ParameterSet&);
  ~	JPsikaskeyPAT	();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
  bool const HasGoodME11(reco::Muon const& muon, double const dxdzCut) const;
  
  /*
  private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;
  
  */



private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

 void CheckL1Triggers(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::string &TrigListNameL1Tmp);
  void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);
 void MatchMuonWithL1L2(const pat::Muon &iMuon, const std::vector<std::string>& TrigListL1L2, std::string &TrigListNameL1L2Tmp);

   
 //int PdgIDatTruthLevel(const reco::Track Track,const edm::Event& iEvent);
  int PdgIDatTruthLevel(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles, int &ParentID, int &PParentID);
  int PdgIDatTruthLevel4(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles, int &ParentID, int &PParentID,int &GPParentID);

  // void ParticleParent(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles);
float Myctau(const RefCountedKinematicParticle &CandMC, const RefCountedKinematicVertex &DecayVertexMC, 
	       const GlobalPoint &PVPtmp, const GlobalError &PVEtmp,float mass_tmp, 
	       float &ctau2Dtmp, float &ctauEtemp, float &ctauE2Dtemp );


  // ----------member data ---------------------------
  /*  
std::string hlTriggerResults_;
  std::string vtxSample;
  std::string genParticles_;
  std::string V0Collection_;
  std::string muonType;
  std::string muonTypeForPAT;
  */
  std::string hltTriggerResults_;
  std::string v0Producer;
  std::string v0Type;
  std::string vtxSample;
  std::string genParticles_;
  std::string v0Collection_;
  std::string muonType;
  std::string muonTypeForPAT;
  bool doMC_;
  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;


  
  std::vector<float>       *priRfVtxX, *priRfVtxY, *priRfVtxZ, *priRfVtxXE, *priRfVtxYE, *priRfVtxZE, *priRfVtxCL;
  std::vector<float>       *priRfVtxXYE, *priRfVtxXZE, *priRfVtxYZE;
  std::vector<int>         *priRfNTrkDif;

  //std::vector<float>       *bctau, *bctau2D, *bctauBS, *bctauBS2D, *bctauRf, *bctau2DRf;
  //std::vector<float>       *bctauE, *bctau2DE, *bctauBSE, *bctauBS2DE, *bctauRfE, *bctau2DRfE;
  //std::vector<float>       *bctau_kaskey, *bctau2D_kaskey;
  //std::vector<float>       *bctauE_kaskey, *bctau2DE_kaskey;
  //std::vector<float>       *bctau_lam, *bctau2D_lam;
  //std::vector<float>       *bctauE_lam, *bctau2DE_lam;

  std::vector<float>       *mumC2;
  std::vector<int>         *mumCat, *mumAngT, *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupCat, *mupAngT, *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;

  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  

  // estas variables son para los vertices y los pt a nivel generacion, es decir, para cuando estamos corriendo el MC.
  // ver funciones fillPsi y fillV0 en el .cc
  //std::vector<float>       *VTrkpMass, *VTrkpPx, *VTrkpPy, *VTrkpPz, *VTrkpD0, *VTrkpD0E;
  //std::vector<float>       *VTrkmMass, *VTrkmPx, *VTrkmPy, *VTrkmPz, *VTrkmD0, *VTrkmD0E;
  
  int                 muAcc, muTrig, weight;

  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;

  float                    priVtxXBS, priVtxYBS, priVtxZBS, priVtxXBSE, priVtxYBSE, priVtxZBSE, priVtxCLBS;
  float                    priVtxXYBSE, priVtxXZBSE, priVtxYZBSE;

  // todos los vertices primarios 
  //std::vector<float>       *pVtxX,  *pVtxY,  *pVtxZ,  *pVtxXE,  *pVtxYE,  *pVtxZE,  *pVtxCL;
  //std::vector<float>       *pVtxXYE, *pVtxXZE, *pVtxYZE;

  // vertice primario CON mejor pointin-angle
   std::vector<float>          *pVtxIPX,  *pVtxIPY, *pVtxIPZ, *pVtxIPXE, *pVtxIPYE, *pVtxIPZE, *pVtxIPCL;
   std::vector<float>          *pVtxIPXYE,  *pVtxIPXZE, *pVtxIPYZE;

  // todos los vertices primarios CON constrain de Beamspot
  //std::vector<float>       *pVtxBSX,  *pVtxBSY,  *pVtxBSZ,  *pVtxBSXE,  *pVtxBSYE,  *pVtxBSZE,  *pVtxBSCL;
  //std::vector<float>       *pVtxBSXYE, *pVtxBSXZE, *pVtxBSYZE;

  // vertice primario CON constrain de Beamspot y mejor pointin-angle
   std::vector<float>        *pVtxBSIPX,  *pVtxBSIPY,  *pVtxBSIPZ, *pVtxBSIPXE, *pVtxBSIPYE, *pVtxBSIPZE, *pVtxBSIPCL;
   std::vector<float>        *pVtxBSIPXYE,  *pVtxBSIPXZE,  *pVtxBSIPYZE;

 // ************ esta es la informacion concerniente al beamspot *************
 //CovarianceMatrix          covarianceBS;
  // float                  PVXBS, PVYBS, PVZBS, PVXBSE, PVYBSE, PVZBSE, PVCLBS;
           
  double                    PVXBS, PVYBS, PVZBS, PVXBSE, PVYBSE, PVZBSE;
  double                    PVXYBSE, PVXZBSE, PVYZBSE;
 
  // ********************************** ************************************************************************

  std::vector<float>       *bDecayVtxX, *bDecayVtxY, *bDecayVtxZ;
  std::vector<double>      *bDecayVtxXE, *bDecayVtxYE, *bDecayVtxZE;
  std::vector<double>      *bDecayVtxXYE, *bDecayVtxXZE, *bDecayVtxYZE;

  std::vector<float>       *VDecayVtxX, *VDecayVtxY, *VDecayVtxZ;
  std::vector<float>       *VDecayVtxXE, *VDecayVtxYE, *VDecayVtxZE;
  std::vector<float>       *VDecayVtxXYE, *VDecayVtxXZE, *VDecayVtxYZE;

  std::vector<float>       *V1DecayVtxX, *V1DecayVtxY, *V1DecayVtxZ;
  std::vector<float>       *V1DecayVtxXE, *V1DecayVtxYE, *V1DecayVtxZE;
  std::vector<float>       *V1DecayVtxXYE,  *V1DecayVtxXZE, *V1DecayVtxYZE;

  std::vector<float>       *JDecayVtxX, *JDecayVtxY, *JDecayVtxZ;
  std::vector<float>       *JDecayVtxXE, *JDecayVtxYE, *JDecayVtxZE;
  std::vector<float>       *JDecayVtxXYE, *JDecayVtxXZE,  *JDecayVtxYZE;

 
  // *************************************

  unsigned int             nB;
  unsigned int             nMu;


  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz;
  
  std::vector<float>       *B_kaskeym_mass, *B_kaskeym_px, *B_kaskeym_py, *B_kaskeym_pz, *B_kaskeym_charge;

  std::vector<float>       *B_kaskeym_pt1, *B_kaskeym_px1, *B_kaskeym_py1, *B_kaskeym_pz1;
  std::vector<float>       *B_kaskeym_pt2, *B_kaskeym_px2, *B_kaskeym_py2, *B_kaskeym_pz2;
  std::vector<int>         *B_kaskeym_charge1, *B_kaskeym_charge2;
  
  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;
  std::vector<float>       *B_J_pt1, *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_pt2, *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

  std::vector<float>       *B_kas_lambda_mass, *B_kas_lambda_px, *B_kas_lambda_py, *B_kas_lambda_pz;
  std::vector<float>       *B_kas_lambda_pt1, *B_kas_lambda_px1, *B_kas_lambda_py1, *B_kas_lambda_pz1;
  std::vector<float>       *B_kas_lambda_pt2, *B_kas_lambda_px2, *B_kas_lambda_py2, *B_kas_lambda_pz2;
  std::vector<int>         *B_kas_lambda_charge1, *B_kas_lambda_charge2;

  std::vector<int>         *B_J_parentId1, *B_J_parentId2;
  std::vector<int>         *B_J_muId1, *B_J_muId2;
  std::vector<int>         *B_lam_parentId1, *B_lam_parentId2;
  std::vector<int>         *B_lam_PId1, *B_lam_piId2, *B_kaskeym_kId3; 
  std::vector<int>         *B_kaskey_parentId1, *B_kaskey_parentId2, *B_kaskey_parentId3;
  std::vector<int>         *B_parentId1, *B_parentId2, *B_parentId3, *B_parentId4, *B_parentId5;
  //std::vector<int>         *B_kaskeym_parentId1, *B_kaskeym_parentId2;
  //std::vector<int>         *B_kaskeym_pId1, *B_kaskeym_pId2;


  //std::vector<float>       *B_kas_lambda_trkP_d0, *B_kas_lambda_trkP_d0E, *B_kas_lambda_trkM_d0, *B_kas_lambda_trkM_d0E;
  //std::vector<float>       *B_kas_pion_d0, *B_kas_pion_d0E;
  
  std::vector<float>       *B_kaskeym_chi2, *B_J_chi2, *B_chi2, *B_kas_lambda_chi2;
  std::vector<float>       *B_Prob, *B_J_Prob, *B_kaskey_Prob, *B_kas_lambda_Prob;


  char triggersL[10000], triggersL1L[10000];

  char triggersMuP[10000],     triggersMuM[10000] ;
  char triggersL1L2_MuP[10000],triggersL1L2_MuM[10000];
  
  char triggersL1[10000];
  int  nTrgL, nTrgL1L,  nMuonTrgL,  nMuonPTrgL,        nMuonMTrgL;
  int  ntriggersL1L2_MuP, ntriggersL1L2_MuM;


  //std::vector<std::string>         *triggersL1L;
  std::vector<std::string>         *triggersMuPL, *triggersMuML;
  std::vector<std::string>         *triggersL1L2_MuPL, *triggersL1L2_MuML;

  int  run, event;
 

};

#endif