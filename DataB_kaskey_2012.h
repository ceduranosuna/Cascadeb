//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 27 17:14:05 2016 by ROOT version 5.34/36
// from TTree ntuple/Bs->J/psi kaskey menos ntuple
// found on file: ../roofilescarpetaC/output_parkedata2012_kaskey_3912.root
//////////////////////////////////////////////////////////

#ifndef DataB_kaskey_2012_h
#define DataB_kaskey_2012_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class DataB_kaskey_2012 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          nB;
   UInt_t          nMu;
   vector<float>   *B_mass;
   vector<float>   *B_px;
   vector<float>   *B_py;
   vector<float>   *B_pz;
   vector<float>   *B_kaskeym_mass;
   vector<float>   *B_kaskeym_px;
   vector<float>   *B_kaskeym_py;
   vector<float>   *B_kaskeym_pz;
   vector<float>   *B_kaskeym_charge;
   vector<float>   *B_J_mass;
   vector<float>   *B_J_px;
   vector<float>   *B_J_py;
   vector<float>   *B_J_pz;
   vector<float>   *B_kas_lambda_mass;
   vector<float>   *B_kas_lambda_px;
   vector<float>   *B_kas_lambda_py;
   vector<float>   *B_kas_lambda_pz;
   vector<float>   *B_kaskeym_pt1;
   vector<float>   *B_kaskeym_px1;
   vector<float>   *B_kaskeym_py1;
   vector<float>   *B_kaskeym_pz1;
   vector<int>     *B_kaskeym_charge1;
   vector<float>   *B_kaskeym_pt2;
   vector<float>   *B_kaskeym_px2;
   vector<float>   *B_kaskeym_py2;
   vector<float>   *B_kaskeym_pz2;
   vector<int>     *B_kaskeym_charge2;
   vector<float>   *B_J_pt1;
   vector<float>   *B_J_px1;
   vector<float>   *B_J_py1;
   vector<float>   *B_J_pz1;
   vector<int>     *B_J_charge1;
   vector<float>   *B_J_pt2;
   vector<float>   *B_J_px2;
   vector<float>   *B_J_py2;
   vector<float>   *B_J_pz2;
   vector<int>     *B_J_charge2;
   vector<float>   *B_kas_lambda_pt1;
   vector<float>   *B_kas_lambda_px1;
   vector<float>   *B_kas_lambda_py1;
   vector<float>   *B_kas_lambda_pz1;
   vector<int>     *B_kas_lambda_charge1;
   vector<float>   *B_kas_lambda_pt2;
   vector<float>   *B_kas_lambda_px2;
   vector<float>   *B_kas_lambda_py2;
   vector<float>   *B_kas_lambda_pz2;
   vector<int>     *B_kas_lambda_charge2;
   vector<int>     *B_J_parentId1;
   vector<int>     *B_J_parentId2;
   vector<int>     *B_J_muId1;
   vector<int>     *B_J_muId2;
   vector<int>     *B_lam_parentId1;
   vector<int>     *B_lam_parentId2;
   vector<int>     *B_lam_PId1;
   vector<int>     *B_lam_piId2;
   vector<int>     *B_kaskeym_kId3;
   vector<int>     *B_kaskey_parentId1;
   vector<int>     *B_kaskey_parentId2;
   vector<int>     *B_kaskey_parentId3;
   vector<int>     *B_parentId1;
   vector<int>     *B_parentId2;
   vector<int>     *B_parentId3;
   vector<int>     *B_parentId4;
   vector<int>     *B_parentId5;
   vector<float>   *B_chi2;
   vector<float>   *B_kaskeym_chi2;
   vector<float>   *B_J_chi2;
   vector<float>   *B_kas_lambda_chi2;
   vector<float>   *B_Prob;
   vector<float>   *B_kaskey_Prob;
   vector<float>   *B_J_Prob;
   vector<float>   *B_kas_lambda_Prob;
   Float_t         priVtxX;
   Float_t         priVtxY;
   Float_t         priVtxZ;
   Float_t         priVtxXE;
   Float_t         priVtxYE;
   Float_t         priVtxZE;
   Float_t         priVtxXYE;
   Float_t         priVtxXZE;
   Float_t         priVtxYZE;
   Float_t         priVtxCL;
   Float_t         priVtxXBS;
   Float_t         priVtxYBS;
   Float_t         priVtxZBS;
   Float_t         priVtxXBSE;
   Float_t         priVtxYBSE;
   Float_t         priVtxZBSE;
   Float_t         priVtxXYBSE;
   Float_t         priVtxXZBSE;
   Float_t         priVtxYZBSE;
   Float_t         priVtxCLBS;
   vector<float>   *pVtxIPX;
   vector<float>   *pVtxIPY;
   vector<float>   *pVtxIPZ;
   vector<float>   *pVtxIPXE;
   vector<float>   *pVtxIPYE;
   vector<float>   *pVtxIPZE;
   vector<float>   *pVtxIPXYE;
   vector<float>   *pVtxIPXZE;
   vector<float>   *pVtxIPYZE;
   vector<float>   *pVtxIPCL;
   vector<float>   *pVtxBSIPX;
   vector<float>   *pVtxBSIPY;
   vector<float>   *pVtxBSIPZ;
   vector<float>   *pVtxBSIPXE;
   vector<float>   *pVtxBSIPYE;
   vector<float>   *pVtxBSIPZE;
   vector<float>   *pVtxBSIPXYE;
   vector<float>   *pVtxBSIPXZE;
   vector<float>   *pVtxBSIPYZE;
   vector<float>   *pVtxBSIPCL;
   UInt_t          nVtx;
   Int_t           run;
   Int_t           event;
   Int_t           nTrgL;
   Int_t           nTrgL1L;
   Char_t          triggersL[136];   //[nTrgL]
   Char_t          triggersL1L[2314];   //[nTrgL1L]
   vector<string>  *triggersMuPL;
   vector<string>  *triggersMuML;
   vector<string>  *triggersL1L2_MuPL;
   vector<string>  *triggersL1L2_MuML;
   vector<float>   *priRfVtxX;
   vector<float>   *priRfVtxY;
   vector<float>   *priRfVtxZ;
   vector<float>   *priRfVtxXE;
   vector<float>   *priRfVtxYE;
   vector<float>   *priRfVtxZE;
   vector<float>   *priRfVtxXYE;
   vector<float>   *priRfVtxXZE;
   vector<float>   *priRfVtxYZE;
   vector<float>   *priRfVtxCL;
   vector<int>     *priRfNTrkDif;
   Double_t        PVXBS;
   Double_t        PVYBS;
   Double_t        PVZBS;
   Double_t        PVXBSE;
   Double_t        PVYBSE;
   Double_t        PVZBSE;
   Double_t        PVXYBSE;
   Double_t        PVXZBSE;
   Double_t        PVYZBSE;
   vector<float>   *bDecayVtxX;
   vector<float>   *bDecayVtxY;
   vector<float>   *bDecayVtxZ;
   vector<double>  *bDecayVtxXE;
   vector<double>  *bDecayVtxYE;
   vector<double>  *bDecayVtxZE;
   vector<double>  *bDecayVtxXYE;
   vector<double>  *bDecayVtxXZE;
   vector<double>  *bDecayVtxYZE;
   vector<float>   *VDecayVtxX;
   vector<float>   *VDecayVtxY;
   vector<float>   *VDecayVtxZ;
   vector<float>   *VDecayVtxXE;
   vector<float>   *VDecayVtxYE;
   vector<float>   *VDecayVtxZE;
   vector<float>   *VDecayVtxXYE;
   vector<float>   *VDecayVtxXZE;
   vector<float>   *VDecayVtxYZE;
   vector<float>   *V1DecayVtxX;
   vector<float>   *V1DecayVtxY;
   vector<float>   *V1DecayVtxZ;
   vector<float>   *V1DecayVtxXE;
   vector<float>   *V1DecayVtxYE;
   vector<float>   *V1DecayVtxZE;
   vector<float>   *V1DecayVtxXYE;
   vector<float>   *V1DecayVtxXZE;
   vector<float>   *V1DecayVtxYZE;
   vector<float>   *JDecayVtxX;
   vector<float>   *JDecayVtxY;
   vector<float>   *JDecayVtxZ;
   vector<float>   *JDecayVtxXE;
   vector<float>   *JDecayVtxYE;
   vector<float>   *JDecayVtxZE;
   vector<float>   *JDecayVtxXYE;
   vector<float>   *JDecayVtxXZE;
   vector<float>   *JDecayVtxYZE;
   vector<float>   *mumC2;
   vector<int>     *mumCat;
   vector<int>     *mumAngT;
   vector<int>     *mumNHits;
   vector<int>     *mumNPHits;
   vector<float>   *mupC2;
   vector<int>     *mupCat;
   vector<int>     *mupAngT;
   vector<int>     *mupNHits;
   vector<int>     *mupNPHits;
   vector<float>   *mumdxy;
   vector<float>   *mupdxy;
   vector<float>   *mumdz;
   vector<float>   *mupdz;
   vector<float>   *muon_dca;
   vector<bool>    *mu1soft;
   vector<bool>    *mu2soft;
   vector<bool>    *mu1tight;
   vector<bool>    *mu2tight;
   vector<bool>    *mu1PF;
   vector<bool>    *mu2PF;
   vector<bool>    *mu1loose;
   vector<bool>    *mu2loose;

   // List of branches
   TBranch        *b_nB;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_B_mass;   //!
   TBranch        *b_B_px;   //!
   TBranch        *b_B_py;   //!
   TBranch        *b_B_pz;   //!
   TBranch        *b_B_kaskeym_mass;   //!
   TBranch        *b_B_kaskeym_px;   //!
   TBranch        *b_B_kaskeym_py;   //!
   TBranch        *b_B_kaskeym_pz;   //!
   TBranch        *b_B_kaskeym_charge;   //!
   TBranch        *b_B_J_mass;   //!
   TBranch        *b_B_J_px;   //!
   TBranch        *b_B_J_py;   //!
   TBranch        *b_B_J_pz;   //!
   TBranch        *b_B_kas_lambda_mass;   //!
   TBranch        *b_B_kas_lambda_px;   //!
   TBranch        *b_B_kas_lambda_py;   //!
   TBranch        *b_B_kas_lambda_pz;   //!
   TBranch        *b_B_kaskeym_pt1;   //!
   TBranch        *b_B_kaskeym_px1;   //!
   TBranch        *b_B_kaskeym_py1;   //!
   TBranch        *b_B_kaskeym_pz1;   //!
   TBranch        *b_B_kaskeym_charge1;   //!
   TBranch        *b_B_kaskeym_pt2;   //!
   TBranch        *b_B_kaskeym_px2;   //!
   TBranch        *b_B_kaskeym_py2;   //!
   TBranch        *b_B_kaskeym_pz2;   //!
   TBranch        *b_B_kaskeym_charge2;   //!
   TBranch        *b_B_J_pt1;   //!
   TBranch        *b_B_J_px1;   //!
   TBranch        *b_B_J_py1;   //!
   TBranch        *b_B_J_pz1;   //!
   TBranch        *b_B_J_charge1;   //!
   TBranch        *b_B_J_pt2;   //!
   TBranch        *b_B_J_px2;   //!
   TBranch        *b_B_J_py2;   //!
   TBranch        *b_B_J_pz2;   //!
   TBranch        *b_B_J_charge2;   //!
   TBranch        *b_B_kas_lambda_pt1;   //!
   TBranch        *b_B_kas_lambda_px1;   //!
   TBranch        *b_B_kas_lambda_py1;   //!
   TBranch        *b_B_kas_lambda_pz1;   //!
   TBranch        *b_B_kas_lambda_charge1;   //!
   TBranch        *b_B_kas_lambda_pt2;   //!
   TBranch        *b_B_kas_lambda_px2;   //!
   TBranch        *b_B_kas_lambda_py2;   //!
   TBranch        *b_B_kas_lambda_pz2;   //!
   TBranch        *b_B_kas_lambda_charge2;   //!
   TBranch        *b_B_J_parentId1;   //!
   TBranch        *b_B_J_parentId2;   //!
   TBranch        *b_B_J_muId1;   //!
   TBranch        *b_B_J_muId2;   //!
   TBranch        *b_B_lam_parentId1;   //!
   TBranch        *b_B_lam_parentId2;   //!
   TBranch        *b_B_lam_PId1;   //!
   TBranch        *b_B_lam_piId2;   //!
   TBranch        *b_B_kaskeym_kId3;   //!
   TBranch        *b_B_kaskey_parentId1;   //!
   TBranch        *b_B_kaskey_parentId2;   //!
   TBranch        *b_B_kaskey_parentId3;   //!
   TBranch        *b_B_parentId1;   //!
   TBranch        *b_B_parentId2;   //!
   TBranch        *b_B_parentId3;   //!
   TBranch        *b_B_parentId4;   //!
   TBranch        *b_B_parentId5;   //!
   TBranch        *b_B_chi2;   //!
   TBranch        *b_B_kaskeym_chi2;   //!
   TBranch        *b_B_J_chi2;   //!
   TBranch        *b_B_kas_lambda_chi2;   //!
   TBranch        *b_B_Prob;   //!
   TBranch        *b_B_kaskey_Prob;   //!
   TBranch        *b_B_J_Prob;   //!
   TBranch        *b_B_kas_lambda_Prob;   //!
   TBranch        *b_priVtxX;   //!
   TBranch        *b_priVtxY;   //!
   TBranch        *b_priVtxZ;   //!
   TBranch        *b_priVtxXE;   //!
   TBranch        *b_priVtxYE;   //!
   TBranch        *b_priVtxZE;   //!
   TBranch        *b_priVtxXYE;   //!
   TBranch        *b_priVtxXZE;   //!
   TBranch        *b_priVtxYZE;   //!
   TBranch        *b_priVtxCL;   //!
   TBranch        *b_priVtxXBS;   //!
   TBranch        *b_priVtxYBS;   //!
   TBranch        *b_priVtxZBS;   //!
   TBranch        *b_priVtxXBSE;   //!
   TBranch        *b_priVtxYBSE;   //!
   TBranch        *b_priVtxZBSE;   //!
   TBranch        *b_priVtxXYBSE;   //!
   TBranch        *b_priVtxXZBSE;   //!
   TBranch        *b_priVtxYZBSE;   //!
   TBranch        *b_priVtxCLBS;   //!
   TBranch        *b_pVtxIPX;   //!
   TBranch        *b_pVtxIPY;   //!
   TBranch        *b_pVtxIPZ;   //!
   TBranch        *b_pVtxIPXE;   //!
   TBranch        *b_pVtxIPYE;   //!
   TBranch        *b_pVtxIPZE;   //!
   TBranch        *b_pVtxIPXYE;   //!
   TBranch        *b_pVtxIPXZE;   //!
   TBranch        *b_pVtxIPYZE;   //!
   TBranch        *b_pVtxIPCL;   //!
   TBranch        *b_pVtxBSIPX;   //!
   TBranch        *b_pVtxBSIPY;   //!
   TBranch        *b_pVtxBSIPZ;   //!
   TBranch        *b_pVtxBSIPXE;   //!
   TBranch        *b_pVtxBSIPYE;   //!
   TBranch        *b_pVtxBSIPZE;   //!
   TBranch        *b_pVtxBSIPXYE;   //!
   TBranch        *b_pVtxBSIPXZE;   //!
   TBranch        *b_pVtxBSIPYZE;   //!
   TBranch        *b_pVtxBSIPCL;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nTrgL;   //!
   TBranch        *b_nTrgL1L;   //!
   TBranch        *b_triggersL;   //!
   TBranch        *b_triggersL1L;   //!
   TBranch        *b_triggersMuPL;   //!
   TBranch        *b_triggersMuML;   //!
   TBranch        *b_triggersL1L2_MuPL;   //!
   TBranch        *b_triggersL1L2_MuML;   //!
   TBranch        *b_priRfVtxX;   //!
   TBranch        *b_priRfVtxY;   //!
   TBranch        *b_priRfVtxZ;   //!
   TBranch        *b_priRfVtxXE;   //!
   TBranch        *b_priRfVtxYE;   //!
   TBranch        *b_priRfVtxZE;   //!
   TBranch        *b_priRfVtxXYE;   //!
   TBranch        *b_priRfVtxXZE;   //!
   TBranch        *b_priRfVtxYZE;   //!
   TBranch        *b_priRfVtxCL;   //!
   TBranch        *b_priRfNTrkDif;   //!
   TBranch        *b_PVXBS;   //!
   TBranch        *b_PVYBS;   //!
   TBranch        *b_PVZBS;   //!
   TBranch        *b_PVXBSE;   //!
   TBranch        *b_PVYBSE;   //!
   TBranch        *b_PVZBSE;   //!
   TBranch        *b_PVXYBSE;   //!
   TBranch        *b_PVXZBSE;   //!
   TBranch        *b_PVYZBSE;   //!
   TBranch        *b_bDecayVtxX;   //!
   TBranch        *b_bDecayVtxY;   //!
   TBranch        *b_bDecayVtxZ;   //!
   TBranch        *b_bDecayVtxXE;   //!
   TBranch        *b_bDecayVtxYE;   //!
   TBranch        *b_bDecayVtxZE;   //!
   TBranch        *b_bDecayVtxXYE;   //!
   TBranch        *b_bDecayVtxXZE;   //!
   TBranch        *b_bDecayVtxYZE;   //!
   TBranch        *b_VDecayVtxX;   //!
   TBranch        *b_VDecayVtxY;   //!
   TBranch        *b_VDecayVtxZ;   //!
   TBranch        *b_VDecayVtxXE;   //!
   TBranch        *b_VDecayVtxYE;   //!
   TBranch        *b_VDecayVtxZE;   //!
   TBranch        *b_VDecayVtxXYE;   //!
   TBranch        *b_VDecayVtxXZE;   //!
   TBranch        *b_VDecayVtxYZE;   //!
   TBranch        *b_V1DecayVtxX;   //!
   TBranch        *b_V1DecayVtxY;   //!
   TBranch        *b_V1DecayVtxZ;   //!
   TBranch        *b_V1DecayVtxXE;   //!
   TBranch        *b_V1DecayVtxYE;   //!
   TBranch        *b_V1DecayVtxZE;   //!
   TBranch        *b_V1DecayVtxXYE;   //!
   TBranch        *b_V1DecayVtxXZE;   //!
   TBranch        *b_V1DecayVtxYZE;   //!
   TBranch        *b_JDecayVtxX;   //!
   TBranch        *b_JDecayVtxY;   //!
   TBranch        *b_JDecayVtxZ;   //!
   TBranch        *b_JDecayVtxXE;   //!
   TBranch        *b_JDecayVtxYE;   //!
   TBranch        *b_JDecayVtxZE;   //!
   TBranch        *b_JDecayVtxXYE;   //!
   TBranch        *b_JDecayVtxXZE;   //!
   TBranch        *b_JDecayVtxYZE;   //!
   TBranch        *b_mumC2;   //!
   TBranch        *b_mumCat;   //!
   TBranch        *b_mumAngT;   //!
   TBranch        *b_mumNHits;   //!
   TBranch        *b_mumNPHits;   //!
   TBranch        *b_mupC2;   //!
   TBranch        *b_mupCat;   //!
   TBranch        *b_mupAngT;   //!
   TBranch        *b_mupNHits;   //!
   TBranch        *b_mupNPHits;   //!
   TBranch        *b_mumdxy;   //!
   TBranch        *b_mupdxy;   //!
   TBranch        *b_mumdz;   //!
   TBranch        *b_mupdz;   //!
   TBranch        *b_muon_dca;   //!
   TBranch        *b_mu1soft;   //!
   TBranch        *b_mu2soft;   //!
   TBranch        *b_mu1tight;   //!
   TBranch        *b_mu2tight;   //!
   TBranch        *b_mu1PF;   //!
   TBranch        *b_mu2PF;   //!
   TBranch        *b_mu1loose;   //!
   TBranch        *b_mu2loose;   //!

   DataB_kaskey_2012(TTree *tree=0);
   virtual ~DataB_kaskey_2012();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DataB_kaskey_2012_cxx
DataB_kaskey_2012::DataB_kaskey_2012(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../roofilescarpetaC/output_parkedata2012_kaskey_3912.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../roofilescarpetaC/output_parkedata2012_kaskey_3912.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../roofilescarpetaC/output_parkedata2012_kaskey_3912.root:/mkcands");
      dir->GetObject("ntuple",tree);

   }
   Init(tree);
}

DataB_kaskey_2012::~DataB_kaskey_2012()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DataB_kaskey_2012::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DataB_kaskey_2012::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DataB_kaskey_2012::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   B_mass = 0;
   B_px = 0;
   B_py = 0;
   B_pz = 0;
   B_kaskeym_mass = 0;
   B_kaskeym_px = 0;
   B_kaskeym_py = 0;
   B_kaskeym_pz = 0;
   B_kaskeym_charge = 0;
   B_J_mass = 0;
   B_J_px = 0;
   B_J_py = 0;
   B_J_pz = 0;
   B_kas_lambda_mass = 0;
   B_kas_lambda_px = 0;
   B_kas_lambda_py = 0;
   B_kas_lambda_pz = 0;
   B_kaskeym_pt1 = 0;
   B_kaskeym_px1 = 0;
   B_kaskeym_py1 = 0;
   B_kaskeym_pz1 = 0;
   B_kaskeym_charge1 = 0;
   B_kaskeym_pt2 = 0;
   B_kaskeym_px2 = 0;
   B_kaskeym_py2 = 0;
   B_kaskeym_pz2 = 0;
   B_kaskeym_charge2 = 0;
   B_J_pt1 = 0;
   B_J_px1 = 0;
   B_J_py1 = 0;
   B_J_pz1 = 0;
   B_J_charge1 = 0;
   B_J_pt2 = 0;
   B_J_px2 = 0;
   B_J_py2 = 0;
   B_J_pz2 = 0;
   B_J_charge2 = 0;
   B_kas_lambda_pt1 = 0;
   B_kas_lambda_px1 = 0;
   B_kas_lambda_py1 = 0;
   B_kas_lambda_pz1 = 0;
   B_kas_lambda_charge1 = 0;
   B_kas_lambda_pt2 = 0;
   B_kas_lambda_px2 = 0;
   B_kas_lambda_py2 = 0;
   B_kas_lambda_pz2 = 0;
   B_kas_lambda_charge2 = 0;
   B_J_parentId1 = 0;
   B_J_parentId2 = 0;
   B_J_muId1 = 0;
   B_J_muId2 = 0;
   B_lam_parentId1 = 0;
   B_lam_parentId2 = 0;
   B_lam_PId1 = 0;
   B_lam_piId2 = 0;
   B_kaskeym_kId3 = 0;
   B_kaskey_parentId1 = 0;
   B_kaskey_parentId2 = 0;
   B_kaskey_parentId3 = 0;
   B_parentId1 = 0;
   B_parentId2 = 0;
   B_parentId3 = 0;
   B_parentId4 = 0;
   B_parentId5 = 0;
   B_chi2 = 0;
   B_kaskeym_chi2 = 0;
   B_J_chi2 = 0;
   B_kas_lambda_chi2 = 0;
   B_Prob = 0;
   B_kaskey_Prob = 0;
   B_J_Prob = 0;
   B_kas_lambda_Prob = 0;
   pVtxIPX = 0;
   pVtxIPY = 0;
   pVtxIPZ = 0;
   pVtxIPXE = 0;
   pVtxIPYE = 0;
   pVtxIPZE = 0;
   pVtxIPXYE = 0;
   pVtxIPXZE = 0;
   pVtxIPYZE = 0;
   pVtxIPCL = 0;
   pVtxBSIPX = 0;
   pVtxBSIPY = 0;
   pVtxBSIPZ = 0;
   pVtxBSIPXE = 0;
   pVtxBSIPYE = 0;
   pVtxBSIPZE = 0;
   pVtxBSIPXYE = 0;
   pVtxBSIPXZE = 0;
   pVtxBSIPYZE = 0;
   pVtxBSIPCL = 0;
   triggersMuPL = 0;
   triggersMuML = 0;
   triggersL1L2_MuPL = 0;
   triggersL1L2_MuML = 0;
   priRfVtxX = 0;
   priRfVtxY = 0;
   priRfVtxZ = 0;
   priRfVtxXE = 0;
   priRfVtxYE = 0;
   priRfVtxZE = 0;
   priRfVtxXYE = 0;
   priRfVtxXZE = 0;
   priRfVtxYZE = 0;
   priRfVtxCL = 0;
   priRfNTrkDif = 0;
   bDecayVtxX = 0;
   bDecayVtxY = 0;
   bDecayVtxZ = 0;
   bDecayVtxXE = 0;
   bDecayVtxYE = 0;
   bDecayVtxZE = 0;
   bDecayVtxXYE = 0;
   bDecayVtxXZE = 0;
   bDecayVtxYZE = 0;
   VDecayVtxX = 0;
   VDecayVtxY = 0;
   VDecayVtxZ = 0;
   VDecayVtxXE = 0;
   VDecayVtxYE = 0;
   VDecayVtxZE = 0;
   VDecayVtxXYE = 0;
   VDecayVtxXZE = 0;
   VDecayVtxYZE = 0;
   V1DecayVtxX = 0;
   V1DecayVtxY = 0;
   V1DecayVtxZ = 0;
   V1DecayVtxXE = 0;
   V1DecayVtxYE = 0;
   V1DecayVtxZE = 0;
   V1DecayVtxXYE = 0;
   V1DecayVtxXZE = 0;
   V1DecayVtxYZE = 0;
   JDecayVtxX = 0;
   JDecayVtxY = 0;
   JDecayVtxZ = 0;
   JDecayVtxXE = 0;
   JDecayVtxYE = 0;
   JDecayVtxZE = 0;
   JDecayVtxXYE = 0;
   JDecayVtxXZE = 0;
   JDecayVtxYZE = 0;
   mumC2 = 0;
   mumCat = 0;
   mumAngT = 0;
   mumNHits = 0;
   mumNPHits = 0;
   mupC2 = 0;
   mupCat = 0;
   mupAngT = 0;
   mupNHits = 0;
   mupNPHits = 0;
   mumdxy = 0;
   mupdxy = 0;
   mumdz = 0;
   mupdz = 0;
   muon_dca = 0;
   mu1soft = 0;
   mu2soft = 0;
   mu1tight = 0;
   mu2tight = 0;
   mu1PF = 0;
   mu2PF = 0;
   mu1loose = 0;
   mu2loose = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nB", &nB, &b_nB);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("B_mass", &B_mass, &b_B_mass);
   fChain->SetBranchAddress("B_px", &B_px, &b_B_px);
   fChain->SetBranchAddress("B_py", &B_py, &b_B_py);
   fChain->SetBranchAddress("B_pz", &B_pz, &b_B_pz);
   fChain->SetBranchAddress("B_kaskeym_mass", &B_kaskeym_mass, &b_B_kaskeym_mass);
   fChain->SetBranchAddress("B_kaskeym_px", &B_kaskeym_px, &b_B_kaskeym_px);
   fChain->SetBranchAddress("B_kaskeym_py", &B_kaskeym_py, &b_B_kaskeym_py);
   fChain->SetBranchAddress("B_kaskeym_pz", &B_kaskeym_pz, &b_B_kaskeym_pz);
   fChain->SetBranchAddress("B_kaskeym_charge", &B_kaskeym_charge, &b_B_kaskeym_charge);
   fChain->SetBranchAddress("B_J_mass", &B_J_mass, &b_B_J_mass);
   fChain->SetBranchAddress("B_J_px", &B_J_px, &b_B_J_px);
   fChain->SetBranchAddress("B_J_py", &B_J_py, &b_B_J_py);
   fChain->SetBranchAddress("B_J_pz", &B_J_pz, &b_B_J_pz);
   fChain->SetBranchAddress("B_kas_lambda_mass", &B_kas_lambda_mass, &b_B_kas_lambda_mass);
   fChain->SetBranchAddress("B_kas_lambda_px", &B_kas_lambda_px, &b_B_kas_lambda_px);
   fChain->SetBranchAddress("B_kas_lambda_py", &B_kas_lambda_py, &b_B_kas_lambda_py);
   fChain->SetBranchAddress("B_kas_lambda_pz", &B_kas_lambda_pz, &b_B_kas_lambda_pz);
   fChain->SetBranchAddress("B_kaskeym_pt1", &B_kaskeym_pt1, &b_B_kaskeym_pt1);
   fChain->SetBranchAddress("B_kaskeym_px1", &B_kaskeym_px1, &b_B_kaskeym_px1);
   fChain->SetBranchAddress("B_kaskeym_py1", &B_kaskeym_py1, &b_B_kaskeym_py1);
   fChain->SetBranchAddress("B_kaskeym_pz1", &B_kaskeym_pz1, &b_B_kaskeym_pz1);
   fChain->SetBranchAddress("B_kaskeym_charge1", &B_kaskeym_charge1, &b_B_kaskeym_charge1);
   fChain->SetBranchAddress("B_kaskeym_pt2", &B_kaskeym_pt2, &b_B_kaskeym_pt2);
   fChain->SetBranchAddress("B_kaskeym_px2", &B_kaskeym_px2, &b_B_kaskeym_px2);
   fChain->SetBranchAddress("B_kaskeym_py2", &B_kaskeym_py2, &b_B_kaskeym_py2);
   fChain->SetBranchAddress("B_kaskeym_pz2", &B_kaskeym_pz2, &b_B_kaskeym_pz2);
   fChain->SetBranchAddress("B_kaskeym_charge2", &B_kaskeym_charge2, &b_B_kaskeym_charge2);
   fChain->SetBranchAddress("B_J_pt1", &B_J_pt1, &b_B_J_pt1);
   fChain->SetBranchAddress("B_J_px1", &B_J_px1, &b_B_J_px1);
   fChain->SetBranchAddress("B_J_py1", &B_J_py1, &b_B_J_py1);
   fChain->SetBranchAddress("B_J_pz1", &B_J_pz1, &b_B_J_pz1);
   fChain->SetBranchAddress("B_J_charge1", &B_J_charge1, &b_B_J_charge1);
   fChain->SetBranchAddress("B_J_pt2", &B_J_pt2, &b_B_J_pt2);
   fChain->SetBranchAddress("B_J_px2", &B_J_px2, &b_B_J_px2);
   fChain->SetBranchAddress("B_J_py2", &B_J_py2, &b_B_J_py2);
   fChain->SetBranchAddress("B_J_pz2", &B_J_pz2, &b_B_J_pz2);
   fChain->SetBranchAddress("B_J_charge2", &B_J_charge2, &b_B_J_charge2);
   fChain->SetBranchAddress("B_kas_lambda_pt1", &B_kas_lambda_pt1, &b_B_kas_lambda_pt1);
   fChain->SetBranchAddress("B_kas_lambda_px1", &B_kas_lambda_px1, &b_B_kas_lambda_px1);
   fChain->SetBranchAddress("B_kas_lambda_py1", &B_kas_lambda_py1, &b_B_kas_lambda_py1);
   fChain->SetBranchAddress("B_kas_lambda_pz1", &B_kas_lambda_pz1, &b_B_kas_lambda_pz1);
   fChain->SetBranchAddress("B_kas_lambda_charge1", &B_kas_lambda_charge1, &b_B_kas_lambda_charge1);
   fChain->SetBranchAddress("B_kas_lambda_pt2", &B_kas_lambda_pt2, &b_B_kas_lambda_pt2);
   fChain->SetBranchAddress("B_kas_lambda_px2", &B_kas_lambda_px2, &b_B_kas_lambda_px2);
   fChain->SetBranchAddress("B_kas_lambda_py2", &B_kas_lambda_py2, &b_B_kas_lambda_py2);
   fChain->SetBranchAddress("B_kas_lambda_pz2", &B_kas_lambda_pz2, &b_B_kas_lambda_pz2);
   fChain->SetBranchAddress("B_kas_lambda_charge2", &B_kas_lambda_charge2, &b_B_kas_lambda_charge2);
   fChain->SetBranchAddress("B_J_parentId1", &B_J_parentId1, &b_B_J_parentId1);
   fChain->SetBranchAddress("B_J_parentId2", &B_J_parentId2, &b_B_J_parentId2);
   fChain->SetBranchAddress("B_J_muId1", &B_J_muId1, &b_B_J_muId1);
   fChain->SetBranchAddress("B_J_muId2", &B_J_muId2, &b_B_J_muId2);
   fChain->SetBranchAddress("B_lam_parentId1", &B_lam_parentId1, &b_B_lam_parentId1);
   fChain->SetBranchAddress("B_lam_parentId2", &B_lam_parentId2, &b_B_lam_parentId2);
   fChain->SetBranchAddress("B_lam_PId1", &B_lam_PId1, &b_B_lam_PId1);
   fChain->SetBranchAddress("B_lam_piId2", &B_lam_piId2, &b_B_lam_piId2);
   fChain->SetBranchAddress("B_kaskeym_kId3", &B_kaskeym_kId3, &b_B_kaskeym_kId3);
   fChain->SetBranchAddress("B_kaskey_parentId1", &B_kaskey_parentId1, &b_B_kaskey_parentId1);
   fChain->SetBranchAddress("B_kaskey_parentId2", &B_kaskey_parentId2, &b_B_kaskey_parentId2);
   fChain->SetBranchAddress("B_kaskey_parentId3", &B_kaskey_parentId3, &b_B_kaskey_parentId3);
   fChain->SetBranchAddress("B_parentId1", &B_parentId1, &b_B_parentId1);
   fChain->SetBranchAddress("B_parentId2", &B_parentId2, &b_B_parentId2);
   fChain->SetBranchAddress("B_parentId3", &B_parentId3, &b_B_parentId3);
   fChain->SetBranchAddress("B_parentId4", &B_parentId4, &b_B_parentId4);
   fChain->SetBranchAddress("B_parentId5", &B_parentId5, &b_B_parentId5);
   fChain->SetBranchAddress("B_chi2", &B_chi2, &b_B_chi2);
   fChain->SetBranchAddress("B_kaskeym_chi2", &B_kaskeym_chi2, &b_B_kaskeym_chi2);
   fChain->SetBranchAddress("B_J_chi2", &B_J_chi2, &b_B_J_chi2);
   fChain->SetBranchAddress("B_kas_lambda_chi2", &B_kas_lambda_chi2, &b_B_kas_lambda_chi2);
   fChain->SetBranchAddress("B_Prob", &B_Prob, &b_B_Prob);
   fChain->SetBranchAddress("B_kaskey_Prob", &B_kaskey_Prob, &b_B_kaskey_Prob);
   fChain->SetBranchAddress("B_J_Prob", &B_J_Prob, &b_B_J_Prob);
   fChain->SetBranchAddress("B_kas_lambda_Prob", &B_kas_lambda_Prob, &b_B_kas_lambda_Prob);
   fChain->SetBranchAddress("priVtxX", &priVtxX, &b_priVtxX);
   fChain->SetBranchAddress("priVtxY", &priVtxY, &b_priVtxY);
   fChain->SetBranchAddress("priVtxZ", &priVtxZ, &b_priVtxZ);
   fChain->SetBranchAddress("priVtxXE", &priVtxXE, &b_priVtxXE);
   fChain->SetBranchAddress("priVtxYE", &priVtxYE, &b_priVtxYE);
   fChain->SetBranchAddress("priVtxZE", &priVtxZE, &b_priVtxZE);
   fChain->SetBranchAddress("priVtxXYE", &priVtxXYE, &b_priVtxXYE);
   fChain->SetBranchAddress("priVtxXZE", &priVtxXZE, &b_priVtxXZE);
   fChain->SetBranchAddress("priVtxYZE", &priVtxYZE, &b_priVtxYZE);
   fChain->SetBranchAddress("priVtxCL", &priVtxCL, &b_priVtxCL);
   fChain->SetBranchAddress("priVtxXBS", &priVtxXBS, &b_priVtxXBS);
   fChain->SetBranchAddress("priVtxYBS", &priVtxYBS, &b_priVtxYBS);
   fChain->SetBranchAddress("priVtxZBS", &priVtxZBS, &b_priVtxZBS);
   fChain->SetBranchAddress("priVtxXBSE", &priVtxXBSE, &b_priVtxXBSE);
   fChain->SetBranchAddress("priVtxYBSE", &priVtxYBSE, &b_priVtxYBSE);
   fChain->SetBranchAddress("priVtxZBSE", &priVtxZBSE, &b_priVtxZBSE);
   fChain->SetBranchAddress("priVtxXYBSE", &priVtxXYBSE, &b_priVtxXYBSE);
   fChain->SetBranchAddress("priVtxXZBSE", &priVtxXZBSE, &b_priVtxXZBSE);
   fChain->SetBranchAddress("priVtxYZBSE", &priVtxYZBSE, &b_priVtxYZBSE);
   fChain->SetBranchAddress("priVtxCLBS", &priVtxCLBS, &b_priVtxCLBS);
   fChain->SetBranchAddress("pVtxIPX", &pVtxIPX, &b_pVtxIPX);
   fChain->SetBranchAddress("pVtxIPY", &pVtxIPY, &b_pVtxIPY);
   fChain->SetBranchAddress("pVtxIPZ", &pVtxIPZ, &b_pVtxIPZ);
   fChain->SetBranchAddress("pVtxIPXE", &pVtxIPXE, &b_pVtxIPXE);
   fChain->SetBranchAddress("pVtxIPYE", &pVtxIPYE, &b_pVtxIPYE);
   fChain->SetBranchAddress("pVtxIPZE", &pVtxIPZE, &b_pVtxIPZE);
   fChain->SetBranchAddress("pVtxIPXYE", &pVtxIPXYE, &b_pVtxIPXYE);
   fChain->SetBranchAddress("pVtxIPXZE", &pVtxIPXZE, &b_pVtxIPXZE);
   fChain->SetBranchAddress("pVtxIPYZE", &pVtxIPYZE, &b_pVtxIPYZE);
   fChain->SetBranchAddress("pVtxIPCL", &pVtxIPCL, &b_pVtxIPCL);
   fChain->SetBranchAddress("pVtxBSIPX", &pVtxBSIPX, &b_pVtxBSIPX);
   fChain->SetBranchAddress("pVtxBSIPY", &pVtxBSIPY, &b_pVtxBSIPY);
   fChain->SetBranchAddress("pVtxBSIPZ", &pVtxBSIPZ, &b_pVtxBSIPZ);
   fChain->SetBranchAddress("pVtxBSIPXE", &pVtxBSIPXE, &b_pVtxBSIPXE);
   fChain->SetBranchAddress("pVtxBSIPYE", &pVtxBSIPYE, &b_pVtxBSIPYE);
   fChain->SetBranchAddress("pVtxBSIPZE", &pVtxBSIPZE, &b_pVtxBSIPZE);
   fChain->SetBranchAddress("pVtxBSIPXYE", &pVtxBSIPXYE, &b_pVtxBSIPXYE);
   fChain->SetBranchAddress("pVtxBSIPXZE", &pVtxBSIPXZE, &b_pVtxBSIPXZE);
   fChain->SetBranchAddress("pVtxBSIPYZE", &pVtxBSIPYZE, &b_pVtxBSIPYZE);
   fChain->SetBranchAddress("pVtxBSIPCL", &pVtxBSIPCL, &b_pVtxBSIPCL);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nTrgL", &nTrgL, &b_nTrgL);
   fChain->SetBranchAddress("nTrgL1L", &nTrgL1L, &b_nTrgL1L);
   fChain->SetBranchAddress("triggersL", triggersL, &b_triggersL);
   fChain->SetBranchAddress("triggersL1L", triggersL1L, &b_triggersL1L);
   fChain->SetBranchAddress("triggersMuPL", &triggersMuPL, &b_triggersMuPL);
   fChain->SetBranchAddress("triggersMuML", &triggersMuML, &b_triggersMuML);
   fChain->SetBranchAddress("triggersL1L2_MuPL", &triggersL1L2_MuPL, &b_triggersL1L2_MuPL);
   fChain->SetBranchAddress("triggersL1L2_MuML", &triggersL1L2_MuML, &b_triggersL1L2_MuML);
   fChain->SetBranchAddress("priRfVtxX", &priRfVtxX, &b_priRfVtxX);
   fChain->SetBranchAddress("priRfVtxY", &priRfVtxY, &b_priRfVtxY);
   fChain->SetBranchAddress("priRfVtxZ", &priRfVtxZ, &b_priRfVtxZ);
   fChain->SetBranchAddress("priRfVtxXE", &priRfVtxXE, &b_priRfVtxXE);
   fChain->SetBranchAddress("priRfVtxYE", &priRfVtxYE, &b_priRfVtxYE);
   fChain->SetBranchAddress("priRfVtxZE", &priRfVtxZE, &b_priRfVtxZE);
   fChain->SetBranchAddress("priRfVtxXYE", &priRfVtxXYE, &b_priRfVtxXYE);
   fChain->SetBranchAddress("priRfVtxXZE", &priRfVtxXZE, &b_priRfVtxXZE);
   fChain->SetBranchAddress("priRfVtxYZE", &priRfVtxYZE, &b_priRfVtxYZE);
   fChain->SetBranchAddress("priRfVtxCL", &priRfVtxCL, &b_priRfVtxCL);
   fChain->SetBranchAddress("priRfNTrkDif", &priRfNTrkDif, &b_priRfNTrkDif);
   fChain->SetBranchAddress("PVXBS", &PVXBS, &b_PVXBS);
   fChain->SetBranchAddress("PVYBS", &PVYBS, &b_PVYBS);
   fChain->SetBranchAddress("PVZBS", &PVZBS, &b_PVZBS);
   fChain->SetBranchAddress("PVXBSE", &PVXBSE, &b_PVXBSE);
   fChain->SetBranchAddress("PVYBSE", &PVYBSE, &b_PVYBSE);
   fChain->SetBranchAddress("PVZBSE", &PVZBSE, &b_PVZBSE);
   fChain->SetBranchAddress("PVXYBSE", &PVXYBSE, &b_PVXYBSE);
   fChain->SetBranchAddress("PVXZBSE", &PVXZBSE, &b_PVXZBSE);
   fChain->SetBranchAddress("PVYZBSE", &PVYZBSE, &b_PVYZBSE);
   fChain->SetBranchAddress("bDecayVtxX", &bDecayVtxX, &b_bDecayVtxX);
   fChain->SetBranchAddress("bDecayVtxY", &bDecayVtxY, &b_bDecayVtxY);
   fChain->SetBranchAddress("bDecayVtxZ", &bDecayVtxZ, &b_bDecayVtxZ);
   fChain->SetBranchAddress("bDecayVtxXE", &bDecayVtxXE, &b_bDecayVtxXE);
   fChain->SetBranchAddress("bDecayVtxYE", &bDecayVtxYE, &b_bDecayVtxYE);
   fChain->SetBranchAddress("bDecayVtxZE", &bDecayVtxZE, &b_bDecayVtxZE);
   fChain->SetBranchAddress("bDecayVtxXYE", &bDecayVtxXYE, &b_bDecayVtxXYE);
   fChain->SetBranchAddress("bDecayVtxXZE", &bDecayVtxXZE, &b_bDecayVtxXZE);
   fChain->SetBranchAddress("bDecayVtxYZE", &bDecayVtxYZE, &b_bDecayVtxYZE);
   fChain->SetBranchAddress("VDecayVtxX", &VDecayVtxX, &b_VDecayVtxX);
   fChain->SetBranchAddress("VDecayVtxY", &VDecayVtxY, &b_VDecayVtxY);
   fChain->SetBranchAddress("VDecayVtxZ", &VDecayVtxZ, &b_VDecayVtxZ);
   fChain->SetBranchAddress("VDecayVtxXE", &VDecayVtxXE, &b_VDecayVtxXE);
   fChain->SetBranchAddress("VDecayVtxYE", &VDecayVtxYE, &b_VDecayVtxYE);
   fChain->SetBranchAddress("VDecayVtxZE", &VDecayVtxZE, &b_VDecayVtxZE);
   fChain->SetBranchAddress("VDecayVtxXYE", &VDecayVtxXYE, &b_VDecayVtxXYE);
   fChain->SetBranchAddress("VDecayVtxXZE", &VDecayVtxXZE, &b_VDecayVtxXZE);
   fChain->SetBranchAddress("VDecayVtxYZE", &VDecayVtxYZE, &b_VDecayVtxYZE);
   fChain->SetBranchAddress("V1DecayVtxX", &V1DecayVtxX, &b_V1DecayVtxX);
   fChain->SetBranchAddress("V1DecayVtxY", &V1DecayVtxY, &b_V1DecayVtxY);
   fChain->SetBranchAddress("V1DecayVtxZ", &V1DecayVtxZ, &b_V1DecayVtxZ);
   fChain->SetBranchAddress("V1DecayVtxXE", &V1DecayVtxXE, &b_V1DecayVtxXE);
   fChain->SetBranchAddress("V1DecayVtxYE", &V1DecayVtxYE, &b_V1DecayVtxYE);
   fChain->SetBranchAddress("V1DecayVtxZE", &V1DecayVtxZE, &b_V1DecayVtxZE);
   fChain->SetBranchAddress("V1DecayVtxXYE", &V1DecayVtxXYE, &b_V1DecayVtxXYE);
   fChain->SetBranchAddress("V1DecayVtxXZE", &V1DecayVtxXZE, &b_V1DecayVtxXZE);
   fChain->SetBranchAddress("V1DecayVtxYZE", &V1DecayVtxYZE, &b_V1DecayVtxYZE);
   fChain->SetBranchAddress("JDecayVtxX", &JDecayVtxX, &b_JDecayVtxX);
   fChain->SetBranchAddress("JDecayVtxY", &JDecayVtxY, &b_JDecayVtxY);
   fChain->SetBranchAddress("JDecayVtxZ", &JDecayVtxZ, &b_JDecayVtxZ);
   fChain->SetBranchAddress("JDecayVtxXE", &JDecayVtxXE, &b_JDecayVtxXE);
   fChain->SetBranchAddress("JDecayVtxYE", &JDecayVtxYE, &b_JDecayVtxYE);
   fChain->SetBranchAddress("JDecayVtxZE", &JDecayVtxZE, &b_JDecayVtxZE);
   fChain->SetBranchAddress("JDecayVtxXYE", &JDecayVtxXYE, &b_JDecayVtxXYE);
   fChain->SetBranchAddress("JDecayVtxXZE", &JDecayVtxXZE, &b_JDecayVtxXZE);
   fChain->SetBranchAddress("JDecayVtxYZE", &JDecayVtxYZE, &b_JDecayVtxYZE);
   fChain->SetBranchAddress("mumC2", &mumC2, &b_mumC2);
   fChain->SetBranchAddress("mumCat", &mumCat, &b_mumCat);
   fChain->SetBranchAddress("mumAngT", &mumAngT, &b_mumAngT);
   fChain->SetBranchAddress("mumNHits", &mumNHits, &b_mumNHits);
   fChain->SetBranchAddress("mumNPHits", &mumNPHits, &b_mumNPHits);
   fChain->SetBranchAddress("mupC2", &mupC2, &b_mupC2);
   fChain->SetBranchAddress("mupCat", &mupCat, &b_mupCat);
   fChain->SetBranchAddress("mupAngT", &mupAngT, &b_mupAngT);
   fChain->SetBranchAddress("mupNHits", &mupNHits, &b_mupNHits);
   fChain->SetBranchAddress("mupNPHits", &mupNPHits, &b_mupNPHits);
   fChain->SetBranchAddress("mumdxy", &mumdxy, &b_mumdxy);
   fChain->SetBranchAddress("mupdxy", &mupdxy, &b_mupdxy);
   fChain->SetBranchAddress("mumdz", &mumdz, &b_mumdz);
   fChain->SetBranchAddress("mupdz", &mupdz, &b_mupdz);
   fChain->SetBranchAddress("muon_dca", &muon_dca, &b_muon_dca);
   fChain->SetBranchAddress("mu1soft", &mu1soft, &b_mu1soft);
   fChain->SetBranchAddress("mu2soft", &mu2soft, &b_mu2soft);
   fChain->SetBranchAddress("mu1tight", &mu1tight, &b_mu1tight);
   fChain->SetBranchAddress("mu2tight", &mu2tight, &b_mu2tight);
   fChain->SetBranchAddress("mu1PF", &mu1PF, &b_mu1PF);
   fChain->SetBranchAddress("mu2PF", &mu2PF, &b_mu2PF);
   fChain->SetBranchAddress("mu1loose", &mu1loose, &b_mu1loose);
   fChain->SetBranchAddress("mu2loose", &mu2loose, &b_mu2loose);
   Notify();
}

Bool_t DataB_kaskey_2012::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DataB_kaskey_2012::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DataB_kaskey_2012::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DataB_kaskey_2012_cxx
