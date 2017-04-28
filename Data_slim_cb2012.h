//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 27 17:30:40 2016 by ROOT version 5.34/36
// from TTree treeS/signal
// found on file: ROOTSB_kaskey_2012parkedata_C.root
//////////////////////////////////////////////////////////

#ifndef Data_slim_cb2012_h
#define Data_slim_cb2012_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Data_slim_cb2012 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        massB;
   Double_t        massJ;
   Double_t        masslamb;
   Double_t        masskm;
   Double_t        Bdl;
   Double_t        BdlBS;
   Double_t        BdlRf;
   Double_t        BdlIP;
   Double_t        BSct;
   Double_t        BdlIPBSc;
   Double_t        BdlE;
   Double_t        BdlBSE;
   Double_t        BdlRfE;
   Double_t        BdlIPE;
   Double_t        BSctE;
   Double_t        BdlIPBScE;
   Double_t        Bdlkminus;
   Double_t        BdlkminusE;
   Double_t        Bdllamb;
   Double_t        BdllambE;
   Double_t        Beta;
   Double_t        Bphi;
   Double_t        sigLxyJ;
   Double_t        cosalfaJ;
   Double_t        Bpt;
   Double_t        lambpt;
   Double_t        Jpsipt;
   Double_t        kmpt;
   Double_t        mu1pt;
   Double_t        mu2pt;
   Double_t        Pi1pt;
   Double_t        Pi2pt;
   Double_t        Pi3pt;
   Double_t        mu1phi;
   Double_t        mu2phi;
   Double_t        Pi1phi;
   Double_t        Pi2phi;
   Double_t        Jeta;
   Double_t        mu1eta;
   Double_t        mu2eta;
   Double_t        Pi1eta;
   Double_t        Pi2eta;
   Double_t        maxPtMu;
   Double_t        maxPtPi;
   Double_t        minPtMu;
   Double_t        minPtPi;
   Double_t        lambchi2;
   Double_t        Jchi2;
   Double_t        Bchi2;
   Double_t        kmchi2;
   Double_t        lambprob;
   Double_t        Jprob;
   Double_t        kmprob;
   Double_t        Bprob;
   Int_t           B_kaskeym_charge1;
   Int_t           B_kas_lambda_charge2;
   Double_t        massOmem;
   Double_t        massKs0;
   Int_t           mupDisplace;
   Int_t           mumDisplace;
   Int_t           runn;
   Int_t           eventn;
   Int_t           mupcat;
   Int_t           mumcat;
   Bool_t          softmu1;
   Bool_t          softmu2;
   Bool_t          tightmu1;
   Bool_t          tightmu2;
   Bool_t          PFmu1;
   Bool_t          PFmu2;
   Int_t           triggerv4;
   Int_t           triggerv5;
   Int_t           triggerv6;
   Int_t           triggerv7;

   // List of branches
   TBranch        *b_massB;   //!
   TBranch        *b_massJ;   //!
   TBranch        *b_masslamb;   //!
   TBranch        *b_masskm;   //!
   TBranch        *b_Bdl;   //!
   TBranch        *b_BdlBS;   //!
   TBranch        *b_BdlRf;   //!
   TBranch        *b_BdlIP;   //!
   TBranch        *b_BSct;   //!
   TBranch        *b_BdlIPBSc;   //!
   TBranch        *b_BdlE;   //!
   TBranch        *b_BdlBSE;   //!
   TBranch        *b_BdlRfE;   //!
   TBranch        *b_BdlIPE;   //!
   TBranch        *b_BSctE;   //!
   TBranch        *b_BdlIPBScE;   //!
   TBranch        *b_Bdlkminus;   //!
   TBranch        *b_BdlkminusE;   //!
   TBranch        *b_Bdllamb;   //!
   TBranch        *b_BdllambE;   //!
   TBranch        *b_Beta;   //!
   TBranch        *b_Bphi;   //!
   TBranch        *b_sigLxyJ;   //!
   TBranch        *b_cosalfaJ;   //!
   TBranch        *b_Bpt;   //!
   TBranch        *b_lambpt;   //!
   TBranch        *b_Jpsipt;   //!
   TBranch        *b_kmpt;   //!
   TBranch        *b_mu1pt;   //!
   TBranch        *b_mu2pt;   //!
   TBranch        *b_Pi1pt;   //!
   TBranch        *b_Pi2pt;   //!
   TBranch        *b_Pi3pt;   //!
   TBranch        *b_mu1phi;   //!
   TBranch        *b_mu2phi;   //!
   TBranch        *b_Pi1phi;   //!
   TBranch        *b_Pi2phi;   //!
   TBranch        *b_Jeta;   //!
   TBranch        *b_mu1eta;   //!
   TBranch        *b_mu2eta;   //!
   TBranch        *b_Pi1eta;   //!
   TBranch        *b_Pi2eta;   //!
   TBranch        *b_maxPtMu;   //!
   TBranch        *b_maxPtPi;   //!
   TBranch        *b_minPtMu;   //!
   TBranch        *b_minPtPi;   //!
   TBranch        *b_lambchi2;   //!
   TBranch        *b_Jchi2;   //!
   TBranch        *b_Bchi2;   //!
   TBranch        *b_kmchi2;   //!
   TBranch        *b_lambprob;   //!
   TBranch        *b_Jprob;   //!
   TBranch        *b_kmprob;   //!
   TBranch        *b_Bprob;   //!
   TBranch        *b_B_kaskeym_charge1;   //!
   TBranch        *b_B_kas_lambda_charge2;   //!
   TBranch        *b_massOmem;   //!
   TBranch        *b_massKs0;   //!
   TBranch        *b_mupDisplace;   //!
   TBranch        *b_mumDisplace;   //!
   TBranch        *b_runn;   //!
   TBranch        *b_eventn;   //!
   TBranch        *b_mupcat;   //!
   TBranch        *b_mumcat;   //!
   TBranch        *b_softmu1;   //!
   TBranch        *b_softmu2;   //!
   TBranch        *b_tightmu1;   //!
   TBranch        *b_tightmu2;   //!
   TBranch        *b_PFmu1;   //!
   TBranch        *b_PFmu2;   //!
   TBranch        *b_triggerv4;   //!
   TBranch        *b_triggerv5;   //!
   TBranch        *b_triggerv6;   //!
   TBranch        *b_triggerv7;   //!

   Data_slim_cb2012(TTree *tree=0);
   virtual ~Data_slim_cb2012();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Data_slim_cb2012_cxx
Data_slim_cb2012::Data_slim_cb2012(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ROOTSB_kaskey_2012parkedata_C.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ROOTSB_kaskey_2012parkedata_C.root");
      }
      f->GetObject("treeS",tree);

   }
   Init(tree);
}

Data_slim_cb2012::~Data_slim_cb2012()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Data_slim_cb2012::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Data_slim_cb2012::LoadTree(Long64_t entry)
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

void Data_slim_cb2012::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("massB", &massB, &b_massB);
   fChain->SetBranchAddress("massJ", &massJ, &b_massJ);
   fChain->SetBranchAddress("masslamb", &masslamb, &b_masslamb);
   fChain->SetBranchAddress("masskm", &masskm, &b_masskm);
   fChain->SetBranchAddress("Bdl", &Bdl, &b_Bdl);
   fChain->SetBranchAddress("BdlBS", &BdlBS, &b_BdlBS);
   fChain->SetBranchAddress("BdlRf", &BdlRf, &b_BdlRf);
   fChain->SetBranchAddress("BdlIP", &BdlIP, &b_BdlIP);
   fChain->SetBranchAddress("BSct", &BSct, &b_BSct);
   fChain->SetBranchAddress("BdlIPBSc", &BdlIPBSc, &b_BdlIPBSc);
   fChain->SetBranchAddress("BdlE", &BdlE, &b_BdlE);
   fChain->SetBranchAddress("BdlBSE", &BdlBSE, &b_BdlBSE);
   fChain->SetBranchAddress("BdlRfE", &BdlRfE, &b_BdlRfE);
   fChain->SetBranchAddress("BdlIPE", &BdlIPE, &b_BdlIPE);
   fChain->SetBranchAddress("BSctE", &BSctE, &b_BSctE);
   fChain->SetBranchAddress("BdlIPBScE", &BdlIPBScE, &b_BdlIPBScE);
   fChain->SetBranchAddress("Bdlkminus", &Bdlkminus, &b_Bdlkminus);
   fChain->SetBranchAddress("BdlkminusE", &BdlkminusE, &b_BdlkminusE);
   fChain->SetBranchAddress("Bdllamb", &Bdllamb, &b_Bdllamb);
   fChain->SetBranchAddress("BdllambE", &BdllambE, &b_BdllambE);
   fChain->SetBranchAddress("Beta", &Beta, &b_Beta);
   fChain->SetBranchAddress("Bphi", &Bphi, &b_Bphi);
   fChain->SetBranchAddress("sigLxyJ", &sigLxyJ, &b_sigLxyJ);
   fChain->SetBranchAddress("cosalfaJ", &cosalfaJ, &b_cosalfaJ);
   fChain->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
   fChain->SetBranchAddress("lambpt", &lambpt, &b_lambpt);
   fChain->SetBranchAddress("Jpsipt", &Jpsipt, &b_Jpsipt);
   fChain->SetBranchAddress("kmpt", &kmpt, &b_kmpt);
   fChain->SetBranchAddress("mu1pt", &mu1pt, &b_mu1pt);
   fChain->SetBranchAddress("mu2pt", &mu2pt, &b_mu2pt);
   fChain->SetBranchAddress("Pi1pt", &Pi1pt, &b_Pi1pt);
   fChain->SetBranchAddress("Pi2pt", &Pi2pt, &b_Pi2pt);
   fChain->SetBranchAddress("Pi3pt", &Pi3pt, &b_Pi3pt);
   fChain->SetBranchAddress("mu1phi", &mu1phi, &b_mu1phi);
   fChain->SetBranchAddress("mu2phi", &mu2phi, &b_mu2phi);
   fChain->SetBranchAddress("Pi1phi", &Pi1phi, &b_Pi1phi);
   fChain->SetBranchAddress("Pi2phi", &Pi2phi, &b_Pi2phi);
   fChain->SetBranchAddress("Jeta", &Jeta, &b_Jeta);
   fChain->SetBranchAddress("mu1eta", &mu1eta, &b_mu1eta);
   fChain->SetBranchAddress("mu2eta", &mu2eta, &b_mu2eta);
   fChain->SetBranchAddress("Pi1eta", &Pi1eta, &b_Pi1eta);
   fChain->SetBranchAddress("Pi2eta", &Pi2eta, &b_Pi2eta);
   fChain->SetBranchAddress("maxPtMu", &maxPtMu, &b_maxPtMu);
   fChain->SetBranchAddress("maxPtPi", &maxPtPi, &b_maxPtPi);
   fChain->SetBranchAddress("minPtMu", &minPtMu, &b_minPtMu);
   fChain->SetBranchAddress("minPtPi", &minPtPi, &b_minPtPi);
   fChain->SetBranchAddress("lambchi2", &lambchi2, &b_lambchi2);
   fChain->SetBranchAddress("Jchi2", &Jchi2, &b_Jchi2);
   fChain->SetBranchAddress("Bchi2", &Bchi2, &b_Bchi2);
   fChain->SetBranchAddress("kmchi2", &kmchi2, &b_kmchi2);
   fChain->SetBranchAddress("lambprob", &lambprob, &b_lambprob);
   fChain->SetBranchAddress("Jprob", &Jprob, &b_Jprob);
   fChain->SetBranchAddress("kmprob", &kmprob, &b_kmprob);
   fChain->SetBranchAddress("Bprob", &Bprob, &b_Bprob);
   fChain->SetBranchAddress("B_kaskeym_charge1", &B_kaskeym_charge1, &b_B_kaskeym_charge1);
   fChain->SetBranchAddress("B_kas_lambda_charge2", &B_kas_lambda_charge2, &b_B_kas_lambda_charge2);
   fChain->SetBranchAddress("massOmem", &massOmem, &b_massOmem);
   fChain->SetBranchAddress("massKs0", &massKs0, &b_massKs0);
   fChain->SetBranchAddress("mupDisplace", &mupDisplace, &b_mupDisplace);
   fChain->SetBranchAddress("mumDisplace", &mumDisplace, &b_mumDisplace);
   fChain->SetBranchAddress("runn", &runn, &b_runn);
   fChain->SetBranchAddress("eventn", &eventn, &b_eventn);
   fChain->SetBranchAddress("mupcat", &mupcat, &b_mupcat);
   fChain->SetBranchAddress("mumcat", &mumcat, &b_mumcat);
   fChain->SetBranchAddress("softmu1", &softmu1, &b_softmu1);
   fChain->SetBranchAddress("softmu2", &softmu2, &b_softmu2);
   fChain->SetBranchAddress("tightmu1", &tightmu1, &b_tightmu1);
   fChain->SetBranchAddress("tightmu2", &tightmu2, &b_tightmu2);
   fChain->SetBranchAddress("PFmu1", &PFmu1, &b_PFmu1);
   fChain->SetBranchAddress("PFmu2", &PFmu2, &b_PFmu2);
   fChain->SetBranchAddress("triggerv4", &triggerv4, &b_triggerv4);
   fChain->SetBranchAddress("triggerv5", &triggerv5, &b_triggerv5);
   fChain->SetBranchAddress("triggerv6", &triggerv6, &b_triggerv6);
   fChain->SetBranchAddress("triggerv7", &triggerv7, &b_triggerv7);
   Notify();
}

Bool_t Data_slim_cb2012::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Data_slim_cb2012::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Data_slim_cb2012::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Data_slim_cb2012_cxx
