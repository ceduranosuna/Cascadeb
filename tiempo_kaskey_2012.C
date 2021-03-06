{
  //------ Plotting
  
  gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white
  gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(1.0);
  //gStyle->SetMarkerSize(0.5);
  gStyle->SetMarkerSize(1.0);
  
  bool showbox =false;
    
  Double_t Mmin = 5.6;//5.45
  Double_t Mmax = 6.0;//5.8
    
  Double_t ctmax = 0.01;//0.01
  Double_t ctmin = 0.00001;//0.0001  
  
 
   //datos a procesar
  //2012
  TChain *ch2012 = new TChain("treeS",""); 
  
  ch2012->Add("ROOTSB_kaskey_2012parkedata_B.root/treeS");
  ch2012->Add("ROOTSB_kaskey_2012parkedata_C.root/treeS");
  ch2012->Add("ROOTSB_kaskey_2012parkedata_D.root/treeS");

    //tomamos las clases para abrir las tples
  gROOT->ProcessLine(".L Data_slim_cb2012.C+");
  
  TTree *tree2012 = (TTree*) ch2012; 
  Data_slim_cb2012 t2012(tree2012);
  
  ////DAta2012
  Long64_t nentries2012 = t2012.fChain->GetEntries();
  cout<<" Entries : "<<nentries2012<<endl;
  
  using namespace RooFit;
  using namespace std;
  
  
  //Definicion de variables
  RooRealVar M("M","Invariant Mass #Xi_{b}^{-}(J/#psi #Xi^{-}) [GeV/c^{2}]",Mmin,Mmax);
  RooRealVar Mlam("Mlam","PDL [cm]",0.02,0.5);//
  RooRealVar MlamE("MlamE","Error de tiempo de vida",ctmin,ctmax);// esta para el error de la funcion de tiempo de vida
  RooDataSet data2012("data2012","data2012",RooArgSet(M,Mlam,MlamE));

  TH1F *hmass2012 = new TH1F("Mass (2012)","hmass2012",40,Mmin,Mmax);

    
 
  /*
    en el for que sigue se hacen los ultimos cortes y ademas se llenan las variables a grafiacar
    con sus respectivos dactos. Donde Btau es la variable de la longitud propia de decaimiento
    y massB es la masa del meson B.
  */
  

  Int_t nTen = nentries2012/10;
  TLorentzVector Jpsi,Xi;
  Int_t k=0;
  Int_t nbytes = 0, nb = 0;
  for(Long64_t jentry=0; jentry<nentries2012;jentry++)
    {
      Long64_t ientry = t2012.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t2012.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
      if(jentry==nentries2012-1) cout<<endl;
      
      //Mass windows cuts
      if(t2012.massB<=Mmin || t2012.massB>=Mmax) continue;
      if(t2012.masskm<=1.311 || t2012.masskm>=1.332) continue;
      if(t2012.masslamb<=1.1096 || t2012.masslamb>=1.1216) continue;
      if(t2012.Pi1pt<1.5) continue;//1.2
      if(t2012.Pi2pt<0.5) continue;
      if(t2012.Pi3pt<0.5) continue;
      if(t2012.Bdlkminus<0.0) continue; // -->IVAN: Estaba en 0.02!
      
      
      if(t2012.Bpt<13.0) continue;
      if(t2012.mu1pt<4.0) continue;
      if(t2012.mu2pt<4.0) continue;
      
      if(t2012.Jprob<0.005)continue;   
      //if(t2012.Bprob<0.1)continue;

      //if(t2012.lambchi2>7.0)continue;//1   //--> Lo quitamos por instrucciones de Eduard: default 26
      //if(t2012.kmchi2>26)continue;
      //if(t2012.Jchi2>7.0)continue;  //  --> Lo quitamos por instrucciones de Eduard: default 26


      //conservacion del numero varionico
      if(t2012.B_kaskeym_charge1!=t2012.B_kas_lambda_charge2)continue;


      if(TMath::IsNaN(t2012.massB)==1) continue;


      /*
      if(t2012.Bdl<0.02) continue;
      if(t2012.Bdl>0.5) continue;  //-->IVAN: Porque estaba este corte
      if(TMath::IsNaN(t2012.Bdl)==1) continue;
      if(TMath::IsNaN(t2012.BdlE)==1) continue;
      if(t2012.BdlE<=ctmin || t2012.BdlE>=ctmax) continue;
      */

	/*
	if(t2012.BdlIP<0.02) continue; // resultado  4.4868e-02 +/-  1.30e-03
	if(t2012.BdlIP>0.5) continue;
	if(TMath::IsNaN(t2012.BdlIP)==1) continue;
	if(TMath::IsNaN(t2012.BdlIPE)==1) continue;
	if(t2012.BdlIPE<=ctmin || t2012.BdlIPE>=ctmax) continue;
	*/

		
	if(t2012.BdlIPBSc<0.02) continue; // resultado 4.4705e-02 +/-  1.29e-03
	if(t2012.BdlIPBSc>0.5) continue;
	if(TMath::IsNaN(t2012.BdlIPBSc)==1) continue;
	if(TMath::IsNaN(t2012.BdlIPBScE)==1) continue;
	if(t2012.BdlIPBScE<=ctmin || t2012.BdlIPBScE>=ctmax) continue;
	

	/*
	if(t2012.BSct<0.02) continue; // resultado 4.4670e-02 +/-  1.29e-03
	if(t2012.BSct>0.5) continue;
	if(TMath::IsNaN(t2012.BSct)==1) continue;
	if(TMath::IsNaN(t2012.BSctE)==1) continue;
	if(t2012.BSctE<=ctmin || t2012.BSctE>=ctmax) continue;
	*/

	/*
	if(t2012.BdlBS<0.02) continue;  // resultados 4.4630e-02 +/-  1.29e-03
	if(t2012.BdlBS>0.5) continue;
	if(TMath::IsNaN(t2012.BdlBS)==1) continue;
	if(TMath::IsNaN(t2012.BdlBSE)==1) continue;
	if(t2012.BdlBSE<=ctmin || t2012.BdlBSE>=ctmax) continue;
	*/

	/*
	if(t2012.BdlRf<0.02) continue;  // resultado 4.4830e-02 +/-  1.30e-03
	if(t2012.BdlRf>0.5) continue;
	if(TMath::IsNaN(t2012.BdlRf)==1) continue;
	if(TMath::IsNaN(t2012.BdlRfE)==1) continue;
	if(t2012.BdlRfE<=ctmin || t2012.BdlRfE>=ctmax) continue;
	*/

      if(t2012.mumDisplace!=1)continue;
      if(t2012.mupDisplace!=1)continue;
      if(t2012.massKs0>0.4876 && t2012.massKs0<0.5076)continue; //massKs0(PDG)=0.4976 pongo un veto a 10 MEG

      //if(t2012.triggerv4!=1)continue;//
      //if(t2012.triggerv5!=1)continue;//
      //if(t2012.triggerv6!=1)continue;// 
      //if(t2012.triggerv7!=1)continue; //
      

      M=t2012.massB;


      //Mlam=t2012.Bdl;
      //MlamE=t2012.BdlE;

      //Mlam=t2012.BdlIP;
      //MlamE=t2012.BdlIPE;
      
      Mlam=t2012.BdlIPBSc;
      MlamE=t2012.BdlIPBScE;
      
      //Mlam=t2012.BSct;
      //MlamE=t2012.BSctE;
      
      //Mlam=t2012.BdlBS;
      //MlamE=t2012.BdlBSE;
      
      //Mlam=t2012.BdlRf;
      //MlamE=t2012.BdlRfE;

      
      data2012.add(RooArgSet(M,Mlam,MlamE));
        
      hmass2012->Fill(t2012.massB);

      
    }
  
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  data2012.Print("v");
  
  //Mass model 2011 y 2012
  Double_t supM = Mmax;
  Double_t infM = Mmin;
  
  Double_t dlbsI = 1.3415e-02;
  Double_t dlblI = 4.6207e-02;
  //Double_t dlbsnI= 2.3516e-02;
  //Double_t fbsI = 0.03;
  //Double_t fblI = 0.03;
  //Double_t fpnI=0.72;
  Double_t fpI = 0.69;
  
  //Construcción de modelo de bkg para masa:
  RooRealVar a0("a0","a0",-20.0,20.0);
  RooRealVar a1("a1","a1",-20.0,20.0);
  RooRealVar a2("a2","a2",-20.0,20.0);
  RooRealVar a3("a3","a3",-20.0,20.0);
  RooRealVar b0("b0","b0",-20.0,20.0);
  RooRealVar b1("b1","b1",-20.0,20.0);
  RooRealVar b2("b2","b2",-20.0,20.0);
  RooRealVar b3("b3","b3",-20.0,20.0);
  RooPolynomial bkg2012("bkg2012","Background",M,RooArgList(b1));
  
  //Construcción de modelo de señal para masa:
  RooRealVar mean("#mu","Mass mean",5.795,5.745,5.835,"GeV");
  //RooRealVar width2012("#sigma (2012)"," Mass width2012",0.014,0.0,0.1,"GeV");
  RooRealVar width2012("width2012","width2012",0.0,0.03,"GeV");
  RooRealVar widthB2012("widthB2012"," Mass width B2012",0.0,0.2,"GeV");

  //Gausiana delgada, 2011 y 2012:
  RooGaussian Sig2012("Sig2012"," Signal PDF2012",M,mean,width2012);

  //Gausiana gruesa, 2011 y 2012:
  RooGaussian SigB2012("SigB2012"," Signal PDF B2012",M,mean,widthB2012);
  
  //********Esta seria la PDF de la suma de las dos gausianas de señal********
  RooRealVar fb2012("fb2012","fb2012",0.,1.);
  RooAddPdf sumgau2012("sumgau2012","suma de gaussianas2012",RooArgList(Sig2012,SigB2012),fb2012);
  
  //Modelo total de masa extendido:
  RooRealVar Ns2012("N_{s}(2012)","Ns2012",0.,data2012->numEntries());
  RooRealVar Nb2012("N_{b}(2012)","Nb2012",0.,data2012->numEntries());
  RooAddPdf MassModel2012("MassModel2012","MassModel2012",RooArgList(Sig2012,bkg2012),RooArgList(Ns2012,Nb2012)); //--> Sólo utilizamos una gaussiana //sumgau2012

  //Modelo total de masa no-extendido:
  RooRealVar fs2012("f_{s}(2012)","fs2012",0.,1.);
  RooAddPdf MassModel2012NE("MassModel2012NE","MassModel2012NE",RooArgList(Sig2012,bkg2012),fs2012); //--> Sólo utilizamos una gaussiana //sumgau2012

    
    
  //----------funcion de resolucion----------------------------
  
  RooRealVar s2012("s(2012)","Scale factor2012",1.0,0.1,8.0);
  RooGaussModel R2012("R2012","Resolution2012",Mlam,RooConst(0),s2012,MlamE);
  
  //Para pruebas, no resolución:
  RooTruthModel noRes2011("noRes2011","Truth Model 2011",Mlam);

  
  //----------- funcion de decaimiento-------------
  
  //señal pdl 2011 y 2012:
  RooRealVar dl("c#tau","c#tau (cm)",0.045, 0.0,0.15,"cm");
  RooDecay dmS2012("dmS2012","decay model2012",Mlam,dl,R2012,RooDecay::SingleSided);
   
  //bkg pdl 2012:
  RooRealVar dlbs2012("#lambda^{ +}_{1} (2012)","Short Proper Decay Length2012",0.015,0.0,1.,"cm");
  RooRealVar dlbl2012("#lambda^{ +}_{2} (2012)","Long Proper Decay Length2012",0.045,0.0,2.,"cm");
  RooDecay dmBl2012("dmBl2012","Sidebands decay model2012",Mlam,dlbl2012,R2012,RooDecay::SingleSided);// esta es exponencial a la derecha	       
  RooDecay dmBs2012("dmBs2012","Sidebands decay model2012",Mlam,dlbs2012,R2012,RooDecay::SingleSided);
  RooRealVar fp2012("h (2012)","Prompt fraction2012",1.0,0.0,1.0);  //fp2011 es la fracción de bkg largo.
  RooAddPdf BkgLFModel2012("BkgLFModel2012","Sideband Shape2012",RooArgList(dmBl2012,dmBs2012),RooArgList(fp2012));
  
  //Modelo total de PDL 2011 y 2012, sin eficiencia.
  RooAddPdf dmModel2012("dmModel2012","Total lf Model2012",RooArgList(dmS2012,dmBl2012),RooArgList(Ns2012,Nb2012)); //--> SOLO UNA COMPONENTE LARGA: dmBl2012 (BkgLFModel2012)
  RooAddPdf dmModel2012NE("dmModel2012NE","Total lf Model2012 NE",RooArgList(dmS2012,dmBl2012),fs2012);  //-->SOLO UNA COMPONENTE LARGA: dmBl2012 (BkgLFModel2012)

  //Correción (eficiencia): solo la tenemos para 2011.
  RooRealVar a("a","a",0.991559); //
  RooRealVar b("b","b",0.091770);
  RooFormulaVar gp("gp","a+b*Mlam",RooArgList(a,b,Mlam)) ;
        
  //Modelo de PDL señal 2012 corregido con eficiencia:
  RooEffProd dmSeff2012("dmseff2012","model with efficiency2012",dmS2012,gp) ;
  
  //Modelo total de PDL 2011 y 2012, CON eficiencia.
  //FALTA PROGRAMAR. No necesario por ahora.
    
  /*
    Modelación de error en PDL.
  */

  
  //señal2012 error en PDL: exp (conv) Gauss
  RooRealVar wE2012("#sigma_{#sigma} (2012)"," Error width2012",1e-4,1e-6,0.01);
  RooRealVar mE2012("#mu_{#sigma} (2012)","error media2012",0.002,0.,0.1);
  RooGaussModel RE2012("RE2012","Resolution2012",MlamE,mE2012,wE2012);
  RooRealVar dlE2012("#gamma (2012)","c#tauE2012 (cm)",1e-4,0.,1e-3);// //1e-8
  RooDecay SigE2012("SigE2012","decay model2012",MlamE,dlE2012,RE2012,RooDecay::SingleSided);
  
  //Ruido2012 error en PDL: expBkg (conv) GaussBkg //Se podría usar: [feb*expShort + (1-feb)*expLong] (conv) GaussBkg
  //RooRealVar feb2012("feb2012","feb2012",0.5,0.,1.0);
  RooRealVar widthEB2012("#sigma^{B}_{#sigma} (2012)"," Error width for bkg2012",1e-5,0.1);
  RooRealVar meanEB2012("#mu^{B}_{#sigma} (2012)","Error mean for bkg2012",1e-4,0.1);
  RooGaussModel REB2012("REB2012","Error resolution for bkg2012",MlamE,meanEB2012,widthEB2012);
  RooRealVar dlEB2012("#gamma^{B} (2012)","Error decay for bkg2012",1e-5,0.,1e-3);
  //RooRealVar dlEB22012("dlEB22012","Error decay for bkg2012",1e-5,0.1);
  RooDecay bkgE2012("bkgE2012","Error decay model for bkg2012",MlamE,dlEB2012,REB2012,RooDecay::SingleSided);
  //RooDecay bkgE22012("bkgE22012","Error decay model for bkg2012",MlamE,dlEB22012,REB2012,RooDecay::SingleSided);
  //RooAddPdf bkgET2012("bkgET2012", "error total2012",RooArgList(bkgE2012,bkgE22012),RooArgList(feb2012) );
 
  //Modelo error total 2012:
  RooAddPdf sumE2012("sumE2012","Total lfE Model2012",RooArgList(SigE2012,bkgE2012),RooArgList(Ns2012,Nb2012));
  
  
  /*
   Se construyen modelos totales
  */
    
  //Señal total 2011 y 2012:
  //RooProdPdf dmModeltotal2012("dmModeltotalm2012","decay model2012",RooArgList(sumgau2012,dmSeff2012,SigE2012));
  RooProdPdf dmModeltotal2012("dmModeltotalm2012","decay model2012",RooArgList(Sig2012,dmS2012,SigE2012));

  
  //Bkg total 2011 y 2012:
  RooProdPdf Background2012("Background2012","Background2012",RooArgList(bkg2012,dmBl2012,bkgE2012));  //---> ANTES: SOLO SE ESTA USANDO COMPONENTE LARGA (dmBl2012,BkgLFModel2012)

  
  //Señal + Bkg total 2011 Y 2012:
  RooAddPdf sum2012("sum2012","Total lf Model2012",RooArgList(dmModeltotal2012,Background2012),RooArgList(Ns2012,Nb2012));
    
  //Señal + Bkg total 2011 Y 2012 no extendido:
  RooAddPdf sumNE2012("sumBE2012","Total lf Model2012 No-Extended",RooArgList(dmModeltotal2012,Background2012),fs2012);
    
  //Guarda set de parametros en:
  //RooArgSet* variables2011 = new RooArgSet(Ns2011,Nb2011,mean,width2011,dl,a1,dlbl2011,s2011);

 
 
  
  /* 
     Ajustando los parametros para el modelo de masa
  */
  
  //Ajuste solo Masa 2011
  //Nb2011.setVal(800);
  mean.setVal(5.795);
  //  width2011.setConstant(kTRUE);
  a1.setVal(0);//15
  a1.setConstant(kTRUE); //Obligamos a una constante (fit consistente a 3 sigma parabólico), Minos tiene problemas para calcular error positivo.
    
    //Ajuste solo Masa 2012:
    Ns2012.setVal(80);// 2049,389
    //Nb2011.setVal(1000);
    mean.setVal(5.795);
    //width2012.setVal(0.014);
    //widthB2012.setVal(0.02);
    b1.setVal(0);//15
    b1.setConstant(kTRUE); //Obligamos a una constante (fit completamente consistente con cero con error asimétrico).

    
  RooFitResult* fitres2012 = MassModel2012.fitTo(data2012,Extended(),Minos(kTRUE),Save(kTRUE), NumCPU(4));
  fitres2012->Print("v");

  
  TCanvas *cm2 = new TCanvas("cm2","",750,550);
  cm2.cd(1);
  RooPlot* Mframeonly2012 = M.frame(Mmin,Mmax,40);
  data2012.plotOn(Mframeonly2012);
  MassModel2012->plotOn(Mframeonly2012);
  MassModel2012->paramOn(Mframeonly2012,Layout(0.55, 0.99,0.99),Format("NELU", AutoPrecision(2)));
  Mframeonly2012->Draw();
  return; 
    
  /*
   Ajuste de distribuciones de error.
  */
    
  //Define sidebands
  Double_t sL = 5.795-3*0.0155;
  Double_t sR = 5.7950+3*0.0155;
  Double_t infM=Mmin;
  Double_t supM=Mmax;
  TString sdCut1 = "M>"; sdCut1+=sR; sdCut1+=" && "; sdCut1+=" M<";sdCut1+=supM; 
  TString sdCut2 = "M>"; sdCut2+=infM; sdCut2+=" && "; sdCut2+=" M<";sdCut2+=sL;

 
  
  //Define 2012 sidebands dataset
  RooDataSet* dataSidebands2012 = (RooDataSet*)data2012.reduce(sdCut1);
  RooDataSet* dataSidebands22012 = (RooDataSet*)data2012.reduce(sdCut2);
  dataSidebands2012->append(*dataSidebands22012);


  //Define 2011 and 2012 signal regions:
  TString sdCut = "M>"; sdCut+=sL; sdCut+=" && "; sdCut+=" M<";sdCut+=sR;
  RooDataSet* dataSignalreg2011 = (RooDataSet*)data2011.reduce(sdCut);
  RooDataSet* dataSignalreg2012 = (RooDataSet*)data2012.reduce(sdCut);
 
  M.setRange("SignalRegion",sL,sR);
  RooAbsReal* fracSigRange2011 = Sig2011.createIntegral(M,M,"SignalRegion");
  RooAbsReal* fracBkgRange2011 = bkg2011.createIntegral(M,M,"SignalRegion");
  RooAbsReal* fracSigRange2012 = Sig2012.createIntegral(M,M,"SignalRegion");
  RooAbsReal* fracBkgRange2012 = bkg2012.createIntegral(M,M,"SignalRegion");
  Double_t S_SignalRegion2011 = Ns2011.getVal() *  fracSigRange2011->getVal();
  Double_t B_Signalregion2011 = Nb2011.getVal() *  fracBkgRange2011->getVal();
  Double_t S_SignalRegion2012 = Ns2012.getVal() *  fracSigRange2012->getVal();
  Double_t B_Signalregion2012 = Nb2012.getVal() *  fracBkgRange2012->getVal();
  cout << "Eventos de señal en region de señal 2011: " << S_SignalRegion2011 << endl;
  cout << "Eventos de background en region de señal 2011: " << B_Signalregion2011 << endl;
  cout << "Eventos de señal en region de señal 2012: " << S_SignalRegion2012 << endl;
  cout << "Eventos de background en region de señal 2012: " << B_Signalregion2012 << endl;

  //Fit to 2011 error distribution for background
  meanEB2011.setVal(0.0015); //si es una sola exponencial es 0.002464
  widthEB2011.setVal(0.0003); //si es una sola exponencial es 0.000490
  dlEB2011.setVal(0.0009); 
  
  RooFitResult* ErrBFit2011 = bkgE2011.fitTo(*dataSidebands2011,/* Extended(),*/Hesse(kTRUE),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  ErrBFit2011->Print("v");

  TCanvas *c1 = new TCanvas("c1","",550,550);
  c1.Divide(1,2);
  c1.cd(1);
  RooPlot* MframebkgE2011 = MlamE.frame(ctmin,ctmax,20);
  dataSidebands2011.plotOn(MframebkgE2011 /*,DataError(RooAbsData::SumW2)*/);
  //data.plotOn(MframebkgE /*,DataError(RooAbsData::SumW2)*/);
  bkgE2011->plotOn(MframebkgE2011);
  Double_t chi2_bkg = MframebkgE2011->chiSquare();
  //bkgE2011->paramOn(MframebkgE2011,Layout(0.60),Format("NELU", AutoPrecision(2)));
  MframebkgE2011->Draw();
  //return;  
  
  //Fit to 2012 error distribution for background
  meanEB2012.setVal(0.0021); //si es una sola exponencial es 0.002464
  widthEB2012.setVal(0.0006); //si es una sola exponencial es 0.000490
  dlEB2012.setVal(0.00049); 
  
  RooFitResult* ErrBFit2012 = bkgE2012.fitTo(*dataSidebands2012,/* Extended(),*/Hesse(kTRUE),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  ErrBFit2012->Print("v");
    
  c1.cd(2);
  RooPlot* MframebkgE2012 = MlamE.frame(ctmin,ctmax,25);
  dataSidebands2012.plotOn(MframebkgE2012 /*,DataError(RooAbsData::SumW2)*/);
  //data.plotOn(MframebkgE /*,DataError(RooAbsData::SumW2)*/);
  bkgE2012->plotOn(MframebkgE2012);
  Double_t chi2_bkg = MframebkgE2012->chiSquare();
  //bkgE2012->paramOn(MframebkgE2012,Layout(0.60),Format("NELU", AutoPrecision(2)));
  MframebkgE2012->Draw();
  //return;  
  //Fit to 2011 error distribution (signal +bkg). Fix bkg pdl error parameters, Ns and Nb, and just let signal params to float.
  meanEB2011.setConstant(kTRUE);
  widthEB2011.setConstant(kTRUE);
  dlEB2011.setConstant(kTRUE);

/* METHOD 1: FIT SIGNAL + BKG IN ALL MASS REGION
  Ns2011.setConstant(kTRUE);
  Nb2011.setConstant(kTRUE);
  
  RooFitResult* ErrSFit2011 = sumE2011.fitTo(data2011, Extended(), Hesse(kTRUE),Minos(kTRUE),Save(kTRUE), NumCPU(4));
  ErrSFit2011->Print("v");
 
  TCanvas *c2 = new TCanvas("c2","",550,550);
  c2.Divide(1,2);
  c2.cd(1);
  RooPlot* MframebkgET2011 = MlamE.frame(ctmin,ctmax,25);
  data2011.plotOn(MframebkgET2011 ); //DataError(RooAbsData::SumW2)
  sumE2011->plotOn(MframebkgET2011);
  Double_t chi2_bkg2011 = MframebkgET2011->chiSquare();
  sumE2011->paramOn(MframebkgET2011,Layout(0.70),Format("NELU", AutoPrecision(2)));
  MframebkgET2011->Draw();
 
 //End method 1.
 */

  //METHOD 2: FIT SIGNAL + BKG IN ONLY SIGNAL REGION in 2011
  Double_t Ns2011tmp = Ns2011.getVal();
  Double_t Nb2011tmp = Nb2011.getVal();
  Ns2011.setVal(S_SignalRegion2011);
  Nb2011.setVal(B_Signalregion2011);
  Ns2011.setConstant(kTRUE);
  Nb2011.setConstant(kTRUE);

  
  wE2011.setVal(0.0013);
  mE2011.setVal(0.0022);
  dlE2011.setVal(0.0007);
  

  RooFitResult* ErrSFit2011 = sumE2011.fitTo(*dataSignalreg2011, Extended(), Hesse(kTRUE),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  ErrSFit2011->Print("v");
    
  TCanvas *c2 = new TCanvas("c2","",550,550);
  c2.Divide(1,2);
  c2.cd(1);
  RooPlot* MframebkgET2011 = MlamE.frame(ctmin,ctmax,25);
  dataSignalreg2011->plotOn(MframebkgET2011 ); //DataError(RooAbsData::SumW2)
  sumE2011->plotOn(MframebkgET2011);
  sumE2011->plotOn(MframebkgET2011,Components(bkgE2011),LineStyle(kDashed), LineColor(kGray));
  sumE2011->plotOn(MframebkgET2011,Components(SigE2011),LineStyle(kDotted), LineColor(kRed));
  Double_t chi2_bkg2011 = MframebkgET2011->chiSquare();
  //sumE2011->paramOn(MframebkgET2011,Layout(0.60),Format("NELU", AutoPrecision(2)));
  MframebkgET2011->Draw();
  //return;
  Ns2011.setConstant(kFALSE);
  Nb2011.setConstant(kFALSE);
  Ns2011.setVal(Ns2011tmp);
  Nb2011.setVal(Nb2011tmp);

  //End Method 2 for 2011.

  mE2011.setConstant(kTRUE);
  wE2011.setConstant(kTRUE);
  dlE2011.setConstant(kTRUE);
  Ns2011.setConstant(kFALSE);
  Nb2011.setConstant(kFALSE);
    
  //Fit to 2012 error distribution (signal +bkg). Fix bkg pdl error parameters, Ns and Nb, and just let signal params to float.
  meanEB2012.setConstant(kTRUE);
  widthEB2012.setConstant(kTRUE);
  dlEB2012.setConstant(kTRUE);

  /* METHOD 1: FIT SIGNAL + BKG IN ALL MASS REGION in 2012
  Ns2012.setConstant(kTRUE);
  Nb2012.setConstant(kTRUE);
  
  RooFitResult* ErrSFit2012 = sumE2012.fitTo(data2012, Extended(), Hesse(kTRUE),Minos(kTRUE),Save(kTRUE), NumCPU(4));
  ErrSFit2012->Print("v");
    
  c2.cd(2);
  RooPlot* MframebkgET2012 = MlamE.frame(ctmin,ctmax,20);
  data2012.plotOn(MframebkgET2012 ); //DataError(RooAbsData::SumW2)
  sumE2012->plotOn(MframebkgET2012);
  Double_t chi2_bkg2012 = MframebkgET2012->chiSquare();
  sumE2012->paramOn(MframebkgET2012,Layout(0.70),Format("NELU", AutoPrecision(2)));
  MframebkgET2012->Draw();
    
  */

  //METHOD 2: FIT SIGNAL + BKG IN ONLY SIGNAL REGION in 2012
  Double_t Ns2012tmp = Ns2012.getVal();
  Double_t Nb2012tmp = Nb2012.getVal();
  Ns2012.setVal(S_SignalRegion2012);
  Nb2012.setVal(B_Signalregion2012);
  Ns2012.setConstant(kTRUE);
  Nb2012.setConstant(kTRUE);

  wE2012.setVal(0.0004);
  mE2012.setVal(0.0017);
  dlE2012.setVal(0.0003);
    
  RooFitResult* ErrSFit2012 = sumE2012.fitTo(*dataSignalreg2012, Extended(), Hesse(kTRUE),Minos(kFALSE),Save(kTRUE), NumCPU(4));
  ErrSFit2012->Print("v");
    
  c2.cd(2);
  RooPlot* MframebkgET2012 = MlamE.frame(ctmin,ctmax,25);
  dataSignalreg2012->plotOn(MframebkgET2012 ); //DataError(RooAbsData::SumW2)
  sumE2012->plotOn(MframebkgET2012);
  sumE2012->plotOn(MframebkgET2012,Components(bkgE2012),LineStyle(kDashed), LineColor(kGray));
  sumE2012->plotOn(MframebkgET2012,Components(SigE2012),LineStyle(kDotted), LineColor(kRed));
  Double_t chi2_bkg2012 = MframebkgET2012->chiSquare();
  //sumE2012->paramOn(MframebkgET2012,Layout(0.60),Format("NELU", AutoPrecision(2)));
  MframebkgET2012->Draw();
  //return;
  Ns2012.setConstant(kFALSE);
  Nb2012.setConstant(kFALSE);
  Ns2012.setVal(Ns2012tmp);
  Nb2012.setVal(Nb2012tmp);
    
  //End Method 2 for 2012.

    
  mE2012.setConstant(kTRUE);
  wE2012.setConstant(kTRUE);
  dlE2012.setConstant(kTRUE);
  Ns2012.setConstant(kFALSE);
  Nb2012.setConstant(kFALSE);
    
  //ajustando los parametros para todas las PDFs juntas.
  
  //Inicialización:
  dl.setVal(0.045);
  //2011
  s2011.setVal(1.0);
  s2011.setConstant(kTRUE);
  dlbl2011.setVal(0.05);//0.1
  //fp2011.setVal(0.5);
  //2012
  //s2012.setVal(1.);
  //s2012.setConstant(kTRUE);
  dlbl2012.setVal(0.05);//0.1
  //fp2012.setVal(0.5);
  fs2011->setVal(Ns2011->getVal()/data2011->numEntries());
  fs2012->setVal(Ns2012->getVal()/data2012->numEntries());
  
  
  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;                                                                          
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;
    
  //  width2011.setConstant(kFALSE);
  
  //Ajustes previos por año:
  //width2011.setConstant(kTRUE); //cannot fit 2011 mass width, so set to combined 2011-2012 (could also fix to MC resolution).
  RooFitResult* fitres1= sumNE2011.fitTo(data2011,ConditionalObservables(MlamE),Minos(kFALSE),Save(kTRUE), NumCPU(4),PrintLevel(2));//,Extended()
  fitres1->Print("v");
  variables2011->writeToFile("fit2011.txt");
    
  //width2011.setConstant(kFALSE); //in 2012 mass width can be determined, so release this variable.
  RooFitResult* fitres2= sumNE2012.fitTo(data2012,ConditionalObservables(MlamE),Minos(kFALSE),Save(kTRUE), NumCPU(4),PrintLevel(2));//,Extended()
  fitres2->Print("v");
  variables2012->writeToFile("fit2012.txt");

  //Ajuste simultaneo:
  RooFitResult* fitres= simPdf.fitTo(combData,ConditionalObservables(MlamE),Hesse(kTRUE),Minos(kFALSE),Save(kTRUE), NumCPU(4),PrintLevel(2));//,Extended()
  fitres->Print("v");
  variablessim->writeToFile("fit2011_2012.txt");

  
  //Plot de masas:
  TCanvas *c3 = new TCanvas("c3","",550,550);
  c3.Divide(1,2);
    
  c3.cd(1);   
  RooPlot* Mframe2011 = M.frame(infM,supM,40);
  combData.plotOn(Mframe2011,Cut("sample==sample::sample2011"));
  simPdf.plotOn(Mframe2011,Slice(sample,"sample2011"),Components("bkg2011"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGray)) ;
  simPdf.plotOn(Mframe2011,Slice(sample,"sample2011"),Components("Sig2011"),ProjWData(sample,combData),LineStyle(kDotted),LineColor(kRed)) ;
  simPdf.plotOn(Mframe2011,Slice(sample,"sample2011"),ProjWData(sample,combData)) ;
  //Mframe2011->SetMinimum(1.0);
  if(showbox)  MassModel2011NE.paramOn(Mframe2011, Format("NELU", AutoPrecision(2)), Layout(0.55, 0.9,0.9) );
  Mframe2011->Draw();
    
  c3.cd(2); 
  RooPlot* Mframe2012 = M.frame(infM,supM,40);
  combData.plotOn(Mframe2012,Cut("sample==sample::sample2012")) ;
  simPdf.plotOn(Mframe2012,Slice(sample,"sample2012"),Components("bkg2012"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kGray)) ;
  simPdf.plotOn(Mframe2012,Slice(sample,"sample2012"),Components("Sig2012"),ProjWData(sample,combData),LineStyle(kDotted),LineColor(kRed)) ;
  simPdf.plotOn(Mframe2012,Slice(sample,"sample2012"),ProjWData(sample,combData)) ;
  //Mframe2012->SetMinimum(1.0);
  if(showbox)  MassModel2012NE.paramOn(Mframe2012, Format("NELU", AutoPrecision(2)), Layout(0.55, 0.9,0.9) );
  Mframe2012->Draw();

  //Plot de PDLs:
  TCanvas *c4 = new TCanvas("c4","",550,550);
  c4.Divide(1,2);
      
  c4.cd(1);
  gPad->SetLogy();
  RooPlot* Mframelam2011 = Mlam.frame(0.02,0.5,60);
  combData.plotOn(Mframelam2011,Cut("sample==sample::sample2011")) ;
  //simPdf.plotOn(Mframelam11,Slice(sample,"sample2011"),ProjWData(sample,combData)) ;
  //simMPdf.plotOn(Mframelam11,Slice(sample,"sample2011"),Components("dms11"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));
  //simPdf.plotOn(Mframelam11,Slice(sample,"sample2011"),ProjWData(RooArgSet(M,MlamE,sample),combData),LineColor(kBlue),LineWidth(3.0));
  dmBl2011.plotOn(Mframelam2011,ProjWData(MlamE,*dataSidebands2011),  Normalization(Nb2011->getVal(),RooAbsReal::NumEvent) ,LineStyle(kDashed),LineColor(kGray));
  dmS2011.plotOn(Mframelam2011,ProjWData(MlamE,*dataSignalreg2011),  Normalization(Ns2011->getVal(),RooAbsReal::NumEvent) ,LineStyle(kDotted),LineColor(kRed));
  dmModel2011.plotOn(Mframelam2011,ProjWData(MlamE,data2011),  Normalization(data2011->numEntries(),RooAbsReal::NumEvent) ,LineColor(kBlue),LineWidth(3.0));
  Double_t chi2tmp2011 = Mframelam2011.chiSquare();
  Mframelam2011->SetMinimum(1.); // esto es para integrar la variable MlamE antes de graficar
  if(showbox)  dmModel2011NE.paramOn(Mframelam2011, Format("NELU", AutoPrecision(2)), Layout(0.55, 0.9,0.9) );
  Mframelam2011->Draw();
  
  c4.cd(2);
  gPad->SetLogy();
  RooPlot* Mframelam2012 = Mlam.frame(0.02,0.5,60);
  combData.plotOn(Mframelam2012,Cut("sample==sample::sample2012")) ;
  //simPdf.plotOn(Mframelam12,Slice(sample,"sample2012"),ProjWData(sample,combData)) ;
  //simMPdf.plotOn(Mframelam12,Slice(sample,"sample2012"),Components("dms12"),ProjWData(sample,combData),LineStyle(kDashed),LineColor(kRed));
  //    dmModel12.plotOn(Mframelam12,Slice(sample,"sample2012"),ProjWData(RooArgSet(Mlam,MlamE,sample),combData),LineColor(kBlue),LineWidth(3.0));
  dmBl2012.plotOn(Mframelam2012,ProjWData(MlamE,*dataSidebands2012),  Normalization(Nb2012->getVal(),RooAbsReal::NumEvent) ,LineStyle(kDashed),LineColor(kGray)); //Solo componente larga.
  //BkgLFModel2012.plotOn(Mframelam2012,ProjWData(MlamE,*dataSidebands2012),  Normalization(Nb2012->getVal(),RooAbsReal::NumEvent) ,LineStyle(kDashed),LineColor(kMagenta));
  dmS2012.plotOn(Mframelam2012,ProjWData(MlamE,*dataSignalreg2012),  Normalization(Ns2012->getVal(),RooAbsReal::NumEvent) ,LineStyle(kDotted),LineColor(kRed));
  dmModel2012.plotOn(Mframelam2012,ProjWData(MlamE,data2012),  Normalization(data2012->numEntries(),RooAbsReal::NumEvent) ,LineColor(kBlue),LineWidth(3.0));
  Double_t chi2tmp2012 = Mframelam2012.chiSquare();
  //simPdf.plotOn(Mframelam12,Slice(sample,"sample2012"),ProjWData(RooArgSet(M,MlamE,sample),combData),LineColor(kBlue),LineWidth(3.0));
  Mframelam2012->SetMinimum(1.); // esto es para integrar la variable MlamE antes de graficar
  if(showbox)  dmModel2012NE.paramOn(Mframelam2012, Format("NELU", AutoPrecision(2)), Layout(0.55, 0.9,0.9) );
  Mframelam2012->Draw();
  
  
  cout<<" Chi2 2011 : "<<chi2tmp2011<< endl;
  cout<<" Chi2 2012 : "<<chi2tmp2012<< endl;
  
    
    TCanvas *ch1 = new TCanvas("ch1","",750,550);
    hmass2011->SetXTitle("J/#psi #Xi^{-} invariant mass (GeV)");
    hmass2011->SetYTitle(Form("Events / %0.2f",(Mmax-Mmin)/40));
    ch1->cd();
    hmass2011->Draw();
    
    TCanvas *ch2 = new TCanvas("ch2","",750,550);
    hmass2012->SetXTitle("J/#psi #Xi^{-} invariant mass (GeV)");
    hmass2012->SetYTitle(Form("Events / %0.2f",(Mmax-Mmin)/40));
    ch2->cd();
    hmass2012->Draw();

    
 // return;
 cout<<"Signal fraction"<<"  &  $f_s$  &  $"<<fs2011.getVal()<<" \\pm "<<fs2011.getError()<<"$  &  $"<<fs2012.getVal()<<" \\pm "<<fs2012.getError()<<"$ \\\\"<<endl;
 cout<<"\\hline"<<endl;

 cout<<"Mass"<<"  &  $M_B$ (GeV) &  $"<<mean.getVal()<<" \\pm "<<mean.getError()<<"$ & simult. to 2011  \\\\" <<endl;
 cout<<"\\hline"<<endl;
 
 cout<<"Width of the  mass Gaussian"<<"  &  $\\sigma_1$ (GeV) &  $"<<width2011.getVal()<<" \\pm "<<width2011.getError()<< "$  & (forced) simult. to 2011 \\\\"<<endl;
 cout<<"\\hline"<<endl;


 cout<<"Bkg. mass slope"<<"  &    $p$ (GeV$^{-1})$     &             0 (fixed)                      & 0 (fixed)     \\\\"<<endl;
 cout<<"\\hline"<<endl;

 cout<<"Lifetime"<<" &   $c\\tau_B$ (cm)  &  $"<<dl.getVal()<<" \\pm "<<dl.getError()<<"$ & simult. to 2011  \\\\" <<endl;
 cout<<"\\hline"<<endl;

cout<<"Bkg. PDL, exp. dec. param. & $\\lambda_1^{+}$ (cm) & $"<<dlbl2011.getVal()<<" \\pm "<<dlbl2011.getError()<<"$  &  $"<<dlbl2012.getVal()<<" \\pm "<<dlbl2012.getError()<<"$ \\\\"<<endl;
 cout<<"\\hline"<<endl;


 //cout<<"Fraction of the 1st PDL bkg. dec. & $\\lambda_1^{+}$ (cm) & $"<<fp1.getVal()<<" \\pm "<<fp1.getError()<<"$  &  $"<<fp.getVal()<<" \\pm "<<fp.getError()<<"$ \\\\"<<endl;
 // cout<<"\\hline"<<endl;

 //cout<<"Resolution scale factor 	  &  $s$ &    $"<<s1.getVal()<<" \\pm "<<s1.getError()<<"$  &  $"<<s.getVal()<<" \\pm "<<s.getError()<<"$ \\\\"<<endl;
 //cout<<"\\hline"<<endl;

    
  cm1->Print("massonly_2011Xib.pdf");
  cm2->Print("massonly_2012Xib.pdf");
  c1->Print("epdlbkg_2011_2012Xib.pdf");
  c2->Print("epdl_signalregion_2011_2012Xib.pdf");
  c3->Print("mass3DfitXib.pdf");
  c4->Print("pdl3DfitXib.pdf");
  ch1->Print("masshisto_2011Xib.pdf");
  ch2->Print("masshisto_2012Xib.pdf");

  /*
    cm1->Print("massonly_2011.eps");
    cm2->Print("massonly_2012.eps");
    c1->Print("epdlbkg_2011_2012.eps");
    c2->Print("epdl_signalregion_2011_2012.eps");
    c3->Print("mass3Dfit.eps");
    c4->Print("pdl3Dfit.eps");
    ch1->Print("masshisto_2011.eps");
    ch2->Print("masshisto_2012.eps");
  */

/*
    cm1->Print("massonly_2011.jpg");
    cm2->Print("massonly_2012.jpg");
    c1->Print("epdlbkg_2011_2012.jpg");
    c2->Print("epdl_signalregion_2011_2012.jpg");
    c3->Print("mass3Dfit.jpg");
    c4->Print("pdl3Dfit.jpg");
    
    cm1->Print("massonly_2011.gif");
    cm2->Print("massonly_2012.gif");
    c1->Print("epdlbkg_2011_2012.gif");
    c2->Print("epdl_signalregion_2011_2012.gif");
    c3->Print("mass3Dfit.gif");
    c4->Print("pdl3Dfit.gif");
    
    cm1->Print("massonly_2011.png");
    cm2->Print("massonly_2012.png");
    c1->Print("epdlbkg_2011_2012.png");
    c2->Print("epdl_signalregion_2011_2012.png");
    c3->Print("mass3Dfit.png");
    c4->Print("pdl3Dfit.png");
 */
    

    
}
