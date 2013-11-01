// Alexandra Oliveira
//
// template by Serguei and Chiara
// 
using namespace RooFit;
using namespace RooStats ;
//
// this code is frozen to use 2 b tag categories with leveled exponential
//
// declare the functions
void AddSigData(RooWorkspace* w, Float_t, const char* filename);
void AddBkgData(RooWorkspace* w, Float_t, const char* filename, bool dobands);
void MakeDataCard(RooWorkspace* sig, RooWorkspace* bkg, const char* fileBaseName, const char* fileBkgName);
void MakeDataCardoneCat(RooWorkspace* sig, RooWorkspace* bkg, const char* fileBaseName, const char* fileBkgName);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);
//
// container for the fit results
RooFitResult* fitresult[2]; //NCAT=2 
//
const int minfit =100, maxfit=180;
RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mtot  = new RooRealVar("mtot","M(#gamma#gamma jj)",200,1200,"GeV");
  RooRealVar* mjj  = new RooRealVar("mjj","M(jj)",90,190,"GeV");
  RooRealVar* mgg  = new RooRealVar("mgg","M(#gamma#gamma)",80,180,"GeV");
  RooRealVar* evWeight   = new RooRealVar("evWeight","HqT x PUwei",0,1000000000,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 2") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  //
  RooArgSet* ntplVars = new RooArgSet(*mtot, *mgg, *cut_based_ct, *evWeight, *mjj);
//  RooArgSet* ntplVars = new RooArgSet(*mtot, *mgg, *cut_based_ct, *evWeight);
  ntplVars->add(*mtot);
  ntplVars->add(*mjj);
  ntplVars->add(*mgg);
  ntplVars->add(*cut_based_ct);
  ntplVars->add(*evWeight);
  return ntplVars;
}
////////////////////////////////////////////////////////////////////////
void runfits(const Float_t mass=120, Int_t mode=1)
{
  style();
  TString fileBaseName(TString::Format("hgg.mH300_8TeV", mass));
  TString fileBkgName(TString::Format("hgg.mH300.inputbkg_8TeV", mass));
  RooFitResult* fitresults;
  // Add the signal and background models to the workspace.
  // the minitree to be addeed
  TString ssignal = "MiniTrees/OlivierOc13/v15_base_mgg_0_massCutVersion0/02013-10-30-Radion_m300_8TeV_nm_m300.root";
  TString ddata   = "MiniTrees/OlivierOc13/v15_base_mgg_0_massCutVersion0/02013-10-30-Data_m300.root";
  //TString ddata   = "MiniTrees/OlivierAug13/v02_regkin_mgg_0/Data_regression-m300_minimal.root";
  //
  cout<<"Signal: "<< ssignal<<endl;
  cout<<"Data: "<< ddata<<endl;
  // declare the worspace outside the function to  construct the datacard
  // for signal
  RooWorkspace *wall = new RooWorkspace("w_all","w_all");
  TString card_nameS("models_mgg_sig.rs"); // put the model parameters here!
  HLFactory hlfS("w_allS", card_nameS, false);
  RooWorkspace* wall = hlfS.GetWs(); // Get the RooWorkspace containing the models and variables
  AddSigData(wall,mass,ssignal,fileBaseName); //OK!
  // for BKG
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  TString card_nameB("models_mgg_bkg.rs"); // put the model parameters here!
  HLFactory hlfB("w_allB", card_nameB, false);
  RooWorkspace* wAll = hlfB.GetWs(); // Get the RooWorkspace containing the models and variables 
  bool dobands=false;
  AddBkgData(wAll,ddata,fileBkgName,dobands);
  //
  MakeDataCard(wall,wAll,fileBaseName,fileBkgName);
  MakeDataCardoneCat(wall,wAll,fileBaseName,fileBkgName);
  cout<< "here"<<endl;
  return;
} // close runfits
////////////////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddSigData(RooWorkspace* wall, Float_t mass, TString signalfile, const char* fileBasename) { 
  // file input
  TFile sigFile(signalfile);  
  TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
  //Float_t Lum = 19785.0; // pb-1
  //RooRealVar lumi("lumi","lumi",Lum);
  //wall->import(lumi); 
  //
  RooArgSet* ntplVars = defineVariables();
  //
  const Int_t ncat = 2;
  Float_t MASS(mass);
  // common preselection cut
  TString mainCut("1");
  //--------------------------------------------
  // re-scale signal to the ratio r to be in fb
  // no cross section added on weights 
  //--------------------------------------------
  // one channel with default weights
  /* RooDataSet sigScaled1(
	"sigScaled1",
	"dataset",
	sigTree, // all variables of RooArgList
	*ntplVars,
	mainCut);
  RooRealVar *k = new RooRealVar("k", "k", 0.0019706 ,"GeV"); // 19.706/10000
  RooFormulaVar *nw = new RooFormulaVar("nw", "nw", "@0*@1", RooArgSet(*evWeight, *k));
  sigScaled1.addColumn(*nw);
  //sigScaled1.setWeightVar("nw"); 
  //sigScaled1->Print();
  RooArgSet* ntplVars1 = sigScaled1.get();  
  //ntplVars->add(nw);
  */
  RooRealVar *evWeight = (RooRealVar*) (*ntplVars)["evWeight"] ;
  RooDataSet sigScaled(
	"sigScaled", 
	"dataset",
	sigTree, 
	*ntplVars,
	mainCut, 
	"evWeight");
  RooDataSet* sigToFit[ncat];
  //----------------------------------------------
  // cuts
  TString cut0 = "&& 1>0";// "&& mtot > 450 && mtot < 550 "; //"&& 1>0";//
  TString cut1 ="&& 1>0";// "&& mtot > 450 && mtot < 550 "; // "&& 1>0";//
  //
  TString cutj0 ="&& 1>0";// "&& mjj > 90 && mjj < 170 "; //"&& 1>0";//
  TString cutj1 ="&& 1>0";// "&& mjj > 100 && mjj < 160 "; // "&& 1>0";//
  //
  // we take only mtot to fit to the workspace, we include the cuts
  sigToFit[0] = (RooDataSet*) sigScaled.reduce(
	*wall->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  wall->import(*sigToFit[0],Rename(TString::Format("Sig_cat%d",0)));
    //
  sigToFit[1] = (RooDataSet*) sigScaled.reduce(
	*wall->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  wall->import(*sigToFit[1],Rename(TString::Format("Sig_cat%d",1)));

  //////////////////////////////////////////////////////////////////////
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel:  " << sigToFit[0]->sumEntries() +sigToFit[1]->sumEntries() << endl; 
  cout << "---- one channel:  " << sigScaled.sumEntries() << endl; 
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = sigToFit[c]->sumEntries();
    cout << TString::Format("nEvt exp.  cat%d : ",c) << nExpEvt 
	 << TString::Format("   eff x Acc  cat%d : ",c) 
	 << "%" 
	 << endl; 
  }
  cout << "======================================================================" << endl;
  sigScaled.Print("v");
  // we end adding signal function
  //////////////////////////////////////////////////////////////////////////
  // we do the model fitting on the same function and save to anew WS
  Float_t MASS(mass);
  // fit range
  Float_t minMassFit(minfit),maxMassFit(maxfit); 
  //
  RooDataSet* sigFit[ncat];
  RooAbsPdf* mggSig[ncat];
  for (int c = 0; c < ncat; ++c) {
    // import sig from workspace
    sigFit[c] = (RooDataSet*) wall->data(TString::Format("Sig_cat%d",c));
    mggSig[c]     = (RooAbsPdf*)  wall->pdf(TString::Format("mggSig_cat%d",c));
    cout << "OK up to now..." <<MASS<< endl;
    ((RooRealVar*) wall->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(MASS);
    // Fit model as M(x|y) to D(x,y)
    mggSig[c]->fitTo(*sigFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE));
    cout << "OK up to now again ..." <<MASS<< endl;
    // IMPORTANT: fix all pdf parameters to constant
    wall->defineSet(TString::Format("SigPdfParam_cat%d",c), 
        RooArgSet(
	 *wall->var(TString::Format("mgg_sig_m0_cat%d",c)),   
	 *wall->var(TString::Format("mgg_sig_alpha_cat%d",c)),
	 *wall->var(TString::Format("mgg_sig_n_cat%d",c)), 
	 *wall->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
	 *wall->var(TString::Format("mgg_sig_sigma_cat%d",c)),
	 *wall->var(TString::Format("mgg_sig_frac_cat%d",c))) 
	);
    SetConstantParams(wall->set(TString::Format("SigPdfParam_cat%d",c)));
  } // close for ncat
  // copy in annother WS before to save
  RooAbsPdf* mggSigPdfTosave[ncat];
  RooWorkspace *wAllSave = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    // import pdfs from previous fit calculation
    mggSigPdfTosave[c] = (RooAbsPdf*)  wall->pdf(TString::Format("mggSig_cat%d",c)); 
    // import datasets from previous fit calculation
    wAllSave->import(*wall->pdf(TString::Format("mggSig_cat%d",c))); 
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAllSave->factory("CMS_hgg_sig_m0_absShift[1,1,1]"); 
  //
  wAllSave->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAllSave->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  // (3) Systematics on resolution
  // 1,1,1 statistical to be treated on the datacard
  wAllSave->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  //
  wAllSave->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAllSave->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  //
  wAllSave->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAllSave->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  // save the other parameters
  for (int c = 0; c < ncat; ++c) {
    wAllSave->factory(
	TString::Format("CMS_hgg_sig_alpha_cat%d[%g,0.5,5]", 
	c, wall->var(TString::Format("mgg_sig_alpha_cat%d",c))->getVal()));
  
    wAllSave->factory(
	TString::Format("CMS_hgg_sig_n_cat%d[%g,0.5,20]", 
	c, wall->var(TString::Format("mgg_sig_n_cat%d",c))->getVal()));
  
    wAllSave->factory(
	TString::Format("CMS_hgg_sig_frac_cat%d[%g,0.0,1.0]", 
	c, wall->var(TString::Format("mgg_sig_frac_cat%d",c))->getVal()));
  }
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) wAllSave->factory(
		  TString::Format("EDIT::CMS_hgg_sig_cat%d(mggSig_cat%d,",c,c) +
		  TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_alpha_cat%d=CMS_hgg_sig_alpha_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_n_cat%d=CMS_hgg_sig_n_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_frac_cat%d=CMS_hgg_sig_frac_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d)", c,c)
  );
  //
  // save to a file
  TString filename(TString(fileBasename)+".inputsig.root");
  wAllSave->writeToFile(filename);
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  cout << "Write signal workspace in: " << filename << " file" << endl;
  //wall->Print();
  //return;
  //
  // we finish fitting, 
  //////////////////////////////////////////////////////////////////////
  // we do plots
  std::vector<TString> catdesc;
  catdesc.push_back("2 btag");
  catdesc.push_back("1 btag");
  // retrieve data sets from the workspace
  // blinded dataset
  RooRealVar* mgg     = wall->var("mgg");  
  mgg->setUnit("GeV");
  Int_t nBinsMass(93); // just need to plot
  RooPlot* plotmtotAll = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
  //signalAll->plotOn(plotmtotAll);
  gStyle->SetOptTitle(0);
  // 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmtot[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmtot[c] = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    sigFit[c]  ->plotOn(plotmtot[c]);  
    mggSig[c]  ->plotOn(plotmtot[c]);
    mggSig[c]  ->plotOn(
	plotmtot[c],
	Components(TString::Format("GaussSig_cat%d",c)),
	LineStyle(kDashed),LineColor(kGreen));
    mggSig[c]  ->plotOn(
	plotmtot[c],
	Components(TString::Format("CBSig_cat%d",c)),
	LineStyle(kDashed),LineColor(kRed));

    mggSig[c]  ->paramOn(plotmtot[c]);
    sigFit[c]  ->plotOn(plotmtot[c]);
    TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
    plotmtot[c]->SetTitle("CMS preliminary 19.62/fb ");      
    plotmtot[c]->SetMinimum(0.0);
    plotmtot[c]->SetMaximum(1.40*plotmtot[c]->GetMaximum());
    plotmtot[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,501,501);
    plotmtot[c]->Draw();  
    plotmtot[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.52,0.65,0.95,0.9);

    legmc->AddEntry(plotmtot[c]->getObject(0),"Simulation","LPE");
    legmc->AddEntry(plotmtot[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmtot[c]->getObject(2),"Crystal Ball component","L");
    legmc->AddEntry(plotmtot[c]->getObject(3),"Gaussian Outliers","L");

    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    

    TLatex *lat  = new TLatex(
	minMassFit+1.5,0.85*plotmtot[c]->GetMaximum(),
	" WP4 300 GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(
	minMassFit+1.5,0.75*plotmtot[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
  } // close categories
  return;
} // end with signal
////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddBkgData(RooWorkspace* wAll, TString datafile, const char* fileBaseName, bool dobands) { 
  std::vector<TString> catdesc;
  catdesc.push_back("2 btag");
  catdesc.push_back("1 btag");
  const Int_t ncat = 2;
  // common preselection cut
  TString mainCut("1");
  // Variables
  RooArgSet* ntplVars = defineVariables();
  TFile dataFile(datafile); // haaaaa  
  TTree* dataTree     = (TTree*) dataFile.Get("TCVARS");
  RooDataSet data("data","dataset",dataTree,*ntplVars,""); 
  cout<<" data events "<< data.sumEntries()<< endl;

  RooDataSet* dataToFit[ncat];
  RooDataSet* dataToPlot[ncat];
  RooRealVar* mgg     = wAll->var("mgg");
  // cuts
  TString cut0 ="&& 1>0";// "&& mtot > 450 && mtot < 550 "; //"&& 1>0";//
  TString cut1 ="&& 1>0";// "&& mtot > 450 && mtot < 550 "; // "&& 1>0";//
  //
  TString cutj0 ="&& 1>0";// "&& mjj > 90 && mjj < 170 "; //"&& 1>0";//
  TString cutj1 ="&& 1>0";// "&& mjj > 100 && mjj < 160 "; // "&& 1>0";//
  //
  // apply the cuts
    dataToFit[0]   = (RooDataSet*) data.reduce(
	*wAll->var("mgg"),
	TString::Format(" cut_based_ct==%d",0)+cut0+cutj0);
    dataToPlot[0]   = (RooDataSet*) data.reduce(
	*wAll->var("mgg"),
	TString::Format(" cut_based_ct==%d",0)+cut0+cutj0
	+TString::Format(" && (mgg < 120 || mgg > 130)")  // blind
    );
    //
    dataToFit[1]   = (RooDataSet*) data.reduce(
	*wAll->var("mgg"),
	TString::Format(" cut_based_ct==%d",1)+cut1+cutj1);
    dataToPlot[1]   = (RooDataSet*) data.reduce(
	*wAll->var("mgg"),
	TString::Format(" cut_based_ct==%d",1)+cut1+cutj1
	+TString::Format(" && (mgg < 120 || mgg > 130)")  // blind
    );
  //
  //char nameToSave[300];
  for (int c = 0; c < ncat; ++c){
    //
  // convert it into a roodatahist - chiara: come capisco i bin?? 
  //RooArgSet setqui2(*mtot);
  //sprintf(nameToSave,"data_obs_cat%d",c);
  //RooDataHist dataBinned(nameToSave, nameToSave, setqui2, data);  
  //wAll->import(dataBinned);
    // unbinned dataset only
    // combine understand data as data_obs
  wAll->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
    //wAll->import(*dataToPlot[c],Rename(TString::Format("Dataplot_cat%d",c)));
  }
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = dataToFit[c]->sumEntries();
    cout << TString::Format("nEvt exp.  cat%d : ",c) << nExpEvt 
	 << TString::Format("   eff x Acc  cat%d : ",c) 
	 << "%" 
	 << endl; 
  }
  cout << "======================================================================" << endl;
  data.Print("v");
  for (int c = 0; c < ncat; ++c) {
    std::cout << "  " << wAll->data(TString::Format("Data_cat%d",c))->sumEntries();
  }
  // we adeed signal
  ///////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  // we do the model fitting on the same function
  //
  // fit range
  Float_t minMassFit(minfit),maxMassFit(maxfit); 
  RooPlot* plotmtotBkg[ncat];
  // dumb dataset
  RooDataSet* datafake[ncat];
  // declare pdf - same for both categories
  //
  Int_t nBinsMass(80);
  for (int c = 0; c < ncat; ++c) {
    // the parameters are set positive, because they are squared!!!!!!!!!!!!!!
    RooFormulaVar *p1mod = new RooFormulaVar(
	TString::Format("p1mod_cat%d",c),
	"","@0*@0",
	*wAll->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(
	TString::Format("p2mod_cat%d",c)
	,"","@0*@0",
	*wAll->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(
	TString::Format("p3mod_cat%d",c),
	"","@0*@0",
	*wAll->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
    ////////////////////////////////////////////////////////////////////
    RooAbsPdf* mggBkgTmp0 = 0; // declare a empty pdf
    // adding pdf's, using the variables
    mggBkgTmp0 = new  RooBernstein( // fill the pdf with the floating parameters
				   TString::Format("mggBkgTmp0_cat%d",c),
				   "", *mgg, 
				   RooArgList(RooConst(1.0),*p1mod, *p2mod, *p3mod)); 
    cout<<"here 6 "<< c<<endl;
    // the normalization variable
  //  wAll->factory(TString::Format("mtot_bkg_8TeV_norm_cat%d[1.0,0.0,100000]",c));
    // we copy the pdf to normalize
    wAll->factory(TString::Format("mgg_bkg_8TeV_cat%d_norm[1.0,0.0,100000]",c));
    RooExtendPdf mggBkgTmp( 
	TString::Format("mggBkg_cat%d",c),
	"",*mggBkgTmp0,
	*wAll->var(TString::Format("mgg_bkg_8TeV_cat%d_norm",c)));
   //
    fitresult[c] = mggBkgTmp.fitTo( // fit with normalized pdf,and return values
	*dataToFit[c], // bkg
	Strategy(1), // MINUIT strategy
	Minos(kFALSE), // interpretation on the errors, nonlinearities
	Range(minMassFit,maxMassFit),
	SumW2Error(kTRUE), 
	Save(kTRUE));
    wAll->import(mggBkgTmp); //store the normalized pdf on wp
    //
    // we plot
  TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,501,501);
  plotmtotBkg[c] = mgg->frame(nBinsMass);  
  //dataToPlot[c]   = (RooDataSet*) wAll->data(TString::Format("Data_cat%d",c));
  dataToPlot[c]->plotOn(plotmtotBkg[c],LineColor(kWhite),MarkerColor(kWhite));  
  mggBkgTmp.plotOn(
	plotmtotBkg[c],
	LineColor(kBlue),
	Range("fitrange"),NormRange("fitrange")); 
  dataToPlot[c]->plotOn(plotmtotBkg[c]); // blind 

  plotmtotBkg[c]->SetTitle("CMS preliminary 19.62/fb");      
  plotmtotBkg[c]->SetMinimum(0.0);
  plotmtotBkg[c]->SetMaximum(30.0*plotmtotBkg[c]->GetMaximum());
  plotmtotBkg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
  plotmtotBkg[c]->Draw();  
  cout << "!!!!!!!!!!!!!!!!!" << endl;

  if (dobands) {
    RooAbsPdf *cpdf; cpdf = mtotBkgTmp0;
    TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();   
    RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
    nlim->removeRange();   
    RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmtotBkg[c]->getObject(1));     
    for (int i=1; i<(plotmtotBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
      double lowedge = plotmtotBkg[c]->GetXaxis()->GetBinLowEdge(i);
      double upedge  = plotmtotBkg[c]->GetXaxis()->GetBinUpEdge(i);
      double center  = plotmtotBkg[c]->GetXaxis()->GetBinCenter(i);
      double nombkg = nomcurve->interpolate(center);
      nlim->setVal(nombkg);
      mgg->setRange("errRange",lowedge,upedge);
      RooAbsPdf *epdf = 0;
      epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
      RooAbsReal *nll = epdf->createNLL(*(dataToFit[c]),Extended());
      RooMinimizer minim(*nll);
      minim.setStrategy(0);
      double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
      double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
      minim.migrad();
      minim.minos(*nlim);   
      onesigma->SetPoint(i-1,center,nombkg);
      onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
      minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); 
      // the 0.5 is because qmu is -2*NLL
      // eventually if cl = 0.95 this is the usual 1.92!     
      minim.migrad();
      minim.minos(*nlim);  
      twosigma->SetPoint(i-1,center,nombkg);
      twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
      delete nll;
      delete epdf;
    }
    mgg->setRange("errRange",minMassFit,maxMassFit);
    twosigma->SetLineColor(kGreen);
    twosigma->SetFillColor(kGreen);
    twosigma->SetMarkerColor(kGreen);
    twosigma->Draw("L3 SAME");
    onesigma->SetLineColor(kYellow);
    onesigma->SetFillColor(kYellow);
    onesigma->SetMarkerColor(kYellow);
    onesigma->Draw("L3 SAME");     
    plotmtotBkg[c]->Draw("SAME");    
    } else {  
	cout << "!!!!!!!!!!!!!!!!!" << endl; // now we fit the gaussian on signal 
	plotmtotBkg[c]->Draw("SAME"); 
    }// close dobands
  plotmtotBkg[c]->GetYaxis()->SetRangeUser(0.0000001,20);
  //plotmtotBkg[c]->Draw("AC");
  ctmp->SetLogy(0);
  ctmp->SetGrid(1);
  cout << "!!!!!!!!!!!!!!!!!" << endl; 
  TLegend *legmc = new TLegend(0.52,0.65,0.92,0.9);
  legmc->AddEntry(plotmtotBkg[c]->getObject(2),"19.706 fb-1 ",""); //"LPE" blind
  legmc->AddEntry(plotmtotBkg[c]->getObject(1),"3* Berntein fit","L");
  if(dobands)legmc->AddEntry(twosigma,"two sigma ","F"); // not...
  if(dobands)legmc->AddEntry(onesigma,"one sigma","F");
  legmc->SetHeader("WP4 300 GeV");
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw();    
  TLatex *lat2 = new TLatex(minMassFit+1.5,0.75*plotmtotBkg[c]->GetMaximum(),catdesc.at(c));
  lat2->Draw();
  ctmp->SaveAs(TString::Format("databkgoversig_cat%d.pdf",c));
  cout<<"here 2 "<< c<<endl;
  ctmp->SaveAs(TString::Format("databkgoversig_cat%d.png",c));
  cout<<"here 3 "<< c<<endl;

} // close to each category
  // create WS before to save it
  RooAbsPdf* mggBkgPdfTosave[ncat];
  RooWorkspace *wAllSave = new RooWorkspace("w_all","w_all");
  RooDataSet* datatosave[ncat];
  for (int c = 0; c < ncat; ++c) {
    // import pdfs from previous fit calculation
    mggBkgPdfTosave[c] = (RooAbsPdf*)  wAll->pdf(TString::Format("mggBkg_cat%d",c)); 
    // import datasets from previous fit calculation
    wAllSave->import(*wAll->pdf(TString::Format("mggBkg_cat%d",c))); 
    // import observed
    //wAllSave->import(*wAll->pdf(TString::Format("data_obs_cat%d",c))); 
    datatosave[c] = (RooDataSet*) wAll->data(TString::Format("Data_cat%d",c));
    wAllSave->import(*datatosave[c], Rename(TString::Format("data_obs_cat%d",c)));
  }
  // (1) import everything functions / depends on quantity of parameters...
  for (int c = 0; c < ncat; ++c) {
    wAllSave->factory(
	TString::Format("CMS_hgg_bkg_8TeV_cat%d_norm[%g,0.0,100000]", 
	c, wAll->var(TString::Format("mgg_bkg_8TeV_cat%d_norm",c))->getVal()));
    wAllSave->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,-10,10]", 
	c, wAll->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c))->getVal()));
    wAllSave->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope2_cat%d[%g,-10,10]", 
	c, wAll->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c))->getVal()));
  }
  // (2) do reparametrization of background
  for (int c = 0; c < ncat; ++c){ 
    wAllSave->factory(
	TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(mggBkg_cat%d,",c,c) +
	TString::Format(" mgg_bkg_8TeV_cat%d_norm=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)+
	TString::Format(" mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d,", c,c)+
	TString::Format(" mgg_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d)", c,c)
  	);
  } // close for cat
  // import also observed
  //
  TString filename(TString(fileBaseName)+".root");
  wAllSave->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;
  // import data to workspace
  std::cout << "observation ";
  for (int c = 0; c < ncat; ++c) {
    std::cout << "  " << wAll->data(TString::Format("Data_cat%d",c))->sumEntries();
  }
  std::cout << std::endl;
  std::cout << "observation saved";
  for (int c = 0; c < ncat; ++c) {
    std::cout << "  " << wAllSave->data(TString::Format("data_obs_cat%d",c))->sumEntries();
  }
  std::cout << std::endl;
  // we finish fitting, 
  ///////////////////////////////////////////////////////////////////
  cout<<"here out loop"<<endl;
  return;
} // close add data .. 


//////////////////////////////////////////////////
void MakeDataCard(RooWorkspace* wall, RooWorkspace* wAll, const char* fileBaseName, const char* fileBkgName) {
  //TString cardDir = "datacards/";
  cout<< " signal " <<endl;
  wall->Print();
  cout<< " BKG " <<endl;
  wAll->Print();
  //
  const Int_t ncat = 2;
  RooDataSet* Data[ncat];
  RooDataSet* sigToFit[ncat];
  cout << "======== Start datacard maker =====================" << endl; 
  for (int c = 0; c < ncat; ++c) {
    Data[c]        = (RooDataSet*) wAll->data(TString::Format("Data_cat%d",c));
    sigToFit[c]    = (RooDataSet*) wall->data(TString::Format("Sig_cat%d",c));
  }

//  RooRealVar*  lumi = wAll->var("lumi"); 
  cout << "======== Expected signal Events Number =====================" << endl;  
  //Float_t siglikeErr[2];
  for (int c = 0; c < 2; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    //siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }

  cout << "======== Expected Events Number =====================" << endl;  
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c)<< Data[c]->sumEntries()  << endl;
  }
  cout << "====================================================" << endl;  

  TString filename(TString(fileBaseName)+".txt");
  ofstream outFile(filename);
  // name of files
  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run: combine -M Asymptotic hgg.mH300.0_8TeV.txt" << endl;
  //outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax "<< ncat << endl;
  outFile << "jmax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
  cout<<"here"<<endl;
  outFile << "# BKG" << endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs  cat1 "<<  TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "# the name after w_all is the name of the rooextpdf we want to use, we have both saved" << endl;
  outFile << "#shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:mggBkg_cat0" << endl;
  outFile << "#shapes mtotBkg   cat1 "<<  TString(fileBkgName)+".root" << " w_all:mggBkg_cat1" << endl;
  outFile << "shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_cat0" << endl;
  outFile << "shapes mtotBkg   cat1 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "#shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_allS_ws:mggSig_cat0" << endl;
  outFile << "#shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_allS_ws:mggSig_cat1" << endl;
  outFile << "shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_allS_ws:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_allS_ws:CMS_hgg_sig_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  // begin declaration
  outFile << "bin          cat0   cat1 " << endl;
  outFile <<  "observation   "  
	<<  Data[0]->sumEntries() << "  " 
	<<  Data[1]->sumEntries() << "  "
	<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0      cat1       cat1" << endl;
  outFile << "process                 mtotSig     mtotBkg     mtotSig    mtotBkg" << endl;
  outFile << "process                    0          1          0         1" << endl;
  outFile <<  "rate                      " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  Data[0]->sumEntries()  
	   << "  " << sigToFit[1]->sumEntries() << "  " <<  Data[1]->sumEntries()
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -     "
	<< "1.022   -     " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -   "
	<< "1.015        -   "
	<<"# JER and JES " << endl; // change latter
  outFile << "btag_eff          lnN " 
	<< "1.06        -  "
	<< "1.03        -  "
	<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g       lnN "
  	<< "1.010        -   "
  	<< "1.010        -   "
  	<< "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
	<< "1.01         -   "
	<< "1.01         -   "
	<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maajj_acceptance       lnN "
  	<< "1.10        -   "
  	<< "1.10        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "#mgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;

  outFile << "#mgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "#mgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "#mgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  /////////////////////////////////////

  /////////////////////////////////////

  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
} // close write datacard

//////////////////////////////////////////////////
void MakeDataCardoneCat(RooWorkspace* wall, RooWorkspace* wAll, const char* fileBaseName, const char* fileBkgName) {
  //TString cardDir = "datacards/";
  cout<< " signal " <<endl;
  wall->Print();
  cout<< " BKG " <<endl;
  wAll->Print();
  //
  const Int_t ncat = 2;
  RooDataSet* Data[ncat];
  RooDataSet* sigToFit[ncat];
  cout << "======== Start datacard maker =====================" << endl; 
  for (int c = 0; c < ncat; ++c) {
    Data[c]        = (RooDataSet*) wAll->data(TString::Format("Data_cat%d",c));
    sigToFit[c]    = (RooDataSet*) wall->data(TString::Format("Sig_cat%d",c));
  }

//  RooRealVar*  lumi = wAll->var("lumi"); 
  cout << "======== Expected signal Events Number =====================" << endl;  
  //Float_t siglikeErr[2];
  for (int c = 0; c < 2; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    //siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }

  cout << "======== Expected Events Number =====================" << endl;  
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c)<< Data[c]->sumEntries()  << endl;
  }
  cout << "====================================================" << endl;  

  TString filename(TString(fileBaseName)+"_onecat.txt");
  ofstream outFile(filename);
  // name of files
  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run: combine -M Asymptotic hgg.mH300.0_8TeV.txt" << endl;
  //outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax 1 "<< endl;
  outFile << "jmax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
  cout<<"here"<<endl;
  outFile << "# BKG" << endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  //outFile << "shapes data_obs  cat1 "<<  TString(fileBkgName)+".root" << " w_allB_ws:data_obs_cat1" << endl;
  outFile << "#shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:mggBkg_cat0" << endl;
   outFile << "shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_cat0" << endl;
  //outFile << "shapes mtotBkg   cat1 "<<  TString(fileBkgName)+".root" << " w_allB_ws:mggBkg_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_allS_ws:mggSig_cat0" << endl;
  //outFile << "shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_allS_ws:mggSig_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  // begin declaration
  outFile << "bin          cat0    " << endl;
  outFile <<  "observation   "  
	<<  Data[0]->sumEntries() << "  " 
	//<<  Data[1]->sumEntries() << "  "
	<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0     " << endl;
  outFile << "process                 mtotSig     mtotBkg     " << endl;
  outFile << "process                    0          1        " << endl;
  outFile <<  "rate                      " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  Data[0]->sumEntries()  
	   //<< "  " << sigToFit[1]->sumEntries() << "  " <<  Data[1]->sumEntries()
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -     "
	//<< "1.022   -     " 
        << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -   "
	//<< "1.015        -   "
	<<"# JER and JES " << endl; // change latter
  outFile << "btag_eff          lnN " 
	<< "1.06        -  "
	//<< "1.03        -  "
	<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g       lnN "
  	<< "1.010        -   "
  	//<< "1.010        -   "
  	<< "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
	<< "1.01         -   "
	//<< "1.01         -   "
	<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maajj_acceptance       lnN "
  	<< "1.10        -   "
  	//<< "1.10        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "mgg_bkg_8TeV_norm_cat0           flatParam  # Normalization uncertainty on background slope" << endl;
  //outFile << "mgg_bkg_8TeV_norm_cat1           flatParam  # Normalization uncertainty on background slope" << endl;

  outFile << "mgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "mgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "mgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "mgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "mgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "mgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  /////////////////////////////////////

  /////////////////////////////////////

  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
} // close write datacard


void SetConstantParams(const RooArgSet* params) { 
  /*
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }
  */
  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }    
} // close set const parameters
//////////////////////////////////////////////////////////////////////////////////////
void style(){
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000); 
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0); 
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0);  // For the axis titles:

    defaultStyle->SetTitleColor(1, "XYZ");
    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.06, "XYZ");
 
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    defaultStyle->SetTitleXOffset(0.9);
    defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
    defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(510, "XYZ");
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
  return;
}
