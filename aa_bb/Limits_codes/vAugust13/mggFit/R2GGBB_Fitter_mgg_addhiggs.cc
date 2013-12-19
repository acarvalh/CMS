/** \macro H2GGFitter.cc
 *
 * $Id$
 *
 * Software developed for the CMS Detector at LHC
 *
 *
 *  Template Serguei Ganjour - CEA/IRFU/SPP, Saclay
 *  
 *
 * Macro is implementing the unbinned maximum-likelihood model for 
 * the Higgs to gamma gamma analysis. PDF model and RooDataSets 
 * are stored in the workspace which is feeded to  HiggsAnalysis/CombinedLimit tools:
 * 
 */
// this one is for mgg fit               

using namespace RooFit;
using namespace RooStats ;

const Int_t NCAT = 2;
bool addHiggs=true;

void AddSigData(RooWorkspace*, Float_t);
void AddHigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*);
void SigModelFit(RooWorkspace*, Float_t);
void MakePlots(RooWorkspace*, Float_t, RooFitResult* );
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeHigWS(RooWorkspace* w, const char* filename);
void MakeBkgWS(RooWorkspace* w, const char* filename);
void MakeDataCard(RooWorkspace* w, const char* filename,  const char* filename1,  const char* filename2);
void MakeDataCardREP(RooWorkspace* w, const char* filename,  const char* filename1);
void MakeDataCardonecat(RooWorkspace* w, const char* filename,  const char* filename1,  const char* filename2);
void MakeDataCardLNU(RooWorkspace* w, const char* filename,  const char* filename1);
void MakeDataCardonecatnohiggs(RooWorkspace* w, const char* filename,  const char* filename1,  const char* filename2);
//  MakeDataCardonecat(w, fileBaseName, fileBkgName);
//  MakeDataCardREP(w, fileBaseName, fileBkgName);
//void MakeDataCardnohiggs(RooWorkspace* w, const char* filename,  const char* filename1,  const char* filename2);
//void MakeDataCardonecat(RooWorkspace* w, const char* filename,  const char* filename1,  const char* filename2);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);

RooFitResult* fitresult[NCAT]; // container for the fit results
RooFitResult*  BkgModelFitBernstein(RooWorkspace*, Bool_t);

RooArgSet* defineVariables()
{
  RooRealVar* mgg  = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  //RooRealVar* mtot  = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  //RooRealVar* mjj  = new RooRealVar("mjj","M(jj)",100,1600,"GeV");
  RooRealVar* evWeight   = new RooRealVar("evWeight","HqT x PUwei",0,100,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  //
  RooArgSet* ntplVars = new RooArgSet(*mgg, * cut_based_ct, *evWeight);
  ntplVars->add(*mgg);
  //ntplVars->add(*mtot);
  //ntplVars->add(*mjj);
  ntplVars->add(*cut_based_ct);
  return ntplVars;
}

void runfits(const Float_t mass=120, Int_t mode=1, Bool_t dobands = false)
{
  style();
  TString fileBaseName(TString::Format("hgg.mH%.1f_8TeV", mass));
  TString fileHiggsName(TString::Format("hgg.hig.mH%.1f_8TeV", mass));
  TString fileBkgName(TString::Format("hgg.inputbkg_8TeV", mass));
  TString card_name("models_test.rs"); // put the model parameters here!
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  RooFitResult* fitresults;
  bool cutbased=true;
  // the minitree to be addeed
  //
  TString hhiggs = "MiniTrees/ChiaraNov13/v20/finalizedTrees_Radion_V07__fitToGG__noKinFit/HToGG125__selez300.root";
  TString ssignal = "MiniTrees/ChiaraNov13/v20/finalizedTrees_Radion_V07__fitToGG__noKinFit/RadionSignal_m300.root";
  TString ddata   = "MiniTrees/ChiaraNov13/v20/finalizedTrees_Radion_V07__fitToGG__noKinFit/Data__selez300.root";
  //
//  TString hhiggs = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Radion_m300_8TeV_nm_m300.root";
//  TString ssignal = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Radion_m300_8TeV_nm_m300.root";
//  TString ddata   = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Data_m300.root";
   //
//  TString hhiggs = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Radion_m300_8TeV_nm_m300.root";
//  TString ssignal = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Radion_m300_8TeV_nm_m300.root";
//  TString ddata   = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Data_m300.root";

  //
//  cout<<"Higgs: "<<hhiggs<<endl;
  cout<<"Signal: "<<ssignal<<endl;
  cout<<"Data: "<<ddata<<endl;
  AddSigData(w, mass,ssignal);
  cout<<"SIGNAL ADDED"<<endl;
  AddHigData(w, mass,hhiggs);
  cout<<"HIGGS ADDED"<<endl;
  AddBkgData(w,ddata);
  cout<<"BKG ADDED"<<endl;
  w->Print("v");
  // construct the models to fit
  SigModelFit(w, mass); // constructing signal pdf
  HigModelFit(w, mass); // constructing higgs pdf

  bool dobands=false;
  cout<<" did dignal WS's"<<endl;
  fitresults = BkgModelFitBernstein(w, dobands);  // this is berestein 3

  // Construct points workspace
  MakeSigWS(w, fileBaseName);
  MakeBkgWS(w, fileBkgName);
  MakeHigWS(w, fileHiggsName);
  MakePlots(w, mass, fitresults);
  MakeDataCardonecat(w, fileBaseName, fileBkgName, fileHiggsName);
  MakeDataCardREP(w, fileBaseName, fileBkgName);
  MakeDataCardLNU(w, fileBaseName, fileBkgName);
  MakeDataCard(w, fileBaseName, fileBkgName, fileHiggsName);//MakeDataCardnohiggs
  MakeDataCardonecatnohiggs(w, fileBaseName, fileBkgName, fileHiggsName);//MakeDataCardnohiggs
 // MakeDataCardonecat(w, fileBaseName, fileBkgName, fileHiggsName);//MakeDataCardnohiggs
  cout<< "here"<<endl;
  return;
} // close runfits
////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddSigData(RooWorkspace* w, Float_t mass, TString signalfile) { 
  cout << "================= Add Signal==============================" << endl;
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  // Luminosity:
  Float_t Lum = 19785.0; // pb-1
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  RooArgSet* ntplVars = defineVariables();
  TFile sigFile(signalfile);  
  TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
  // common preselection cut
  TString mainCut("1"); 
  RooDataSet sigScaled(
	"sigScaled",
	"dataset",
	sigTree, 
	*ntplVars,
	mainCut,
	"evWeight");
  cout << "======================================================================" << endl;
  RooDataSet* sigToFit[ncat]; 
  TString cut0 = " && 1>0";// "&& mtot > 955 && mtot < 1150 "; //
  TString cut1 = " && 1>0";//"&& mtot > 955 && mtot < 1150 "; // "&& 1>0";//
  //
  TString cutj0 = " && 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 170 "; //"&& 1>0";//
  TString cutj1 = " && 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 160 "; // "&& 1>0";//
  //
  // we take only mtot to fit to the workspace, we include the cuts
  sigToFit[0] = (RooDataSet*) sigScaled.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*sigToFit[0],Rename(TString::Format("Sig_cat%d",0)));
    //
  sigToFit[1] = (RooDataSet*) sigScaled.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*sigToFit[1],Rename(TString::Format("Sig_cat%d",1)));
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll  = (RooDataSet*) sigScaled.reduce(*w->var("mgg"),mainCut);
  cout << "======================================================================" << endl;
  w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel:  " << sigScaled.sumEntries() << endl; 
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = sigToFitAll[c]->sumEntries();
    cout << TString::Format("nEvt exp.  cat%d : ",c) << nExpEvt 
	 << TString::Format("   eff x Acc  cat%d : ",c) 
	 << "%" 
	 << endl; 
  } // close ncat
  cout << "======================================================================" << endl;
  sigScaled.Print("v");
  return;
} // end add signal function
///////////////////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddBkgData(RooWorkspace* w, TString datafile) { 
  const Int_t ncat = NCAT;
  // common preselection cut
  TString mainCut("1");
  RooArgSet* ntplVars = defineVariables();
  RooRealVar weightVar("weightVar","",1,0,1000);
  weightVar.setVal(1.);
  // no common preselection cut applied yet; 
  TFile dataFile(datafile); 
  TTree* dataTree     = (TTree*) dataFile.Get("TCVARS");
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","weightVar"); 
  // evweight is 1 anyway...
  RooDataSet* dataToFit[ncat];
  RooDataSet* dataToPlot[ncat];
  TString cut0 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; //"&& 1>0";//
  TString cut1 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; //"&& 1>0";//
  //
  TString cutj0 = "&& 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 170 "; //"&& 1>0";//
  TString cutj1 = "&& 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 160 "; // "&& 1>0";//
  //
  cout<<" HERE TAKE DATASET"<<endl;

	dataToFit[0]   = (RooDataSet*) Data.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d",0)+cut0+cutj0);
    dataToPlot[0]   = (RooDataSet*) Data.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d",0)
	+TString::Format(" && (mgg > 130 || mgg < 120)")+cut0+cutj0); // blind
   
	dataToFit[1]   = (RooDataSet*) Data.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d",1)+cut1);
    dataToPlot[1]   = (RooDataSet*) Data.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d",1)
	+TString::Format(" && (mgg > 130 || mgg < 120)")+cut1); // blind

  for (int c = 0; c < ncat; ++c) {   
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
    w->import(*dataToPlot[c],Rename(TString::Format("Dataplot_cat%d",c)));
  }
  // Create full data set without categorization
  RooDataSet* data    = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut);
  w->import(*data, Rename("Data"));
  data->Print("v");
  return;
} // close add data .. 
/////////////////////////////////////////////////////////////////////
// we make the fit model
void SigModelFit(RooWorkspace* w, Float_t mass) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  //******************************************//
  // Fit signal with model pdfs
  //******************************************//
  // four categories to fit
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggSig[ncat];
  // fit range
  Float_t minMassFit(100),maxMassFit(180); 
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    sigToFit[c]   = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("mggSig_cat%d",c));
      ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(MASS);
    //RooRealVar* peak = w->var(TString::Format("mgg_sig_m0_cat%d",c));
    //peak->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)
    mggSig[c]->fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE));
    // IMPORTANT: fix all pdf parameters to constant, why?
    // parameters for signal model, different rages, why?
    w->defineSet(TString::Format("SigPdfParam_cat%d",c), 
        RooArgSet(
	 *w->var(TString::Format("mgg_sig_m0_cat%d",c)),   
	 *w->var(TString::Format("mgg_sig_sigma_cat%d",c)),
	 *w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
	 *w->var(TString::Format("mgg_sig_n_cat%d",c)), 
	 *w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
	 *w->var(TString::Format("mgg_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
  } // close for ncat
} // close signal model fit
/////////////////////////////////////////
// we make the higgs model
void HigModelFit(RooWorkspace* w, Float_t mass) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  // four categories to fit
  RooDataSet* higToFit[ncat];
  RooAbsPdf* mggHig[ncat];
  // fit range
  Float_t minMassFit(100),maxMassFit(180); 
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    higToFit[c]   = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c));
    mggHig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("mggHig_cat%d",c));
    //RooRealVar* peak = w->var(TString::Format("mgg_hig_m0_cat%d",c));
    //peak->setVal(MASS);
    ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_cat%d",c)))->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)
    mggHig[c]->fitTo(*higToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE));
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("HigPdfParam_cat%d",c), 
        RooArgSet(
	 *w->var(TString::Format("mgg_hig_m0_cat%d",c)),   
	 *w->var(TString::Format("mgg_hig_sigma_cat%d",c)),
	 *w->var(TString::Format("mgg_hig_alpha_cat%d",c)),
	 *w->var(TString::Format("mgg_hig_n_cat%d",c)), 
	 *w->var(TString::Format("mgg_hig_gsigma_cat%d",c)),
	 *w->var(TString::Format("mgg_hig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("HigPdfParam_cat%d",c)));
  } // close for ncat
} // close higgs model fit
////////////////////////////////////////////////////////////
// BKG model berestein 3
RooFitResult* BkgModelFitBernstein(RooWorkspace* w, Bool_t dobands) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("2 btag");
  catdesc.push_back("1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  //******************************************//
  // Fit background with model pdfs
  //******************************************//
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[ncat];
  RooDataSet* dataplot[ncat]; // the data
  RooBernstein* mggBkg[ncat];// the polinomial of 4* order
   RooGaussian* Higgs[ncat]; // the higgs to sum
  RooPlot* plotmggBkg[ncat];
  RooDataSet* sigToFit[ncat]; 
  RooAbsPdf*  mggSig[ncat]; 
  Float_t minMassFit(100),maxMassFit(180); 
  // Fit data with background pdf for data limit
  RooRealVar* mgg     = w->var("mgg");
  mgg->setUnit("GeV");
  //
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //
  for (int c = 0; c < ncat; ++c) { // to each category  
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    cout << "!!!!!!!!!!!!!" << endl;
    ////////////////////////////////////
    // these are the parameters for the bkg polinomial
    // one slope by category - float from -10 > 10
    // the parameters are squared
    RooFormulaVar *p1mod = new RooFormulaVar(
	TString::Format("p1mod_cat%d",c),
	"","@0*@0",
	*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
    if(c==1){RooFormulaVar *p2mod = new RooFormulaVar(
	TString::Format("p2mod_cat%d",c)
	,"","@0*@0",
	*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(
	TString::Format("p3mod_cat%d",c)
	,"","@0*@0",
	*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
    /*
    RooFormulaVar *p4mod = new RooFormulaVar(
	TString::Format("p4mod_cat%d",c)
	,"","@0*@0",
	*w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",c)));
    RooFormulaVar *p5mod = new RooFormulaVar(
	TString::Format("p5mod_cat%d",c)
	,"","@0*@0",
	*w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",c)));
    */
    }
    ////////////////////////////////////////////////////////////////////
    RooAbsPdf* mggBkgTmp0 = 0; // declare a empty pdf
    // adding pdf's, using the variables
    if(c==0){
    mggBkgTmp0 = new  RooBernstein( // fill the pdf with the floating parameters
				   TString::Format("mggBkgTmp0_cat%d",c),
				   "", *mgg, 
				   RooArgList(RooConst(1.0),*p1mod)); 
   }
    if(c==1){
    mggBkgTmp0 = new  RooBernstein( // fill the pdf with the floating parameters
				   TString::Format("mggBkgTmp0_cat%d",c),
				   "", *mgg, 
				   RooArgList(RooConst(1.0),*p1mod, *p2mod, *p3mod));//, *p4mod, *p5mod)); 
   }
   // we first wrap the normalization of mggBkgTmp0
    w->factory(TString::Format("mgg_bkg_8TeV_norm_cat%d[1.0,0.0,100000]",c));
    RooExtendPdf mggBkgTmp( // we copy the pdf? normalized
	TString::Format("mggBkg_cat%d",c),
	"",*mggBkgTmp0,
	*w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)) // normalization only on full bkg
    );
    fitresult[c] = mggBkgTmp.fitTo( // fit with normalized pdf,and return values
	*data[c], // bkg
	Strategy(1), // MINUIT strategy
	Minos(kFALSE), // interpretation on the errors, nonlinearities
	Range(minMassFit,maxMassFit),
	SumW2Error(kTRUE), 
	Save(kTRUE));
    w->import(mggBkgTmp); 
   //************************************************//
   // Plot mgg background fit results per categories 
   //************************************************//
   TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,500,500);
   Int_t nBinsMass(80);
   plotmggBkg[c] = mgg->frame(nBinsMass);
    cout<<" here 1"<<endl;
   dataplot[c]   = (RooDataSet*) w->data(TString::Format("Dataplot_cat%d",c));
   cout<<" here 1"<<endl;
   data[c]->plotOn(plotmggBkg[c],LineColor(kWhite),MarkerColor(kWhite));  
   mggBkgTmp.plotOn(
	plotmggBkg[c],
	LineColor(kBlue),
	Range("fitrange"),NormRange("fitrange")); 
    dataplot[c]->plotOn(plotmggBkg[c]); 
    plotmggBkg[c]->Draw();  
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!" << endl; // now we fit the gaussian on signal
    // plot signal also
    /*
    sigToFit[c]       = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    double norm     = 1.0*sigToFit[c]->sumEntries(); // 
    mggSig[c]       = (RooAbsPdf*)  w->pdf(TString::Format("mggSig_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig[c]       ->plotOn(
	plotmggBkg[c],
	Normalization(norm,RooAbsPdf::NumEvent),
	DrawOption("F"),
	LineColor(kRed),FillStyle(1001),FillColor(19));
    mggSig[c]->plotOn(
	plotmggBkg[c],
	Normalization(norm,RooAbsPdf::NumEvent),LineColor(kRed)); 
    */
    plotmggBkg[c]->SetTitle("CMS preliminary 19.702/fb");      
    plotmggBkg[c]->SetMinimum(0.0);
    plotmggBkg[c]->SetMaximum(1.40*plotmggBkg[c]->GetMaximum());
    plotmggBkg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    //double test = sigToFit[c]->sumEntries();
    //cout<<"number of events on dataset "<<test<<endl;
  if (dobands) {
      RooAbsPdf *cpdf; cpdf = mggBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
      for (int i=1; i<(plotmggBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
        double lowedge = plotmggBkg[c]->GetXaxis()->GetBinLowEdge(i);
        double upedge  = plotmggBkg[c]->GetXaxis()->GetBinUpEdge(i);
        double center  = plotmggBkg[c]->GetXaxis()->GetBinCenter(i); 
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        mgg->setRange("errRange",lowedge,upedge);
        RooAbsPdf *epdf = 0;
        epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange"); 
        RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
        RooMinimizer minim(*nll);
        minim.setStrategy(0);
        double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
        double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        minim.migrad();
        minim.minos(*nlim);
        // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
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
         } // close for bin
     mgg->setRange("errRange",minMassFit,maxMassFit);     
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");    
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      plotmggBkg[c]->Draw("SAME");
    } else plotmggBkg[c]->Draw("SAME"); // close dobands
    cout << "!!!!!!!!!!!!!!!!!" << endl; 
    TLegend *legmc = new TLegend(0.60,0.72,0.92,0.9);
    legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","LPE"); // not...
    legmc->AddEntry(plotmggBkg[c]->getObject(1),"Fit","L");
    if(dobands)legmc->AddEntry(twosigma,"two sigma ","F"); // not...
    if(dobands)legmc->AddEntry(onesigma,"one sigma","F");
    legmc->SetHeader(" 300 GeV");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    TLatex *lat2 = new TLatex(minMassFit+1.5,0.75*plotmggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    //
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d.png",c));
  // ctmp->SaveAs(TString::Format("databkgoversig_cat%d.C",c));
  } // close to each category
  RooBernstein mggBkgAll("mggBkgAll", "", *mgg, 
	RooArgList(RooConst(1.0), 
	 *w->var("mgg_bkg_8TeV_slope1"), 
	 *w->var("mgg_bkg_8TeV_slope2"), 
	 *w->var("mgg_bkg_8TeV_slope3")));
  w->import(mggBkgAll);
  RooFitResult* fitresults;
  fitresults =  w->pdf("mggBkgAll")->fitTo( // save results to workspace
	*w->data("Data"), 
	Range(minMassFit,maxMassFit),
	SumW2Error(kTRUE), Save(kTRUE));
  fitresults->Print();
  return fitresults;
} // close berestein 3
///////////////////////////////////////////////////////////////
void MakeSigWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir   = "workspaces/";
  const Int_t ncat = NCAT;
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save
  // for statistical tests. 
  //**********************************************************************//
  RooAbsPdf* mggSigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    mggSigPdf[c] = (RooAbsPdf*)  w->pdf(TString::Format("mggSig_cat%d",c)); 
    wAll->import(*w->pdf(TString::Format("mggSig_cat%d",c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]"); 
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  //wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");

  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  // save the other parameters
/*  for (int c = 0; c < ncat; ++c) {
    wAll->factory(
	TString::Format("CMS_hgg_sig_alpha_cat%d[%g,0.5,5]", 
	c, wAll->var(TString::Format("mgg_sig_alpha_cat%d",c))->getVal()));
  
    wAll->factory(
	TString::Format("CMS_hgg_sig_n_cat%d[%g,0.5,20]", 
	c, wAll->var(TString::Format("mgg_sig_n_cat%d",c))->getVal()));
  
    wAll->factory(
	TString::Format("CMS_hgg_sig_frac_cat%d[%g,0.0,1.0]", 
	c, wAll->var(TString::Format("mgg_sig_frac_cat%d",c))->getVal()));
  }
*/
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) wAll->factory(
		  TString::Format("EDIT::CMS_hgg_sig_cat%d(mggSig_cat%d,",c,c) +
		  TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d, ", c,c) +
//		  TString::Format(" mgg_sig_alpha_cat%d=CMS_hgg_sig_alpha_cat%d, ", c,c) +
//		  TString::Format(" mgg_sig_n_cat%d=CMS_hgg_sig_n_cat%d, ", c,c) +
//		  TString::Format(" mgg_sig_frac_cat%d=CMS_hgg_sig_frac_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d)", c,c)
  );
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make signal WP
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir   = "workspaces/";
  const Int_t ncat = NCAT;  
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests. 
  //**********************************************************************//
  RooDataSet* data[ncat];
  RooAbsPdf* mggBkgPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    data[c]      = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    mggBkgPdf[c] = (RooAbsPdf*)  w->pdf(TString::Format("mggBkg_cat%d",c));
    wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));
    wAll->import(*w->pdf(TString::Format("mggBkg_cat%d",c)));
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_cat%d_norm[%g,0.0,100000.0]", 
	c, w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c))->getVal()));
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,-10,10]", 
	c, w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c))->getVal()));
    if(c==1){ wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope2_cat%d[%g,-10,10]", 
	c, w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c))->getVal()));
  //
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope3_cat%d[%g,-10,10]", c, 
        w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c))->getVal()));
    
/*
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope4_cat%d[%g,-10,10]", 
	c, w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",c))->getVal()));
  //
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope5_cat%d[%g,-10,10]", c, 
        w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",c))->getVal()));
*/
    }
  } // close ncat
  // (2) do reparametrization of background
  for (int c = 0; c < ncat; ++c){ 
	if(c==0) wAll->factory(
	TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(mggBkg_cat%d,",c,c) +
	TString::Format(" mgg_bkg_8TeV_norm_cat%d=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)+
	TString::Format(" mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d)", c,c)
	//TString::Format(" mgg_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d,", c,c)+
	//TString::Format(" mgg_bkg_8TeV_slope3_cat%d=CMS_hgg_bkg_8TeV_slope3_cat%d)", c,c)
  	);
	if(c==1) wAll->factory(
	TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(mggBkg_cat%d,",c,c) +
	TString::Format(" mgg_bkg_8TeV_norm_cat%d=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)+
	TString::Format(" mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d,", c,c)+
	TString::Format(" mgg_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d,", c,c)+
	TString::Format(" mgg_bkg_8TeV_slope3_cat%d=CMS_hgg_bkg_8TeV_slope3_cat%d)", c,c)
//	TString::Format(" mgg_bkg_8TeV_slope4_cat%d=CMS_hgg_bkg_8TeV_slope4_cat%d,", c,c)+
//	TString::Format(" mgg_bkg_8TeV_slope5_cat%d=CMS_hgg_bkg_8TeV_slope5_cat%d)", c,c)
  	);
 } // close for cat

  TString filename(wsDir+TString(fileBaseName)+".root");
  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;
  std::cout << "observation ";
  for (int c = 0; c < ncat; ++c) {
    std::cout << "  " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries();
  }
  std::cout << std::endl;
  return;
} // close make BKG workspace
////////////////////////////////////////////////////////////////////
void SetConstantParams(const RooArgSet* params) { 
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
} // close set const parameters
////////////////////////////////////////////////////////////////////////
void MakePlots(RooWorkspace* w, Float_t Mass, RooFitResult* fitresults) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("   2 btag");
  catdesc.push_back("   1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  // retrieve data sets from the workspace
//  RooDataSet* dataAll         = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll       = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll       = (RooDataSet*) w->data("Hig");
  // blinded dataset
//  RooDataSet* data[ncat];  
  RooDataSet* sigToFit[ncat];
  RooAbsPdf*  mggGaussSig[ncat];
  RooAbsPdf*  mggCBSig[ncat];
  RooAbsPdf*  mggSig[ncat];
  //
  RooAbsPdf*  mggBkg[ncat];  
  for (int c = 0; c < ncat; ++c) {
  //  data[c]         = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]     = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c)); 
    mggGaussSig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("mggGaussSig_cat%d",c));
    mggCBSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("mggCBSig_cat%d",c));
    mggSig[c]       = (RooAbsPdf*)  w->pdf(TString::Format("mggSig_cat%d",c));
    mggBkg[c]       = (RooAbsPdf*)  w->pdf(TString::Format("mggBkg_cat%d",c));
  } // close categories
  RooRealVar* mgg     = w->var("mgg");  
  mgg->setUnit("GeV");
  RooAbsPdf* mggGaussSigAll  = w->pdf("mggGaussSig");
  RooAbsPdf* mggCBSigAll     = w->pdf("mggCBSig");
  RooAbsPdf* mggSigAll       = w->pdf("mggSig");
  //RooAbsPdf* mggBkgAll       = w->pdf("mggBkg_cat1");
  //
  //****************************//
  // Plot mgg Fit results
  //****************************//
  // Set P.D.F. parameter names
  // WARNING: Do not use it if Workspaces are created
  //  SetParamNames(w);
  Float_t minMassFit(100),maxMassFit(180); 
  Float_t MASS(Mass);  
  Int_t nBinsMass(50); // just need to plot
  RooPlot* plotmggAll = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
  signalAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas("c1","mgg",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories 
  //****************************//
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmgg[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmgg[c] = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmgg[c],LineColor(kWhite),MarkerColor(kWhite));    
    mggSig[c]  ->plotOn(plotmgg[c]);
    mggSig[c]  ->plotOn(
	plotmgg[c],
	Components(TString::Format("GaussSig_cat%d",c)),
	LineStyle(kDashed),LineColor(kGreen));
    mggSig[c]  ->plotOn(
	plotmgg[c],
	Components(TString::Format("CBSig_cat%d",c)),
	LineStyle(kDashed),LineColor(kRed));
    mggSig[c]  ->paramOn(plotmgg[c]);
    sigToFit[c]  ->plotOn(plotmgg[c]);
//    TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
    plotmgg[c]->SetTitle("CMS preliminary 19.702/fb ");      
    plotmgg[c]->SetMinimum(0.0);
    plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
    plotmgg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,500,500);
    plotmgg[c]->Draw();  
    plotmgg[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.95,0.9);
    legmc->AddEntry(plotmgg[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmgg[c]->getObject(2),"Crystal Ball component","L");
    legmc->AddEntry(plotmgg[c]->getObject(3),"Gaussian Outliers","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    //    float effS = effSigma(hist);
    TLatex *lat  = new TLatex(
	minMassFit+1.5,0.85*plotmgg[c]->GetMaximum(),
	" Resonance - 300 GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(
	minMassFit+1.5,0.75*plotmgg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
  } // close categories
    return;
} // close makeplots signal
////////////////////////////////////////////////////////////////////////
void MakePlotsHiggs(RooWorkspace* w, Float_t Mass, RooFitResult* fitresults) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("   2 btag");
  catdesc.push_back("   1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  //RooDataSet* higToFit[ncat];
  RooAbsPdf*  mggGaussHig[ncat];
  RooAbsPdf*  mggCBHig[ncat];
  RooAbsPdf*  mggHig[ncat];
  //
  for (int c = 0; c < ncat; ++c) {
    higToFit[c]     = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c)); 
    mggGaussHig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("mggGaussHig_cat%d",c));
    mggCBHig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("mggCBHig_cat%d",c));
    mggHig[c]       = (RooAbsPdf*)  w->pdf(TString::Format("mggHig_cat%d",c));
  } // close categories
  RooRealVar* mgg     = w->var("mgg");  
  mgg->setUnit("GeV");
  RooAbsPdf* mggGaussHigAll  = w->pdf("mggGaussHig");
  RooAbsPdf* mggCBHigAll     = w->pdf("mggCBHig");
  RooAbsPdf* mggHigAll       = w->pdf("mggHig");
  //
  Float_t minMassFit(100),maxMassFit(180); 
  Float_t MASS(Mass);  
  Int_t nBinsMass(50); // just need to plot
  RooPlot* plotmggAll = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
  signalAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  higgsAll->plotOn(plotmggAll);
  mggHigAll->plotOn(plotmggAll);
  mggHigAll->plotOn(
	plotmggAll,Components("mggGaussHig"),
	LineStyle(kDashed),LineColor(kGreen));
  mggHigAll->plotOn(
	plotmggAll,Components("mggCBHig"),
	LineStyle(kDashed),LineColor(kRed));
  mggHigAll->paramOn(
	plotmggAll, 
	ShowConstants(true), 
	Layout(0.15,0.55,0.9), 
	Format("NEU",AutoPrecision(2)));
  //
  plotmggAll->getAttText()->SetTextSize(0.03);
  TCanvas* c1 = new TCanvas("c1","mgg",0,0,500,500);
  c1->cd(1);
  //****************************//
  // Plot higgs Background 
  //****************************//
  RooPlot* plotmggh[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmggh[c] = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    higToFit[c]->plotOn(plotmggh[c],LineColor(kWhite),MarkerColor(kWhite));    
    mggHig[c]  ->plotOn(plotmggh[c]);
    mggHig[c]  ->plotOn(
	plotmggh[c],
	Components(TString::Format("GaussHig_cat%d",c)),
	LineStyle(kDashed),LineColor(kGreen));
    mggHig[c]  ->plotOn(
	plotmggh[c],
	Components(TString::Format("CBHig_cat%d",c)),
	LineStyle(kDashed),LineColor(kRed));
    mggHig[c]  ->paramOn(plotmggh[c]);
    higToFit[c]  ->plotOn(plotmggh[c]);
    //TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    //TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
    plotmggh[c]->SetTitle("CMS preliminary 19.702/fb ");      
    plotmggh[c]->SetMinimum(0.0);
    plotmggh[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
    plotmggh[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    TCanvas* ctmp2 = new TCanvas("ctmp2","Background Categories",0,0,500,500);
    plotmggh[c]->Draw();  
    plotmggh[c]->Draw("SAME");  
    legmc->Draw();    
    //    float effS = effSigma(hist);
    TLatex *lath  = new TLatex(
	minMassFit+1.5,0.85*plotmggh[c]->GetMaximum(),
	"Higgs - 300 GeV");
    lath->Draw();
    TLatex *lat2h = new TLatex(
	minMassFit+1.5,0.75*plotmggh[c]->GetMaximum(),catdesc.at(c));
    lat2h->Draw();
    ctmp2->SaveAs(TString::Format("higmodel_cat%d.pdf",c));
    ctmp2->SaveAs(TString::Format("higmodel_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
  } // close categories

    return;
} // close makeplots signal
////////////////////////////////////////////////////////////////////
// we add the higgs to the workspace in categories
void AddHigData(RooWorkspace* w, Float_t mass, TString signalfile) { 
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  RooArgSet* ntplVars = defineVariables();
  TFile higFile(signalfile);  
  TTree* higTree = (TTree*) higFile.Get("TCVARS");
  //higTree->AddVar(scaleWeightVar);
  // common preselection cut
  TString mainCut("1"); 
  // one channel with right weights
  //RooDataSet higScaled1(
  RooDataSet higScaled(
	"higScaled1",
	"dataset",
	higTree, // all variables of RooArgList
	*ntplVars,
	mainCut,
	"evWeight");
  /*
  RooRealVar *evWeight = (RooRealVar*) (*ntplVars)["evWeight"] ;
  RooRealVar *k = new RooRealVar("k", "k", 0.0006424);
  RooFormulaVar *nw = new RooFormulaVar("nw", "nw", "@1", RooArgSet(*evWeight, *k));
  higScaled1.addColumn(*nw);
  RooArgSet *ntplVars1 = higScaled1.get();  
  RooDataSet *higScaled = new RooDataSet("higScaled", "dataset",higTree, *ntplVars1,"", "nw");
  */
  //
  RooDataSet* higToFit[ncat]; 
  TString cut0 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; // 
  TString cut1 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; // "&& 1>0";//
  //
  TString cutj0 = "&& 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 160 "; //"&& 1>0";//
  TString cutj1 = "&& 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 170 "; // "&& 1>0";//
  //
  // we take only mtot to fit to the workspace, we include the cuts
  higToFit[0] = (RooDataSet*) higScaled.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*higToFit[0],Rename(TString::Format("Hig_cat%d",0)));
    //
  higToFit[1] = (RooDataSet*) higScaled.reduce(
	*w->var("mgg"),
	mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*higToFit[1],Rename(TString::Format("Hig_cat%d",1)));  // Create full signal data set without categorization
  RooDataSet* higToFitAll  = (RooDataSet*) higScaled->reduce(*w->var("mgg"),mainCut);
  w->import(*higToFitAll,Rename("Hig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel:  " << higScaled->sumEntries() << endl; 
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = higToFit[c]->sumEntries();
    cout << TString::Format("nEvt exp.  cat%d : ",c) << nExpEvt 
	 << TString::Format("   eff x Acc  cat%d : ",c) 
	 << "%" 
	 << endl; 
  }
  cout << "======================================================================" << endl;
  higScaled.Print("v");
  return;
} // end add higgs function
///////////////////////////////////////////////////////////////
void MakeHigWS(RooWorkspace* w, const char* fileHiggsName) {
  TString wsDir   = "workspaces/";
  const Int_t ncat = NCAT;
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests. 
  //**********************************************************************//
  RooAbsPdf* mggHigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    mggHigPdf[c] = (RooAbsPdf*)  w->pdf(TString::Format("mggHig_cat%d",c)); 
    wAll->import(*w->pdf(TString::Format("mggHig_cat%d",c))); 
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_hig_m0_absShift[1,1,1]"); 
  wAll->factory("prod::CMS_hgg_hig_m0_cat0(mgg_hig_m0_cat0, CMS_hgg_hig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_hig_m0_cat1(mgg_hig_m0_cat1, CMS_hgg_hig_m0_absShift)");
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_hig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_hig_sigma_cat0(mgg_hig_sigma_cat0, CMS_hgg_hig_sigmaScale)");

  wAll->factory("prod::CMS_hgg_hig_sigma_cat1(mgg_hig_sigma_cat1, CMS_hgg_hig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_hig_gsigma_cat0(mgg_hig_gsigma_cat0, CMS_hgg_hig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_hig_gsigma_cat1(mgg_hig_gsigma_cat1, CMS_hgg_hig_sigmaScale)");
  // save the other parameters
/*  for (int c = 0; c < ncat; ++c) {
    wAll->factory(
	TString::Format("CMS_hgg_sig_alpha_cat%d[%g,0.5,5]", 
	c, wAll->var(TString::Format("mgg_sig_alpha_cat%d",c))->getVal()));
  
    wAll->factory(
	TString::Format("CMS_hgg_sig_n_cat%d[%g,0.5,20]", 
	c, wAll->var(TString::Format("mgg_sig_n_cat%d",c))->getVal()));
  
    wAll->factory(
	TString::Format("CMS_hgg_sig_frac_cat%d[%g,0.0,1.0]", 
	c, wAll->var(TString::Format("mgg_sig_frac_cat%d",c))->getVal()));
  }
*/
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) wAll->factory(
		  TString::Format("EDIT::CMS_hgg_hig_cat%d(mggHig_cat%d,",c,c) +
		  TString::Format(" mgg_hig_m0_cat%d=CMS_hgg_hig_m0_cat%d, ", c,c) +
		  TString::Format(" mgg_hig_sigma_cat%d=CMS_hgg_hig_sigma_cat%d, ", c,c) +
//		  TString::Format(" mgg_sig_alpha_cat%d=CMS_hgg_sig_alpha_cat%d, ", c,c) +
//		  TString::Format(" mgg_sig_n_cat%d=CMS_hgg_sig_n_cat%d, ", c,c) +
//		  TString::Format(" mgg_sig_frac_cat%d=CMS_hgg_sig_frac_cat%d, ", c,c) +
		  TString::Format(" mgg_hig_gsigma_cat%d=CMS_hgg_hig_gsigma_cat%d)", c,c)
  );
  TString filename(wsDir+TString(fileHiggsName)+".inputhig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make higgs WP
///////////////////////////////////////////////////////////
// declare histos or what -> NOT USED
Double_t effSigma(TH1 *hist) {
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();
  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...
  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {jbm++; xj+=bwid; bin=hist->GetBinContent(jbm); total+=bin; if(total > rlim) break;} else ierr=1;
      if(kbm > 0) {kbm--; xk-=bwid; bin=hist->GetBinContent(kbm); total+=bin; if(total > rlim) break; } else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) { widmin=wid; ismin=iscan; }
  } // Scan window centre
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;
  return widmin;
} // close effSigma
//////////////////////////////////////////////////
// with higgs
void MakeDataCard(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , const char* fileHiggsName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooDataSet* higToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]      = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    higToFit[c]      = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c));
  }
  RooRealVar*  lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << "19785" << " pb-1 ............................" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
//  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << ".........Expected Signal for L = " << "19785" << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("Data")->sumEntries()  << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;  
  TString filename(cardDir+TString(fileBaseName)+".txt");
  ofstream outFile(filename);

  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH300.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
//  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi =  " << "19785" << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 2" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

cout<<"here"<<endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs  cat1 "<<  TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mggBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "shapes mggBkg   cat1 "<<  TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mggSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggHig cat0 " << TString(fileHiggsName)+".inputhig.root" << " w_all:CMS_hgg_hig_cat0" << endl;
  outFile << "shapes mggHig cat1 " << TString(fileHiggsName)+".inputhig.root" << " w_all:CMS_hgg_hig_cat1" << endl;

  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(addHiggs) { //
  outFile << "bin          cat0   cat1 " << endl;
  outFile <<  "observation   "  
	<<  data[0]->sumEntries() << "  " 
	<<  data[1]->sumEntries() << "  "     << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0       cat0      cat1       cat1       cat1" << endl;
  outFile << "process                 mggSig     mggBkg     mggHig     mggSig    mggBkg     mggHig" << endl;
  outFile << "process                    0          1          2          0         1          2" << endl;
  outFile <<  "rate                     " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  1 << "  " << higToFit[0]->sumEntries() << "  " 
	   << "  " << sigToFit[1]->sumEntries() << "  " <<  1 << "  " << higToFit[1]->sumEntries()
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -    1.022 "
	<< "1.022   -    1.022 " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -  1.015 "
	<< "1.015        -  1.015 "
	<<"# JER and JES " << endl;
  outFile << "btag_eff          lnN " 
	<< "1.06        -  1.06 "
	<< "1.03        -  1.03 "
	<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g       lnN "
  	<< "1.010        -   1.010  "
  	<< "1.010        -   1.010 "
  	<< "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
	<< "1.01         -   1.010 "
	<< "1.01         -   1.010 "
	<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maa_acceptance       lnN "
  	<< "1.10        -   1.10 "
  	<< "1.10        -   1.10 "
  	<< "# photon energy resolution" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on higgs " << endl;
  outFile << "CMS_hgg_hig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_hig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mgg fit - slopes" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "#CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "#CMS_hgg_bkg_8TeV_slope3_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope3_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "#CMS_hgg_bkg_8TeV_slope4_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  outFile << "#CMS_hgg_bkg_8TeV_slope5_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  } // if ncat ==2
  /////////////////////////////////////

  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write full datacard
//////////////////////////////////////////////////
///// datacards witout higgs
//////////////////////////////////////////////////
// with reparametrization of BKG
void MakeDataCardREP(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]      = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
  }
  RooRealVar*  lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("Data")->sumEntries()  << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;  
  TString filename(cardDir+TString(fileBaseName)+"rep.txt");
  ofstream outFile(filename);
  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH300.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
cout<<"here"<<endl;
  outFile << "# the name after w_all is the name of the rooextpdf we want to use, we have both saved" << endl;
  outFile << "# BKG" << endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs  cat1 "<<  TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "shapes mtotBkg   cat1 "<<  TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  /////////////////////////////////////
  outFile << "bin          cat0   cat1 " << endl;
  outFile <<  "observation   "  
	<<  data[0]->sumEntries() << "  " 
	<<  data[1]->sumEntries() << "  "
	<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0      cat1       cat1" << endl;
  outFile << "process                 mtotSig     mtotBkg     mtotSig    mtotBkg" << endl;
  outFile << "process                    0          1          0         1" << endl;
  outFile <<  "rate                      " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  1  
	   << "  " << sigToFit[1]->sumEntries() << "  " <<  1
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -     "
	<< "1.022   -     " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -   "
	<< "1.015        -   "
	<<"# JER and JES " << endl;
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
  outFile << "maa_acceptance       lnN "
  	<< "1.10        -   "
  	<< "1.10        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## with reparametrization" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#CMS_hgg_bkg_8TeV_slope3_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope3_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#CMS_hgg_bkg_8TeV_slope4_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#CMS_hgg_bkg_8TeV_slope5_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard with rep
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////
// withou reparametrization of BKG
void MakeDataCardLNU(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]      = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
  }
  RooRealVar*  lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("Data")->sumEntries()  << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;  
  TString filename(cardDir+TString(fileBaseName)+"lnu.txt");
  ofstream outFile(filename);
  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH300.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
cout<<"here"<<endl;
  outFile << "# the name after w_all is the name of the rooextpdf we want to use, we have both saved" << endl;
  outFile << "# BKG" << endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs  cat1 "<<  TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:mggBkg_cat0" << endl;
  outFile << "shapes mtotBkg   cat1 "<<  TString(fileBkgName)+".root" << " w_all:mggBkg_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  /////////////////////////////////////
  outFile << "bin          cat0   cat1 " << endl;
  outFile <<  "observation   "  
	<<  data[0]->sumEntries() << "  " 
	<<  data[1]->sumEntries() << "  "
	<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0      cat1       cat1" << endl;
  outFile << "process                 mtotSig     mtotBkg     mtotSig    mtotBkg" << endl;
  outFile << "process                    0          1          0         1" << endl;
  outFile <<  "rate                      " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  data[0]->sumEntries()  
	   << "  " << sigToFit[1]->sumEntries() << "  " <<  data[1]->sumEntries()
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -     "
	<< "1.022   -     " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -   "
	<< "1.015        -   "
	<<"# JER and JES " << endl;
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
  outFile << "maa_acceptance       lnN "
  	<< "1.10        -   "
  	<< "1.10        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "mggBkg       lnU "
  	<< " -        2   "
  	<< " -        2   "
  	<< "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## with reparametrization" << endl;
  outFile << "mgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope3_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mgg_bkg_8TeV_slope3_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope4_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "#mgg_bkg_8TeV_slope5_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard without rep
////////////////////////////////////////////////////////////////////////////////
void MakeDataCardonecat(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , const char* fileHiggsName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooDataSet* higToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]      = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    higToFit[c]      = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c));
  }
  RooRealVar*  lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << "19785" << " pb-1 ............................" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
//  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << ".........Expected Signal for L = " << "19785" << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("Data")->sumEntries()  << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;  
  TString filename(cardDir+TString(fileBaseName)+"onecat.txt");
  ofstream outFile(filename);

  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH300.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
//  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi =  " << "19785" << " pb-1" << endl;
  outFile << "imax 1" << endl;
  outFile << "jmax 2" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

cout<<"here"<<endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mggBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggHig cat0 " << TString(fileHiggsName)+".inputhig.root" << " w_all:CMS_hgg_hig_cat0" << endl;

  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(addHiggs) { //
  outFile << "bin          cat0  " << endl;
  outFile <<  "observation   "  
	<<  data[0]->sumEntries() << "  "   << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0       cat0   " << endl;
  outFile << "process                 mggSig     mggBkg     mggHig  " << endl;
  outFile << "process                    0          1          2    " << endl;
  outFile <<  "rate                     " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  1 << "  " << higToFit[0]->sumEntries() 
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -    1.022 " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -  1.015 "
	<<"# JER and JES " << endl;
  outFile << "btag_eff          lnN " 
	<< "1.06        -  1.06 "
	<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g       lnN "
  	<< "1.010        -   1.010  "
  	<< "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
	<< "1.01         -   1.010 "
	<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maa_acceptance       lnN "
  	<< "1.10        -   1.10 "
  	<< "# photon energy resolution" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on higgs " << endl;
  outFile << "CMS_hgg_hig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_hig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mgg fit - slopes" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;


  outFile << "#CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;


  outFile << "#CMS_hgg_bkg_8TeV_slope3_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  } // if ncat ==2
  /////////////////////////////////////

  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard one cat
////////////////////////////////////////////////////////////////////////////////
void MakeDataCardonecatnohiggs(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , const char* fileHiggsName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooDataSet* higToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]      = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    higToFit[c]      = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c));
  }
  RooRealVar*  lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << "19785" << " pb-1 ............................" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
//  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << ".........Expected Signal for L = " << "19785" << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("Data")->sumEntries()  << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;  
  TString filename(cardDir+TString(fileBaseName)+"onecatnohiggs.txt");
  ofstream outFile(filename);

  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH300.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
//  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi =  " << "19785" << " pb-1" << endl;
  outFile << "imax 1" << endl;
  outFile << "jmax 1" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

cout<<"here"<<endl;
  outFile << "shapes data_obs  cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mggBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggHig cat0 " << TString(fileHiggsName)+".inputhig.root" << " w_all:CMS_hgg_hig_cat0" << endl;

  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(addHiggs) { //
  outFile << "bin          cat0  " << endl;
  outFile <<  "observation   "  
	<<  data[0]->sumEntries() << "  "   << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0      " << endl;
  outFile << "process                 mggSig     mggBkg    " << endl;
  outFile << "process                    0          1       " << endl;
  outFile <<  "rate                     " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  1  
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -    " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015   -   "
	<<"# JER and JES " << endl;
  outFile << "btag_eff          lnN " 
	<< "1.06        -  "
	<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g       lnN "
  	<< "1.010        -     "
  	<< "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
	<< "1.01         -    "
	<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maa_acceptance       lnN "
  	<< "1.10        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on higgs " << endl;
  outFile << "CMS_hgg_hig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_hig_sigmaScale     param   1   0.30   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mgg fit - slopes" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;


  outFile << "#CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;


  outFile << "#CMS_hgg_bkg_8TeV_slope3_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;

  } // if ncat ==2
  /////////////////////////////////////

  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard one cat



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

void SetParamNames(RooWorkspace* w) { // not used it if Workspaces are created => float fit
  const Int_t ncat = NCAT;
  //****************************//
  // mgg signal all categories
  //****************************//
  RooRealVar* mgg_sig_m0     = w->var("mgg_sig_m0");  
  RooRealVar* mgg_sig_sigma  = w->var("mgg_sig_sigma");
  RooRealVar* mgg_sig_alpha  = w->var("mgg_sig_alpha"); 
  RooRealVar* mgg_sig_n      = w->var("mgg_sig_n"); 
  RooRealVar* mgg_sig_gsigma = w->var("mgg_sig_gsigma");
  RooRealVar* mgg_sig_frac   = w->var("mgg_sig_frac");
  mgg_sig_m0    ->SetName("m_{0}");
  mgg_sig_sigma ->SetName("#sigma_{CB}");
  mgg_sig_alpha ->SetName("#alpha");
  mgg_sig_n     ->SetName("n");
  mgg_sig_gsigma->SetName("#sigma_G");  
  mgg_sig_frac  ->SetName("f_G");  
  mgg_sig_m0    ->setUnit("GeV");
  mgg_sig_sigma ->setUnit("GeV");
  mgg_sig_gsigma->setUnit("GeV"); 
  //****************************//
  // mgg background  
  //****************************//
  RooRealVar* mgg_bkg_8TeV_slope1  = w->var("mgg_bkg_8TeV_slope1");
  mgg_bkg_8TeV_slope1              ->SetName("a_{B}");
  mgg_bkg_8TeV_slope1              ->setUnit("1/GeV");
  RooRealVar* mgg_bkg_8TeV_slope2  = w->var("mgg_bkg_8TeV_slope2");
  mgg_bkg_8TeV_slope2              ->SetName("a_{B}");
  mgg_bkg_8TeV_slope2              ->setUnit("1/GeV");
  //****************************//
  // mgg per category  
  //****************************//
  for (int c = 0; c < ncat; ++c) {
    mgg_sig_m0     = (RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c));
    mgg_sig_sigma  = (RooRealVar*) w->var(TString::Format("mgg_sig_sigma_cat%d",c));
    mgg_sig_alpha  = (RooRealVar*) w->var(TString::Format("mgg_sig_alpha_cat%d",c));
    mgg_sig_n      = (RooRealVar*) w->var(TString::Format("mgg_sig_n_cat%d",c));
    mgg_sig_gsigma = (RooRealVar*) w->var(TString::Format("mgg_sig_gsigma_cat%d",c));
    mgg_sig_frac   = (RooRealVar*) w->var(TString::Format("mgg_sig_frac_cat%d",c));
    mgg_sig_m0     ->SetName("m_{0}"); 
    mgg_sig_sigma  ->SetName("#sigma_{CB}");
    mgg_sig_alpha  ->SetName("#alpha");
    mgg_sig_n      ->SetName("n"); 
    mgg_sig_gsigma ->SetName("#sigma_{G}");
    mgg_sig_frac   ->SetName("f_{G}");
    mgg_sig_m0     ->setUnit("GeV");
    mgg_sig_sigma  ->setUnit("GeV");
    mgg_sig_gsigma ->setUnit("GeV");
    mgg_bkg_8TeV_slope1 = w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c));
    mgg_bkg_8TeV_slope1 ->SetName("p_{B}^{1}");
    mgg_bkg_8TeV_slope1 ->setUnit("1/GeV");
    mgg_bkg_8TeV_slope2 = w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c));
    mgg_bkg_8TeV_slope2 ->SetName("p_{B}^{2}");
    mgg_bkg_8TeV_slope2 ->setUnit("1/GeV^{2}");
//    RooRealVar* mgg_bkg_8TeV_frac = w->var(TString::Format("mgg_bkg_8TeV_frac_cat%d",c));
//    mgg_bkg_8TeV_frac ->SetName("f");
  }
} // close setparameters
