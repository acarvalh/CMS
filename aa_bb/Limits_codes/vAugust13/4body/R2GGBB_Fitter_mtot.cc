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
  // this one is for 4 body fit               

#include "GaussExp.h" 
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

//using namespace RooFit;
//using namespace RooStats ;

const Int_t NCAT = 2;

// declare the functions
void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*);
void SigModelFit(RooWorkspace*, Float_t);
void MakePlots(RooWorkspace*, Float_t, RooFitResult* );
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeBkgWS(RooWorkspace* w, const char* filename);
void MakeDataCardonecat(RooWorkspace* w, const char* filename,  const char* filename1);
void MakeDataCardREP(RooWorkspace* w, const char* filename,  const char* filename1);
void MakeDataCardLnU(RooWorkspace* w, const char* filename,  const char* filename1);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);

RooFitResult* fitresult[NCAT]; // container for the fit results
RooFitResult*  BkgModelFitBernstein(RooWorkspace*, Bool_t);

const int minfit =150, maxfit=1200;

RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mtot  = new RooRealVar("mtot","M(#gamma#gamma jj)",150,1200,"GeV");
  RooRealVar* mgg  = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  RooRealVar* evWeight   = new RooRealVar("evWeight","HqT x PUwei",0,100,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 2") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  //
  RooArgSet* ntplVars = new RooArgSet(*mtot, *mgg, * cut_based_ct, *evWeight);
  ntplVars->add(*cut_based_ct);
  ntplVars->add(*mtot);
  ntplVars->add(*mgg);
  ntplVars->add(*evWeight);
  return ntplVars;
}

void runfits(const Float_t mass=120, Int_t mode=1, Bool_t dobands = false)
{
  style();
  TString fileBaseName(TString::Format("hgg.mH%.1f_8TeV", mass));
  TString fileBkgName(TString::Format("hgg.inputbkg_8TeV", mass));
  TString card_name("models_mtot_exp.rs"); // fit model parameters
  // declare a  first WS
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs(); // Get models and variables
  RooFitResult* fitresults;
//  TString ssignal = "MiniTrees/OlivierAug13/v02_regkin_mggjj_0/Radion_m500_regression-m500_minimal.root";
//  TString ddata   = "MiniTrees/OlivierAug13/v02_regkin_mggjj_0/Data_regression-m500_minimal.root";
//  TString ssignal = "MiniTrees/OlivierOc13/v15_base_mggjj_0/02013-10-30-Radion_m550_8TeV_nm_m550.root";
//  TString ddata   = "MiniTrees/OlivierOc13/v15_base_mggjj_0/02013-10-30-Data_m550.root";
  //
  TString ssignal = "MiniTrees/ChiaraNov13/v20/finalizedTrees_Radion_V07__fitToGGJJ__withKinFit/RadionSignal_m550.root";
  TString ddata   = "MiniTrees/ChiaraNov13/v20/finalizedTrees_Radion_V07__fitToGGJJ__withKinFit/Data.root";
  //
  cout<<"Signal: "<< ssignal<<endl;
  cout<<"Data: "<< ddata<<endl;
  AddSigData(w, mass,ssignal);
  AddBkgData(w,ddata);
  w->Print("v");
  // construct the models to fit
  SigModelFit(w, mass); 
  bool dobands=false;
  fitresults = BkgModelFitBernstein(w, dobands);
  // Construct points workspace
  MakeSigWS(w, fileBaseName);
  MakeBkgWS(w, fileBkgName);
  MakePlots(w, mass, fitresults);

  MakeDataCardonecat(w, fileBaseName, fileBkgName);
  MakeDataCardREP(w, fileBaseName, fileBkgName);
  MakeDataCardLnU(w, fileBaseName, fileBkgName);
  cout<< "here"<<endl;
  return;
} // close runfits
////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddSigData(RooWorkspace* w, Float_t mass, TString signalfile) { 
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  // Luminosity:
  Float_t Lum = 19620.0; // pb-1
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  RooArgSet* ntplVars = defineVariables();
  // file input
  TFile sigFile(signalfile);  
  TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
  // common preselection cut
  TString mainCut("1"); 
  // one channel with right weights
  RooDataSet sigScaled(
	"sigScaled",
	"dataset",
	sigTree, 
	*ntplVars,
	mainCut,
	"evWeight");
  RooDataSet* sigToFit[ncat]; 
  TString cut0 = "&& mgg > 120 && mgg < 130 "; // " && 1>0";//
  for (int c = 0; c < ncat; ++c) {
    sigToFit[c] = (RooDataSet*) sigScaled.reduce(
	*w->var("mtot"),
	mainCut+TString::Format(" && cut_based_ct==%d",c)+cut0);
    w->import(*sigToFit[c],Rename(TString::Format("Sig_cat%d",c)));
  } // close ncat
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll  = (RooDataSet*) sigScaled.reduce(*w->var("mtot"),mainCut);
  w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
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
  return;
} // end add signal function
///////////////////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddBkgData(RooWorkspace* w, TString datafile) { 
  const Int_t ncat = NCAT;
  // common preselection cut
  TString mainCut("1");
  RooArgSet* ntplVars = defineVariables();
  // weith data as with 1
  RooRealVar weightVar("weightVar","",1,0,1000);
  weightVar.setVal(1.);
  TFile dataFile(datafile);
  TTree* dataTree     = (TTree*) dataFile.Get("TCVARS");
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","weightVar"); 

  RooDataSet* dataToFit[ncat];
  RooDataSet* dataToPlot[ncat];
  TString cut0 = "&& mgg > 120 && mgg < 130 "; // " && 1>0";//
  for (int c = 0; c < ncat; ++c) {
    dataToFit[c]   = (RooDataSet*) Data.reduce(
	*w->var("mtot"),
	TString::Format(" cut_based_ct==%d",c)+cut0);
    dataToPlot[c]   = (RooDataSet*) Data.reduce(
	*w->var("mtot"),
	TString::Format(" cut_based_ct==%d",c) 
	+TString::Format(" && (mtot > 2050)") + cut0 // blind
    );
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
    w->import(*dataToPlot[c],Rename(TString::Format("Dataplot_cat%d",c)));
  }
  // Create full data set without categorization
  RooDataSet* data    = (RooDataSet*) Data.reduce(*w->var("mtot"));
  w->import(*data, Rename("Data"));
  data->Print("v");
  return;
} // close add data .. 
/////////////////////////////////////////////////////////////////////
// we make the fit model
void SigModelFit(RooWorkspace* w, Float_t mass) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);

  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mtotSig[ncat];
  // fit range
  Float_t minMassFit(minfit),maxMassFit(maxfit); 
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    sigToFit[c]   = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mtotSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("mtotSig_cat%d",c));
    RooRealVar* peak = w->var(TString::Format("mtot_sig_m0_cat%d",c));
    peak->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)
    mtotSig[c]->fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE));
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("SigPdfParam_cat%d",c), 
        RooArgSet(
	 *w->var(TString::Format("mtot_sig_m0_cat%d",c)),   
	 *w->var(TString::Format("mtot_sig_alpha_cat%d",c)),
	 *w->var(TString::Format("mtot_sig_n_cat%d",c)), 
	 *w->var(TString::Format("mtot_sig_gsigma_cat%d",c)),
	 *w->var(TString::Format("mtot_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
  } // close for ncat
} // close signal model fit
////////////////////////////////////////////////////////////
// BKG model berestein 3
RooFitResult* BkgModelFitBernstein(RooWorkspace* w, Bool_t dobands) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("2 btag");
  catdesc.push_back("1 btag");
  //******************************************//
  // Fit background with model pdfs
  //******************************************//
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[ncat];
  RooDataSet* dataplot[ncat]; 
  //RooBernstein* mtotBkg[ncat];
  RooPlot* plotmtotBkg[ncat];
  RooDataSet* sigToFit[ncat]; // the signal
  RooAbsPdf*  mtotSig[ncat]; // polinomial of 4* order for signal
  Float_t minMassFit(minfit),maxMassFit(maxfit); 
  // Fit data with background pdf for data limit
  RooRealVar* mtot     = w->var("mtot");
  mtot->setUnit("GeV");
  //
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  // 
  int order=-2;
  //
  for (int c = 0; c < ncat; ++c) { 
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    cout << "!!!!!!!!!!!!!" << endl;
    /*
    w->factory(TString::Format("mtot_bkg_8TeV_norm_cat%d[1.0,0.0,100000]",c)); // is is on the .rs
    */
    ////////////////////////////////////
    // these are the parameters for exponential fit
    // one slope by category -  not squared slopes
    //((RooRealVar*) w->var(TString::Format("mtot_bkg_8TeV_slope3_cat%d",c)))->setConstant(true);
    //cout << "---------------- Parameter 3 set to const" << endl;
  cout<<"here 0 "<< c<<endl;
  RooFormulaVar *p1mod = new RooFormulaVar( // parameter 0 is the mean of the Gaussian core,
	TString::Format("p1mod_cat%d",c),"","@0",*w->var(TString::Format("mtot_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p2mod = new RooFormulaVar( // parameter 1 is the std dev of the Gaussian
	TString::Format("p2mod_cat%d",c),"","@0",*w->var(TString::Format("mtot_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p3mod = new RooFormulaVar( // parameter 2 is the number of std dev from the mean where the exponential fall off starts on the high side of the Gaussian. Parameter number 2, incidentally, is identical to theexponent of the exponential, by requiring analytical continuity at the joining point
	TString::Format("p3mod_cat%d",c),"","@0",*w->var(TString::Format("mtot_bkg_8TeV_slope3_cat%d",c)));

    RooAbsPdf* mtotBkgTmp0 = 0; // declare a empty pdf
    // adding pdf's, using the variables
    mtotBkgTmp0 = new  GaussExp( // fill the pdf with the floating parameters
				   TString::Format("mggBkgTmp0_cat%d",c),
				   "", *mtot, 
				   *p1mod, *p2mod, *p3mod); 

/*
  RooAbsPdf* mtotBkgTmp0 = new RooGenericPdf(  // if exp function
		TString::Format("DijetBackground_%d",c), 
		"exp(-@0/(@1*@1+@2*@2*@0))", 
		RooArgList(*mtot, *p1mod, *p2mod));
*/
  w->factory(TString::Format("mtot_bkg_8TeV_cat%d_norm[1.0,0.0,100000]",c));

  RooExtendPdf mtotBkgTmp( // we copy the pdf? normalized
	TString::Format("mtotBkg_cat%d",c),
	"",*mtotBkgTmp0,
	*w->var(TString::Format("mtot_bkg_8TeV_cat%d_norm",c)));
  fitresult[c] = mtotBkgTmp.fitTo( // fit with normalized pdf,and return values
	*data[c], // bkg
	Strategy(1), // MINUIT strategy
	Minos(kFALSE), // interpretation on the errors, nonlinearities
	Range(minMassFit,maxMassFit),
	SumW2Error(kTRUE), 
	Save(kTRUE));
  w->import(mtotBkgTmp); //store the normalized pdf on wp
   //************************************************//
   // Plot mtot background fit results per categories 
   //************************************************//
   TCanvas* ctmp = new TCanvas("ctmp","mtot Background Categories",0,0,501,501);
   Int_t nBinsMass(80);
   plotmtotBkg[c] = mtot->frame(nBinsMass);
   
   dataplot[c]   = (RooDataSet*) w->data(TString::Format("Dataplot_cat%d",c));
   data[c]->plotOn(plotmtotBkg[c],LineColor(kWhite),MarkerColor(kWhite));  
   mtotBkgTmp.plotOn(
	plotmtotBkg[c],
	LineColor(kBlue),
	Range("fitrange"),NormRange("fitrange")); 
//    dataplot[c]->plotOn(plotmtotBkg[c]); // blind 
    plotmtotBkg[c]->Draw();  
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!" << endl; // now we fit the gaussian on signal
    // we draw signal on the same
    /*
    sigToFit[c]       = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    double norm     = 0.01*sigToFit[c]->sumEntries(); // 
    mtotSig[c]       = (RooAbsPdf*)  w->pdf(TString::Format("mtotSig_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mtotSig[c]       ->plotOn(
	plotmtotBkg[c],
	Normalization(norm,RooAbsPdf::NumEvent),
	DrawOption("F"),
	LineColor(kRed),FillStyle(1001),FillColor(19));
    mtotSig[c]->plotOn(
	plotmtotBkg[c],
	Normalization(norm,RooAbsPdf::NumEvent),LineColor(kRed)); 
    */ 
    plotmtotBkg[c]->SetTitle("CMS preliminary 19.702/fb");      
    plotmtotBkg[c]->SetMinimum(0.0);
    plotmtotBkg[c]->SetMaximum(20*plotmtotBkg[c]->GetMaximum());
    plotmtotBkg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma jj} (GeV)");
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
        mtot->setRange("errRange",lowedge,upedge);
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
        minim.migrad();
        minim.minos(*nlim); 
        twosigma->SetPoint(i-1,center,nombkg);
        twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
       delete nll;
       delete epdf;
      } // close for each bin
      mtot->setRange("errRange",minMassFit,maxMassFit);     
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      plotmtotBkg[c]->Draw("SAME");
    } else plotmtotBkg[c]->Draw("SAME"); // close dobands
   plotmtotBkg[c]->GetYaxis()->SetRangeUser(0.0000001,100); // this one works
   //plotmtotBkg[c]->Draw("AC");
   ctmp->SetLogy(1);
   ctmp->SetGrid(1);
   cout << "!!!!!!!!!!!!!!!!!" << endl; 

    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
//    legmc->AddEntry(plotmtotBkg[c]->getObject(2),"Data ",""); //"LPE" blind
    legmc->AddEntry(plotmtotBkg[c]->getObject(1),"Exponential fit","L");
    if(dobands)legmc->AddEntry(twosigma,"two sigma ","F"); 
    if(dobands)legmc->AddEntry(onesigma,"one sigma","F");
    legmc->SetHeader("WP4 550 GeV");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    TLatex *lat2 = new TLatex(363.0,0.85*plotmtotBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d.pdf",c));
  cout<<"here 2 "<< c<<endl;
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d.png",c));
  cout<<"here 3 "<< c<<endl;
    //ctmp->SaveAs(TString::Format("databkgoversig_cat%d.C",c));
  } // close to each category
  cout<<"here out loop"<<endl;

  return fitresult;
} // close berestein 3
///////////////////////////////////////////////////////////////
void MakeSigWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir   = "workspaces/";
  const Int_t ncat = NCAT;
  //**********************************************************************//
  // Write pdfs and datasets into the workspace 
  // for statistical tests. 
  //**********************************************************************//
  RooAbsPdf* mtotSigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all"); 
  // import pdfs from previous fit calculation before to save
  for (int c = 0; c < ncat; ++c) {
    mtotSigPdf[c] = (RooAbsPdf*)  w->pdf(TString::Format("mtotSig_cat%d",c)); 
    wAll->import(*w->pdf(TString::Format("mtotSig_cat%d",c))); // import datasets from previous fit calculation
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]"); 
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat1[1,1,1]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat2[1,1,1]"); //CHANGE
  //
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mtot_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mtot_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  // (3) Systematics on resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  //
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mtot_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mtot_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  //
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mtot_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mtot_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) wAll->factory(
		  TString::Format("EDIT::CMS_hgg_sig_cat%d(mtotSig_cat%d,",c,c) +
		  TString::Format(" mtot_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d, ", c,c) +
		  TString::Format(" mtot_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d, ", c,c) +
		  TString::Format(" mtot_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d)", c,c)
  );
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make signal WP
////////////////////////////////////////////////////////////////
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir   = "workspaces/";
  const Int_t ncat = NCAT;  
  //**********************************************************************//
  // Write pdfs and datasets into the workspace 
  // for statistical tests. 
  //**********************************************************************//
  RooDataSet* data[ncat];
  RooAbsPdf* mtotBkgPdf[ncat];
  // (1) import everything functions before to save to a file
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    data[c]      = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    mtotBkgPdf[c] = (RooAbsPdf*)  w->pdf(TString::Format("mtotBkg_cat%d",c));
    wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));
    wAll->import(*w->pdf(TString::Format("mtotBkg_cat%d",c)));
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_cat%d_norm[%g,0.0,100000.0]", 
	c, wAll->var(TString::Format("mtot_bkg_8TeV_cat%d_norm",c))->getVal())); 
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,-500,500]", 
	c, wAll->var(TString::Format("mtot_bkg_8TeV_slope1_cat%d",c))->getVal()));
    wAll->factory(
	TString::Format("CMS_hgg_bkg_8TeV_slope2_cat%d[%g,-100,100]", 
	c, wAll->var(TString::Format("mtot_bkg_8TeV_slope2_cat%d",c))->getVal()));
//    wAll->factory(
//	TString::Format("CMS_hgg_bkg_8TeV_slope3_cat%d[%g,-100,100]", 
//	c, wAll->var(TString::Format("mtot_bkg_8TeV_slope3_cat%d",c))->getVal()));
  }
  // (2) do reparametrization of background
  for (int c = 0; c < ncat; ++c){ 
    wAll->factory(
	TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(mtotBkg_cat%d,",c,c) +
	TString::Format(" mtot_bkg_8TeV_cat%d_norm=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)+
	TString::Format(" mtot_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d,", c,c)+
	TString::Format(" mtot_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d)", c,c)//+
//	TString::Format(" mtot_bkg_8TeV_slope3_cat%d=CMS_hgg_bkg_8TeV_slope3_cat%d)", c,c)
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
  catdesc.push_back("2 btag");
  catdesc.push_back("1 btag");
  RooDataSet* signalAll       = (RooDataSet*) w->data("Sig");
  // blinded dataset 
  RooDataSet* sigToFit[ncat];
  RooAbsPdf*  mtotGaussSig[ncat];
  RooAbsPdf*  mtotCBSig[ncat];
  RooAbsPdf*  mtotSig[ncat]; 
  for (int c = 0; c < ncat; ++c) {
    sigToFit[c]     = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c)); 
    mtotGaussSig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("mtotGaussSig_cat%d",c));
    mtotCBSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("mtotCBSig_cat%d",c));
    mtotSig[c]       = (RooAbsPdf*)  w->pdf(TString::Format("mtotSig_cat%d",c));
  } // close categories
  RooRealVar* mtot     = w->var("mtot");  
  mtot->setUnit("GeV");
  RooAbsPdf* mtotGaussSigAll  = w->pdf("mtotGaussSig");
  RooAbsPdf* mtotCBSigAll     = w->pdf("mtotCBSig");
  RooAbsPdf* mtotSigAll       = w->pdf("mtotSig");
  RooAbsPdf* mtotBkgAll       = w->pdf("mtotBkg_cat1");
  //****************************//
  // Plot mtot Fit results
  //****************************//
  Float_t minMassFit(minfit),maxMassFit(maxfit); 
  Float_t MASS(Mass);  
  Int_t nBinsMass(93); 
  RooPlot* plotmtotAll = mtot->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
/*  signalAll->plotOn(plotmtotAll);
  gStyle->SetOptTitle(0);
  mtotSigAll->plotOn(plotmtotAll);
  mtotSigAll->plotOn(
	plotmtotAll,Components("mtotGaussSig"),
	LineStyle(kDashed),LineColor(kGreen));
  //mtotSigAll->plotOn(
	//plotmtotAll,Components("mtotCBSig"),
	//LineStyle(kDashed),LineColor(kRed));
  mtotSigAll->paramOn(
	plotmtotAll, 
	ShowConstants(true), 
	Layout(0.15,0.55,0.9), 
	Format("NEU",AutoPrecision(2)));
  plotmtotAll->getAttText()->SetTextSize(0.03);
*/
  TCanvas* c1 = new TCanvas("c1","mtot",0,0,501,501);
  c1->cd(1);
  plotmtotAll->Draw();  
  //****************************//
  // Plot Signal Categories 
  //****************************//
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmtot[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmtot[c] = mtot->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmtot[c],LineColor(kWhite),MarkerColor(kWhite));    
    mtotSig[c]  ->plotOn(plotmtot[c]);
    mtotSig[c]  ->plotOn(
	plotmtot[c],
	Components(TString::Format("mtotGaussSig_cat%d",c)),
	LineStyle(kDashed),LineColor(kGreen));
    mtotSig[c]  ->plotOn(
	plotmtot[c],
	Components(TString::Format("mtotCBSig_cat%d",c)),
	LineStyle(kDashed),LineColor(kRed));
    mtotSig[c]  ->paramOn(plotmtot[c]);
    sigToFit[c]  ->plotOn(plotmtot[c]);
    TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
    plotmtot[c]->SetTitle("CMS preliminary 19.702/fb ");      
    plotmtot[c]->SetMinimum(0.0);
    plotmtot[c]->SetMaximum(1.40*plotmtot[c]->GetMaximum());
    plotmtot[c]->GetXaxis()->SetTitle("M_{#gamma#gamma jj} (GeV)");
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,501,501);
    plotmtot[c]->Draw();  
    plotmtot[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.95,0.9);
    legmc->AddEntry(plotmtot[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotmtot[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmtot[c]->getObject(3),"Crystal Ball component","L");
    legmc->AddEntry(plotmtot[c]->getObject(2),"Gaussian Outliers","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    //    float effS = effSigma(hist);
    TLatex *lat  = new TLatex(
	minMassFit+1.5,0.85*plotmtot[c]->GetMaximum(),
	" WP4 550 GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(
	minMassFit+1.5,0.75*plotmtot[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
  } // close categories
    return;
} // close makeplots
///////////////////////////////////////////////////////////
// declare histos or what -> NOT  USED
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
  outFile << "#Run with: combine -d hgg.mH130.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
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
  	<< "1.005        -   "
  	<< "1.005        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.11   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## with reparametrization" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard with rep
//////////////////////////////////////////////////
// with reparametrization of BKG
void MakeDataCardLnU(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* datashape[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c]      = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    datashape[c] = (RooAbsPdf*) w->data(TString::Format("mtotBkg_ca%d",c));
  }
  RooRealVar*  lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
//  cout << "#Events expected:        "  << endl;
//  for (int c = 0; c < ncat; ++c) {
//    cout << TString::Format("#Events expected cat%d:   ",c) << datashape[c]->getAnalyticalIntegral()
// << endl;
//  }
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
  outFile << "#Run with: combine hgg.mH500.0_8TeV.txt -M Asymptotic" << endl;
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
  outFile << "############## without reparametrization" << endl;
  outFile << "shapes mtotBkg   cat0 " << TString(fileBkgName)+".root" << " w_all:mtotBkg_cat0" << endl;
  outFile << "shapes mtotBkg   cat1 "<<  TString(fileBkgName)+".root" << " w_all:mtotBkg_cat1" << endl;
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
  	<< "1.005        -   "
  	<< "1.005        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  outFile << "mtotBkg       lnU "
  	<< "-        2   "
  	<< "-        2   "
  	<< "# k means vary between 1/k and k times the value in the workspace" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.11   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## without reparametrization" << endl;
  outFile << "mtot_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mtot_bkg_8TeV_slope1_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mtot_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mtot_bkg_8TeV_slope2_cat1         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard with LnU

void MakeDataCardonecat(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
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
  TString filename(cardDir+TString(fileBaseName)+"onecat.txt");
  ofstream outFile(filename);
  outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH130.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax 1" << endl;
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
  outFile << "bin          cat0 " << endl;
  outFile <<  "observation   "  
	<<  data[0]->sumEntries() << "  " 
	<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                      cat0       cat0   " << endl;
  outFile << "process                 mtotSig     mtotBkg" << endl;
  outFile << "process                    0          1    " << endl;
  outFile <<  "rate                      " 
	   << "  " << sigToFit[0]->sumEntries() << "  " <<  1  
	   << "  " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV           lnN "
	<< "1.022   -     "
	<< endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance              lnN " 
	<< "1.015        -   "
	<<"# JER and JES " << endl;
  outFile << "btag_eff          lnN " 
	<< "1.06        -  "
	<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g       lnN "
  	<< "1.010        -   "
  	<< "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
	<< "1.01         -   "
	<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maa_acceptance       lnN "
  	<< "1.005        -   "
  	<< "1.005        -   "
  	<< "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift    param   1   0.006   # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale     param   1   0.11   # optimistic estimative of resolution uncertainty  " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## with reparametrization" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm           flatParam  # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope2_cat0         flatParam  # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard

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
  // mtot signal all categories
  //****************************//
  RooRealVar* mtot_sig_m0     = w->var("mtot_sig_m0");  
  RooRealVar* mtot_sig_sigma  = w->var("mtot_sig_sigma");
  RooRealVar* mtot_sig_alpha  = w->var("mtot_sig_alpha"); 
  RooRealVar* mtot_sig_n      = w->var("mtot_sig_n"); 
  RooRealVar* mtot_sig_gsigma = w->var("mtot_sig_gsigma");
  RooRealVar* mtot_sig_frac   = w->var("mtot_sig_frac");
  mtot_sig_m0    ->SetName("m_{0}");
  mtot_sig_sigma ->SetName("#sigma_{CB}");
  mtot_sig_alpha ->SetName("#alpha");
  mtot_sig_n     ->SetName("n");
  mtot_sig_gsigma->SetName("#sigma_G");  
  mtot_sig_frac  ->SetName("f_G");  
  mtot_sig_m0    ->setUnit("GeV");
  mtot_sig_sigma ->setUnit("GeV");
  mtot_sig_gsigma->setUnit("GeV"); 
  //****************************//
  // mtot background  
  //****************************//
  RooRealVar* mtot_bkg_8TeV_slope1  = w->var("mtot_bkg_8TeV_slope1");
  mtot_bkg_8TeV_slope1              ->SetName("a_{B}");
  mtot_bkg_8TeV_slope1              ->setUnit("1/GeV");
  RooRealVar* mtot_bkg_8TeV_slope2  = w->var("mtot_bkg_8TeV_slope2");
  mtot_bkg_8TeV_slope2              ->SetName("a_{B}");
  mtot_bkg_8TeV_slope2              ->setUnit("1/GeV");
  //****************************//
  // mtot per category  
  //****************************//
  for (int c = 0; c < ncat; ++c) {
    mtot_sig_m0     = (RooRealVar*) w->var(TString::Format("mtot_sig_m0_cat%d",c));
    mtot_sig_sigma  = (RooRealVar*) w->var(TString::Format("mtot_sig_sigma_cat%d",c));
    mtot_sig_alpha  = (RooRealVar*) w->var(TString::Format("mtot_sig_alpha_cat%d",c));
    mtot_sig_n      = (RooRealVar*) w->var(TString::Format("mtot_sig_n_cat%d",c));
    mtot_sig_gsigma = (RooRealVar*) w->var(TString::Format("mtot_sig_gsigma_cat%d",c));
    mtot_sig_frac   = (RooRealVar*) w->var(TString::Format("mtot_sig_frac_cat%d",c));
    mtot_sig_m0     ->SetName("m_{0}"); 
    mtot_sig_sigma  ->SetName("#sigma_{CB}");
    mtot_sig_alpha  ->SetName("#alpha");
    mtot_sig_n      ->SetName("n"); 
    mtot_sig_gsigma ->SetName("#sigma_{G}");
    mtot_sig_frac   ->SetName("f_{G}");
    mtot_sig_m0     ->setUnit("GeV");
    mtot_sig_sigma  ->setUnit("GeV");
    mtot_sig_gsigma ->setUnit("GeV");
    mtot_bkg_8TeV_slope1 = w->var(TString::Format("mtot_bkg_8TeV_slope1_cat%d",c));
    mtot_bkg_8TeV_slope1 ->SetName("p_{B}^{1}");
    mtot_bkg_8TeV_slope1 ->setUnit("1/GeV");
    mtot_bkg_8TeV_slope2 = w->var(TString::Format("mtot_bkg_8TeV_slope2_cat%d",c));
    mtot_bkg_8TeV_slope2 ->SetName("p_{B}^{2}");
    mtot_bkg_8TeV_slope2 ->setUnit("1/GeV^{2}");
//    RooRealVar* mtot_bkg_8TeV_frac = w->var(TString::Format("mtot_bkg_8TeV_frac_cat%d",c));
//    mtot_bkg_8TeV_frac ->SetName("f");
  }
} // close setparameters
