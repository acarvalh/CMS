{
  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I/afs/cern.ch/cms/slc5_amd64_gcc472/lcg/roofit/5.34.04-cms2/include");
  incpath.Append(" -L/afs/cern.ch/cms/slc5_amd64_gcc472/lcg/roofit/5.34.04-cms2/lib");
  gSystem->SetIncludePath(incpath.Data());
  gSystem->Load("libRooFitCore") ;
  gSystem->Load("libRooStats") ;
  using namespace RooFit;
  using namespace RooStats ;
  gStyle->SetStripDecimals(kFALSE);
  gROOT->ProcessLine(".L GaussExp.cxx+");
  //gROOT->LoadMacro("/afs/cern.ch/work/a/acarvalh/CMSSW_5_3_7/ggfits/R2GGBBFitter_mtot_side.cc");
  gROOT->LoadMacro("R2GGBBFitter_mtot.cc");
  runfits(500,1);
}
