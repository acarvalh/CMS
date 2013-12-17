mtot[150,1200];

mtot_sig_m0_cat0[1200, 1100, 1300];
mtot_sig_sigma_cat0[20, 5., 40.0];
mtot_sig_alpha_cat0[1.0, 0.5, 3]; 
mtot_sig_n_cat0[4.0, 0.5, 10]; 
mtot_sig_gsigma_cat0[20, 5., 40.0];
mtot_sig_frac_cat0[0.5, 0, 1.0];

mtot_sig_m0_cat1[1200, 1100, 1300];
mtot_sig_sigma_cat1[20, 5., 40.0];
mtot_sig_alpha_cat1[1.0, 0.5, 3]; 
mtot_sig_n_cat1[4.0, 0.5, 10]; 
mtot_sig_gsigma_cat1[20.0, 10.0, 50.0];
mtot_sig_frac_cat1[0.5, 0, 1.0];

mtotGaussSig_cat0 = Gaussian(mtot, mtot_sig_m0_cat0, mtot_sig_gsigma_cat0);
mtotCBSig_cat0    = CBShape(mtot, mtot_sig_m0_cat0, mtot_sig_sigma_cat0, mtot_sig_alpha_cat0, mtot_sig_n_cat0);
mtotSig_cat0      = AddPdf(mtotGaussSig_cat0, mtotCBSig_cat0, mtot_sig_frac_cat0);

mtotGaussSig_cat1 = Gaussian(mtot, mtot_sig_m0_cat1, mtot_sig_gsigma_cat1);
mtotCBSig_cat1    = CBShape(mtot, mtot_sig_m0_cat1, mtot_sig_sigma_cat1, mtot_sig_alpha_cat1, mtot_sig_n_cat1);
mtotSig_cat1      = AddPdf(mtotGaussSig_cat1, mtotCBSig_cat1, mtot_sig_frac_cat1);

mtot_bkg_8TeV_slope1_cat0[300.,200, 400.0];
mtot_bkg_8TeV_slope2_cat0[1., 0.0, 100.0];
mtot_bkg_8TeV_slope3_cat0[0.5,0, 10.0];

mtot_bkg_8TeV_slope1_cat1[270.,200, 400.0];
mtot_bkg_8TeV_slope2_cat1[1., 0.0, 100.0];
mtot_bkg_8TeV_slope3_cat1[0.5,0.0, 10.0];

wei[1,0,10];
sqrtS[8000., 8000., 8000., 8000.]
