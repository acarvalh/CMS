mgg[100,180];

mgg_bkg_8TeV_slope_cat0[1.0,0, 1];
mgg_bkg_8TeV_slope1_cat0[0.9,-10.0, 10.0];
mgg_bkg_8TeV_slope2_cat0[0.6,-10.0, 10.0];
mgg_bkg_8TeV_slope3_cat0[0.9,-10.0, 10.0];

mgg_bkg_8TeV_slope_cat1[1.0,0, 1];
mgg_bkg_8TeV_slope1_cat1[0.9,-10.0, 10.0];
mgg_bkg_8TeV_slope2_cat1[0.6,-10.0, 10.0];
mgg_bkg_8TeV_slope3_cat1[0.3,-10.0, 10.0];

mgg_bkg_8TeV_slope_cat2[1.0,0, 1];
mgg_bkg_8TeV_slope1_cat2[0.9,-10.0, 10.0];
mgg_bkg_8TeV_slope2_cat2[0.6,-10.0, 10.0];
mgg_bkg_8TeV_slope3_cat2[0.3,-10.0, 10.0];

mgg_sig_m0_cat0[110.0, 105, 155];
mgg_sig_sigma_cat0[1.2, 0.8, 1.6];
mgg_sig_alpha_cat0[2.0, 1.0, 2.5]; 
mgg_sig_n_cat0[2.0, 1.0, 5.0]; 
mgg_sig_gsigma_cat0[5.0, 3.0, 8.0];
mgg_sig_frac_cat0[0.1, 0, 0.2];

mgg_sig_m0_cat1[110.0, 70, 160];
mgg_sig_sigma_cat1[1.8, 1.3, 2.5];
mgg_sig_alpha_cat1[2.0, 1.2, 5]; 
mgg_sig_n_cat1[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat1[5.0, 3.0, 8.0];
mgg_sig_frac_cat1[0.1, 0, 0.3];

GaussSig_cat0 = Gaussian(mgg, mgg_sig_m0_cat0, mgg_sig_gsigma_cat0);
CBSig_cat0    = CBShape(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha_cat0, mgg_sig_n_cat0);
mggSig_cat0      = AddPdf(GaussSig_cat0, CBSig_cat0, mgg_sig_frac_cat0);

GaussSig_cat1 = Gaussian(mgg, mgg_sig_m0_cat1, mgg_sig_gsigma_cat1);
CBSig_cat1    = CBShape(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha_cat1, mgg_sig_n_cat1);
mggSig_cat1      = AddPdf(GaussSig_cat1, CBSig_cat1, mgg_sig_frac_cat1);

mgg_hig_m0_cat0[110.0, 105, 155];
mgg_hig_sigma_cat0[1.2, 0.8, 1.6];
mgg_hig_alpha_cat0[2.0, 1.0, 2.5]; 
mgg_hig_n_cat0[2.0, 1.0, 5.0]; 
mgg_hig_gsigma_cat0[5.0, 3.0, 8.0];
mgg_hig_frac_cat0[0.1, 0, 0.2];

mgg_hig_m0_cat1[110.0, 70, 160];
mgg_hig_sigma_cat1[1.8, 1.3, 2.5];
mgg_hig_alpha_cat1[2.0, 1.2, 5]; 
mgg_hig_n_cat1[2.0, 1.5, 10]; 
mgg_hig_gsigma_cat1[5.0, 3.0, 8.0];
mgg_hig_frac_cat1[0.1, 0, 0.3];

GaussHig_cat0 = Gaussian(mgg, mgg_hig_m0_cat0, mgg_hig_gsigma_cat0);
CBHig_cat0    = CBShape(mgg, mgg_hig_m0_cat0, mgg_hig_sigma_cat0, mgg_hig_alpha_cat0, mgg_hig_n_cat0);
mggHig_cat0      = AddPdf(GaussHig_cat0, CBHig_cat0, mgg_hig_frac_cat0);

GaussHig_cat1 = Gaussian(mgg, mgg_hig_m0_cat1, mgg_hig_gsigma_cat1);
CBHig_cat1    = CBShape(mgg, mgg_hig_m0_cat1, mgg_hig_sigma_cat1, mgg_hig_alpha_cat1, mgg_hig_n_cat1);
mggHig_cat1      = AddPdf(GaussHig_cat1, CBHig_cat1, mgg_hig_frac_cat1);
