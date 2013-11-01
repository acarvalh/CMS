mtot[320,1200];
mtot_sig_m0_cat0[500, 450, 550];
mtot_sig_sigma_cat0[10, 3., 40.0];
mtot_sig_alpha_cat0[1.0, 0.5, 3]; 
mtot_sig_n_cat0[10.0, 0.5, 20]; 
mtot_sig_gsigma_cat0[20.0, 10.0, 180.0];
mtot_sig_frac_cat0[0.5, 0, 1.0];

mtot_sig_m0_cat1[500, 450, 550];
mtot_sig_sigma_cat1[10, 3., 40.0];
mtot_sig_alpha_cat1[1.0, 0.5, 5]; 
mtot_sig_n_cat1[10.0, 0.5, 20]; 
mtot_sig_gsigma_cat1[20.0, 10.0, 180.0];
mtot_sig_frac_cat1[0.5, 0, 1.0];

GaussSig_cat0 = Gaussian(mtot, mtot_sig_m0_cat0, mtot_sig_gsigma_cat0);
CBSig_cat0    = CBShape(mtot, mtot_sig_m0_cat0, mtot_sig_sigma_cat0, mtot_sig_alpha_cat0, mtot_sig_n_cat0);
mtotSig_cat0      = AddPdf(GaussSig_cat0, CBSig_cat0, mtot_sig_frac_cat0);

GaussSig_cat1 = Gaussian(mtot, mtot_sig_m0_cat1, mtot_sig_gsigma_cat1);
CBSig_cat1    = CBShape(mtot, mtot_sig_m0_cat1, mtot_sig_sigma_cat1, mtot_sig_alpha_cat1, mtot_sig_n_cat1);
mtotSig_cat1      = AddPdf(GaussSig_cat1, CBSig_cat1, mtot_sig_frac_cat1);
