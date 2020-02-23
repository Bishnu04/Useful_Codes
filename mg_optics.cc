// copied from HT_optics.cc
// like to include the TT dta for tune. January 02, 2020 
// On jan 06, 2020, I am adding some flags HHRflag for example to include the Al events for tune. A copy of the code is saved on backup
// We would like to add all of the Al events to see 
// Jan 07 I am adding some more variables with _4 eg npeak_4 to include the  Al/T data for tune

extern double calcf2t_th(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_ph(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_mom(double* P, 
			  double xf, double xpf,
			  double yf, double ypf,double);

const double  XFPm=-0.7,  XpFPm=-0.15; // m is the mean from the old definition
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; // tm = target offset.. MOmm is the momentum offset
const double  XFPr=1.3,   XpFPr=0.27; // r is the scaling factor or range
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; // tr is the target range
const double  PLm = 25.4, PLr=0.7; // m is the offset and PLr is the path laegth range
const double  Ztm = -0.15,Ztr=0.35; //Ztm  z position at target  point offset
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double tune(double* pa, int j);

// Lorenzian Peak function
Double_t cauchy(Double_t *x, Double_t *par){
  return (0.5*par[0]*par[1]/TMath::Pi()) / 
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2]) 
		+ 0.25*par[1]*par[1]);
}

const double N_b = 6.671;
Double_t user_gaus(double*x,double*par){  // gaussian function for nnl peak 
  double g1;
  g1 = par[0] * exp(-0.5*pow((x[0]-par[1])/par[2],2.0));
  
  return g1 + N_b;
}

const int nmax = 4000; // was 3000 before Dec5
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;

const int npeak = 2;
double Lambda_width[npeak] = {2.52, 2.48}; //6.5,8.8}
double Lambda_cent[npeak] ={1115.68,1192.30};
double Lambda_real[npeak] ={1115.683,1192.642}; // Mev



double p10[nmax],p11[nmax],p12[nmax];
double p13[nmax],p14[nmax],p15[nmax];
double p16[nmax],p17[nmax],p18[nmax],p19[nmax];
double phir[nmax];
double phil[nmax];
// ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
const int nmax_2 = 3500;   // 2400 before Dec 5
double x_2[nmax_2], y_2[nmax_2]; 
double xp_2[nmax_2], yp_2[nmax_2];
double z_recon_2[nmax_2];
int foil_flag_2[nmax_2];

const int npeak_2 = 1;
double Lambda_width_2[npeak_2] = {2.42}; //6.5,8.8}
double Lambda_cent_2[npeak_2] ={1115.68};
double Lambda_real_2[npeak_2] ={1115.683}; // Mev

double p10_2[nmax_2],p11_2[nmax_2],p12_2[nmax_2];
double p13_2[nmax_2],p14_2[nmax_2],p15_2[nmax_2];
double p16_2[nmax_2];
double phir_2[nmax_2];
double phil_2[nmax_2];
int ntune_event_2 = 0;

//===============================================Jan 07 2020==========================================
const int nmax_4 = 3500;   
double x_4[nmax_4], y_4[nmax_4]; 
double xp_4[nmax_4], yp_4[nmax_4];
double z_recon_4[nmax_4];
int foil_flag_4[nmax_4];

const int npeak_4 = 4; // JAn 07
double Lambda_width_4[npeak_4] ={0.94,0.80,1.15,1.15};
double Lambda_cent_4[npeak_4] ={-3.907,5.918,20.108,30.226};
double Lambda_real_4[npeak_4] ={-3.907,5.918,20.108,30.226};

double p10_4[nmax_4],p11_4[nmax_4],p12_4[nmax_4];
double p13_4[nmax_4],p14_4[nmax_4],p15_4[nmax_4];
double p16_4[nmax_4];
double phir_4[nmax_4];
double phil_4[nmax_4];
int ntune_event_4 = 0;

//))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
// ((((((((((((((((((((((((((((((((   t3   (((((((((((((((((((((((((((((((((((((((((((((((((((   // Jan_02
const int nmax_3 =3500;
double x_3[nmax_3], y_3[nmax_3];
double xp_3[nmax_3], yp_3[nmax_3];
double z_recon_3[nmax_3];
int foil_flag_3[nmax_3];

const int npeak_3 = 4;
double Lambda_width_3[npeak_3] ={0.94,0.80,1.15,1.15};
double Lambda_cent_3[npeak_3] ={-3.907,5.918,20.108,30.226};
double Lambda_real_3[npeak_3] ={-3.907,5.918,20.108,30.226};

double p10_3[nmax_3],p11_3[nmax_3],p12_3[nmax_3];
double p13_3[nmax_3],p14_3[nmax_3],p15_3[nmax_3];
double p16_3[nmax_3];
double phir_3[nmax_3];
double phil_3[nmax_3];
int ntune_event_3 = 0;

const double m_Al =25.1267; // Al target mass // BY Dr. Tang on Dec 19 2019
const double m_T = 2.808921; // for tritium target by Gogami Tritium target mass

//))))))))))))))))))))))))))))))))   t3  ))))))))))))))))))))))))))))))))))))))))))))))))))))))
//========================================

const int Total_Par = 126;
double thetaL_opt[nmax];
double phiL_opt[nmax];
double thetaR_opt[nmax];
double phiR_opt[nmax];
double momL_opt[nmax];
double momR_opt[nmax];
const int Mom_Par = 252;
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
const double me = 0.000511;
const double mk = 0.493677;
const double mp = 0.938272;
const double mL = 1.115683;
//extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);// before dec3 
extern double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt);
void mg_optics(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ======================================== 
 
  TChain * t1 = new TChain("T");  
  TChain * t2 = new TChain("T"); 
  TChain * t3 = new TChain("T");// Jan_02
  //// TChain * t5 = new TChain("T");// Jan_02

  t1->Add("./Rootfiles/DEC17_Rootfiles/DEC17_H149_542.root");// replayed on Dec 17, 2019 replay by ole
  t2->Add("./Rootfiles/DEC17_Rootfiles/DEC17_HT_552_716.root");
  t3->Add("./Rootfiles/DEC17_Rootfiles/DEC23_T221_830.root");  //// slightly better,// decoded on Dec23/ 2019 // Jan_02
  //// t5->Add("./Rootfiles/DEC17_Rootfiles/He_Jan12.root");  //// slightly better,// decoded on Dec23/ 2019 // Jan_02
 
  double ent = t1->GetEntries();
  double ent_2 = t2->GetEntries();
  double ent_3 = t3->GetEntries(); // Jan_02
  ////  double ent_5 = t5->GetEntries(); // Jan_12
  
  // ent_3 =5000;
  // ent = 5000;
  cout<<"entry in the t1=="<<ent<<endl;
  cout<<"entry in the t2=="<<ent_2<<endl;
  cout<<"entry in the t3=="<<ent_3<<endl; // Jan_02
  //// cout<<"entry in the t5=="<<ent_5<<endl; // Jan_02
  
  const int max = 100;
  Double_t trig5[max]; // JUly 01, 2019 
  double momL[max];
  double momR[max]; 
  
  double lvz[max],rvz[max];// raster corrected 
  double th1[max], ph1[max];// RHRS angle 
  double th2[max], ph2[max];  
  double delta_pep[max];     // target straggling
  double pep_real[max]; 
  double delta_pk[max];
  double pk_real[max];
  double par_ep[3];
  double par_k[3];
  double mm; 
  double hallap;
  
  double l_th_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double l_y_fp[max];

  double r_th_fp[max];
  double r_ph_fp[max];
  double r_x_fp[max];
  double r_y_fp[max];
  const int n = 16; 
  double ctime; 
 
 
  double z_av[nmax];
  double z_av_1[nmax];
  double a1, a2;

  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  Double_t trig5_2[max]; // JUly 01, 2019 
  double momL_2[max];
  double momR_2[max]; 
  
  double lvz_2[max],rvz_2[max];// raster corrected 
  double th1_2[max], ph1_2[max];// RHRS angle 
  double th2_2[max], ph2_2[max];  
  double delta_pep_2[max];     // target straggling
  double pep_real_2[max]; 
  double delta_pk_2[max];
  double pk_real_2[max];
  double par_ep_2[3];
  double par_k_2[3];
  double mm_2;
  double mm_4;
  double mm_ht;
  double hallap_2;
  
  double l_th_fp_2[max];
  double l_ph_fp_2[max];
  double l_x_fp_2[max];
  double l_y_fp_2[max];

  double r_th_fp_2[max];
  double r_ph_fp_2[max];
  double r_x_fp_2[max];
  double r_y_fp_2[max];
  double ctime_2; 
 
  double z_av_2[nmax];
  double z_av_1_2[nmax];
  double a1_2, a2_2;

  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t3 variables starts here    ((((((((((((((((((((((((((((((((((((((((((((((((((( Jan_02
  Double_t trig5_3[max]; // JUly 01, 2019 
  double momL_3[max];
  double momR_3[max]; 
  double lvz_3[max],rvz_3[max];// raster corrected 
  double th1_3[max], ph1_3[max];// RHRS angle 
  double th2_3[max], ph2_3[max];  
  double delta_pep_3[max];     // target straggling
  double pep_real_3[max]; 
  double delta_pk_3[max];
  double pk_real_3[max];
  double par_ep_3[3];
  double par_k_3[3];
  double mm_3;
  double mm_t;  
  double mm_Al;
  double mm_Al1;
  double a1_3, a2_3;
  double mm_h;

  double hallap_3;

  double l_th_fp_3[max];
  double l_ph_fp_3[max];
  double l_x_fp_3[max];
  double l_y_fp_3[max];

  double r_th_fp_3[max];
  double r_ph_fp_3[max];
  double r_x_fp_3[max];
  double r_y_fp_3[max];
  double ctime_3; 
  
  double z_av_3[nmax];
  double z_av_1_3[nmax];
  //))))))))))))))))))))))))))))))))   t3 variables up to here   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t1 branch address  (((((((((((((((((((((((((((((((((((((((((((((((((((
  t1->SetBranchAddress("HALLA_p", &hallap);  
  t1->SetBranchAddress("DR.T5", &trig5);  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);   
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp); 
  t1->SetBranchAddress("coin_time",  &ctime); 
  t1->SetBranchAddress("ztR_wRC",  &rvz);
  t1->SetBranchAddress("ztL_wRC",  &lvz);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);

  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  t2->SetBranchAddress("HALLA_p", &hallap_2);  
  t2->SetBranchAddress("DR.T5", &trig5_2);  
  t2->SetBranchAddress("L.tr.x",   &l_x_fp_2);
  t2->SetBranchAddress("L.tr.y",   &l_y_fp_2);
  t2->SetBranchAddress("L.tr.th",  &l_th_fp_2);
  t2->SetBranchAddress("L.tr.ph",  &l_ph_fp_2);   
  t2->SetBranchAddress("R.tr.x",   &r_x_fp_2);
  t2->SetBranchAddress("R.tr.y",   &r_y_fp_2);
  t2->SetBranchAddress("R.tr.th",  &r_th_fp_2);
  t2->SetBranchAddress("R.tr.ph",  &r_ph_fp_2); 
  t2->SetBranchAddress("coin_time",  &ctime_2); 
  t2->SetBranchAddress("ztR_wRC",  &rvz_2);
  t2->SetBranchAddress("ztL_wRC",  &lvz_2);
  t2->SetBranchAddress("R.a1.asum_c", &a1_2);
  t2->SetBranchAddress("R.a2.asum_c", &a2_2);

  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t3   ((((((((((((((((((((((((((((((((((((((((((((((((((( Jan_02
  t3->SetBranchAddress("HALLA_p", &hallap_3);  
  t3->SetBranchAddress("DR.T5", &trig5_3);   
  t3->SetBranchAddress("L.tr.x",   &l_x_fp_3);
  t3->SetBranchAddress("L.tr.y",   &l_y_fp_3);
  t3->SetBranchAddress("L.tr.th",  &l_th_fp_3);
  t3->SetBranchAddress("L.tr.ph",  &l_ph_fp_3);
  t3->SetBranchAddress("R.tr.x",   &r_x_fp_3);
  t3->SetBranchAddress("R.tr.y",   &r_y_fp_3);
  t3->SetBranchAddress("R.tr.th",  &r_th_fp_3);
  t3->SetBranchAddress("R.tr.ph",  &r_ph_fp_3); 
  t3->SetBranchAddress("coin_time",  &ctime_3); 
  t3->SetBranchAddress("ztR_wRC",  &rvz_3);
  t3->SetBranchAddress("ztL_wRC",  &lvz_3);
  t3->SetBranchAddress("R.a1.asum_c", &a1_3);
  t3->SetBranchAddress("R.a2.asum_c", &a2_3);

 //))))))))))))))))))))))))))))))))   t3 barnch address up to here    ))))))))))))))))))))))))))))))))))))))))))))))))))))))  
  TFile* fnew = new TFile("./output_root/angle_lhrs.root","recreate"); 
  TTree* tnew = new TTree("tree","For z calibration (LHRS)");
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D");
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D");
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D");
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D");
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D");  
  tnew->Branch("L.tr.tg_th_TH2", &th2, "L.tr.tg_th_TH2[100]/D");
  tnew->Branch("L.tr.tg_ph_PH2", &ph2, "L.tr.tg_ph_PH2[100]/D");
  
 
  double XFP, XpFP;
  double YFP, YpFP;
  double R_XFP, R_XpFP; 
  double R_YFP, R_YpFP;

  // ((((((((((((((((((((((((((((((((((((((((((( for t2 ((((((((((
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double R_XFP_2, R_XpFP_2; 
  double R_YFP_2, R_YpFP_2;

  //)))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((((( for t3 (((((((((( Jan_02
  double XFP_3, XpFP_3;
  double YFP_3, YpFP_3;
  double R_XFP_3, R_XpFP_3; 
  double R_YFP_3, R_YpFP_3; 

  //)))))))))))))))))))))))))))))))))))))))))))))
 // ((((((((((((((((((((((((((((((((((((((((((( for t2 ((((((((((
  double XFP_4, XpFP_4;
  double YFP_4, YpFP_4;  // jan 07 2020
  double R_XFP_4, R_XpFP_4; 
  double R_YFP_4, R_YpFP_4;
  // ===============or LHRS  theta information input==========   
  ntune_event = 0;
  for(int i=0 ; i<Total_Par; i++){
    thetaL_opt[i] = -2222.0;
  }
  
  char name_Angle_L[500]; 
  // sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");//theta_3rd_LHRS_Opt_7.dat
  sprintf(name_Angle_L,"./matrices/theta_L4th_4th_6.dat");// optimized on OCT 23, 2019 with SS data
  ///// sprintf(name_Angle_L,"./matrices/theta_L4th_6th_0.dat"); //optimized on dec 6 ,2019 with HH HT data
  ifstream Angle_L(name_Angle_L);
  double Theta_L[Total_Par];    
  for(int i =0; i<Total_Par;i++){
    double par1 =0.0;
    int p1 =0;
    Angle_L>>par1>>p1>>p1>>p1>>p1>>p1;
    Theta_L[i]=par1;
    thetaL_opt[i] = Theta_L[i];
  }
  Angle_L.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // for LHRS  phi information input
  //--------------------------------------------------------
  
  ntune_event = 0;
  for(int i =0;i<Total_Par;i++){
    phiL_opt[i] = -2222.0;
  }
  char name_angle_phi[500];
  //  sprintf(name_angle_phi,"./matrices/phi_LHRS_3rd_Opt_9.dat");
  sprintf(name_angle_phi,"./matrices/phi_L4th_5th_5.dat");// optimized on OCT 23, 2019  
  //// sprintf(name_angle_phi,"./matrices/phi_L4th_6th_0.dat");//optimized on dec 6 ,2019 with HH HT data

  ifstream angle_phi(name_angle_phi);
  double PHI_L[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par2 =0.0;
    int p2 =0;
    angle_phi>>par2>>p2>>p2>>p2>>p2>>p2;
    PHI_L[i] = par2;
    phiL_opt[i]= PHI_L[i];
  }
  angle_phi.close();
  // LHRS momentum information========================July 20, 2019
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momL_opt[i] = -2222.0;
  }
  char name_Mom_lhrs[500];
  
 
  // sprintf(name_Mom_lhrs,"./All_Matrices/mom_LHRS_5_upto3.dat"); 
  // sprintf(name_Mom_lhrs,"./All_Matrices/mom30_L_6th_1.dat");
  // sprintf(name_Mom_lhrs,"./All_Matrices/momLH_3rd_1.dat");
  //  sprintf(name_Mom_lhrs,"./MOM_MATRICES/mom_LHRS_5_0.dat");// orignal first 5th order matrix Nov 15, 2019
  // sprintf(name_Mom_lhrs,"./MOM_MATRICES/mom_L5_4th_2.dat");
  sprintf(name_Mom_lhrs,"./MOM_MATRICES/LMOM5_36th_4.dat");
  ifstream Mom_lhrs(name_Mom_lhrs);
  double mom_L[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par5 = 0.0;
    int p5 =0;
    Mom_lhrs>>par5>>p5>>p5>>p5>>p5>>p5;
    mom_L[i]= par5;
    momL_opt[i] = mom_L[i];
  }
  Mom_lhrs.close();
  // up to here thelhrs momentum matrix======================


  
  // =======RHRS theta input information
  ntune_event =0;
  for(int i =0;i<Total_Par;i++){
    thetaR_opt[i] = -2222.0;
  }
  char name_Angle_R[500]; 
  sprintf(name_Angle_R,"./All_Matrices/xpt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  // sprintf(name_Angle_R,"./matrices/thetaR4_20th_0.dat"); // tuned with H/H and H/T data on Feb 01, 2020
  ifstream Angle_R(name_Angle_R);
  double Theta_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par3 =0.0;
    int p3 = 0;
    Angle_R>>par3>>p3>>p3>>p3>>p3>>p3;
    Theta_R[i]=par3;
    thetaR_opt[i] = Theta_R[i];
  }
  Angle_R.close();
  //====================================================
  //=======RHRS phi input information===============
  ntune_event = 0;
  for(int i =0;i<Total_Par;i++){
    phiR_opt[i] = -2222.0;
  }
  char name_phi_Rhrs[500];
 
  sprintf(name_phi_Rhrs,"./All_Matrices/ypt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  //  sprintf(name_phi_Rhrs,"./matrices/phiR4_20th_0.dat"); // tuned with H/H and H/T data on Feb 01, 2020
  ifstream phi_Rhrs(name_phi_Rhrs);
  double PHI_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par4 =0.0;
    int p4 =0;
    phi_Rhrs>>par4>>p4>>p4>>p4>>p4>>p4;
    PHI_R[i] = par4;
    phiR_opt[i]= PHI_R[i];
  }
  phi_Rhrs.close();
  //==================================================
  // =====RHRS momentum recon==========================6
 
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momR_opt[i] = -2222.0;
  }
  char name_Mom_rhrs[500];
 
  //sprintf(name_Mom_rhrs,"./All_Matrices/MOMR_30th_5.dat");
  //  sprintf(name_Mom_rhrs,"./All_Matrices/mom_RHRS_5_upto3.dat");
  //  sprintf(name_Mom_rhrs,"./All_Matrices/mom30_R_2nd_0.dat");
  // sprintf(name_Mom_rhrs,"./All_Matrices/momRH_3rd_4.dat");
  // sprintf(name_Mom_rhrs,"./MOM_MATRICES/mom_RHRS_5_0.dat"); // matrix prodeced on the Nov 15, 2019
  //  sprintf(name_Mom_rhrs,"./MOM_MATRICES/mom_R5_8th_2.dat"); // matrix prodeced on the Nov 15, 2019
  sprintf(name_Mom_rhrs,"./MOM_MATRICES/RMOM5_36th_6.dat"); // DEc3 started
  ifstream Mom_rhrs(name_Mom_rhrs);
  double mom_R[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par6 = 0.0;
    int p6 =0;
    Mom_rhrs>>par6>>p6>>p6>>p6>>p6>>p6;
    mom_R[i]= par6;
    momR_opt[i] = mom_R[i];
  }
  Mom_rhrs.close();
  // =====RHRS momentum recon up to here=============== 
  
  TH1F *h = new TH1F("h"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); //
  TH1F *h_lo = new TH1F("h_lo"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); //
  TH1F *h_2 = new TH1F("h_2"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); // fot lambda and sigma 
  TH1F *h_2lo = new TH1F("h_2lo"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); //cauchy function
  
  TH1F *hh = new TH1F("hh","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/ 0.5 MeV ",600,-100,200);
  TH1F *ht = new TH1F("ht","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/ 0.5 MeV ",600,-100,200);
  TH1F *htt = new TH1F("htt","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/ 0.5 MeV ",600,-100,200);

  TH1F *hb = new TH1F("hb","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/ 0.5 MeV ",600,-100,200);
  TH1F *hb1 = new TH1F("hb1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/ 0.5 MeV ",600,-100,200);
  TH1F *hb2 = new TH1F("hb2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/ 0.5 MeV ",600,-100,200);
  
  TH1F *h20 = new TH1F("h20","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/MeV ",310,-100,210);
  TH1F *h21 = new TH1F("h21","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/MeV ",310,-100,210 ); 
  TH1F *HT = new TH1F("HT","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/MeV ",310,-100,210);

 
  TH1F *h20_b = new TH1F("h20_b","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/MeV ",310,-100,210);
  TH1F *h21_b = new TH1F("h21_b","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/MeV ",310,-100,210 ); 
  TH1F *HT_b = new TH1F("HT_b","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/MeV ",310,-100,210);

  TH1F *hd = new TH1F("hd","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/0.75 MeV ",416,-100,212);
  TH1F *he = new TH1F("he","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/0.75 MeV ",416,-100,212); 
  TH1F *hf = new TH1F("hf","Al Spectrum(H/T+ TT); -B_{#Lambda}(MeV);Counts/0.75 MeV ",416,-100,212);

  TH1F *hb5 = new TH1F("hb5","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/0.75 MeV ",416,-100,212);
  TH1F *hb6 = new TH1F("hb6","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/0.75 MeV ",416,-100,212); 
  TH1F *hb7 = new TH1F("hb7","Al Spectrum(H/T+ TT); -B_{#Lambda}(MeV);Counts/0.75 MeV ",416,-100,212);

  TH1F *hl = new TH1F("hl","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",248,-100,210);
  TH1F *hm = new TH1F("hm","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",248,-100,210 ); 
  TH1F *hn = new TH1F("hn","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/1.25 MeV ",248,-100,210);
 
  TH1F *hb10 = new TH1F("hb10","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",248,-100,210);
  TH1F *hb11 = new TH1F("hb11","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",248,-100,210 ); 
  TH1F *hb12 = new TH1F("h12","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/1.25 MeV ",248,-100,210);

  TH1F *h50 = new TH1F("h50","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",200,-100,150);
  TH1F *h50b = new TH1F("h50b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",200,-100,150);
  
  TH1F *h52 = new TH1F("h52","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ MeV ",250,-100,150);
  TH1F *h53 = new TH1F("h53","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100,152);
  TH1F *h54 = new TH1F("h54","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",125,-100,150);
  
  TH1F *hb53 = new TH1F("hb53","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100,152);
  TH1F *hb54 = new TH1F("hb54","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",125,-100,150);

  TH1F *h1_2 = new TH1F("h1_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  TH1F *h1_2b = new TH1F("h1_2b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);

  TH1F *h1_qu = new TH1F("h1_qu","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  TH1F *h1_qub = new TH1F("h1_qub","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  // TH1F *h1_3 = new TH1F("h1_3","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  // TH1F *h1_35 = new TH1F("h1_35","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.53 MeV ",164,-100,151);
  // TH1F *h1_4 = new TH1F("h1_4","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.54 MeV ",172,-100,165);
  // TH1F *h1_45 = new TH1F("h1_45","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.56 MeV ",166,-100,159);
  
  // TH1F *h1_55 = new TH1F("h1_55","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.41 MeV ",190,-110,158);
  // TH1F *h1_6 = new TH1F("h1_6","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.42 MeV ",176,-100,150);
  // TH1F *h1_65 = new TH1F("h1_65","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.43 MeV ",179,-100,156);
  // TH1F *h1_7 = new TH1F("h1_7","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.7 MeV ",150,-100,155);
  // TH1F *h1_75 = new TH1F("h1_75","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.75 MeV ",144,-100,152);
  
  
  TH1F *h_h = new TH1F("h_h","H in  T/T data  ;-B_{#Lambda}(MeV);Counts/ 1.5MeV ",148,-100,122);
  
  TH2F *h5 = new TH2F("h5",";Al_Mg B_{#Lambda}(MeV); Tritium Kinematics in the Al region B_{#Lambda}(MeV)",416,-100,212,250,-100,150);
  TH1F *hz = new TH1F("hz"," Tritum Data; Z-average(m);Counts/ mm",400,-0.2,0.2);

  TH1F *h53_t = new TH1F("h53_t","nnL Spectrum, T/T data(-100.5, 151.5) ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100.5,151.5);
  TH1F *h53_tb = new TH1F("h53_tb","nnL Spectrum, T/T data(-100.5, 151.5) ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100.5,151.5);
 

  // TH1F *ha = new TH1F("ha","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.75 MeV ",229,-150,250);
  // TH1F *hb = new TH1F("hb","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/1.75 MeV ",229,-150,250 ); 
  // TH1F *hc = new TH1F("hc","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/1.75 MeV ",229,-150,250);
  // TH1F *He = new TH1F("He","Al Spectrum, Al/He data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",267,-150,250);
  // TH1F *He1 = new TH1F("He1","Al Spectrum, Al/He data  ; -B_{#Lambda}(MeV);Counts/1.25 MeV ",320,-150,250);

  // TH1F *H50 = new TH1F("H50","Al Spectrum(Al/HT + Al/T+ Al/He data); -B_{#Lambda}(MeV);Counts/1.5 MeV ",267,-150,250);
  TH1F *H25 = new TH1F("H25","H/H data;Coincidence Time (ns);Counts/0.15 ns  ",400,-30,30);


  gStyle->SetOptStat(111111);
  
  TH1F *h_1 = new TH1F("h_1"," ",120,1000, 1300);
  TH2F * h10 = new TH2F("h10","H data with H kinematics;Missing Mass(MeV/c^{2});L.tr.th",200,1000,1300, 200,1.2, 2.5);
  h->GetXaxis()->CenterTitle(); 
  h->GetYaxis()->CenterTitle(); 
  h_2->GetXaxis()->CenterTitle();
  h_2->GetYaxis()->CenterTitle();
 
  TH1F *h6 = new TH1F("h6",";RHRS reconstructed Momentum;Counts/ 14.4 mev",250,1.7,2.0); 
  char tempc[500];
  // ======================================================
 
  //===================================================================  
 
  
  bool rtrig = false; 
  for(int i=0; i<nmax; i++){
    x[i]    = -2222.0; 
    y[i]    = -2222.0; 
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    z_av[i] = -2222.0;
    z_av_1[i] = -2222.0;
    phir[i] = -2222.0;
    phil[i] = -2222.0;
    z_recon[i] = -2222.0;  /// Jan 04, 2019
    foil_flag[i] = -1;
  }

  // ((((((((((((((((((((((((((((((((((((((((((((

  bool rtrig_2 = false; 
  for(int i=0; i<nmax_2; i++){
    x_2[i]    = -2222.0; 
    y_2[i]    = -2222.0; 
    xp_2[i]   = -2222.0;
    yp_2[i]   = -2222.0;
    z_av_2[i] = -2222.0;
    z_av_1_2[i] = -2222.0;
    phir_2[i] = -2222.0;
    phil_2[i] = -2222.0;
    z_recon_2[i] = -2222.0; ///Jan 04, 2019
    foil_flag_2[i] = -1;
    //   for the _4 variables for Al/HT data in tune jan 07, 2020

    x_4[i]    = -2222.0; 
    y_4[i]    = -2222.0; 
    xp_4[i]   = -2222.0;
    yp_4[i]   = -2222.0; 
   
    phir_4[i] = -2222.0;
    phil_4[i] = -2222.0;
    z_recon_4[i] = -2222.0; ///Jan 04, 2019
    foil_flag_4[i] = -1;
  }
  
  //)))))))))))))))))))))))))))))))))))))))))))))

  // +++++++++++++++++++++++++ for t1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int i=0; i< ent; i++){
    for(int j=0; j<max; j++){
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0; 
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      th1[j] = -2222.0;
      th2[j] = -2222.0;
      ph1[j] =-2222.0;
      ph2[j] =-2222.0;    
     
      delta_pep[j]= -2222.0;
      pep_real[j] =-2222.0;
      delta_pk[j]= -2222.0;
      pk_real[j] = -2222.0;
     
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
    
      
      
      trig5[j] = 0.0;
      rtrig = false;
    }
   
    trig5[0] = 0.0;
    rtrig = false;
   
    t1->GetEntry(i);
   
   
    if(trig5[0]>1.0) rtrig = true; //JUly 01, 2019
    else rtrig = false;

    z_av[0] = (lvz[0] + rvz[0])/2.0;
    z_av_1[0] =  z_av[0];
   
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    R_XFP   = r_x_fp[0];
    R_XpFP  = r_th_fp[0];
    R_YFP   = r_y_fp[0];
    R_YpFP  = r_ph_fp[0];
 
    if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.040  && fabs(z_av_1[0])<0.10){ // need change del_z cut  
       ///// && a1<120.0 && a2>1650.0 && a2<6800.0
      
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      
      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;

      z_av[0] =(z_av[0]- Ztm)/Ztr;
      
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av_1[0]);
      th2[0] = th2[0]*Xptr + Xptm; 
     
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP, z_av_1[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;
          
     
      momL[0] =  calcf2t_mom(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]); 
      momL[0] = momL[0]*Momr + Momm;
    
      

      // par_ep[0] = momL[0] ; // Dec4 3 lines
      par_ep[1] = -th2[0]; // right handed system
      par_ep[2] = -ph2[0]; // right handed system
      double holiang;


      
      // Target struggling LHRS step #7
      if( z_av_1[0]<8.0e-2){	
	
	holiang =  par_ep[2] + hrs_ang;
	holiang=-holiang;
	delta_pep[0] = -1.35758 * sin(-4.59571 * holiang) + 2.09093;

      } 
      else{
	holiang =  par_ep[2] + hrs_ang;
	holiang=-holiang;	
	delta_pep[0] = 6.23409e-3 * holiang + 4.03363e-1;
      }       
      pep_real[0] = momL[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV
      
      par_ep[0] = pep_real[0] ; 
    



      // RHRS angle and momentum calculation      
      th1[0] = calcf2t_th(Theta_R, R_XFP, R_XpFP, R_YFP, R_YpFP,   z_av[0]);
      th1[0] = th1[0]*Xptr + Xptm;
     
      ph1[0] = calcf2t_ph(PHI_R, R_XFP, R_XpFP, R_YFP, R_YpFP, z_av[0]);
      ph1[0] = ph1[0]*Yptr + Yptm;
     
      
      momR[0] =  calcf2t_mom(mom_R, R_XFP, R_XpFP, R_YFP, R_YpFP,  z_av[0]);
      momR[0] = momR[0]*Momr+Momm;
     
      
      par_k[1] = -th1[0]; /// DEc4 2019 2 lines
      par_k[2] = -ph1[0];
      double holiang1;
 
      // target struggling step #11	
      if(z_av_1[0]<8.0e-2){ 
	holiang1=  par_k[2] - hrs_ang;
	delta_pk[0] =-1.31749 * sin(-4.61513* holiang1) + 2.03687;
	
      } 
      else{
 	holiang1=  par_k[2] - hrs_ang;
	delta_pk[0] = 3.158e-2 * holiang1 + 4.05819e-1;	
      }
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; // kaon momentum at the reaction point
      
      par_k[0] = pk_real[0];     
           
      
      // missing mass calculation==============================      
     
      hallap = hallap - 0.1843 ;// must be -ve
      hallap = hallap/1000.0; // MeV-->GeV
     
      mm = CalcMM(hallap, par_ep, par_k, mp);     
      mm = (mm)*1000.; // MeV--->GeV
      h->Fill(mm);
      h_lo->Fill(mm);

      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      
      
      R_XFP  = R_XFP*XFPr +XFPm ; 
      R_XpFP = R_XpFP*XpFPr+XpFPm;
      R_YFP  = R_YFP*YFPr+ YFPm;
      R_YpFP = R_YpFP*YpFPr +YpFPm; 
      
      z_av[0] =z_av[0]*Ztr + Ztm;
    
      tnew->Fill();
      
      bool lambdaflag=false;  
      int peak_with_hit= -1; 
      for(int j=0; j<npeak; j++){
	if(Lambda_cent[j]-Lambda_width[j]<mm
	   &&mm < Lambda_cent[j]+Lambda_width[j]){	 
	  
	  lambdaflag=true;
	  peak_with_hit=j; 
	  h_1 ->Fill(mm);
	  h_1 ->SetLineColor(j+2);
	  
	}
	else lambdaflag=false;  

	if(ntune_event<nmax && lambdaflag==true){
	  foil_flag[ntune_event] = peak_with_hit;
	  
	  p10[ntune_event]  = par_ep[0];
	  p11[ntune_event]  = par_ep[1];
	  p12[ntune_event]  = par_ep[2];
	  p13[ntune_event]  = par_k[0];
	  p14[ntune_event]  = par_k[1];
	  p15[ntune_event]  = par_k[2];
	  p16[ntune_event]  = hallap;
	    
	  // x[ntune_event]  = R_XFP; ////RHRS
	  // y[ntune_event]  = R_YFP;
	  // xp[ntune_event] = R_XpFP;
	  // yp[ntune_event] = R_YpFP;

	  x[ntune_event]  = XFP; ////LHRS
	  y[ntune_event]  = YFP;
	  xp[ntune_event] = XpFP;
	  yp[ntune_event] = YpFP;

	  z_recon[ntune_event] = z_av_1[0];
	  phir[ntune_event] =ph1[0];
	  phil[ntune_event] =ph2[0];
	 
	  ntune_event++;
	    
	}
	  
      }//int j	
    }
    
  }
  
  tnew->Write();
  // ((((((((((((((((((((((((((((((((((((((((( t2 ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  for (int i=0 ; i< ent_2 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_2[j]  = -2222.0;
      l_th_fp_2[j] = -2222.0; 
      l_y_fp_2[j]  = -2222.0;
      l_ph_fp_2[j] = -2222.0;
      th1_2[j] = -2222.0;
      th2_2[j] = -2222.0;
      ph1_2[j] =-2222.0;
      ph2_2[j] =-2222.0;    
     
      delta_pep_2[j]= -2222.0;
      pep_real_2[j] =-2222.0;
      delta_pk_2[j]= -2222.0;
      pk_real_2[j] = -2222.0;
     
      r_x_fp_2[j]  = -2222.0;
      r_th_fp_2[j] = -2222.0;
      r_y_fp_2[j]  = -2222.0;
      r_ph_fp_2[j] = -2222.0;
    
      
      
      trig5_2[j] = 0.0;
      rtrig_2 = false;
    }
   
    trig5_2[0] = 0.0;
    rtrig_2 = false;
   
    t2->GetEntry(i);
   
   
    if(trig5_2[0]>1.0) rtrig_2 = true; //JUly 01, 2019
    else rtrig_2 = false;

    z_av_2[0] = (lvz_2[0] + rvz_2[0])/2.0;
    z_av_1_2[0] =  z_av_2[0];
   
    XFP_2   = l_x_fp_2[0];
    XpFP_2  = l_th_fp_2[0];
    YFP_2   = l_y_fp_2[0];
    YpFP_2  = l_ph_fp_2[0];
    
    R_XFP_2   = r_x_fp_2[0];
    R_XpFP_2  = r_th_fp_2[0];
    R_YFP_2   = r_y_fp_2[0];
    R_YpFP_2  = r_ph_fp_2[0];
    
       
    XFP_2  = (XFP_2-XFPm)/XFPr;
    XpFP_2 = (XpFP_2-XpFPm)/XpFPr;
    YFP_2  = (YFP_2-YFPm)/YFPr;
    YpFP_2 = (YpFP_2-YpFPm)/YpFPr;
    
    R_XFP_2  = (R_XFP_2-XFPm)/XFPr; 
    R_XpFP_2 = (R_XpFP_2-XpFPm)/XpFPr;
    R_YFP_2  = (R_YFP_2-YFPm)/YFPr;
    R_YpFP_2 = (R_YpFP_2-YpFPm)/YpFPr;

    /// ================= for Al event in tune ==== jan 07, 2020

    XFP_4  = XFP_2; 
    XpFP_4 = XpFP_2;
    YFP_4  = YFP_2; 
    YpFP_4 = YpFP_2;

    R_XFP_4  =R_XFP_2;
    R_XpFP_4 = R_XpFP_2;
    R_YFP_4  =R_YFP_2;
    R_YpFP_4 =R_YpFP_2;


    z_av_2[0] =(z_av_2[0]- Ztm)/Ztr;
    
    th2_2[0] = calcf2t_th(Theta_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_1_2[0]);
    th2_2[0] = th2_2[0]*Xptr + Xptm; 
    
    ph2_2[0] = calcf2t_ph(PHI_L, XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_1_2[0] );
    ph2_2[0] = ph2_2[0]*Yptr + Yptm;   
    
    
    
    momL_2[0] =  calcf2t_mom(mom_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2[0]); 
    momL_2[0] =momL_2[0]*Momr + Momm; 
    
    momL_2[0] = momL_2[0]*2.21819/2.10; //  Original for H/T and tritium only 
    
    
    par_ep_2[1] = -th2_2[0]; // Dec4, 2019
    par_ep_2[2] = -ph2_2[0];
    
    double holiang2;
    // Target struggling LHRS step #7
    if( z_av_1_2[0]<8.0e-2){
      holiang2 = par_ep_2[2] + hrs_ang;	
      holiang2 = - holiang2;
      
      delta_pep_2[0] = -1.35758*sin(-4.59571* holiang2) + 2.09093;
      
    } 
    else{
      holiang2 = par_ep_2[2] + hrs_ang;
      holiang2 = - holiang2;
      delta_pep_2[0] = 6.23409e-3* holiang2 + 4.03363e-1;
      
    } 
    
    pep_real_2[0] = momL_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
    par_ep_2[0] = pep_real_2[0] ; 
        
    // RHRS angle and momentum calculation      
    th1_2[0] = calcf2t_th(Theta_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,   z_av_2[0]);
    th1_2[0] = th1_2[0]*Xptr + Xptm;
    
    ph1_2[0] = calcf2t_ph(PHI_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2, z_av_2[0]);
    ph1_2[0] = ph1_2[0]*Yptr + Yptm;    
    
    momR_2[0] =  calcf2t_mom(mom_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,  z_av_2[0]);
    momR_2[0] = momR_2[0]*Momr+Momm;    
    
    par_k_2[1] = -th1_2[0]; // Dec2, 2019
    par_k_2[2] = -ph1_2[0];  
    
    double holiang3;
    // target struggling step #11	
    if(z_av_1_2[0]<8.0e-2){ 
      holiang3 = par_k_2[2] - hrs_ang;
      holiang3 = holiang3;
      delta_pk_2[0] =-1.31749*sin(-4.61513*holiang3) + 2.03687;      
    } 
    else{
      holiang3 = par_k_2[2] - hrs_ang;
      holiang3 = holiang3;
      delta_pk_2[0] = 3.158e-2*holiang3 + 4.05819e-1;
    }
    pk_real_2[0] = momR_2[0] + delta_pk_2[0]/1000.0; // kaon momentum at the reaction point     
    par_k_2[0] = pk_real_2[0];    
    
    // missing mass calculation==============================    
    hallap_2 = hallap_2-0.1843 ;// must be -ve
    hallap_2 = hallap_2/1000.0; // MeV-->GeV
    
    mm_2 = CalcMM(hallap_2, par_ep_2, par_k_2, mp);     
    mm_2 = (mm_2)*1000.; // MeV--->GeV

    bool HTflag = false;   
    if(rtrig_2==true &&  fabs(ctime_2)<1.0  && fabs(lvz_2[0]-rvz_2[0])<0.040){ 
      HTflag = true;
    }
    else  HTflag = false;    
       
    if(HTflag == true  && fabs(z_av_1_2[0])<0.10){
      
      h_2->Fill(mm_2);
      h_2lo->Fill(mm_2);
     
      XFP_2 = XFP_2 * XFPr + XFPm;
      XpFP_2 = XpFP_2 * XpFPr + XpFPm;
      YFP_2 = YFP_2 * YFPr + YFPm;
      YpFP_2 = YpFP_2 * YpFPr + YpFPm;
      
      
      R_XFP_2  = R_XFP_2*XFPr +XFPm ; 
      R_XpFP_2 = R_XpFP_2*XpFPr+XpFPm;
      R_YFP_2  = R_YFP_2*YFPr+ YFPm;
      R_YpFP_2 = R_YpFP_2*YpFPr +YpFPm; 
      
      z_av_2[0] =z_av_2[0]*Ztr + Ztm;
      
      tnew->Fill();
      
      bool lambdaflag_2=false;  
      int peak_with_hit_2= -1; 
      for(int j=0; j<npeak_2; j++){
	if(Lambda_cent_2[j]-Lambda_width_2[j]<mm_2
	   &&mm_2 < Lambda_cent_2[j]+Lambda_width_2[j]){	 
	  
	  lambdaflag_2=true;
	  peak_with_hit_2=j; 
	 
	  
	}
	else lambdaflag_2=false;  

	if(ntune_event_2<nmax_2 && lambdaflag_2==true){
	  foil_flag_2[ntune_event_2] = peak_with_hit_2;
	  
	  p10_2[ntune_event_2]  = par_ep_2[0];
	  p11_2[ntune_event_2]  = par_ep_2[1];
	  p12_2[ntune_event_2]  = par_ep_2[2];
	  p13_2[ntune_event_2]  = par_k_2[0];
	  p14_2[ntune_event_2]  = par_k_2[1];
	  p15_2[ntune_event_2]  = par_k_2[2];
	  p16_2[ntune_event_2]  = hallap_2;
	    
	  // x_2[ntune_event_2]  = R_XFP_2; ////RHRS
	  // y_2[ntune_event_2]  = R_YFP_2;
	  // xp_2[ntune_event_2] = R_XpFP_2;
	  // yp_2[ntune_event_2] = R_YpFP_2;

	  x_2[ntune_event_2]  = XFP_2; ////LHRS
	  y_2[ntune_event_2]  = YFP_2;
	  xp_2[ntune_event_2] = XpFP_2;
	  yp_2[ntune_event_2] = YpFP_2;

	  z_recon_2[ntune_event_2] = z_av_1_2[0];
	  phir_2[ntune_event_2] =ph1_2[0];
	  phil_2[ntune_event_2] =ph2_2[0];
	 
	  ntune_event_2++;
	    
	}
	  
      }//int j	
    }
    /// ==================== to include the Al/HT evenys for tune========== jan 07 2020    
    mm_4 = CalcMM(hallap_2, par_ep_2, par_k_2, m_Al); 
    mm_4 = (mm_4)*1000.0;
    mm_ht = mm_4 - 25.3123*1000;  

    bool htflag1 = false;   
    if(rtrig_2==true &&  fabs(ctime_2)<1.0  && fabs(lvz_2[0]-rvz_2[0])<0.040){ 
      htflag1 = true;
    }
    else  htflag1 = false;  
  
   
    bool HTallflag = false;
    if(a1_2<120.0 && a2_2>1650.0 && a2_2<6800.0 &&  //chnged  Jan 11 2020   
       ((z_av_1_2[0] > -0.14 && z_av_1_2[0]< -0.11)||(z_av_1_2[0] > 0.11 && z_av_1_2[0]< 0.14))){
      HTallflag = true;
    }
    else HTallflag = false; 

    ///// AL background analysis
    bool hbgflag = false;
    if(rtrig_2==true &&((ctime_2>-49.39 && ctime_2 < -9.06)||(ctime_2>13.18 && ctime_2 < 48.6))&& fabs(lvz_2[0]-rvz_2[0])<0.040){ 
      hbgflag = true;
    }
    else  hbgflag = false;
    if(HTallflag == true && hbgflag == true)
      {
	hb-> Fill(mm_ht);
	hb5-> Fill(mm_ht);
	hb10-> Fill(mm_ht);
	h20_b-> Fill(mm_ht);
      } 

    ///// Real spectrum
    if(htflag1 == true  && HTallflag == true /*&& ((mm_ht> -5.3 && mm_ht <-1.1)||(mm_ht> 19.5 && mm_ht <32.4))*/){
     
      h21->Fill(mm_ht);      
      hh->Fill(mm_ht);   
      he->Fill(mm_ht); 
      hl->Fill(mm_ht); 

      XFP_4 = XFP_4 * XFPr + XFPm;
      XpFP_4 = XpFP_4 * XpFPr + XpFPm;
      YFP_4 = YFP_4 * YFPr + YFPm;
      YpFP_4 = YpFP_4 * YpFPr + YpFPm;
      
      
      R_XFP_4  = R_XFP_4*XFPr +XFPm ; 
      R_XpFP_4 = R_XpFP_4*XpFPr+XpFPm;
      R_YFP_4  = R_YFP_4*YFPr+ YFPm;
      R_YpFP_4 = R_YpFP_4*YpFPr +YpFPm;    
      
      tnew->Fill(); 

      bool lambdaflag_4=false;  
      int peak_with_hit_4= -1; 
      for(int j=0; j<npeak_4; j++){
	if(Lambda_cent_4[j]-Lambda_width_4[j]<mm_ht
	   &&mm_ht < Lambda_cent_4[j]+Lambda_width_4[j]){	 
	  
	  lambdaflag_4=true;
	  peak_with_hit_4=j; 
	}
	else lambdaflag_4=false;  

	if(ntune_event_4<nmax_4 && lambdaflag_4==true){
	  foil_flag_4[ntune_event_4] = peak_with_hit_4;
	  
	  p10_4[ntune_event_4]  = par_ep_2[0]; // right side should be _2
	  p11_4[ntune_event_4]  = par_ep_2[1];
	  p12_4[ntune_event_4]  = par_ep_2[2];
	  p13_4[ntune_event_4]  = par_k_2[0];
	  p14_4[ntune_event_4]  = par_k_2[1];
	  p15_4[ntune_event_4]  = par_k_2[2];
	  p16_4[ntune_event_4]  = hallap_2;
	    
	  // x_4[ntune_event_4]  = R_XFP_4; ////RHRS
	  // y_4[ntune_event_4]  = R_YFP_4;
	  // xp_4[ntune_event_4] = R_XpFP_4;
	  // yp_4[ntune_event_4] = R_YpFP_4;

	  x_4[ntune_event_4]  = XFP_4; ////LHRS
	  y_4[ntune_event_4]  = YFP_4;
	  xp_4[ntune_event_4] = XpFP_4;
	  yp_4[ntune_event_4] = YpFP_4;

	  z_recon_4[ntune_event_4] = z_av_1_2[0];
	  phir_4[ntune_event_4] =ph1_2[0];
	  phil_4[ntune_event_4] =ph2_2[0];
	  ntune_event_4++;
	   
	}
	
	
	
      }
      
      
    }//  j=0; j< npeak_4
  }
  
  tnew->Write();
  
  // ))))))))))))))))))))))))))))))))))))))))) t2 ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((( t3 ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  bool rtrig_3 = false; 
  for(int i=0; i<nmax_3; i++){
    x_3[i]    = -2222.0; 
    y_3[i]    = -2222.0; 
    xp_3[i]   = -2222.0;
    yp_3[i]   = -2222.0;
    z_av_3[i] = -2222.0;
    z_av_1_3[i] = -2222.0;
    phir_3[i] = -2222.0;
    phil_3[i] = -2222.0;
    z_recon_3[i] = -2222.0;  //// Jan 04, 2019  
    
    foil_flag_3[i] = -1;
  }
  
  // ///////////////////////////////////////////////////////////////////////////////////
  for (int i=0 ; i< ent_3 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_3[j]  = -2222.0;
      l_th_fp_3[j] = -2222.0; 
      l_y_fp_3[j]  = -2222.0;
      l_ph_fp_3[j] = -2222.0;
      th1_3[j] = -2222.0;
      th2_3[j] = -2222.0;
      ph1_3[j] =-2222.0;
      ph2_3[j] =-2222.0;    
      
      delta_pep_3[j]= -2222.0;
      pep_real_3[j] =-2222.0;
      delta_pk_3[j]= -2222.0;
      pk_real_3[j] = -2222.0;
      
      r_x_fp_3[j]  = -2222.0;
      r_th_fp_3[j] = -2222.0;
      r_y_fp_3[j]  = -2222.0;
      r_ph_fp_3[j] = -2222.0;
      
      
      
      trig5_3[j] = 0.0;
      rtrig_3 = false;
    }
    
    trig5_3[0] = 0.0;
    rtrig_3 = false;
    
    t3->GetEntry(i);

    if(trig5_3[0]>1.0) rtrig_3 = true; //JUly 01, 2019
    else rtrig_3 = false;
    
    z_av_3[0] = (lvz_3[0] + rvz_3[0])/2.0;
    z_av_1_3[0] =  z_av_3[0];
    
    XFP_3   = l_x_fp_3[0];
    XpFP_3  = l_th_fp_3[0];
    YFP_3   = l_y_fp_3[0];
    YpFP_3  = l_ph_fp_3[0];
    
    R_XFP_3   = r_x_fp_3[0];
    R_XpFP_3  = r_th_fp_3[0];
    R_YFP_3   = r_y_fp_3[0];
    R_YpFP_3  = r_ph_fp_3[0];
    
   
    XFP_3  = (XFP_3-XFPm)/XFPr;
    XpFP_3 = (XpFP_3-XpFPm)/XpFPr;
    YFP_3  = (YFP_3-YFPm)/YFPr;
    YpFP_3 = (YpFP_3-YpFPm)/YpFPr;
      
    R_XFP_3  = (R_XFP_3-XFPm)/XFPr; 
    R_XpFP_3 = (R_XpFP_3-XpFPm)/XpFPr;
    R_YFP_3  = (R_YFP_3-YFPm)/YFPr;
    R_YpFP_3 = (R_YpFP_3-YpFPm)/YpFPr;
      
    z_av_3[0] =(z_av_3[0]- Ztm)/Ztr;
      
    th2_3[0] = calcf2t_th(Theta_L, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_1_3[0]);
    th2_3[0] = th2_3[0]*Xptr + Xptm;     
    ph2_3[0] = calcf2t_ph(PHI_L, XFP_3, XpFP_3, YFP_3, YpFP_3, z_av_1_3[0] );
    ph2_3[0] = ph2_3[0]*Yptr + Yptm; 
     
    momL_3[0] =  calcf2t_mom(mom_L, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_3[0]); 
    momL_3[0] =momL_3[0]*Momr + Momm;     
    momL_3[0] = momL_3[0]*2.21819/2.10; //  Original for H/T and tritium only 
     

    par_ep_3[1] = -th2_3[0]; // Dec4, 2019
    par_ep_3[2] = -ph2_3[0];

    double holiang5;
    // Target struggling LHRS step #7
    if( z_av_1_3[0]<8.0e-2){
      holiang5 = par_ep_3[2] + hrs_ang;	
      holiang5 = - holiang5;

      delta_pep_3[0] = -1.35758*sin(-4.59571* holiang5) + 2.09093;

    } 
    else{
      holiang5 = par_ep_3[2] + hrs_ang;
      holiang5 = - holiang5;
      delta_pep_3[0] = 6.23409e-3* holiang5 + 4.03363e-1;
	
    } 
      
    pep_real_3[0] = momL_3[0] + delta_pep_3[0]/1000.0; //LHRS  momentum at the reaction point in GeV
      
    par_ep_3[0] = pep_real_3[0] ; 
      

    // RHRS angle and momentum calculation      
    th1_3[0] = calcf2t_th(Theta_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3,   z_av_3[0]);
    th1_3[0] = th1_3[0]*Xptr + Xptm;
    
    ph1_3[0] = calcf2t_ph(PHI_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3, z_av_3[0]);
    ph1_3[0] = ph1_3[0]*Yptr + Yptm;   
           
    momR_3[0] =  calcf2t_mom(mom_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3,  z_av_3[0]);
    momR_3[0] = momR_3[0]*Momr+Momm;
      

    par_k_3[1] = -th1_3[0]; // Dec2, 2019
    par_k_3[2] = -ph1_3[0];
     
    double holiang6;
    // target struggling step #11	
    if(z_av_1_3[0]<8.0e-2){ 
      holiang6 = par_k_3[2] - hrs_ang;
      holiang6 = holiang6;
      delta_pk_3[0] =-1.31749*sin(-4.61513*holiang6) + 2.03687;
	
    } 
    else{
      holiang6 = par_k_3[2] - hrs_ang;
      holiang6 = holiang6;
      delta_pk_3[0] = 3.158e-2*holiang6 + 4.05819e-1;
    }
    pk_real_3[0] = momR_3[0] + delta_pk_3[0]/1000.0; // kaon momentum at the reaction point       
     
    par_k_3[0] = pk_real_3[0];

    // missing mass calculation==============================
    hallap_3 = hallap_3-0.1843 ;// must be -ve
    hallap_3 = hallap_3/1000.0; // MeV-->GeV

    mm_h = CalcMM(hallap_3, par_ep_3, par_k_3, mp); //// to see hydrogen in tritium data
    mm_h = (mm_h)*1000.;
    mm_h =  mm_h -1115.683;

    mm_3 = CalcMM(hallap_3, par_ep_3, par_k_3, m_T); 
    mm_3 = (mm_3)*1000.; // GeV--->MeV
    mm_t = mm_3 -2994.814; // for tritium target only By TOSHI when consider the tritium mass

    mm_Al = CalcMM(hallap_3, par_ep_3, par_k_3, m_Al);
    mm_Al = (mm_Al)*1000.0; 
    mm_Al1 = mm_Al -25.3123*1000; // for Al  Al kinematics only. when consider Al as target
   
  
    //// ========================= from here is tritium data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    bool Tritium_flag = false;
    if(rtrig_3==true && fabs(lvz_3[0]-rvz_3[0])<0.053 &&
		       a1_3<160.0 && a2_3>1685.0 && a2_3<8000.0 && fabs(z_av_1_3[0])<0.10){//chnged  Jan 11 2020  
      Tritium_flag = true;
    }
    else Tritium_flag = false;

    //// for nnl real analysis
    if(Tritium_flag == true && fabs(ctime_3)<1.0){
      h50->Fill(mm_t); // For nnL spectrum
      h52->Fill(mm_t);
      h53->Fill(mm_t); // For nnL spectrum
      h54->Fill(mm_t);
      h_h->Fill(mm_h);

      h53_t->Fill(mm_t); 
   
      h1_2->Fill(mm_t);
      h1_qu->Fill(mm_t); 
 
      // h1_3->Fill(mm_t); 
      // h1_35->Fill(mm_t);
      // h1_4->Fill(mm_t);
      // h1_45->Fill(mm_t);

      // h1_55->Fill(mm_t); 
      // h1_6->Fill(mm_t); 
      // h1_65->Fill(mm_t);
      // h1_7->Fill(mm_t);
      // h1_75->Fill(mm_t);

    }
    
    //// for nnl backkground analysis
    bool bg_nnlflag = false;
    if((ctime_3>-49.39 && ctime_3 < -9.06)||(ctime_3> 13.18 && ctime_3 < 48.6)){//chnged  Jan 11 2020  
      bg_nnlflag = true;
    }
    else bg_nnlflag = false;
   
    if(Tritium_flag == true && bg_nnlflag ==true){
      hb53->Fill(mm_t);
      hb54->Fill(mm_t);
      h50b->Fill(mm_t);
     
      h53_tb->Fill(mm_t);
      h1_2b->Fill(mm_t);
      h1_qub->Fill(mm_t);
    }

    ////// ========================= upt o here is tritium data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    // //// ========================= from here is Al data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    ////// for tune a1 nd a2 = 120 1650, 6800& 0.04 and to see the spectrum a1 , a2 = 160, 1585, 8000 & 0.53
   
    bool TTallflag = false;
    if(rtrig_3==true && fabs(lvz_3[0]-rvz_3[0])<0.040 && a1_3<120.0
       && a2_3>1650.0 && a2_3<6800.0 &&
       ((z_av_1_3[0]>-0.14 && z_av_1_3[0]<-0.11) ||(z_av_1_3[0]>0.11  && z_av_1_3[0]<0.14))){//chnged  Jan 11 2020  
      TTallflag = true;
    }
    else TTallflag = false;

    ///// for Al background analysis ======================================================================
    bool bg1_flag = false;
    if((ctime_3> -49.39 && ctime_3 < -9.06) ||(ctime_3> 13.18 && ctime_3 < 48.6)){
      bg1_flag = true;
    }
    else bg1_flag = false;
   
    if(TTallflag == true && bg1_flag  == true)
      {
	hb1->Fill(mm_Al1);
	hb6->Fill(mm_Al1);
	hb11->Fill(mm_Al1);
	h21_b->Fill(mm_Al1);
      } 


    /// real AL spectrum 
    if(TTallflag == true && fabs(ctime_3) < 1.0 /*&&((mm_Al1 > -5.3 &&mm_Al1<-1.1)||(mm_Al1 >19.5 &&mm_Al1<32.4))*/){// to plot  Al  spectrum
    
      
      //   h50->Fill(mm_t);
      ///  h5->Fill(mm_t,mm_Al1);
      h20->Fill(mm_Al1);
      ht->Fill(mm_Al1);
      hd->Fill(mm_Al1);
      hm->Fill(mm_Al1);
      // h51->Fill(mm_t);

      // h5 -> Fill(mm_Al1, mm_t);  // 2 D
      
      XFP_3 = XFP_3 * XFPr + XFPm;
      XpFP_3 = XpFP_3 * XpFPr + XpFPm;
      YFP_3 = YFP_3 * YFPr + YFPm;
      YpFP_3 = YpFP_3 * YpFPr + YpFPm;
      
      
      R_XFP_3  = R_XFP_3*XFPr +XFPm ; 
      R_XpFP_3 = R_XpFP_3*XpFPr+XpFPm;
      R_YFP_3  = R_YFP_3*YFPr+ YFPm;
      R_YpFP_3 = R_YpFP_3*YpFPr +YpFPm; 
      
      z_av_3[0] =z_av_3[0]*Ztr + Ztm;
      
      tnew->Fill();
         
      bool lambdaflag_3=false;  // need adjustment for tune Jan_02
      int peak_with_hit_3= -1; 
      for(int j=0; j<npeak_3; j++){
  	if(Lambda_cent_3[j]-Lambda_width_3[j]<mm_Al1  // mm_3 need to be adjusted for event selection
  	   &&mm_Al1 < Lambda_cent_3[j]+Lambda_width_3[j]){	 
	  
  	  lambdaflag_3=true;
  	  peak_with_hit_3=j; 
	  
  	}
  	else lambdaflag_3=false;  
	
  	if(ntune_event_3<nmax_3 && lambdaflag_3==true){
  	  foil_flag_3[ntune_event_3] = peak_with_hit_3;
	  
  	  p10_3[ntune_event_3]  = par_ep_3[0];
  	  p11_3[ntune_event_3]  = par_ep_3[1];
  	  p12_3[ntune_event_3]  = par_ep_3[2];
  	  p13_3[ntune_event_3]  = par_k_3[0];
  	  p14_3[ntune_event_3]  = par_k_3[1];
  	  p15_3[ntune_event_3]  = par_k_3[2];
  	  p16_3[ntune_event_3]  = hallap_3;
	    
  	  // x_3[ntune_event_3]  = R_XFP_3; ////RHRS
  	  // y_3[ntune_event_3]  = R_YFP_3;
  	  // xp_3[ntune_event_3] = R_XpFP_3;
  	  // yp_3[ntune_event_3] = R_YpFP_3;

  	  x_3[ntune_event_3]  = XFP_3; ////LHRS
  	  y_3[ntune_event_3]  = YFP_3;
  	  xp_3[ntune_event_3] = XpFP_3;
  	  yp_3[ntune_event_3] = YpFP_3;

  	  z_recon_3[ntune_event_3] = z_av_1_3[0];
  	  phir_3[ntune_event_3] =ph1_3[0];
  	  phil_3[ntune_event_3] =ph2_3[0];
	 
  	  ntune_event_3++;
	    
  	}
	  
      }//int j	
    }
    
  }
  tnew->Write();
  
  // ))))))))))))))))))))))))))))))))))))))))) t3 ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
 
  // fnew->Close(); 
  
  // =================================== 
  // ======== Draw histograms ========== 
  // ===================================
  
  // // HH data
  TF1 *f1 = new TF1("f1","gaus",1112.38,1118.93);
  TF1 *f2 = new TF1("f2","gaus",1189.71,1195.42);
  /*
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h->Draw();
  f1->SetLineWidth(1);
  f2->SetLineWidth(1);
  h->Fit("f1","","",1112.38,1118.93);
  h->Fit("f2","","",1189.71,1195.42);
  f1->Draw("same"); 
    
  TLatex l;
  l.SetTextSize(0.025);
  l.DrawLatex(1125,60,Form("#Lambda"));
 
  l.DrawLatex(1125,70,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  l.DrawLatex(1125,80,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  l.DrawLatex(1200,20,Form("#Sigma^{0}"));
  l.DrawLatex(1200,30,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  l.DrawLatex(1200,40,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));

  /// lorenzian function to the H/H data
  TF1 *f_lo = new TF1("f_lo", cauchy,1110.0,1152.0,3);
  TF1 *f_lo1 = new TF1("f_lo1", cauchy,1188.55,1207.0,3);
  f_lo->SetLineColor(kRed);
  f_lo->SetParameter(0,220); // These are the parameters given  may not be the exact but approx
  f_lo->SetParameter(1,1115.68);
  f_lo->SetParameter(2,1.5);
  f_lo->SetLineWidth(1);
  
  f_lo1->SetLineColor(kRed);
  f_lo1->SetParameter(0,70);
  f_lo1->SetParameter(1,1192.64);
  f_lo1->SetParameter(2,1.5);
  f_lo1->SetLineWidth(1);



  //  f_lo->SetParameter(2,1115.68);
  /// HH data with the Cauchy or the Lorenzian function
  TCanvas* c_lo = new TCanvas("c_lo","c_lo",600,600);
  c_lo->cd();
  h_lo->Draw();
  h_lo->Fit("f_lo","R+","ep");
  h_lo->Fit("f_lo1","R+","ep");
  f_lo->Draw("same"); 

  TLatex l_lo;
  l_lo.SetTextSize(0.025);
  l_lo.DrawLatex(1125,70,Form("#Lambda"));
  //l_lo.DrawLatex(1125,70,Form("#color[2]{#sigma = %.6g}",f_lo->GetParameter(0)));
  l_lo.DrawLatex(1125,80,Form("#color[2]{#sigma = %.6g}",f_lo->GetParameter(1)));
  l_lo.DrawLatex(1125,90,Form("#color[2]{mean = %.6g}",f_lo->GetParameter(2)));
  l_lo.DrawLatex(1200,20,Form("#Sigma^{0}"));
  l_lo.DrawLatex(1200,30,Form("#color[2]{#sigma = %.6g}",f_lo1->GetParameter(1)));
  l_lo.DrawLatex(1200,40,Form("#color[2]{mean = %.6g}",f_lo1->GetParameter(2)));
  
  ////  For H data with T kinematics 
  TF1 *f1_2 = new TF1("f1_2","gaus",1112.2,1119.14);
  TCanvas* c2_2 = new TCanvas("c2_2","c2_2",600,600);
  c2_2->cd();
  h_2->Draw();  
  f1_2->SetLineWidth(1);
  h_2->Fit("f1_2","","",1112.2,1119.14);

  TLatex l2;
  l2.SetTextSize(0.025);
  l2.DrawLatex(1140,50,Form("#Lambda"));  
  l2.DrawLatex(1140,55,Form("#color[2]{#sigma = %.6g}",f1_2->GetParameter(2)));
  l2.DrawLatex(1140,60,Form("#color[2]{mean = %.6g}",f1_2->GetParameter(1)));
  
  /// lorenzian function to the H/T data
  TF1 *f_2lo = new TF1("f_2lo", cauchy,1111.5,1140.5,3);
  f_2lo->SetLineColor(kRed);
  f_2lo->SetParameter(0,100);
  f_2lo->SetParameter(1,1115.68);
  f_2lo->SetParameter(2,1.5);
  f_2lo->SetLineWidth(1);

  TCanvas* c_2lo = new TCanvas("c_2lo","c_2lo",600,600);
  c_2lo->cd();
  h_2lo->Draw();
  h_2lo->Fit("f_2lo","R+","ep");

  TLatex l_2lo;
  l_2lo.SetTextSize(0.025);
  l_2lo.DrawLatex(1140,50,Form("#Lambda"));  
  l_2lo.DrawLatex(1140,55,Form("#color[2]{#sigma = %.6g}",f_2lo->GetParameter(1)));
  l_2lo.DrawLatex(1140,60,Form("#color[2]{mean = %.6g}",f_2lo->GetParameter(2)));

  // TCanvas* chh = new TCanvas("chh","chh",600,600);
  // chh->cd();
  // hz->Draw(); 
 
  // TCanvas* cht = new TCanvas("cht","cht",600,600);
  // cht->cd();
  // ht->Draw(); 
  // TH1F *ht1 = (TH1F*)ht->Clone();
  
  // TCanvas* c20 = new TCanvas("c20","c20",600,600);
  // c20->cd();
  // h20->Draw(); 
  // TH1F *htt = (TH1F*)h20->Clone();
  
  // TCanvas* c21 = new TCanvas("c21","c21",600,600);
  // c21->cd();
  // h21->Draw(); 
  // TH1F *ht1 = (TH1F*)h21->Clone();
  // TF1 *ft = new TF1("ft","gaus",-9.69,-5.65);
  TF1 *ft2 = new TF1("ft2","gaus",-4.5686,-2.1805);
  TF1 *f3 = new TF1("f3","gaus",19.4155,22.883);
  TF1 *f4 = new TF1("f4","gaus",29.4463,31.8078);

  HT_b ->Add(h20_b,h21_b,1.0,1.0);
  TCanvas* cT = new TCanvas("cT","cT",600,600);
  cT->cd();
  HT->Add(h20,h21,1.0,1.0);
  HT->Draw();
  
  HT_b->Draw("E2 same"); // background
  HT_b->Scale(1.0/38.0);
  HT_b->SetFillStyle(3002);
  HT_b->SetMarkerStyle(28);
  HT_b->SetMarkerColor(kGreen);
  
  // ft->SetLineWidth(1);
  // ft2->SetLineWidth(1);
  // f3->SetLineWidth(1);
  // f4->SetLineWidth(1);

  // //  HT->Fit("ft","","",-9.69,-5.65);
  // HT->Fit("ft2","","",-4.5686,-2.1805);
  // HT->Fit("f3","","",19.4155,22.883);
  // HT->Fit("f4","","",29.4463,31.8078);
  // // ft->Draw("same"); 
  // ft2->Draw("same"); 
  // f3->Draw("same"); 

  // TLatex lt;
  // lt.SetTextSize(0.025);
 
  // // lt.DrawLatex(-80,14,Form("#color[2]{#sigma = %.6g}",ft->GetParameter(2)));
  // // lt.DrawLatex(-80,16,Form("#color[2]{I. mean = %.6g}",ft->GetParameter(1)));

  // lt.DrawLatex(-80,22,Form("#color[2]{#sigma = %.6g}",ft2->GetParameter(2)));
  // lt.DrawLatex(-80,24,Form("#color[2]{I. mean = %.6g}",ft2->GetParameter(1)));

  // lt.DrawLatex(-80,32,Form("#color[2]{#sigma = %.6g}",f3->GetParameter(2)));
  // lt.DrawLatex(-80,34,Form("#color[2]{III. mean = %.6g}",f3->GetParameter(1)));
 
  // lt.DrawLatex(-80,40,Form("#color[2]{#sigma = %.6g}",f4->GetParameter(2)));
  // lt.DrawLatex(-80,42,Form("#color[2]{IV. mean = %.6g}",f4->GetParameter(1)));

  // TCanvas* ca = new TCanvas("ca","ca",600,600);
  // ca->cd();
  // ha->Draw(); 
  // TH1F *haa = (TH1F*)ha->Clone();
  
  TCanvas* cb2 = new TCanvas("cb2","cb2",600,600);
  cb2->cd();
  hb2->Add(hb,hb1,1.0,1.0);
  hb2->Draw();
  hb2->Scale(1.0/38.0);
  TH1F *hb3 = (TH1F*)hb2->Clone(); 
 
  
  // counts /0.5
  TF1 *f10 = new TF1("f10","gaus",-4.1589,-2.4788);
  // TF1 *fn = new TF1("fn","gaus",2.428,4.843);
  TF1 *f15 = new TF1("f15","gaus",5.1629,7.281);
  TF1 *f20 = new TF1("f20","gaus",19.7162,21.559);
  TF1 *f30 = new TF1("f30","gaus",29.757,31.769);

  TCanvas* c11 = new TCanvas("c11","c11",600,600);
  c11->cd();
  htt->Add(hh,ht,1.0,1.0);
  htt->Draw();
  f10->SetLineWidth(1);
  //  fn->SetLineWidth(1);
  f15->SetLineWidth(1);
  f20->SetLineWidth(1);
  f30->SetLineWidth(1);

  htt->Fit("f10","","",-4.1589,-2.4788);
  //  htt->Fit("fn","","",2.428,4.843);
  htt->Fit("f15","","",5.1629,7.281);
  htt->Fit("f20","","",19.7162,21.559);
  htt->Fit("f30","","",29.6761,31.769);
  f10->Draw("same"); 
  //  fn->Draw("same"); 
  f15->Draw("same"); 
  f20->Draw("same"); 
 
  hb3->Draw("E2 same"); // background
  hb3->SetFillStyle(3002);
  hb3->SetMarkerStyle(28);
  hb3->SetMarkerColor(kGreen);
 
  TLatex lt2;
  lt2.SetTextSize(0.025);
  
  lt2.DrawLatex(-80,20,Form("#color[2]{#sigma = %.6g}",f10->GetParameter(2)));
  lt2.DrawLatex(-80,21.5,Form("#color[2]{I. mean = %.6g}",f10->GetParameter(1)));
  // lt2.DrawLatex(-80,12,Form("#color[2]{#sigma = %.6g}",fn->GetParameter(2)));
  // lt2.DrawLatex(-80,13.5,Form("#color[2]{II. mean = %.6g}",fn->GetParameter(1)));

  lt2.DrawLatex(-80,23,Form("#color[2]{#sigma = %.6g}",f15->GetParameter(2)));
  lt2.DrawLatex(-80,25.5,Form("#color[2]{II. mean = %.6g}",f15->GetParameter(1)));
  
  lt2.DrawLatex(-80,27,Form("#color[2]{#sigma = %.6g}",f20->GetParameter(2)));
  lt2.DrawLatex(-80,29.5,Form("#color[2]{III. mean = %.6g}",f20->GetParameter(1)));
  
  lt2.DrawLatex(-80,31,Form("#color[2]{#sigma = %.6g}",f30->GetParameter(2)));
  lt2.DrawLatex(-80,32.5,Form("#color[2]{IV. mean = %.6g}",f30->GetParameter(1)));
  
  
  // TCanvas* c5 = new TCanvas("c5","c5",600,600);
  // c5->cd();
  // h5->Draw("colz"); 
 
   
  TCanvas* cb6 = new TCanvas("cb6","cb6",600,600);
  cb6->cd();
  hb7->Add(hb5,hb6,1.0,1.0);
  hb7->Draw();
  hb7->Scale(1.0/38.0);
  TH1F *hb8 = (TH1F*)hb7->Clone(); 
 

  TF1 *f5 = new TF1("f5","gaus",5.0445,7.043);
  TF1 *f6 = new TF1("f6","gaus",-4.906,-2.383);
  TF1 *f7 = new TF1("f7","gaus",19.117,21.855);
  TF1 *f8 = new TF1("f8","gaus",29.6706,31.148);
  
  TCanvas* cf = new TCanvas("cf","cf",600,600);
  cf->cd();
  hf->Add(hd,he,1.0,1.0);
  hf->Draw();
  // f5->SetLineWidth(1);
  // f6->SetLineWidth(1);
  // f7->SetLineWidth(1);
  // f8->SetLineWidth(1);
  // hf->Fit("f5","","",5.0445,7.043);
  // hf->Fit("f6","","",-4.906,-2.383);
  // hf->Fit("f7","","",19.117,21.855);
  // hf->Fit("f8","","",29.6706,31.148);
  // f5->Draw("same"); 
  // f6->Draw("same"); 
  // f7->Draw("same"); 

  hb8->Draw("E2 same"); // background
  hb8->SetFillStyle(3002);
  hb8->SetMarkerStyle(28);
  hb8->SetMarkerColor(kGreen);

  // TLatex lt1;
  // lt1.SetTextSize(0.025);
  
  // lt1.DrawLatex(-80,22,Form("#color[2]{#sigma = %.6g}",f5->GetParameter(2)));
  // lt1.DrawLatex(-80,24,Form("#color[2]{II. mean = %.6g}",f5->GetParameter(1)));
  
  // lt1.DrawLatex(-80,18,Form("#color[2]{#sigma = %.6g}",f6->GetParameter(2)));
  // lt1.DrawLatex(-80,20,Form("#color[2]{I. mean = %.6g}",f6->GetParameter(1)));
  
  // lt1.DrawLatex(-80,27,Form("#color[2]{#sigma = %.6g}",f7->GetParameter(2)));
  // lt1.DrawLatex(-80,29,Form("#color[2]{III. mean = %.6g}",f7->GetParameter(1)));
  
  // lt1.DrawLatex(-80,32,Form("#color[2]{#sigma = %.6g}",f8->GetParameter(2)));
  // lt1.DrawLatex(-80,34,Form("#color[2]{IV. mean = %.6g}",f8->GetParameter(1)));





  TCanvas* cb10 = new TCanvas("cb10","cb10",600,600);
  cb10->cd();
  hb12->Add(hb10,hb11,1.0,1.0);
  hb12->Draw();
  hb12->Scale(1.0/38.0);
  TH1F *hb13 = (TH1F*)hb12->Clone(); 


  TCanvas* cn = new TCanvas("cn","cn",600,600);
  cn->cd();
  hn->Add(hl,hm,1.0,1.0);
  hn->Draw();

  hb13->Draw("E2 same");
  hb13->SetFillStyle(3002);
  hb13->SetMarkerStyle(28);
  hb13->SetMarkerColor(kGreen);


  TCanvas* c50b = new TCanvas("c50b","c50b",600,600);
  c50b->cd();
  h50b->Draw(); 
  h50b ->Scale(1.0/38.0);
  TH1F * h50b1 = (TH1F*)h50b->Clone();
  
 
  TCanvas* c50 = new TCanvas("c50","c50",600,600);
  c50->cd();
  h50->Draw();
  h50b1->Draw("E2 same");
  h50b1->SetFillStyle(3002);
  h50b1->SetMarkerStyle(28);
  h50b1->SetMarkerColor(kGreen);

  
  // TCanvas* c51 = new TCanvas("c51","c51",600,600);
  // c51->cd();
  // h51->Draw(); 
  
  TCanvas* c52 = new TCanvas("c52","c52",600,600);
  c52->cd();
  h52->Draw(); 
  //  auto g1 = new TF1("g1","pol1",4.7434,10.4674);

  TF1 * f53 = new TF1("f53","gaus",4.5207,10.379);
  TF1 * f53_l = new TF1("f53_l","gaus",-1.999,3.1918);
  TCanvas* cb53 = new TCanvas("cb53","cb53",600,600);
  cb53->cd();
  hb53->Draw(); 
  hb53 ->Scale(1.0/38.0);
  TH1F * hb55 = (TH1F*)hb53->Clone();
  

  TCanvas* c53 = new TCanvas("c53","c53",600,600);
  c53->cd();
  h53->Draw();
  hb55->Draw("E2 same");
  hb55->SetFillStyle(3002);
  hb55->SetMarkerStyle(28);
  hb55->SetMarkerColor(kGreen);
  f53->SetLineWidth(1);
  f53_l->SetLineWidth(1);
  h53->Fit("f53","","",4.5207,10.379); 
  h53->Fit("f53_l","","",-1.999,3.1918);
  f53->Draw("same");


  TLatex l53;
  l53.SetTextSize(0.025);
  l53.DrawLatex(-80,45,Form("#color[2]{#sigma = %.6g}",f53_l->GetParameter(2)));
  l53.DrawLatex(-80,47,Form("#color[2]{I. mean = %.6g}",f53_l->GetParameter(1)));
  l53.DrawLatex(-80,49,Form("#color[2]{#sigma = %.6g}",f53->GetParameter(2)));
  l53.DrawLatex(-80,51,Form("#color[2]{II.mean = %.6g}",f53->GetParameter(1)));
 */


  /// following 2 are the  temporary
  // h53_tb ->Scale(1.0/38.0);
  // TCanvas* cb53_t = new TCanvas("cb53_t","cb53_t",600,600);
  // cb53_t->cd();
  // h53_t->Draw();
  // h53_tb->Draw("E2 same");
  // h53_tb->SetFillStyle(3002);
  // h53_tb->SetMarkerStyle(28);
  // h53_tb->SetMarkerColor(kGreen);

  TCanvas* cb54 = new TCanvas("cb54","cb54",600,600);
  cb54->cd();
  hb54->Draw();
  hb54->Scale(1.0/38.0);
  TH1F * hb56 = (TH1F*)hb54->Clone(); 
  
  
  TCanvas* c54 = new TCanvas("c54","c54",600,600);
  c54->cd();
  h54->Draw(); 
  hb56->Draw("E2 same");
  hb56->SetFillStyle(3002);
  hb56->SetMarkerStyle(28);
  hb56->SetMarkerColor(kGreen);

  TCanvas *c_h = new TCanvas("c_h","c_h", 600,600);
  c_h->cd();
  h_h->Draw();

  /// to estimate the quasi free shape//////////////////////////////////////////////////////////////////////
  // double x11[67] = {-2.0, 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,  // to estima the quais free shape
  // 		    10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,
  // 		    20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,
  // 		    30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,
  // 		    40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,
  // 		    50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,

  // 		    60.0,65.0,70.0,75.0,80.0,85.0};
  
  // double y11[67] = {0.0,28.0,50.0,75.0,103.0,126.0,152.0,187.0,215.0,
  // 		    246.0,280.0,314.0,354.0,386.0,425.0,462.0,498.0,
  // 		    547.0,590.0,635.0,675.0,716.0,756.0,805.0,852.0,
  // 		    900.0,952.0,1010.0,1051.0,1104.0,1154.0,1206.0,
  // 		    1268.0,1318.0,1369.0,1435.0,1480.0,1534.0,1576.0,
  // 		    1628.0,1672.0,1725.0,1765.0,1810.0,1857.0,1904.0,
  // 		    1940.0,1972.0,2007.0,2038.0,2072.0,2102.0,2125.0,
  // 		    2154.0,2181.0,2204.0,2222.0,2240.0,2250.0,2262.0,
  // 		    2271.0,2280.0,2256.0,2210.0,2120.0,2016.0,1888.0};
  
  //  TF1 *f_q = new TF1("f_q", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) +[6]*pow(x,6) +[7]*pow(x,7)",-2.0,75.0);
 
  //  for(int q=0;q<67;q++) // loop is to scale the y axis by 48 time
  //   {
  //     x11[q] = q;
  //     y11[q]= (1.0/55.0)*y11[q];
  //     //  y11[q]= y11[q] + N_b;
  //   }


  // TGraphErrors *gr = new TGraphErrors(67,x11,y11,0,0);
  // gr->GetXaxis()->SetLimits(-100.,157.);
  // //  gr->GetYaxis()->SetLimits(0.0,55.0);

  TF1 *f_qq = new TF1("f_qq", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) +[6]*pow(x,6) +[7]*pow(x,7)",-2.0,75.0);
  f_qq->SetParameters(31.3472,20.1666,1.11717,-0.0477299,0.00220349, -5.2519e-05,5.42044e-07,-2.03287e-09);
  
  TF1 *f_qq2 = new TF1("f_qq2", "f_qq+7",-2.0,75.0);
  f_qq2->SetParameters(31.3472,20.1666,1.11717,-0.0477299,0.00220349, -5.2519e-05,5.42044e-07,-2.03287e-09);
  //////// up to here is the quasi frre shape calculation /////////////////////////////////////////////////////////// 
  
  TF1 *fit_fun = new TF1("fit_fun",user_gaus,-3.056,2.7702,3);  // -2 to 2 is  fit range and 3 is the # of parameter
  fit_fun->SetParameters(15,0.0,0.9);
  fit_fun->SetNpx(1000);
  
  h1_2b->Scale(1.0/38.0);
  TCanvas *c1_2 = new TCanvas("c1_2","c1_2", 600,600);
  c1_2->cd();
  h1_2->Draw();
  h1_2b->Draw("E2 same");
  //gr->Draw("Ap*");
  // gr->Fit("f_q","R+");
  h1_2b->SetFillStyle(3002);
  h1_2b->SetMarkerStyle(28);
  h1_2b->SetMarkerColor(kGreen);
  fit_fun->SetLineWidth(1);
  h1_2->Fit("fit_fun","","",-3.056,2.7702);
  fit_fun->Draw("same");
  h1_2->Fit("f_qq2","","",-2.0,75.0);

  // gr->Draw("Ap*");
  // gr->Draw("same");
  //  gr->Fit("f_q","R+");
 
  TLatex l52;
  l52.SetTextSize(0.025);
  l52.DrawLatex(-80,48,Form("#color[2]{#sigma = %.6g}",fit_fun->GetParameter(2)));
  l52.DrawLatex(-80,50,Form("#color[2]{mean = %.6g}",fit_fun->GetParameter(1)));
  
 TH1F *h1_22 = new TH1F("h1_22","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);

  TFile* f_new = new TFile("./output_root/quasi_free.root","recreate");
  //  TObjArray h1_22(0);
  h1_22->Add(h1_2);
  h1_22->Add(h1_2b);
  h1_2->Write();
  h1_2b->Write();
  h1_22->Write();
  // f_new->Write();
  h1_qu->Write();
  h1_qub->Write();
  f_new->Close();
  
  c1_2->SaveAs("h1_2.pdf");

 //  TCanvas *c1_3 = new TCanvas("c1_3","c1_3", 600,600);
 //  c1_3->cd();
 //  h1_3->Draw();

 //  TCanvas *c1_35 = new TCanvas("c1_35","c1_35", 600,600);
 //  c1_35->cd();
 //  h1_35->Draw();

 //  TCanvas *c1_4 = new TCanvas("c1_4","c1_4", 600,600);
 //  c1_4->cd();
 //  h1_4->Draw();

 //  TCanvas *c1_45 = new TCanvas("c1_45","c1_45", 600,600);
 //  c1_45->cd();
 //  h1_45->Draw();


 // TCanvas *c1_55 = new TCanvas("c1_55","c1_55", 600,600);
 //  c1_55->cd();
 //  h1_55->Draw();

 //  TCanvas *c1_6 = new TCanvas("c1_6","c1_6", 600,600);
 //  c1_6->cd();
 //  h1_6->Draw();

 //  TCanvas *c1_65 = new TCanvas("c1_65","c1_65", 600,600);
 //  c1_65->cd();
 //  h1_65->Draw();

 //  TCanvas *c1_7 = new TCanvas("c1_7","c1_7", 600,600);
 //  c1_7->cd();
 //  h1_7->Draw();

 //  TCanvas *c1_75 = new TCanvas("c1_75","c1_75", 600,600);
 //  c1_75->cd();
 //  h1_75->Draw();


  ///To save the histNograms
  // c2->SaveAs("T1.pdf");
  // c_lo->SaveAs("T2.pdf");
  // c2_2->SaveAs("T3.pdf");
  // c_2lo->SaveAs("T4.pdf");
  // c11->SaveAs("T5.pdf");
  // cf->SaveAs("T6.pdf");
  // cT->SaveAs("T7.pdf");
  // c_h->SaveAs("T8.pdf");
  // c50->SaveAs("T9.pdf");
  // c53->SaveAs("T10.pdf");
  // cb53_t->SaveAs("T11.pdf");
  // cb53_t1->SaveAs("T12.pdf");
  // c54->SaveAs("T13.pdf");


  ////  vertex
  // TF1 *f26 = new TF1("f26","gaus",-0.13539,-0.11584);
  // TF1 *f27 = new TF1("f27","gaus",0.1144,0.1338);
  // TCanvas* c25 = new TCanvas("c25","c52",600,600);
  // c25->cd();
  // H25->Draw();
  // f26->SetLineWidth(1);
  // f27->SetLineWidth(1);
  // H25->Fit("f26","","",-0.13539,-0.11584);
  // H25->Fit("f27","","",0.1144,0.1338);
  // f26->Draw("same"); 
  ///// TLatex l;
  // l.SetTextSize(0.025);
  // l.DrawLatex(1125,80,Form("#Lambda"));
 
  // l.DrawLatex(1125,85,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  // l.DrawLatex(1125,90,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  // l.DrawLatex(1200,32,Form("#Sigma^{0}"));
  // l.DrawLatex(1200,37,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  // l.DrawLatex(1200,42,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));

   

  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  const int nite =0;
  double temp[nite]; 
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  for(int i=0 ; i<nite ; i++){
    x[i] = i+1;
    //  temp[i] = tune(phiR_opt,i); // jan 31
    temp[i] = tune(momL_opt,i);
    // temp[i] = tune(momR_opt,i);   //jan 31 
    
    sprintf(tempc, "./MOM_MATRICES/LMOM5_51th_%d.dat",i); // jan 31
    //   sprintf(tempc, "./matrices/phiR4_20th_%d.dat",i);// jan 31
    ofstream * ofs = new ofstream(tempc); 
    int nppp = 0;
    const int nn = 5;  // 5 for 5th order 4 for 4th order jan 31
    for(int i=0; i<nn+1; i++){
      for(int e=0; e<nn+1; e++){
	for(int d=0; d<nn+1; d++){ 
	  for(int c=0; c<nn+1; c++){
	    for(int b=0; b<nn+1; b++){
	      for(int a=0; a<nn+1; a++){  
		if(a+b+c+d+e==i){
		  *ofs <<momL_opt[nppp] // jan 31
		       << " " << a 
		       << " " << b
		       << " " << c
		       << " " << d
		       << " " << e << endl;
		  nppp++; 
		  
		}
	      }
	    }
	  }
	}
      }
    }
    ofs->close();
    ofs->clear();
    
    cout << temp[i]<<endl; 
  }    
  
  if(nite>0){
    TGraph * gr = new TGraph(nite,x,temp);  
    TCanvas * c4 = new TCanvas("c4","",600,600); 
    gr->Draw("*la"); 
  }
} //end of  main function

//////////////////////////////////////////////////
double calcf2t_th(double* P, double xf, double xpf, 
		  double yf, double ypf,double zt)
//////////////////////////////////////////////////
{
  // -----4th order -----   
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
  
  return Y; 
}
// ////////////////////////////////////////////////
//////////////////////////////////////////////////
double calcf2t_ph(double* P, double xf, double xpf, 
		  double yf, double ypf, double zt)
//////////////////////////////////////////////////
{
  // -----4th order -----   
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
  
  return Y; 
}


//////////////////////////////////////////////////
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
//////////////////////////////////////////////////
{
  // -----4th order -----   
  const int nMatT=5;  
  const int nXf=5;
  const int nXpf=5;
  const int nYf=5;
  const int nYpf=5;
  const int nZt=5;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
  
  return Y; 
}

// missing mass function definition====================
double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt){ // Dec 3,2019
  
  double pe = ee;
  double Ee = sqrt(me*me + pe*pe);
  Ee = Ee - 0.0003; // GeV
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = pvec_ep[0];
  double xpep = pvec_ep[1];
  double ypep = pvec_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  vec_ep.RotateX(hrs_ang);
  //double Eep = sqrt(vec_ep * vec_ep);
  double Eep = sqrt(pep*pep + me*me);
  
 
  double pk  = pvec_k[0];
  double xpk = pvec_k[1];
  double ypk = pvec_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
  vec_k.RotateX(-hrs_ang);
  //double Ek = sqrt(vec_k * vec_k);
  double Ek = sqrt(pk*pk + mk*mk);
  
 
  double missingE2 = 0.0, missingP2 = 0.0, missingM2 = 0.0;
  missingE2 = pow(Ee + mt - Ek - Eep, 2.0);
  missingP2 = (vec_e - vec_ep - vec_k) * (vec_e - vec_ep - vec_k);
  missingM2 = missingE2 - missingP2;
  
  double MissingMass = 0.0;
  MissingMass = sqrt(missingM2);

  return MissingMass;
  
}

//############### up to hear missing mass #####################

// #############################################################
double tune(double* pa, int j) // tune fun defn
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam =Mom_Par; // for momentum tune jan 31
  //  int allparam = Total_Par; //  for the angle tune  jan 31

  TMinuit* minuit = new TMinuit(allparam); 
  minuit->SetFCN(fcn); // very imp function setying for chi square 
    
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1); 
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
 
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    //  step[i] = 1.0e-3; /// Original 
    step[i] = 2.0*0.5;  // for rough matrix, when matrix is far from reality
    
    LLim[i] = pa[i] -10; // pa[i]*0.8; // KI
    ULim[i] = pa[i] + 10; //pa[i]*0.8; // KI
    // LLim[i] = pa[i] - pa[i]*0.8; // KI
    // ULim[i] = pa[i] + pa[i]*0.8; // KI


    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  //  arglist[0] = 2.0; // was active before
  arglist[0] = 1.0;  // KI
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ one of the way to get optimized parameter
  arglist[0] = 20000;
  arglist[1] = 0.01; // To make more presise
  minuit -> mnexcm("MINImize",arglist,2,ierflg); 
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat); 
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){  
     
    minuit -> GetParameter(i,momL_opt[i],er); // need change while tuning matrices LHRS
    // minuit -> GetParameter(i,momR_opt[i],er); // need change while tuning matrices  RHRS
    // //  minuit -> GetParameter(i,phiR_opt[i],er); // for angle optimization  jan 31
  }
  
  return chi2; 
}


// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
  
{
  double chi2 = 0.0;
  double chi_12 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  const double sigma = 0.0045; 
  double ref_mm = 0.0; 
  double residual = 0.0;

  double par_ep[3];
  double par_k[3];
  double halla_p;
  double momr[100];
  double moml[100];
  double z_av;
  double z_av_sc;
  double rvz;
  double ph1;
  double th1;
  double ph2;
  double th2;
  double  delta_pk[100];
  double delta_pep[100];
  double pk_real[100];
  double pep_real[100];
  double MM;
  double THL;
  double PHL;
  double THR;
  double PHR;
  // (((((((((((((((((((((((((((((((((((((((( t2 (((((((((((((((((((((((((((((((((
 
  double chi_22 = 0.0;
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double ref_mm_2 = 0.0; 
  double residual_2 = 0.0;

  double par_ep_2[3];
  double par_k_2[3];
  double halla_p_2;
  double momr_2[100];
  double moml_2[100];
  double z_av_2;
  double z_av_sc_2;
  double ph1_2;
  double th1_2;
  double ph2_2;
  double th2_2;
  double delta_pk_2[100];
  double delta_pep_2[100];
  double pk_real_2[100];
  double pep_real_2[100];
  double MM_2;
  double THL_2;
  double PHL_2;
  double THR_2;
  double PHR_2;


  //)))))))))))))))))))))))))))))))))))))))))))))   t2  ))))))))))))))))))))))))))))  
  // (((((((((((((((((((((((((((((((((((((((( t3 (((((((((((((((((((((((((((((((((
  double chi_33 = 0.0;
  double XFP_3, XpFP_3;
  double YFP_3, YpFP_3;
  double ref_mm_3 = 0.0; 
  double residual_3 = 0.0;

  double par_ep_3[3];
  double par_k_3[3];
  double halla_p_3;
  double momr_3[100];
  double moml_3[100];
  double z_av_3;
  double z_av_sc_3;

  double ph1_3;
  double th1_3;
  double ph2_3;
  double th2_3;
  double delta_pk_3[100];
  double delta_pep_3[100];
  double pk_real_3[100];
  double pep_real_3[100];
  double MM_3;
  double THL_3;
  double PHL_3;
  double THR_3;
  double PHR_3;

  //)))))))))))))))))))))))))))))))))))))))))))))   t3  ))))))))))))))))))))))))))))
 // (((((((((((((((((((((((((((((((((((((((( to include the Al/HT data for tune Jan 07, 2020 (((((((((((((((((((((((((((((((((
  double chi_44 = 0.0;
  double XFP_4, XpFP_4;
  double YFP_4, YpFP_4;
  double ref_mm_4 = 0.0; 
  double residual_4 = 0.0;

  double par_ep_4[3];
  double par_k_4[3];
  double halla_p_4;
  double momr_4[100];
  double moml_4[100];
  double z_av_4;
  double z_av_sc_4;

  double ph1_4;
  double th1_4;
  double ph2_4;
  double th2_4;
  double delta_pk_4[100];
  double delta_pep_4[100];
  double pk_real_4[100];
  double pep_real_4[100];
  double MM_4;
  double THL_4;
  double PHL_4;
  double THR_4;
  double PHR_4;

  //)))))))))))))))))))))))))))))))))))))))))))))   up here is _4  ))))))))))))))))))))))))))))

 //)))))))))))))))))))))))))))))))))))))))))))))  ))))))))))))))))))))))))))))
  for(int i=0; i<ntune_event; i++){ 
    residual = 0.0;
    ref_mm = 0.0; 
    ref_mm  = Lambda_real[foil_flag[i]];    
    ref_mm = ref_mm/1000.0;
    
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    z_av = z_recon[i]; 
    ph1 = phir[i];  // open when calibrate the Momentum 
    ph1= -ph1;
    ph2 = phil[i];    
    ph2 = -ph2; 

    XFP   =(XFP -XFPm)/XFPr;  
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    z_av_sc = (z_av - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml[0] =  calcf2t_mom(param, XFP, XpFP, YFP, YpFP,  z_av_sc);
    moml[0] = moml[0]*Momr+Momm; 
    double hang;
    if( z_av<8.0e-2){
      hang = ph2 + hrs_ang;
      hang = -hang;
      delta_pep[0] = -1.35758*sin(-4.59571* hang) + 2.09093;
      
    } 
    else{
      hang = ph2 + hrs_ang;
      hang = -hang;
      delta_pep[0] = 6.23409e-3*hang + 4.03363e-1;
    } 
    
    pep_real[0] = moml[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
   
    // THL =  calcf2t_th(param, XFP, XpFP, YFP, YpFP,  z_av);
    // THL = THL*Xptr + Xptm;

    // PHL =  calcf2t_ph(param, XFP, XpFP, YFP, YpFP,  z_av);
    // PHL = PHL*Yptr + Yptm;
    
    par_ep[0] = pep_real[0];// Wwhen LHRS momentum  tuned
    //  par_ep[0] = p10[i];
    par_ep[1] = p11[i]; 
    par_ep[2] = p12[i];
   
    //// par_ep[1] = -THL;  // for theta optimization   jan 31
    ////  par_ep[2] = - PHL;
   
   
    //  for the RHRS momentum tunning
    momr[0] =  calcf2t_mom(param, XFP, XpFP, YFP, YpFP,  z_av_sc);
    momr[0] = momr[0]*Momr+Momm;
    double hang1;    

    if(z_av<8.0e-2){
      hang1= ph1 - hrs_ang;
      delta_pk[0] =-1.31749*sin(-4.61513* hang1) + 2.03687;
      
    } 
    else{
      hang1= ph1 - hrs_ang;
      delta_pk[0] = 3.158e-2*hang1 + 4.05819e-1; 
    }
    pk_real[0] = momr[0] + delta_pk[0]/1000.0; 


    // // THR =  calcf2t_th(param, XFP, XpFP, YFP, YpFP, z_av_sc);
    // // THR = THR*Xptr + Xptm; 
    // // jan 31    
    // // PHR =  calcf2t_th(param, XFP, XpFP, YFP, YpFP, z_av_sc);
    // // PHR = PHR*Yptr + Yptm;  
    
    //  par_k[0] = pk_real[0];// when RHRS matrix tuned
    par_k[0] = p13[i];    
    par_k[1] = p14[i]; // jan 31    
    par_k[2] = p15[i];

    ///   par_k[1] = -THR;// jan 31 
    ///   par_k[2] = -PHR;
    
    halla_p = p16[i];
    MM = CalcMM(halla_p, par_ep, par_k, mp);    
    residual = MM-ref_mm;
   
    //   chi_12 = chi_12 + pow(residual,2.0);
    ///////  if need to use the sigma statistical weigh

    if(foil_flag[i] ==0)
      {chi_12 = chi_12 + pow(residual,2.0);}
    else
      {chi_12 = chi_12 +pow(residual,2.0);}
    
  }
  // (((((((((((((((((((((((((((((((((((((((( t2 ((((((((((((((((((((((((((((((((( 
  for(int i=0; i<ntune_event_2; i++){ 
    residual_2 = 0.0;
    ref_mm_2 = 0.0; 
    ref_mm_2  = Lambda_real_2[foil_flag_2[i]];    
    ref_mm_2 = ref_mm_2/1000.0;
    
    XFP_2   = x_2[i];
    XpFP_2  = xp_2[i];
    YFP_2   = y_2[i];
    YpFP_2  = yp_2[i];
    z_av_2 = z_recon_2[i]; 
    ph1_2 = phir_2[i];  // open when calibrate the Momentum 
    ph1_2 = - ph1_2;
    ph2_2 = phil_2[i];    
    ph2_2 = - ph2_2;

    XFP_2   =(XFP_2 -XFPm)/XFPr;  
    XpFP_2  =(XpFP_2-XpFPm)/XpFPr;
    YFP_2   =(YFP_2 -YFPm)/YFPr;
    YpFP_2  =(YpFP_2-YpFPm)/YpFPr;
    z_av_sc_2 = (z_av_2 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_2[0] =  calcf2t_mom(param, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_sc_2);
    moml_2[0] = moml_2[0]*Momr+Momm;

    moml_2[0] = moml_2[0]*2.21819/2.1; //for H/T and tritium only
    double hang2;   
    if( z_av_2<8.0e-2){

      hang2 =  ph2_2 + hrs_ang;
      hang2 = -hang2;
      delta_pep_2[0] = -1.35758*sin(-4.59571* hang2) + 2.09093;
     
    } 
    else{
      hang2 =  ph2_2 + hrs_ang;
      hang2 = -hang2;
      delta_pep_2[0] = 6.23409e-3*hang2 + 4.03363e-1;
    }    
    pep_real_2[0] = moml_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV 
    
    // THL_2 =  calcf2t_th(param,XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2);
    // THL_2 = THL_2*Xptr + Xptm;
    // PHL_2 =  calcf2t_ph(param,XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2);
    // PHL_2 = PHL_2*Yptr + Yptm;
    
    par_ep_2[0] = pep_real_2[0];// Wwhen LHRS momentum  tuned
    //  par_ep_2[0] = p10_2[i];
    par_ep_2[1] = p11_2[i];
    par_ep_2[2] = p12_2[i];
    
    ////par_ep_2[1] = - THL_2;
    ////par_ep_2[2] = - PHL_2; 
   
    //  for the RHRS momentum tunning
    momr_2[0] =  calcf2t_mom(param, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_sc_2);
    momr_2[0] = momr_2[0]*Momr+Momm;
    double hang3;    
    if(z_av_2<8.0e-2){

      hang3 = ph1_2 - hrs_ang;
      delta_pk_2[0] =-1.31749*sin(-4.61513* hang3) + 2.03687;
     
    } 
    else{
      hang3 = ph1_2 - hrs_ang;
      delta_pk_2[0] = 3.158e-2*hang3 + 4.05819e-1; 
    }
    pk_real_2[0] = momr_2[0] + delta_pk_2[0]/1000.0;

    //// THR_2 =  calcf2t_th(param,XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_sc_2);
    //// THR_2 = THR_2*Xptr + Xptm;   // jan 31
    //// PHR_2 =  calcf2t_th(param,XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_sc_2);
    //// PHR_2 = PHR_2*Yptr + Yptm;
    
    
    //   par_k_2[0] = pk_real_2[0];// when RHRS matrix tuned
    par_k_2[0] = p13_2[i];    
    par_k_2[1] = p14_2[i];   // jan 31 
    par_k_2[2] = p15_2[i];
 
    ////  par_k_2[1] = -THR_2; // jan 31
    ////  par_k_2[2] = -PHR_2; 
    
    halla_p_2 = p16_2[i];
    MM_2 = CalcMM(halla_p_2, par_ep_2, par_k_2, mp);    
    residual_2 = MM_2-ref_mm_2;
   
    chi_22 = chi_22 +pow(residual_2,2.0);
    // chi_22 = chi_22 +pow(residual_2,2.0);

  }
  //)))))))))))))))))))))))))))))))))))))))))))))   t2  ))))))))))))))))))))))))))))
  // (((((((((((((((((((((((((((((((((((((((( t3 ((((((((((((((((((((((((((((((

  for(int i=0; i<ntune_event_3; i++){ 
    residual_3 = 0.0;
    ref_mm_3 = 0.0; 
    ref_mm_3  = Lambda_real_3[foil_flag_3[i]];    
    ref_mm_3 = ref_mm_3/1000.0;
    
    XFP_3   = x_3[i];
    XpFP_3  = xp_3[i];
    YFP_3   = y_3[i];
    YpFP_3  = yp_3[i];
    z_av_3 = z_recon_3[i]; 
    ph1_3 = phir_3[i];  // open when calibrate the Momentum 
    ph1_3 = - ph1_3;
    ph2_3 = phil_3[i];    
    ph2_3 = - ph2_3;

    XFP_3   =(XFP_3 -XFPm)/XFPr;  
    XpFP_3  =(XpFP_3-XpFPm)/XpFPr;
    YFP_3   =(YFP_3 -YFPm)/YFPr;
    YpFP_3  =(YpFP_3-YpFPm)/YpFPr;
    z_av_sc_3 = (z_av_3 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_3[0] =  calcf2t_mom(param, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_sc_3);
    moml_3[0] = moml_3[0]*Momr+Momm;

    moml_3[0] = moml_3[0]*2.21819/2.1; //for H/T and tritium only
    double hang9;   
    if( z_av_3<8.0e-2){

      hang9 =  ph2_3 + hrs_ang;
      hang9 = -hang9;
      delta_pep_3[0] = -1.35758*sin(-4.59571* hang9) + 2.09093;
     
    } 
    else{
      hang9 =  ph2_3 + hrs_ang;
      hang9 = -hang9;
      delta_pep_3[0] = 6.23409e-3*hang9 + 4.03363e-1;
    }    
    pep_real_3[0] = moml_3[0] + delta_pep_3[0]/1000.0; //LHRS  momentum at the reaction point in GeV 

    
    par_ep_3[0] = pep_real_3[0];// Wwhen LHRS momentum  tuned
    // par_ep_3[0] = p10_3[i];
    par_ep_3[1] = p11_3[i];
    par_ep_3[2] = p12_3[i];
   
   
    //  for the RHRS momentum tunning
    momr_3[0] =  calcf2t_mom(param, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_sc_3);
    momr_3[0] = momr_3[0]*Momr+Momm;
    double hang10;    
    if(z_av_3<8.0e-2){
      hang10 = ph1_3 - hrs_ang;
      delta_pk_3[0] =-1.31749*sin(-4.61513* hang10) + 2.03687;
      
    } 
    else{
      hang10 = ph1_3 - hrs_ang;
      delta_pk_3[0] = 3.158e-2*hang10 + 4.05819e-1; 
    }
    pk_real_3[0] = momr_3[0] + delta_pk_3[0]/1000.0;    
    

    // // THR_3 =  calcf2t_th(param,XFP_3, XpFP_3, YFP_3, YpFP_3, z_av_sc_3);
    // // THR_3 = THR_3*Xptr + Xptm;   // jan 31
    // // PHR_3 =  calcf2t_th(param,XFP_3, XpFP_3, YFP_3, YpFP_3, z_av_sc_3);
    // // PHR_3 = PHR_3*Yptr + Yptm;
    


    //   par_k_3[0] = pk_real_3[0];// when RHRS matrix tuned
    par_k_3[0] = p13_3[i];    
    par_k_3[1] = p14_3[i];  // jan 31  
    par_k_3[2] = p15_3[i];
  
    //// par_k_3[1] = - THR_3 ;// jan 31 
    ////par_k_3[2] = -PHR_3; 
    
    halla_p_3 = p16_3[i];
    MM_3 = CalcMM(halla_p_3, par_ep_3, par_k_3, m_Al);  //need to adjust  
    MM_3 = MM_3 -25.3123;

    residual_3 = MM_3-ref_mm_3;
    //  cout<<" TT data = " << MM_3 << " ref_mm_3 = " << ref_mm_3 << " residual_3 = " << residual_3<<endl;  

    chi_33 = chi_33 +3*pow(residual_3,2.0);
  }
  //)))))))))))))))))))))))))))))))))))))))))))))   t3  ))))))))))))))))))))))))))))

 // ((((((((((((((((((((((((4(((((((((((((((( _4   Jan 07, 2020((((((((((((((((((((((((((((((
 for(int i=0; i<ntune_event_4; i++){ 
    residual_4 = 0.0;
    ref_mm_4 = 0.0; 
    ref_mm_4  = Lambda_real_4[foil_flag_4[i]];    
    ref_mm_4 = ref_mm_4/1000.0;
    
    XFP_4   = x_4[i];
    XpFP_4  = xp_4[i];
    YFP_4   = y_4[i];
    YpFP_4  = yp_4[i];
    z_av_4 = z_recon_4[i]; 
    ph1_4 = phir_4[i];  // open when calibrate the Momentum 
    ph1_4 = - ph1_4;
    ph2_4 = phil_4[i];    
    ph2_4 = - ph2_4;

    XFP_4   =(XFP_4 -XFPm)/XFPr;  
    XpFP_4  =(XpFP_4-XpFPm)/XpFPr;
    YFP_4   =(YFP_4 -YFPm)/YFPr;
    YpFP_4  =(YpFP_4-YpFPm)/YpFPr;
    z_av_sc_4 = (z_av_4 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_4[0] =  calcf2t_mom(param, XFP_4, XpFP_4, YFP_4, YpFP_4,  z_av_sc_4);
    moml_4[0] = moml_4[0]*Momr+Momm;

    moml_4[0] = moml_4[0]*2.21819/2.1; //for H/T and tritium only
    double hang9;   
    if(z_av_4<8.0e-2){
      hang9 =  ph2_4 + hrs_ang;
      hang9 = -hang9;
      delta_pep_4[0] = -1.35758*sin(-4.59571* hang9) + 2.09093;
    } 
    else{
      hang9 =  ph2_4 + hrs_ang;
      hang9 = -hang9;
      delta_pep_4[0] = 6.23409e-3*hang9 + 4.03363e-1;
    }    
    pep_real_4[0] = moml_4[0] + delta_pep_4[0]/1000.0; //LHRS  momentum at the reaction point in GeV 

    
    par_ep_4[0] = pep_real_4[0];// Wwhen LHRS momentum  tuned
    //  par_ep_4[0] = p10_4[i];
    par_ep_4[1] = p11_4[i];
    par_ep_4[2] = p12_4[i];
   
   
    ////  for the RHRS momentum tunning
    momr_4[0] =  calcf2t_mom(param, XFP_4, XpFP_4, YFP_4, YpFP_4,  z_av_sc_4);
    momr_4[0] = momr_4[0]*Momr+Momm;
    double hang10;    
    if(z_av_4<8.0e-2){
      hang10 = ph1_4 - hrs_ang;
      delta_pk_4[0] =-1.31749*sin(-4.61513* hang10) + 2.03687;
    } 
    else{
      hang10 = ph1_4 - hrs_ang;
      delta_pk_4[0] = 3.158e-2*hang10 + 4.05819e-1; 
    }
    pk_real_4[0] = momr_4[0] + delta_pk_4[0]/1000.0;    
    
   
    //// THR_4 =  calcf2t_th(param,XFP_4, XpFP_4, YFP_4, YpFP_4, z_av_sc_4);
    //// THR_4 = THR_4*Xptr + Xptm;   // jan 31

    //// PHR_4 =  calcf2t_th(param,XFP_4, XpFP_4, YFP_4, YpFP_4, z_av_sc_4);
    //// PHR_4 = PHR_4*Yptr + Yptm;

    //   par_k_4[0] = pk_real_4[0];// when RHRS matrix tuned
    par_k_4[0] = p13_4[i];    
    par_k_4[1] = p14_4[i];  // jan 31  
    par_k_4[2] = p15_4[i];
   
    //// par_k_4[1] = - THR_4; 
    //// par_k_4[2] = - PHR_4;

    halla_p_4 = p16_4[i];
    MM_4 = CalcMM(halla_p_4, par_ep_4, par_k_4, m_Al);  //need to adjust  
    MM_4 = MM_4 -25.3123;

    residual_4 = MM_4-ref_mm_4;
    //  cout<<" TT data = " << MM_4 << " ref_mm_4 = " << ref_mm_4 << " residual_4 = " << residual_4<<endl;  

    chi_44 = chi_44 +3*pow(residual_4,2.0);
  }
 //)))))))))))))))))))))))))))))))))))))))))))))   _4 Jan 07, 2020  ))))))))))))))))))))))))))))
 //)))))))))))))))))))))))))))))))))))))))))))))  Jan 12, 2020  ))))))))))))))))))))))))))))
  chi2 = chi_12 +chi_22 + chi_33 + chi_44;
  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  fval = chi2;
}

//cout<<"the value of mm_Al1 is = " << MM_3 << "and thet if al_3 is = "<<al_3<< endl;
