// Feb 23, 2020, studyiing the quasiree shape.
// The data in the QF.dat  is given by Dr. Tang from Atti function Hall C Thesis
// Author Bishnu pandey
// Note on Februray 28, 2020, Sho Nagao updated this code for the Gaussian fitting.
// He added two gaussian and that of the 7th order polynomial  together.
void quasi_free()
{
  TFile * f_qu = new TFile("./root/quasi_free_march2.root"); // root file has just one nnl Hisatogram
  TH1F *h_22 = (TH1F*)f_qu->Get("h1_qu");
  TH1F *h_22b = (TH1F*)f_qu->Get("h1_qub");
  
  TH1F *HT = (TH1F*)f_qu->Get("H_T");
  TH1F *HT_B = (TH1F*)f_qu->Get("H_TB");
  TH1F *HH = (TH1F*)HT->Clone();
  TH1F *H_add = new TH1F("H_add","H inT/T data(BG subtracted);-B_{#Lambda}(MeV);Counts/ 1.52MeV ",150,-100,128);
 
  
  
  // // TH1F *H1_0 = (TH1F*)f_qu->Get("h_75_1");
  // // TH1F *H1_0b = (TH1F*)f_qu->Get("h_75_1b");
  TH1F *H_75_1 = (TH1F*)f_qu->Get("h_75_1");
  TH1F *H_75_1b = (TH1F*)f_qu->Get("h_75_1b");
 

  /////For LAmbda Quasi frre
  const int no = 152;
  char name_quasi[500];
  //  sprintf(name_quasi,"./QF_FEB_27.dat"); 
  sprintf(name_quasi,"./QF_MARCH_3.dat"); // to input the quasifree shape data 
  ifstream quasi(name_quasi);
  double x[no];
  double y[no];
  for (int i=0;i<no;i++){
    double x1=0.0;
    double y1 = 0.0;
    quasi >> x1 >>y1;
    x[i] = x1;
    y[i] = y1; 
    // cout<<"the value of x1 = " <<x1 << " and that of y1 = "<< y1<<endl;   
  }
  quasi.close();

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  TF1 *f1 = new TF1("f1", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f1->SetNpx(1000); // divides the fit interval in 1000 parts in stead the default value 100

  TGraph *gr = new TGraphErrors(no,x,y,0,0);
  gr->SetTitle("up to 7th order polynomial ");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP*"); 
  gr->Fit("f1","R+");

  // // TCanvas *cc = new TCanvas("cc","cc",800,600);
  // // cc->cd(1);
  // // hh->Draw();
  // // gr->Draw("sameP*"); // sho nagao used to fit the function from -2 to 65
  // // gr->Fit("f1","","",-2.0,65.0); 
  
  TF1 *f2 = new TF1("f2", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f2->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
  f2->SetNpx(1000);
  // TCanvas *c2 = new TCanvas("c2","c2",700,500);
  // c2->cd();
  // f2->Draw();
 
  ////// 0.016 is the sclae down factor and 6.651 is the background height
  TF1 *f3 = new TF1("f3","0.016*f2 + 6.651",-2.0, 60.0); // How to scale down the 7th order polynomial
  f3->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
  f3->SetNpx(1000);
  // TCanvas *c3 = new TCanvas("c3","c3",700,700);
  // c3->cd();
  // f3->Draw();
 
  
  // // // 2nd half of the quasi free
  TCanvas *c4 = new TCanvas("c4","c4",700,700);
  c4->cd();
  TF1 *f4 = new TF1("f4","[0]+[1]*x +[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f4->SetNpx(1000);
  TGraph *gr2 = new TGraphErrors(no,x,y,0,0);
  gr2->SetTitle("up to 7th order polynomial");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);
  gr2->Draw("AP*"); 
  gr2->Fit("f4","R+");
  
  TF1 *f5 = new TF1("f5","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f5->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f5->SetNpx(1000);
  // TCanvas *c5 = new TCanvas("c5","c5",700,500);
  // c5->cd();
  // f5->Draw();
  
  TF1 *f6 = new TF1("f6","0.01584*f5 +6.58449",60.0,75.0);
  f6->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f6->SetNpx(1000);
  // TCanvas *c6 = new TCanvas("c6","c6",700,500);
  // c6->cd();
  // f6->Draw();


 ///////+++++++++++++  Sigma quasifree shape   ++++++++++++++++++++++++++++++++++++++
 const int n1 = 16;
 char name_sigma_quasi[500];
 sprintf(name_sigma_quasi,"./SIGMA_March_3.dat");
 ifstream sigma_quasi(name_sigma_quasi);
 double A1[n1];
 double B1[n1];
 for(int i=0;i<n1;i++){
   double a1=0.0;
   double b1=0.0;
   sigma_quasi >>a1>>b1;
   A1[i] = a1;
   B1[i] = b1;
   //  cout<<"the value of a1 = " <<a1 << " and that of b1 = "<< b1<<endl;
 }
 sigma_quasi.close();


 // TCanvas *ca = new TCanvas("ca","ca", 600, 600);
 // ca->cd();
 TF1 *fa = new TF1("fa", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",75.0,150.0);
 fa->SetNpx(1000); // 
 TGraphErrors *ga = new TGraphErrors(n1,A1,B1,0,0);
 ga->SetTitle("up to 7th order polynomial ");
 ga->SetMarkerColor(4);
 ga->SetMarkerStyle(21);
 // ga->Draw("AP*"); 
 // ga->Fit("fa","R+");


 TF1 *fb = new TF1("fb", "0.6*fa",75.0,150.0);
 fb->SetParameters(53.7528,-1.44004,0.0042699,8.9401e-05,1.77886e-07,-4.10384e-09,-3.44382e-11,2.21558e-13);
 fb->SetNpx(1000);
 // TCanvas *cb = new TCanvas("cb","cb", 600, 600);
 // cb->cd();
 // fb->Draw();

///////+++++++++++++ up to here  Sigma quasifree shape   ++++++++++++++++++++++++++++++++++++++

  double pp[8];
  f6->GetParameters(pp);
  for(int i=0;i<8;i++){
    pp[i] = pp[i]*0.01584;
  }
  pp[0] =pp[0]+ 6.58449;

  ////======================+++++++++++
  double pp1[8];
  fb->GetParameters(pp1);
  for(int i=0;i<8;i++){
    pp1[i] = 0.3*pp1[i];
  }
 
  
  TF1 *f50 = new TF1("f50"," 0.01584*f5 +6.58449",75.0,150.0);
  f50->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f50->SetNpx(1000);
  
 TF1 *fc = new TF1("fc", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7) +[8]+[9]*x + [10]*pow(x,2) + [11]*pow(x,3) + [12]*pow(x,4) + [13]*pow(x,5)+[14]*pow(x,6)+[15]*pow(x,7)  ",75.0,150.0);
  fc->SetNpx(1000);
  
  fc->FixParameter(0,pp[0]);
  fc->FixParameter(1,pp[1]);
  fc->FixParameter(2,pp[2]);
  fc->FixParameter(3,pp[3]);
  fc->FixParameter(4,pp[4]);
  fc->FixParameter(5,pp[5]);
  fc->FixParameter(6,pp[6]);
  fc->FixParameter(7,pp[7]);
  
  fc->FixParameter(8,pp1[0]);
  fc->FixParameter(9,pp1[1]);
  fc->FixParameter(10,pp1[2]);
  fc->FixParameter(11,pp1[3]);
  fc->FixParameter(12,pp1[4]);
  fc->FixParameter(13,pp1[5]);
  fc->FixParameter(14,pp1[6]);
  fc->FixParameter(15,pp1[7]);
  
  
  h_22b->Scale(1.0/38.0);
  TCanvas *cc = new TCanvas("cc","cc", 600,600);
  cc->cd();
  h_22->Draw();
  h_22b->Draw("E2 same");
  h_22b->SetFillStyle(3002);
  h_22b->SetMarkerStyle(28);
  h_22b->SetMarkerColor(kGreen);
  //  h_22->Fit("f50","","",75,150.0);
  f3->Draw("same");
  f6->Draw("same");
  f6->SetLineColor(kGreen);
  // f50->Draw("same");
  // f50->SetLineColor(kBlue);
  // fb->Draw("same");
  fc->Draw("same");
  fc->SetLineColor(kBlue);
  /////++++++++++++++++++++++u to here is the addition of two 7th order pol +++++++++++++++++++++++++++++++++++++++++++++++++++++
  /////++++++++++++++++++++++u to here is the addition of two 7th order pol +++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  /*
  TF1 * f7 = new TF1("f7", "gaus", -2.0,1.83);
  TF1 * f8 = new TF1("f8", "gaus", 5.3,11.15);
  //  TF1 *f9 = new TF1("f9","pol7(0)+gaus(8)+gaus(11)", -2.17,60);
  TF1 *f9 = new TF1("f9","gaus(0)+gaus(3)+[6]+[7]*x+[8]*pow(x,2)+[9]*pow(x,3)+[10]*pow(x,4)+[11]*pow(x,5)+[12]*pow(x,6)+[13]*pow(x,7)", -2.0,60);
  f9->SetNpx(1000);
for(int i=0;i<8;i++){
  //  cout<<"the value of pp is = "<<pp[i]<<endl;
}
  f9->FixParameter(6,pp[0]);
  f9->FixParameter(7,pp[1]);
  f9->FixParameter(8,pp[2]);
  f9->FixParameter(9,pp[3]);
  f9->FixParameter(10,pp[4]);
  f9->FixParameter(11,pp[5]);
  f9->FixParameter(12,pp[6]);
  f9->FixParameter(13,pp[7]);
  f9->SetParameter(0,20);
  f9->SetParameter(1,0);
  f9->SetParameter(2,1.0);
  f9->SetParameter(3,25);
  f9->SetParameter(4,8.0);
  f9->SetParameter(5,1.0);
  f9->SetParLimits(0,1,20);
  f9->SetParLimits(1,-1,1);
  f9->SetParLimits(2,0.5,2.0);
  f9->SetParLimits(3,1,25);
  f9->SetParLimits(4,7,9);
  f9->SetParLimits(5,0.5,1.4); // use a smal;l number it displays the same number
  f9->SetLineWidth(1);
  
  h_22b->Scale(1.0/38.0);
  TCanvas *c7 = new TCanvas("c7","c7", 600,600);
  c7->cd();
  h_22->Draw();
  h_22b->Draw("E2 same");
  h_22b->SetFillStyle(3002);
  h_22b->SetMarkerStyle(28);
  h_22b->SetMarkerColor(kGreen);
  f3->Draw("same"); // 1st half og QF
  f6->Draw("same"); //2nd  half og QF
  fb->Draw("same");
  h_22->Fit("f9","","",-2.0,11.15);
  ///// h_22->Fit("f9","R+");
  f9->Draw("same");

  TLatex l;
  l.SetTextSize(0.025);
  l.DrawLatex(-80,42,Form("#color[2]{I.  mean = %.6g +/- %.6g}",f9->GetParameter(1),f9->GetParError(1)));
  l.DrawLatex(-80,40,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f9->GetParameter(2),f9->GetParError(2)));
  l.DrawLatex(-80,48,Form("#color[2]{II. mean = %.6g +/- %.6g}",f9->GetParameter(4),f9->GetParError(4)));
  l.DrawLatex(-80,46,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f9->GetParameter(5),f9->GetParError(5)));
  
  */

  // h_22->Fit("f8","R+");
  // f7->Draw("same");
  // h_22->Fit("fit_fun","","",-2.0,10.5);
  // h_22->Fit("f_it","","",-2.0,10.5);
  //  f_it->Draw("same");

  // // HT_B->Scale(1.0/38.0);
  // // TCanvas *c8 = new TCanvas("c8","c8", 600,600);
  // // c8->cd();
  // // HT->Draw();
  // // HT_B->Draw("E2 same");
  // // HT_B->SetFillStyle(3002);
  // // HT_B->SetMarkerStyle(28);
  // // HT_B->SetMarkerColor(kGreen);
  // // TH1F *HT_B1 = (TH1F*)HT_B->Clone();

  /////// The following lines asre for the different bin size of the nnL spectrum
 ////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++======
  TF1 *f14 = new TF1("f14","0.008*f5 + 3.3225",60.0,150.0);
  f14->SetParameters(8486.83,-185.703,1.27841,0.00263458,-5.96243e-05,-1.35606e-07,3.37053e-09,-9.68808e-12);
  f14->SetNpx(1000);

  TF1 *f10 = new TF1("f10","0.008*f2 + 3.3255",-2.0, 60.0); // How to scale down the 7th order polynomial
  f10->SetParameters(30.546,22.4772,-0.210242,0.167295,-0.0119721,0.000383648,-5.6675e-06,3.08163e-08);
  f10->SetNpx(1000);
  
  double pp2[8];
  f10->GetParameters(pp2);
  for(int i=0;i<8;i++){
    pp2[i] *= 0.008;
  }
  pp2[0] += 3.3255;
  
  // TF1 * f11 = new TF1("f11", "gaus", -2.0,1.83);
  TF1 * f12 = new TF1("f12", "gaus", 5.3,11.15);
  //  TF1 *f13 = new TF1("f13","gaus(0)+gaus(3)+[6]+[7]*x+[8]*pow(x,2)+[9]*pow(x,3)+[10]*pow(x,4)+[11]*pow(x,5)+[12]*pow(x,6)+[13]*pow(x,7)", -2.0,60);
  TF1 *f13 = new TF1("f13","gaus(0)+[3]+[4]*x+[5]*pow(x,2)+[6]*pow(x,3)+[7]*pow(x,4)+[8]*pow(x,5)+[9]*pow(x,6)+[10]*pow(x,7)", -2.0,60);
  
  for(int i=0;i<8;i++){
    //  cout<<"the value of pp is = "<<pp[i]<<endl;
  }
  f13->FixParameter(3,pp2[0]);
  f13->FixParameter(4,pp2[1]);
  f13->FixParameter(5,pp2[2]);
  f13->FixParameter(6,pp2[3]);
  f13->FixParameter(7,pp2[4]);
  f13->FixParameter(8,pp2[5]);
  f13->FixParameter(9,pp2[6]);
  f13->FixParameter(10,pp2[7]);
  // // f13->SetParameter(0,20);
  // // f13->SetParameter(1,0);
  // // f13->SetParameter(2,1.0);
  
  f13->SetParameter(0,25);
  f13->SetParameter(1,8.0);
  f13->SetParameter(2,1.0);
  // f13->SetParLimits(0,1,20);
  // f13->SetParLimits(1,-1,1);
  // f13->SetParLimits(2,0.5,2.0);
  f13->SetParLimits(0,1,25);
  f13->SetParLimits(1,7,9);
  f13->SetParLimits(2,0.5,1.75); // use a smal;l number it displays the same number
  f13->SetLineWidth(1);
  gStyle->SetOptFit(1);  
  
  H_75_1b->Scale(1.0/38.0);
  TCanvas *c10 = new TCanvas("c10","c10", 600,600);
  c10->cd();
  H_75_1->Draw();
  H_75_1b->Draw("E2 same");
  H_75_1b->SetFillStyle(3002);
  H_75_1b->SetMarkerStyle(28);
  H_75_1b->SetMarkerColor(kGreen);
  //  f10->Draw("same");
  H_75_1->Fit("f13","","",5.3,11.15);
  f13->Draw("same");
  f14->Draw("same");
  
  TLatex l1;
  l1.SetTextSize(0.025);
  // l1.DrawLatex(-80,22,Form("#color[2]{I.  mean = %.6g +/- %.6g}",f13->GetParameter(1),f13->GetParError(1)));
  // l1.DrawLatex(-80,20,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f13->GetParameter(2),f13->GetParError(2)));
  l1.DrawLatex(-80,28,Form("#color[2]{II. mean = %.6g +/- %.6g}",f13->GetParameter(1),f13->GetParError(1)));
  l1.DrawLatex(-80,26,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f13->GetParameter(2),f13->GetParError(2))); 
  ////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
}
