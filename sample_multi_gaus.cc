double three_gaus(double *x,double *par){
  double g1,g2,g3;
 g1 = par[0] * exp(-0.5*pow((x[0]-par[1])/par[2],2.0)); // 3 parameters
 g2 = par[3] * exp(-0.5*pow((x[0]-par[4])/par[5],2.0));
 g3 = par[6] * exp(-0.5*pow((x[0]-par[7])/par[8],2.0));
 return g1+g2+g3+25;
}


void sample_multi_gaus(){

  TH1F *h1 = new TH1F("h1","h1",200,-20,20);// histogram defn

  TF1 *f1 = new TF1("f1","gaus"); // First peak in the form og gaussian
  f1->SetParameters(1,-3,0.5);/// amplitude, mean,sigma

  TF1 *f2 = new TF1("f2","gaus"); // 2nd peak in the form og gaussian
  f2->SetParameters(1,0,0.5);

  TF1 *f3 = new TF1("f3","gaus");
  f3->SetParameters(1,7,0.5);

  
  h1->FillRandom("f1",5000); // filling the events in the gaussinan form
  h1->FillRandom("f2",5000);
  h1->FillRandom("f3",5000);

  ///now lets define the fitting function
 TF1* fitfunc = new TF1("fitfunc",three_gaus,-10,15,9); // here 9 is no of parameter, and -10 to 20 is th fitting range
 fitfunc->SetParameter(0,1);// amplitude
 fitfunc->SetParameter(1,-3); // mean 			
 fitfunc->SetParameter(2,0.5); // sigma

 fitfunc->SetParameter(3,1);// amplitude
 fitfunc->SetParameter(4,0); // mean 			
 fitfunc->SetParameter(5,0.5); // sigma

 fitfunc->SetParameter(6,1);// amplitude
 fitfunc->SetParameter(7,7); // mean 			
 fitfunc->SetParameter(8,0.5); // sigma

 // fitfunc->SetParLimits(0,
 //		       -5.0,5.0); // these controls the range of the gaussian

 // fitfunc->SetParLimits(4,
 //			4.0,18.0);
  fitfunc->SetNpx(1000);   
  
  TCanvas *c = new TCanvas("c","c",600,600);
  h1->Draw(); //hist draw
  h1->Fit("fitfunc","Nq","",-10.0,12.0); // Fitting but not be shown
  fitfunc->Draw("same");

 


}
 
