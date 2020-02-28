/*
  test.cc
  
  Toshi, Mar 17, 2019
*/


// ---- Definition of the function you use for fitting ------ //
double twog(double*x, double*par){
  double g1,g2;
   g1 = par[0] * exp(-0.5*pow((x[0]-par[1])/par[2],2.0)); // 3 parameters
   g2 = par[3] * exp(-0.5*pow((x[0]-par[4])/par[5],2.0)); // 3 parameters
  
  return g1 + g2;
}


void GausFitSample(){
  // ----- Sample histograms ------------------- //
  TH1F* h1 = new TH1F("h1","h1",200,-10,20); // making the histograms
  TF1* f1 = new TF1("f1","gaus"); // just a function 
  f1->SetParameters(1.0,0.0,1.0);  // amplitude, mean  sigma or the widthxxxxxx
  TF1* f2 = new TF1("f2","gaus");
  f2->SetParameters(1.0,6.0,1.0);
  TF1* f3 = new TF1("f3","gaus");
  f3->SetParameters(1.0,15.0,2.0);
  
  // ---- Filling histograms with random data -----
  h1->FillRandom("f1",5000); // making the histograms by using the random method
  h1->FillRandom("f2",5000);
  h1->FillRandom("f3",5000);

  // ----- Defining the fitting function ------- //
  TF1* fitfunc = new TF1("fitfunc",twog,-10,20,6); // here 6 is no of parameter, and -10 to 20 is th fitting range
  fitfunc->SetParameters(1,   // 0 (amp)
			 0.0, // 1 (mean)
			 1.0, // 2 (width)
			 1,   // 3 (amp)  
			 7.0, // 4 (mean) 
			 1.0); //5 (width)
  
  fitfunc->SetParLimits(1,
			-3.0,3.0); // these controls the range of the gaussian
  fitfunc->SetParLimits(4,
			4.0,18.0);
  
  //  fitfunc->FixParameter(2,1.0); // CONTROLS THE WIDTH OF THE GAUSSIAN OR THE CONTROLS THE SIGMA VALUE
  //fitfunc->FixParameter(5,1.0); // width
  
  fitfunc->SetNpx(1000); // to make the function smooth BY DIVIDING THE WHOLE RANGE IN TO 1000 EQUAL PARTS IN STEAD OF DEFAULT VALUE 100
  
  
  // ---- fitting -----
  h1->Draw();
  h1->Fit("fitfunc","Nq","",-5.0,10.0); // Fitting but not be shown
  fitfunc->Draw("same");
}
