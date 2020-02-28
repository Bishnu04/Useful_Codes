// Feb 23, 2020, studyiing the quasiree shape.
// The data in the QF.dat  is given by Dr. Tang from Atti function Hall C Thesis
// Author Bishnu pandey
// Tgaraphical cut for the one dimensional hsitogram
void tau_quasi_free()
{
  TFile * f_qu = new TFile("./root/quasi_free_feb25.root");
  TH1F *h_22 = (TH1F*)f_qu->Get("h1_qu");
  TH1F *T_T = (TH1F*)h_22->Clone("T_T");
 
  TH1F *h_22b = (TH1F*)f_qu->Get("h1_qub");
  TH1F *TT_final = new TH1F("TT_final","nnL Spectrum,T/T data(Bg subtracted)  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);


  TH1F *HT = (TH1F*)f_qu->Get("H_T");
  TH1F *HH = (TH1F*)HT->Clone();

  TH1F *HT_B = (TH1F*)f_qu->Get("H_TB"); 
  TH1F *H_add = new TH1F("H_add","H inT/T data(BG subtracted);-B_{#Lambda}(MeV);Counts/ 1.52MeV ",150,-100,128);


  const int no = 152;
  char name_quasi[500];
  sprintf(name_quasi,"./QF_FEB23.dat");  
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


  // TCanvas *c1 = new TCanvas("c1","c1",700,700);
  // // c1->SetFillColor(9);
  // // c1->SetGrid();
  // TF1 *f1 = new TF1("f1", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  // TGraph *gr = new TGraphErrors(no,x,y,0,0);
  // gr->SetTitle("F1 = 0.016*f1 + 6.651 ");
  // // gr->SetMarkerColor(4);
  // gr->SetMarkerStyle(21);
  // gr->Draw("AP*"); 
  // gr->Fit("f1","R+");
 
  TF1 *f2 = new TF1("f2", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f2->SetParameters(34.6742,22.777,-0.590337,0.211099,-0.0140323,0.000431205,-6.20371e-06,3.3182e-08);
  // TCanvas *c2 = new TCanvas("c2","c2",700,500);
  // c2->cd();
  // f2->Draw();
 
  // TCanvas *c3 = new TCanvas("c3","c3",700,700);
  // c3->cd();
  TF1 *f3 = new TF1("f3","0.016*f2 + 6.651",-2.0, 60.0);
  f3->SetParameters(34.6742,22.777,-0.590337,0.211099,-0.0140323,0.000431205,-6.20371e-06,3.3182e-08);
  // f3->Draw();


  // // 2nd half of the quasi free
  // TCanvas *c4 = new TCanvas("c4","c4",700,700);
  // c4->cd();
  // TF1 *f4 = new TF1("f4", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  // TGraph *gr2 = new TGraphErrors(no,x,y,0,0);
  // gr2->SetTitle("F2 = 0.016*f2 + 6.651");
  // //  gr2->SetMarkerColor(3);
  // gr2->SetMarkerStyle(21);
  // gr2->Draw("AP*"); 
  // gr2->Fit("f4","R+");
 
  TF1 *f5 = new TF1("f5", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f5->SetParameters(8486.83,-185.703,1.27841,0.00263458,-5.96243e-05,-1.35606e-07,3.37053e-09,-9.68808e-12);
  // TCanvas *c5 = new TCanvas("c5","c5",700,500);
  // c5->cd();
  // f5->Draw();
 
  TF1 *f6 = new TF1("f6","0.016*f5 + 6.651",60.0,150.0);
  f6->SetParameters(8486.83,-185.703,1.27841,0.00263458,-5.96243e-05,-1.35606e-07,3.37053e-09,-9.68808e-12);
  // TCanvas *c6 = new TCanvas("c6","c6",700,500);
  // c6->cd();
  // f6->Draw();

  h_22b->Scale(1.0/38.0);
  TCanvas *c7 = new TCanvas("c7","c7", 600,600);
  c7->cd();
  h_22->Draw();
  h_22b->Draw("E2 same");
  h_22b->SetFillStyle(3002);
  h_22b->SetMarkerStyle(28);
  h_22b->SetMarkerColor(kGreen);
  TH1F *bg_cpy = (TH1F*)h_22b->Clone();
  f3->Draw("same"); // 1st half og QF
  f6->Draw("same"); //2nd  half og QF

  HT_B->Scale(1.0/38.0);
  TCanvas *c8 = new TCanvas("c8","c8", 600,600);
  c8->cd();
  HT->Draw();
  HT_B->Draw("E2 same");
  HT_B->SetFillStyle(3002);
  HT_B->SetMarkerStyle(28);
  HT_B->SetMarkerColor(kGreen);
  TH1F *HT_B1 = (TH1F*)HT_B->Clone();


  TCanvas *c9 = new TCanvas("c9","c9",600,600);
  // H_add->Add(HH,HT_B1,1.0,-1.0);

  TH1F *H_add1 = (TH1F*)HH->Clone();
  HH->Draw();

  //// using the TGraphical cut to calculate the actual entries on the H peak
  //// for the tritium analyzed as H target
  TH1F * HT_final = new TH1F ("HT_final","H in  T/T data(BG not subtracted)  ;-B_{#Lambda}(MeV);Counts/ 1.52MeV ",150,-100,128);

  TCutG *cut50 = new TCutG("cut50",6); // 6 is the no of points
  cut50->SetPoint(0,-8.93103,92.9007 ); // thses points gives a closed region
  cut50->SetPoint(1,-8.93103,52.0843 ); 
  cut50->SetPoint(2,-8.93103,52.0843 );
  cut50->SetPoint(3,8.81322,52.0843 );
  cut50->SetPoint(4,8.81322,92.9007 );    
  cut50->SetPoint(5,-8.93103,92.9007 );

  Double_t xVal,yVal;
  for(Int_t i=1;i<=H_add1->GetNbinsX();i++){   
    xVal = H_add1->GetXaxis()->GetBinCenter(i);// i is the bin number goes from 1 to 150(bin No used) xVal goes from -100 to 128
      yVal = H_add1->GetBinContent(i);// Bin height
      
      if(cut50->IsInside(xVal,yVal)){
  	cout<<" the value of xval = "<< xVal << " and yVal is = "<<yVal<<endl;
  	HT_final->SetBinContent(i,H_add1->GetBinContent(i) -52.0843);//filling histogram.here 49.978 is lowest y value in cut50 
      }     /// need to change the abobe - number 
  }
  
  TCanvas *c10 = new TCanvas("c10","c10",600,600);
  HT_final->Draw();
  cout<<"The total no of H event = "<< HT_final->Integral(1,HT_final->GetNbinsX())<<endl;
  cut50->Print();
  TLatex lte;
  lte.SetTextSize(0.025);
  lte.DrawLatex(10.0,30.0,Form("#color[4]{The total no of H event = %.6g}",HT_final->Integral(1,HT_final->GetNbinsX())));


  //// for the nnl spectrum
  // TCanvas *c_T = new TCanvas("c_T", "c_T", 600,600);
  // c_T->cd();
  // // TT_final ->Add(T_T,bg_cpy,1.0,-1.0);
  // T_T-> Draw();
  // TH1F * H1 = (TH1F*)T_T->Clone();
 

  // //// for the nnl spectrum
  // TH1F *H2 = new TH1F("H2","nnL Spectrum (Bg not subtracted)  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
 
  // TCutG *cut30 = new TCutG("cut30",6); // 6 is the no of points
  // cut30->SetPoint(0,50.9146,52.8375 ); // thses points gives a closed region
  // cut30->SetPoint(1,40.1679,52.8375 ); 
  // cut30->SetPoint(2,40.1679,40.7522 );// need to change
  // cut30->SetPoint(3,60.4331,40.7522 );
  // cut30->SetPoint(4,60.4331,52.8375 );    
  // cut30->SetPoint(5,50.9146,52.8375);

  // Double_t xVal,yVal;
  // for(Int_t i=1;i<=H1->GetNbinsX();i++){   
  //   xVal = H1->GetXaxis()->GetBinCenter(i);// i is the bin number goes from 1 to 150(bin No used) xVal goes from -100 to 128
  //     yVal = H1->GetBinContent(i);// Bin height
      
  //     if(cut30->IsInside(xVal,yVal)){
  // 	cout<<" the value of xval = "<< xVal << " and yVal is = "<<yVal<<endl;
  // 	H2->SetBinContent(i,H1->GetBinContent(i) -40.7522);//filling histogram.here 49.978 is lowest y value in cut30 
  //     }     /// need to change the abobe - number 
  // }
  
  gStyle->SetOptStat(0);
  // TCanvas *c20 = new TCanvas("c20","c20",600,600);
  // c20->cd();
  // H2->Draw();
  // cout<<"The total no of H event = "<< H2->Integral(1,H2->GetNbinsX())<<endl;
  // cut30->Print();
  // TLatex lte;
  // lte.SetTextSize(0.025);
  // lte.DrawLatex(-80.0,10.0,Form("#color[4]{The total no of H event = %.6g}",H2->Integral(1,H2->GetNbinsX())));
}
