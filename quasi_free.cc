// Feb 23, 2020, studyiing the quasiree shape.
// The data in the QF.dat  is given by Dr. Tang from Atti function Hall C Thesis
// Author Bishnu pandey

void quasi_free()
{
  TFile * f_qu = new TFile("./root/quasi_free.root");
  TH1F *h_22 = (TH1F*)f_qu->Get("h1_qu");
  TH1F *h_22b = (TH1F*)f_qu->Get("h1_qub");
  
  const int no = 152;
  char name_quasi[500];
  sprintf(name_quasi,"./QF.dat");  
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
  // c1->SetFillColor(9);
  // c1->SetGrid();
  TF1 *f1 = new TF1("f1", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  auto gr = new TGraphErrors(no,x,y,0,0);
  gr->SetTitle("#Lambda Q.F. Shape");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP*"); 
  gr->Fit("f1","R+");
 
  TF1 *f2 = new TF1("f2", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f2->SetParameters(35.6934,23.3324,-0.909125,0.253986,-0.0164945,0.000500056,-7.12444e-06,3.78673e-08);
  TCanvas *c2 = new TCanvas("c2","c2",700,500);
  c2->cd();
  f2->Draw();
 
  TCanvas *c3 = new TCanvas("c3","c3",700,700);
  c3->cd();
  TF1 *f3 = new TF1("f3","0.016*f2 + 6.651",-2.0, 60.0);
  f3->SetParameters(35.6934,23.3324,-0.909125,0.253986,-0.0164945,0.000500056,-7.12444e-06,3.78673e-08);
  f3->Draw();


  // 2nd half of the quasi free
  TCanvas *c4 = new TCanvas("c4","c4",700,700);
  c4->cd();
  TF1 *f4 = new TF1("f4", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  auto gr2 = new TGraphErrors(no,x,y,0,0);
  gr2->SetTitle("#Lambda Q.F. Shape");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);
  gr2->Draw("AP*"); 
  gr2->Fit("f4","R+");
 
  TF1 *f5 = new TF1("f5", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f5->SetParameters(19563.2,-759.342,12.4578,-0.0886574,-4.42872e-05,4.89369e-06,-2.92745e-08,5.70998e-11);
  TCanvas *c5 = new TCanvas("c5","c5",700,500);
  c5->cd();
  f5->Draw();
 
  TF1 *f6 = new TF1("f6","0.016*f5 + 6.651",60.0,150.0);
  f6->SetParameters(19563.2,-759.342,12.4578,-0.0886574,-4.42872e-05,4.89369e-06,-2.92745e-08,5.70998e-11);
  TCanvas *c6 = new TCanvas("c6","c6",700,500);
  c6->cd();
  f6->Draw();

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
  



}
