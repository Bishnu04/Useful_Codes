// 03/ 15/ 2020
// In this code(error_study ) I am studying the error by defining the Gaussian test as explained by Dr. Tang

//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TRandom.h>
#include <TFile.h>
#include <stdio.h>
#include <TROOT.h>
#include <TTree.h>

using namespace std;
double RndUni(double x1,double x2)
{
  double Uni=gRandom->Uniform(x1,x2);
  return Uni;
}

void error_study()
{

  double delta_E;// will be the error for the energy
  double delta_E_min;
  double delta_E_max;
  double delta_Ee;

  double test;// is the test between 0 and that of 1
  double test_min;
  double test_max;

  double g1; // Introducing the gaussian type of the variable
  double sigma_e;// one third od the delta_E_min ordelta_E_max
  double pi;
  double A;// To store the variables
  double B;// To store the variables

  double delta_pep_min;// to see the errors in the LHRS  momentum
  double delta_pep_max;
  double delta_pep;
  double sigma_ep;// one third of delta_pep

  double delta_pk_min;// to see the errors in the RHRS momentum
  double delta_pk_max;
  double delta_pk;
  double sigma_pk;// one third of delta_pk

  double delta_thetaep_min;// to see the errors in the LHRS Scattering Angle
  double delta_thetaep_max;
  double delta_thetaep;
  double sigma_thetaep;// one third of delta_thetaep

  double delta_thetak_min;// to see the errors in the RHRS Scattering Angle
  double delta_thetak_max;
  double delta_thetak;
  double sigma_thetak;// one third of delta_thetaep

  double delta_ph12_min;// to see the errors in the phi_1 - phi_2
  double delta_ph12_max;
  double delta_ph12;
  double sigma_ph12;// one third of phi_1 - phi_2


  int loop_num = 10000;
  delta_E_min = -0.0012;// delta E
  delta_E_max = 0.0012; //GeV
  test_min = 0.;
  test_max= 1.0;

  sigma_e= 0.0004;
  pi = 3.14159;
 
  delta_pep_min = -0.000665;
  delta_pep_max = 0.000665;
  sigma_ep= 0.0002218;

  delta_pk_min = -0.000547;
  delta_pk_max = 0.000547;
  sigma_pk= 0.0001823;

  delta_thetaep_min = -0.0105; // need the exact value
  delta_thetaep_max =0.0105;// need the exact value
  sigma_thetaep=0.0035;// need the exact value
   
  delta_thetak_min = -0.0105; // need the exact value
  delta_thetak_max =0.0105;// need the exact value
  sigma_thetak=0.0035;// need the exact value

  delta_ph12_min = -0.0105; // need the exact value
  delta_ph12_max =0.0105;// need the exact value
  sigma_ph12=0.0035;// need the exact value
  
  

  TFile *f2 = new TFile("error_study.root","recreate");
  f2->cd();
  TTree *ktree = new TTree("ktree","generated data");
  ktree->Branch("delta_E", &delta_E,"delta_E/D");
  ktree->Branch("test", &test,"test/D");
  ktree->Branch("g1", &g1,"g1/D");
  ktree->Branch("delta_Ee", &delta_Ee,"delta_Ee/D");
  ktree->Branch("sigma_e", &sigma_e,"sigma_e/D");
  ktree->Branch("A", &A,"A/D");
  ktree->Branch("B", &B,"B/D");

  ktree->Branch("delta_pep", &delta_pep,"delta_pep/D");
  ktree->Branch("delta_pk", &delta_pk,"delta_pk/D");
  ktree->Branch("delta_thetaep", &delta_thetaep,"delta_thetaep/D");
  ktree->Branch("delta_thetak", &delta_thetak,"delta_thetak/D");
  ktree->Branch("delta_ph12", &delta_ph12,"delta_ph12/D");

  for(int i=0; i<loop_num;i++)
    {
      delta_E= RndUni(delta_E_min,delta_E_max);
      test=RndUni(test_min,test_max);
      //A = exp(-delta_E*delta_E/(2*sigma_e*sigma_e));
      A = delta_E*delta_E/(2*sigma_e*sigma_e);
      B=sigma_e*sqrt(2*pi);
      g1 =exp(-A)/B;

      // if(g1<test)
      // 	{delta_Ee=delta_E;}     

      delta_pep = RndUni(delta_pep_min,delta_pep_max);
      delta_pk = RndUni(delta_pk_min,delta_pk_max);
      delta_thetaep = RndUni(delta_thetaep_min,delta_thetaep_max);
      delta_thetak = RndUni(delta_thetak_min,delta_thetak_max);
      delta_ph12 = RndUni(delta_ph12_min,delta_ph12_max);

      ktree->Fill();
    }
  ktree->Write();
  f2->Close();

  cout<<" hello world!"<< endl;
}
