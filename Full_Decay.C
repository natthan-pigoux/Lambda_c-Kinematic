#include <cmath>
#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TImage.h"
#include "TStyle.h"
#include "TRandom.h"
#include "THistPainter.h"
#include "TLegend.h"
#include <tuple>
using namespace std;

/* do : root "name_of_the_macro".C to run the code (or .x "name_of _the_macro".C if you are already in root)*/
/*-> then it will create a pdf named "Full-kinematic.pdf" in the directory with all the plots */

                                                           /*###### Functions #####*/

/* ###### Compute the quadri_impulsions of all final particles from Lambda_b decay in its rest frame ##### */

/*  ########  Define Boost and Rotation ####### */

/* Compute boost from a rest frame to the resonnance frame */

tuple<double, double , double , double> Boost(double  p[4], double beta){

  double gamma = 1 / sqrt( 1 - pow( beta ,2 ) );


  double E,Px,Py,Pz;

  E = gamma*p[0] + gamma*beta*p[3];
  Px = p[1];
  Py = p[2];
  Pz = gamma*p[3] + gamma*beta*p[0];

  return make_tuple(E,Px,Py,Pz); /*make_tuple is simply to return a tuple. C'est la seule façon que j'ai trouver de renvoyer un genre de vecteur*/

}


/* Compute Rotation to the helicity frame of the resonnance */

tuple<double, double, double, double > Rotation(double phi, double theta, double p[4]){

  double E, Px, Py, Pz;

  E = p[0];
  Px = cos(phi)*cos(theta)*p[1] - sin(phi)*p[2] + cos(phi)*sin(theta)*p[3];
  Py = sin(phi)*cos(theta)*p[1] + cos(phi)*p[2] + sin(phi)*sin(theta)*p[3];
  Pz = -sin(theta)*p[1] + cos(theta)*p[3];

  return make_tuple(E,Px,Py,Pz);

}



  tuple <double , double , double > Impulsions(double P, double phi, double theta){

    double p[3];

    p[0] = P*sin(theta)*cos(phi);
    p[1] = P*sin(theta)*sin(phi);
    p[2] = P*cos(theta);

    return make_tuple(p[0],p[1],p[2]);

  }

  /* ####### MAIN ####### */


  void Full_Decay(){

  
  /* Set variables and parameters */
  
  TRandom3 *r1 = new TRandom3();
  TRandom3 *r2 = new TRandom3();
  TRandom3 *r3 = new TRandom3();
  TRandom3 *r4 = new TRandom3();
  TRandom3 *r5 = new TRandom3();
  TRandom3 *r6 = new TRandom3();
  
  r1->SetSeed(0);
  r2->SetSeed(0);
  r3->SetSeed(0);
  r4->SetSeed(0);
  r5->SetSeed(0);
  r6->SetSeed(0);
  
  double p1 = 0.85, p2 = 0.02, p3 = 0.13, prand, Pi = acos(-1.0);

  double theta_lc, phi_lc, cosT_lc, theta_lb, phi_lb,cosT_lb, theta_12, phi_12,cosT_12, theta_32, phi_32,cosT_32, cosT_f, phi_f;

  double a1_lc = 0.094, a2_lc = 0.221, a_lcs[3];

  a_lcs[0] = -1.0, a_lcs[1] = 0.0, a_lcs[2] = 1.0;
  
  double Mb = 5618.6, Mc12s = 2595.2, Mc32s = 2628.1, Mmu = 105.65, Mnu= 0.000001,  Mc = 2286.6, Ml = 1115.68, Mpip = 139.57, Mprot = 938.272, Mpim = 139.57;
  
  double q2, Plcs, Plc, plc[4], plcs[4], beta_lcs, norm =1;

  TF1 *f_cosT_32 = new TF1("f_cosT_32"," (1+ [0]*pow(x,2) )/[1]",-1.0,1.0);
  TF1 *f_cosT_lb = new TF1("f_cosT_lb", "(1+[0]*x)/2",-1.0,1.0);
 f_cosT_lb->SetParameter(0, a1_lc);

 TF1 *f_phi_lb = new TF1("f_phi_lb", "(1 + [0]*cos(x))/(2*[1])",0.0,2*Pi);
 f_phi_lb->SetParameter(0,a2_lc);
 f_phi_lb->SetParameter(1, Pi);


  int N =100000;

  
  /* Define Canvas and Histograms */

    
  for (int j=0; j<3; j++){


    TCanvas *Can = new TCanvas("Can","");
    Can->Divide(2,1);

    TH1D *hcosT = new TH1D("hcosT","cos(#theta) #Lambda_{c}^{+} in #Lambda_{b} RF",10,-1.0,1.0);
    TH1D *hphi = new TH1D("hphi", "#phi #Lambda_{c}^{+} in #Lambda_{b} RF",10, 0.0,2*acos(-1.0));

    f_cosT_32->SetParameter(0,a_lcs[j]);
    f_cosT_32->SetParameter(1, (2 + (2/3)*a_lcs[j]) );
   
    for (int i=0; i<N; i++){
      
      prand = r1->Uniform(0.0,1.0);

      if (prand < p1 ){

	/* Lambda_b -> Lambda_c (1/2) + mu nu */
	
	cosT_lc = f_cosT_lb->GetRandom();
	phi_lc = f_phi_lb->GetRandom();

	theta_lc = acos( cosT_lc);
       	
	q2 = r2->Uniform(pow(Mmu,2)+Mnu,pow(Mb - Mc,2));

	Plc = sqrt(pow( ( pow(Mb,2) + pow(Mc,2) - q2 )/(2*Mb) ,2 ) - pow(Mc,2) );

	phi_f = phi_lc;
	cosT_f = cos(theta_lc);

	hcosT->Fill(cosT_f);
	hphi->Fill(phi_f);
      
      }
   
      if (p1 <= prand && prand < p1 + p2 ){

	/* Lambda_b -> Lambda_c* (1/2) + mu nu */
      
	cosT_lb = f_cosT_lb->GetRandom();
	phi_lb = f_phi_lb->GetRandom();;

	theta_lb = acos( cosT_lb);
	
	q2 = r2->Uniform(pow(Mmu,2)+Mnu,pow(Mb - Mc12s,2));

	Plcs = sqrt(pow( ( pow(Mb,2) + pow(Mc12s,2) - q2 )/(2*Mb) ,2 ) - pow(Mc12s,2) );

	plcs[0] = sqrt( pow(Plcs,2) + pow(Mc12s,2) );

	tie(plcs[1], plcs[2], plcs[3]) = Impulsions(Plcs,phi_lb, theta_lb);

	beta_lcs = Plcs / plcs[0];

      
	/* Lambda_c* (1/2) -> Lambda_c (1/2) + (pi+ pi-) */

	
	cosT_12 = r5->Uniform(-1.0,1.0);
	phi_12 = r6->Uniform(0.0,2*Pi);

	theta_12 = acos( cosT_12);
	
	Plc = sqrt(pow( ( pow(Mc12s,2) + pow(Mc,2) - pow(Mpim+Mpip,2)  )/(2*Mc12s) ,2 ) - pow(Mc,2) );
	
	plc[0] = sqrt( pow(Plc,2) + pow(Mc,2) );

	tie(plc[1], plc[2], plc[3]) = Impulsions(Plc, phi_12, theta_12);
    
	/* Boost */
	
	tie(plc[0], plc[1], plc[2], plc[3]) = Boost(plc, beta_lcs);
      
	/* Rotation */
      
	tie(plc[0], plc[1], plc[2], plc[3]) = Rotation(phi_lb, theta_lb, plc);

	phi_f = atan2(plc[2],plc[1]);

	if(phi_f < 0){

	  phi_f = phi_f + 2*Pi;
	}
	
	cosT_f = plc[3]/( sqrt( pow(plc[1],2) + pow(plc[2],2) + pow(plc[3],2)) );

	hcosT->Fill(cosT_f);
	hphi->Fill(phi_f);
      
      }
      
      if ( p1 + p2 <= prand  ){

	/* Lambda_b -> Lambda_c* (3/2) + mu nu */
      
	cosT_lb =f_cosT_lb->GetRandom();
	phi_lb = f_phi_lb->GetRandom();;

	theta_lb = acos( cosT_lb);
	
	q2 = r2->Uniform(pow(Mmu,2)+Mnu,pow(Mb - Mc32s,2));

	Plcs = sqrt(pow( ( pow(Mb,2) + pow(Mc32s,2) - q2 )/(2*Mb) ,2 ) - pow(Mc32s,2) );

	plcs[0] = sqrt( pow(Plcs,2) + pow(Mc32s,2) );

	tie(plcs[1], plcs[2], plcs[3]) = Impulsions(Plcs,phi_lb, theta_lb);

	beta_lcs = Plcs / plcs[0];

	/* Lambda_c* (3/2) -> Lambda_c (1/2) + (pi+ pi-) */

	cosT_32 = f_cosT_32->GetRandom();
	phi_32 = r5->Uniform(0.0,2*Pi);
	
	Plc = sqrt(pow( ( pow(Mc32s,2) + pow(Mc,2) - pow(Mpim+Mpip,2)  )/(2*Mc32s) ,2 ) - pow(Mc,2) );
	
	plc[0] = sqrt( pow(Plc,2) + pow(Mc,2) );

	tie(plc[1], plc[2], plc[3]) = Impulsions(Plc, phi_32, theta_32);

	
	/* Boost */
	
	tie(plc[0], plc[1], plc[2], plc[3]) = Boost(plc, beta_lcs);
      
	/* ROtation */

	tie(plc[0], plc[1], plc[2], plc[3]) = Rotation(phi_lb, theta_lb, plc);

	phi_f = atan2(plc[2],plc[1]);

	if(phi_f < 0){

          phi_f = phi_f	+ 2*Pi;
	}
	
	cosT_f = plc[3]/( sqrt( pow(plc[1],2) + pow(plc[2],2) + pow(plc[3],2)) );

	hcosT->Fill(cosT_f);
	hphi->Fill(phi_f);

	
      
      }
      

    }/* End for(i... */
    
    /* Select alpha = -1, 0 or 1 and create the pdfs   */
    
    hcosT->SetAxisRange(0.0,12000,"Y");
    hphi->SetAxisRange(0.0,12000,"Y");
    /*
    hcosT->Scale(norm/hcosT->Integral(),"width");
    hphi->Scale(norm/hphi->Integral(),"width");
    
    hcosT->SetAxisRange(0.0,0.6,"Y");
    hphi->SetAxisRange(0.0,0.2,"Y");
    
    TF1 *fit_cos = new TF1("fit_cos","(1 +[0]*x)/(2)",-1.0,1.0);
    TF1 *fit_phi = new TF1("fit_phi","(1 + [0]*cos(x))/(2*acos(-1.0))",0.0,2*Pi);
    */    
/* to write the parameters of the fit in a file */
    //fit_cos->SetLineWidth(2);
    //fit_phi->SetLineWidth(2);
 
    if (j == 0){

      
      Can->cd(1);
      //hcosT->Fit("fit_cos");
      hcosT->Draw("HIST");
      
      Can->cd(2);
      //hphi->Fit("fit_phi");
      hphi->Draw("HIST");
      Can->Print("Weighted_Decays_-1_unnorm_unif.pdf");  
    }

    if (j == 1){
      

     Can->cd(1);
     //hcosT->Fit("fit_cos");
      hcosT->Draw("HIST");

      Can->cd(2);
      //hphi->Fit("fit_phi");
      hphi->Draw("HIST");
      Can->Print("Weighted_Decays_0_unnorm_unif.pdf");
    }

    if (j == 2){
      
     Can->cd(1);
     //hcosT->Fit("fit_cos");
      hcosT->Draw("HIST");

      Can->cd(2);
      // hphi->Fit("fit_phi");
      hphi->Draw("HIST");
      Can->Print("Weighted_Decays_1_unnorm_unif.pdf");
    }

    delete hcosT;
    delete hphi;
    delete Can;
    
  } /* End for( j... */

    


  }
