/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */
  
#define CONSTRAINTS 

//#define MONOJET
//#define HIGGSBOUNDS 
//#define HIGGSSIGNALS
//#define LILITH
//#define SMODELS
  
#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */
   
#define INDIRECT_DETECTION  
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */
      
#define RESET_FORMFACTORS
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

#define CDM_NUCLEUS     
  /* Calculate number of events for 1kg*day and recoil energy distibution
      for various nuclei
  */
#define NEUTRINO    
 /*  Neutrino signal of DM annihilation in Sun and Earth */

#define DECAYS

#define CROSS_SECTIONS 
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"
#include <string>

using namespace std;


int main(int argc,char** argv)
{  int err, i;
   char cdmName[10], lspname[10], nlspname[10];
   int spin2, charge3,cdim;
   double w;

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */

  VZdecay=0; VWdecay=0;  

/*  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
  err=readVar(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
*/           


  FILE *omega = fopen("Block_micrOMEGAs.out","w");
  fprintf(omega,"Block micrOMEGAs\n");

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
  if(CDM1) 
  { 
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2) 
  { 
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

  
#ifdef MASSES_INFO
{
   char*  pname = pdg2name(32);   /* MZp Mass */
   double mass;
     
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);

  printf("\n=== MASSES OF EXTRA GAUGE BOSONS: ===\n");
   if(pname)
   {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   401   %2.6E   #MZp\n",mass);
   }

 
    pname = pdg2name(34);      /* MWp Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   402   %2.6E   #MWp\n",mass);
  }


    pname = pdg2name(6000012);  /* Mne Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   403   %2.6E   #Mne\n",mass);
  }


    pname = pdg2name(6000014);  /* Mnm Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   404   %2.6E   #Mnm\n",mass);
  }


    pname = pdg2name(6000016);  /* Mnt Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   405   %2.6E   #Mnt\n",mass);
  }


    pname = pdg2name(6000001);  /* MDP Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   406   %2.6E   #MDD\n",mass);
  }


    pname = pdg2name(6000003);  /* MSP Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   407   %2.6E   #MDS\n",mass);
  }

    pname = pdg2name(6000005);  /* MBP Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   408   %2.6E   #MDB\n",mass);
  }

    pname = pdg2name(35);    /* MH1 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   409   %2.6E   #mh1\n",mass);
  }

    pname = pdg2name(45);    /* MH2 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   410   %2.6E   #mh2\n",mass);
  }


    pname = pdg2name(55);    /* MH3 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   411   %2.6E   #mh3\n",mass);
  }


    pname = pdg2name(36);    /* MA1 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   412   %2.6E   #a1\n",mass);
  }


    pname = pdg2name(46);    /* MA2 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   413   %2.6E   #a2\n",mass);
  }


    pname = pdg2name(37);    /* MHp1 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   414   %2.6E   #h1+\n",mass);
  }

    pname = pdg2name(47);    /* MHp2 Mass */
  if(pname)
  {
    mass=pMass(pname);
    printf("\n%s : %E \n",pname,mass);
    fprintf(omega,"   415   %2.6E   #h2-\n",mass);
  }


  printMasses(stdout,1);   


}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invizible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q q-\\bar  Cross Section= %.2E pb\n",csLim);
}
#endif

#ifdef MONOJET
{ double CL=monoJet();
  printf(" Monojet signal exclusion CL is %.3e\n", CL);   
}  
#endif

#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0,NHch;  // number of neutral and charged Higgs particles.
    int HB_id[3],HB_result[3];
   double  HB_obsratio[3],HS_observ,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch);
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif 
#ifdef HIGGSSIGNALS
   printf("HiggsSignals(%s)\n",HS_version); 
   printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
#endif   
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char call_lilith[100], Lilith_version[20];
   if(LilithMO("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n",
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   } else printf("LILITH: there is no Higgs candidate\n");
}     
#endif


#ifdef SMODELS
{  int result=0;
   double Rvalue=0;
   char analysis[30]={},topology[30]={};
   int LHCrun=LHC8|LHC13;  // LHC8  - 8TeV; LHC13 - 13TeV; LHC8|LHC13 - both  
#include "../include/SMODELS.inc" 
}   
#endif 


#ifdef OMEGA
{ int fast=1;        /* 0 = best accuracy, 1 = "fast option" accuracy ~1%  */
  double Beps=1.E-4, cut=0.01;
  double Omega;  
  int i,err; 
  printf("\n==== Calculation of relic density =====\n");   

  if(CDM1 && CDM2) 
  {
  
    Omega= darkOmega2(fast,Beps);
  /*
    displayPlot("vs1120F","T", Tend, Tstart,1,1,"vs1120F", 0,vs1120F,NULL);
    displayPlot("vs2200F","T", Tend, Tstart,1,1,"vs2200F", 0,vs2200F,NULL);
    displayPlot("vs1100F","T", Tend, Tstart,1,1,"vs1100F", 0,vs1100F,NULL);
    displayPlot("vs1210F","T", Tend, Tstart,1,1,"vs1210F", 0,vs1210F,NULL);
    displayPlot("vs1122F","T", Tend, Tstart,1,1,"vs1122F", 0,vs1122F,NULL);
    displayPlot("vs2211F","T", Tend, Tstart,1,1,"vs2211F", 0,vs2211F,NULL);

    displayPlot("vs1110F","T", Tend, Tstart,1,1,"vs1110F", 0,vs1110F,NULL);
    displayPlot("vs2220F","T", Tend, Tstart,1,1,"vs2220F", 0,vs2220F,NULL);
    displayPlot("vs1110F","T", Tend, Tstart,1,1,"vs1110F", 0,vs1112F,NULL);
    displayPlot("vs1222F","T", Tend, Tstart,1,1,"vs1222F", 0,vs1222F,NULL);
    displayPlot("vs1220F","T", Tend, Tstart,1,1,"vs1220F", 0,vs1220F,NULL);
    displayPlot("vs2210F","T", Tend, Tstart,1,1,"vs2210F", 0,vs2210F,NULL);
    displayPlot("vs2221F","T", Tend, Tstart,1,1,"vs2221F", 0,vs2221F,NULL);
    displayPlot("vs1211F","T", Tend, Tstart,1,1,"vs1211F", 0,vs1211F,NULL);
  */
  
    printf("Omega_1h^2=%.2E\n", Omega*(1-fracCDM2));
    printf("Omega_2h^2=%.2E\n", Omega*fracCDM2);
  } else
  {  double Xf;
     Omega=darkOmega(&Xf,fast,Beps,&err);
     printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
     if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);

  FILE *channels = fopen("Channels.out","w");
  if(Omega>0)printChannels(Xf,cut,Beps,1,channels);

  }





//  fprintf(omega,"Block micrOMEGAs\n");
                        fprintf(omega,"   %i   %2.6E   # relic density \n",700,Omega);
                        w = 1.;
                        i = 0;
                        while (w>cut)
                        {
                            fprintf(omega,"   %i   %2.6E   # %s %s -> %s %s\n",100+i,omegaCh[i].weight,omegaCh[i].prtcl[0],omegaCh[i].prtcl[1],omegaCh[i].prtcl[2],omegaCh[i].prtcl[3]);
                            i++;
                            w = omegaCh[i].weight;
                        }




}

#endif

#ifdef FREEZEIN
{
  double TR=1E10;
  double omegaFi;  
  toFeebleList(CDM1);
  VWdecay=0; VZdecay=0;
  
  omegaFi=darkOmegaFi(TR,&err);
  printf("omega freeze-in=%.3E\n", omegaFi);
  printChannelsFi(0,0,stdout);
}
#endif



#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */

  printf("sigmaV=%.2E [cm^3/s]\n", sigmaV);
  fprintf(omega,"   306   %2.6E   #sigmaV [cm^3/s]\n",sigmaV);





  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest); 
     fprintf(omega,"   307   %2.6E   #Photon flux [cm^2 s GeV]^{-1}\n",SpectdNdE(Etest, FluxA));

  }

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);
    fprintf(omega,"   308   %2.6E   #Positron flux [cm^2 sr s GeV]^{-1}\n",SpectdNdE(Etest, FluxE));
    
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);
    fprintf(omega,"   309   %2.6E   #Antiproton flux [cm^2 sr s GeV]^{-1}\n",SpectdNdE(Etest, FluxP));

  }
}  
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");
 
//  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
//  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);
  
//                    To restore default form factors of  version 2  call 
//     calcScalarQuarkFF(0.553,18.9,55.,243.5);


  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call 
	  calcScalarQuarkFF(0.56,20.2,34,42);


}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

  if(CDM1)
  {  
    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM1);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);

     fprintf(omega,"   201   %2.6E   #P SI [pb]\n",SCcoeff*pA0[0]*pA0[0]);
     fprintf(omega,"   202   %2.6E   #P SD [pb]\n",3*SCcoeff*pA5[0]*pA5[0]);
     fprintf(omega,"   203   %2.6E   #N SI [pb]\n",SCcoeff*nA0[0]*nA0[0]);
     fprintf(omega,"   204   %2.6E   #N SD [pb]\n",3*SCcoeff*nA5[0]*nA5[0]);

  }
  if(CDM2)
  {
    nucleonAmplitudes(CDM2, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM2);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);

  }     
}
#endif
  
#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;
  double nEventsCut;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n", nEventsCut);
  fprintf(omega,"   301   %2.6E   #73Ge: Total number of events [day/kg]\n",nEvents);
                                                                               #ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 73Ge","E[KeV]",0,200,0,1,"dN/dE",0,dNdERecoil,dNdE);
#endif





  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n", nEventsCut);
  fprintf(omega,"   302   %2.6E   #131Xe: Total number of events [day/kg]\n", nEvents);

#ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 131Xe","E[KeV]",0,200,0,1,"dN/dE",0,dNdERecoil,dNdE);
#endif


  



  
  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n", nEventsCut);
  fprintf(omega,"   303   %2.6E   #23Na: Total number of events [day/kg]\n",nEvents);

#ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 23Na","E[KeV]",0,200,0,1,"dN/dE",0,dNdERecoil,dNdE);
#endif






  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n", nEventsCut);
  fprintf(omega,"   304   %2.6E   #I127: Total number of events [day/kg]\n",nEvents);


#ifdef SHOWPLOTS
  displayPlot("Distribution of recoil energy of 127I","E[KeV]",0,200,0,1,"dN/dE",0,dNdERecoil,dNdE);
#endif
  


}
#endif 

#ifdef NEUTRINO
if(!CDM1 || !CDM2)
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;


  neutrinoFlux(Maxwell,forSun,nu,nu_bar);
  double CLs = 100*exLevIC22(nu,nu_bar,NULL);
  printf("IceCube22 exclusion confidence level = %.2E%%\n",CLs);
  fprintf(omega,"   305   %2.6E   #IceCube22 exclusion confidence level\n",CLs);



  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0, 2,"dnu/dE",0,SpectdNdE,nu,"dnu_bar/dE",0,SpectdNdE,nu_bar);
#endif
{ 
    printf(" E>%.1E GeV neutrino flux       %.2E [1/Year/km^2] \n",Emin,spectrInfo(Emin,nu,NULL));
    printf(" E>%.1E GeV anti-neutrino flux  %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,nu_bar,NULL));  
} 
  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS
  displayPlot("Upward muons[1/Year/km^2/GeV]","E",Emin,Mcdm/2, 0,1,"mu",0,SpectdNdE,mu);  
#endif
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS 
  displayPlot("Contained  muons[1/Year/km^3/GeV]","E",Emin,Mcdm,0,1,"",0,SpectdNdE,mu); 
#endif
  printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin/Mcdm,mu,NULL));
}        
#endif 


#ifdef DECAYS
{ char*  pname = pdg2name(25);
  txtList L;
  double width;

  FILE *TotalWidths_FILE = fopen("TotalWidths.out","w");
  FILE *DECAY_FILE = fopen("DECAYS.out","w");
	
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);
   } 


  pname = pdg2name(6);    /* Top Quark Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);
  }


  pname = pdg2name(23);    /* MZ Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);
  }

  pname = pdg2name(24);    /* MW Mass */  
  if(pname)
  { 
    width=pWidth(pname,&L);  
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  } 


      pname = pdg2name(32);    /* MZp Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);    
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

    pname = pdg2name(34);      /* MWp Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);    
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

    pname = pdg2name(35);     /* MH1 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);    
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(36);   /* MA1 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);    
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);
  }

      pname = pdg2name(37);   /* MHp1 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);    
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(45);   /* MH2 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(46);   /* MA2 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }
  
      pname = pdg2name(47);   /* MHp2 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(55);   /* MH3 Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);
  }

      pname = pdg2name(6000001);   /* MDP Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(6000003);   /* MSP Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(6000005);   /* MBP Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

      pname = pdg2name(6000012);   /* Mne Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }

        pname = pdg2name(6000014);   /* Mnm Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);
  }

        pname = pdg2name(6000016);   /* Mnt Mass */
  if(pname)
  {
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname,width);
    fprintf(TotalWidths_FILE, "%E #total width of %s\n",width,pname);
/*    printTxtList(L,DECAY_FILE); */
    printPartialWidth(width,L,DECAY_FILE);    
  }  

    fclose(TotalWidths_FILE);
    fclose(DECAY_FILE);
}            
#endif


/* double ZpLogic = Zprimelimits(); */



#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;
    
  next=nextOdd(1,&nextM); 
/*  if(next && nextM<1000)  
  { 
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  } 
*/


/* Cross section of W + jet */ 
  char*  outgoing1 = pdg2name(24); /* W boson      */
  char*  outgoing2 = pdg2name(1);  /* Down Quark    */  
  char*  outgoing3 = pdg2name(2);  /* Up Quark      */
  char*  outgoing4 = pdg2name(3);  /* Strange Quark */
  char*  outgoing5 = pdg2name(4);  /* Charm Quark   */  
  char*  outgoing6 = pdg2name(5);  /* Bottom Quark  */
  char*  outgoing7 = pdg2name(6);  /* Top Quark     */  
  char*  outgoing8 = pdg2name(21); /* Gluon         */

  char*  outgoing9 = pdg2name(34);  /* WR boson */
  char*  outgoing10 = pdg2name(6000001); /* dd boson */
  char*  outgoing11 = pdg2name(32); /* Zp boson */


/*  if(outgoing1 && nextM<1000)
  {
     double cs_tot, cs1, cs2, cs3, cs4, cs5, cs6, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);
     Qren=Qfact;

     cs1=hCollider(Pcm,1,nf,Qren, Qfact, outgoing1,outgoing2,pTmin,1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing1,outgoing2, cs1);

     cs2 = hCollider(Pcm,1,nf,Qren, Qfact, outgoing1,antiParticle(outgoing3),pTmin,1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing1,antiParticle(outgoing3), cs2);

     cs3 = hCollider(Pcm,1,nf,Qren, Qfact, outgoing1,outgoing4,pTmin,1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing1,outgoing4, cs3);

     cs4 = hCollider(Pcm,1,nf,Qren, Qfact, outgoing1, antiParticle(outgoing5), pTmin,1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing1, antiParticle(outgoing5), cs4);

     cs5 = hCollider(Pcm,1,nf,Qren, Qfact, outgoing1, outgoing6, pTmin,1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing1, outgoing6, cs5);

     cs6 = hCollider(Pcm,1,nf,Qren, Qfact, outgoing1, outgoing8, pTmin,1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing1, outgoing8, cs6);

     cs_tot = cs1 + cs2 + cs3 + cs4 + cs5 + cs6;
     printf("cs(pp-> %s, jet(%s,%s,%s,%s,%s,%s))=%.2E[pb]\n",outgoing1,outgoing2,outgoing3,outgoing4,outgoing5,outgoing6,outgoing8,cs_tot);     
  }
*/

  /* Cross section of WR WR */
  if(outgoing1)
  {
     double cs_tot, cs1, cs2, cs3, cs4, cs5, cs6, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);
     Qren=Qfact;

     cs1=hCollider(Pcm,1,nf,Qren, Qfact, outgoing9, antiParticle(outgoing9), pTmin, 1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing9, antiParticle(outgoing9), cs1);
     fprintf(omega, "   900   %E   #cs(pp-> %s,%s)\n", cs1, outgoing9, antiParticle(outgoing9));
  }


  
   /* Cross section of WR dd */
  if(outgoing1)
  {
     double cs_tot, cs1, cs2, cs3, cs4, cs5, cs6, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);
     Qren=Qfact;

     cs1=hCollider(Pcm,1,nf,Qren, Qfact, outgoing9, outgoing10, pTmin, 1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing9, outgoing10, cs1);
     fprintf(omega, "   901   %E   #cs(pp-> %s,%s)\n", cs1, outgoing9, outgoing10);
  }



  /* Cross section of dd dd */
  if(outgoing1)
  {
     double cs_tot, cs1, cs2, cs3, cs4, cs5, cs6, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);
     Qren=Qfact;

     cs1=hCollider(Pcm,1,nf,Qren, Qfact, outgoing10, antiParticle(outgoing10), pTmin, 1);
     printf("cs(pp-> %s,%s)=%.2E[pb]\n",outgoing10, antiParticle(outgoing10), cs1);
     fprintf(omega, "   902   %E   #cs(pp-> %s,%s)\n", cs1, outgoing10, antiParticle(outgoing10));
  }



  /* Cross section of Zp */
  if(outgoing1)
  {
     double cs_tot, cs1, cs2, cs3, cs4, cs5, cs6, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     Qren=Qfact=nextM;

     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);
     Qren=Qfact;

     cs1=hCollider(Pcm,1,nf,Qren, Qfact, outgoing11, NULL, pTmin, 1);
     printf("cs(pp-> %s)=%.2E[pb]\n",outgoing11, cs1);
     fprintf(omega, "   903   %E   #cs(pp-> %s)\n", cs1, outgoing11);
  }


}
 
#endif 

#ifdef CLEAN
  system("rm -f HB.* HB.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.in  smodels.log  smodels.out  summary.*");  
#endif 


  killPlots();

  fclose(omega);

  return 0;
}
