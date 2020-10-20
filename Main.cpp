// C++ simulation code for presynatic NMDA receptors and short-term 
// plasticity at hippocampal SC-CA1 synapse. 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <assert.h>
#include <sys/stat.h>

#include <windows.h>

using namespace std;

// g++ -DNDEBG ... to turn off assertions
// g++ -w Main.cpp -DHill -I include/
// g++ -w Main.cpp lib/*.cpp -DHill -I include/  if putting implementation in lib

static const double PI=M_PI;  // use Pi defined in math.h

#include "utilities.h";

#ifdef Hill
#include "vesicles_hill.h";  
#endif 

#ifdef Markov
#include "vesicles_markov.h"  // original matrix version
#endif 

#ifdef Markov6
#include "vesicles_markov_6.h"
#endif 

#ifdef Allosteric
#include "vesicles_allosteric.h"
#endif 

#include "er_ip3_receptor.h";
#include "er_ryr_receptor.h";
#include "er.h";

#include "bouton_receptors.h";
#include "bouton.h";

#include "save.h";
#include "score.h";

#include "simulation.h";

#include "fit.h";


int main(int argc, char* argv[])
{
//   clock_t t1,t2;
//   t1=clock();
   
   // struct Ex contains simulation parameters
   EX ex, optimal;  // experimental setup
   
   double  isi, seconds, trials;
   int AP5, RyR;
   
   double * pr_ACSF_barChart; 
   double * pr_BLOCKER_barChart;
   
   if (argc < 7) 
   {
     fprintf(stderr, "Usage: %s  isi seconds trials AP5={0,1} RyR={0,1} astro={0,1}} \n", argv[0]); 
     exit(1);
   }
   else
   {  
     isi = atof(argv[1]);
     seconds = atof(argv[2]);
     trials  = atof(argv[3]);
     
     AP5 = atoi(argv[4]);  // 0 or 1,  off or on
     RyR = atoi(argv[5]);
   } 
   
   // for unix OS
   /**
   struct stat sb;  
  
   if (stat("csv", &sb) != 0) 
   {
       printf("mkdir csv \n");
       mkdir("csv", 0700);
   }
  **/
   
   // for Windows OS
   if (CreateDirectory("csv", NULL) || ERROR_ALREADY_EXISTS == GetLastError())
    {
        printf("Creat Dir csv failed/n");
    }
   
   if (CreateDirectory("png", NULL) || ERROR_ALREADY_EXISTS == GetLastError())
    {
               printf("Creat Dir png failed/n");
    }

  // for unix OS
  /**
   if (stat("png", &sb) != 0)
   {
       printf("mkdir csv \n");
       mkdir("csv", 0700);
   }
   **/
   
   double deltaT=0.05;
   
   ex.isi=isi, ex.seconds=seconds, ex.trials=trials, ex.deltaT=deltaT, ex.astro=astro;
     
   printf("\n\n =======Sim for isi=%0.0f \n =========", isi);
   
   buildTrain(ex);   // builds spike train
   
   pr_ACSF_barChart     = init_double(ex.bins); 
   pr_BLOCKER_barChart  = init_double(ex.bins); 
   
   if(AP5 == 1) 
   {
    ex.AP5_exp=1;
    ex.RY_exp =0;
    printf(" ==> sim ACSF, AP5: isi=%0.0f \n\n", isi);
   } 
   else if (RyR == 1) 
   {
    ex.AP5_exp=0;
    ex.RY_exp = 1;
    printf(" ==> sim ACSF, RyR: isi=%0.0f \n\n", isi);
   }
   else 
   {
    printf("No blocker was specified.\n");
    return 0;
   }
   
   int fit_hill=0;
   
   if( ! fit_hill)
   {
      sim(pr_ACSF_barChart, pr_BLOCKER_barChart, ex, 1);
      if (ex.isi != 75) {
         score(pr_ACSF_barChart, pr_BLOCKER_barChart, ex);
      }
   }
   
   if(fit_hill) 
   { 
      fit(pr_ACSF_barChart, pr_BLOCKER_barChart, ex);
   }
     
   printf("----------------------------------------\n");
   printf("\n         Simulation done.\n");
//   t2=clock();
//   float diff ((float) t2-(float) t1);
//   printf("time= %f \n", diff/CLOCKS_PER_SEC);
   printf("----------------------------------------\n");
            
   return 0;
}

      