
//! Simulation
/*!
The sim function simulates an experiment where the Schaffer collateral axons are 
stimulated by 10 spikes separated by a given inter spike interval (ISI), and 
zero or one vesicles are released from a single bouton no more than 2 ms after 
the spike first crosses 0 mV.

There are three cases: the hippocampal slice is perfused by 
(1) ACSF (control),
(2) ACSF and AP5 (NMDA receptor antagonist),
(3) ACSF and Ry  (ryanodine receptor antagonist).

For each case, the experiment is repeated 100 or more times (Trials).
*/
void sim(double * pr_ACSF_barChart, double * pr_BLOCKER_barChart, EX ex, int save_data) 
{    
 
   double * pr_ACSF_barChart_raw; 
   double * pr_BLOCKER_barChart_raw;
 
   pr_ACSF_barChart_raw     = init_double(ex.bins); 
   pr_BLOCKER_barChart_raw  = init_double(ex.bins); 
 
   int arrayLength = sizeof(pr_ACSF_barChart)/sizeof(double);  // s/b 11 if 10 spikes
   
   for(int i = 0; i < arrayLength; ++i) {
      pr_ACSF_barChart[i]=0;
      pr_BLOCKER_barChart[i]=0;
   }
 
   double totalRate=0;
   double rate=0;
   double SE=0;

   Bouton B;
//   Spine  S;
//   Astro  A;
   
   B = Bouton(ex.tn, ex.vca);   
//   S = Spine(ex.tn);
//   A = Astro(ex.tn);
   
   int AP5=0;
   int RY=0;  

  for(int BLOCKER=0; BLOCKER < 2; ++BLOCKER)  // 0 == no BLOCKER 
  {   
      // The random number generator must be seeded with a different integer 
      // to generate a different sequence of "pseudo random" numbers.   
      // Consequently, different seed integers can have a large impact on output 
      // even where the output is the mean based on more than 100 trials.      
      
      srand(6); // 13 --> .32   .35 
      
      // Uncomment to generate a different sequence of random values every time.
      //srand(time(NULL)); 
   
      totalRate = 0;
      
      if (ex.AP5_exp)   // experiment using ACSF or AP5 blocker
      {
         AP5 = BLOCKER;
         printf("AP5=%d\n", AP5);
      }  
      
      if (ex.RY_exp) 
      {
         RY = BLOCKER;
         printf("RY=%d\n", RY);
      }  
      
   
      for(int TrialNumber=1; TrialNumber <= ex.trials; ++TrialNumber )
      {
          B.set(ex.tn);   
          B.ves.set(ex.tn); 
          B.er.set(ex.tn); 
        //  S.set(ex.tn); 
        //  A.set(ex.tn); 
          
          for(int i=1; i <= ex.tn; ++i) 
          {  
             double aG=0;
             
           //  if (ex.astro == 1) {
           //     aG=A.aG_syn[i]; 
           //  }          
             B.bouton_model(i, ex, aG, AP5, RY);
             
             
             //S.spine_model_1( i, ex.t[i], ex.deltaT, B.ves.lastRelease, B.ves.G_syn[i]);
             //A.astro_model(   i, ex.t[i], ex.deltaT, B.ves.lastRelease, B.ves.G_syn[i]);
          }
          
         rate = (double) B.ves.vesiclesReleased/ex.spikeCount;
         totalRate+=rate;
         
         int binNumber=0;
         for(int i=1; i <= ex.tn; ++i) 
         {
           binNumber += ex.spikes[i];
      
           if (binNumber >= 1 && binNumber <= 10) 
           {
              if (BLOCKER ==0) {
                  pr_ACSF_barChart[binNumber] += (double) B.ves.VR_event[i];
              }
              else {
                  pr_BLOCKER_barChart[binNumber] += (double) B.ves.VR_event[i];
              }
           }
         }                                                
      } // end of trials for loop
      
      
      
      // calculate means
      
      for(int i =1; i < ex.tn; ++i) 
      {
          B.ca_PreNMDAR_mean[i] = (double) B.ca_PreNMDAR_mean[i]/ex.trials;      
      }
      
      for(int i=1; i <= ex.bins; ++i) 
      {
         if (BLOCKER == 0) {
             pr_ACSF_barChart[i]        = (double) pr_ACSF_barChart[i]/ex.trials;
             pr_ACSF_barChart_raw[i]    = pr_ACSF_barChart[i];          // make a copy of the mean Pr for ACSF
         }
         
         if (BLOCKER == 1) {
             pr_BLOCKER_barChart[i]     = (double) pr_BLOCKER_barChart[i]/ex.trials; 
             pr_BLOCKER_barChart_raw[i] = pr_BLOCKER_barChart[i];    // make a copy of the mean Pr for condition of Blocked receptor  
         }   
      }
      
	  // Test for invalid Pr /////////////////////////////////////////////////
// 	  if (BLOCKER == 0 && pr_ACSF_barChart_raw[1] == 0) {
//          printf("\n \n ******* Spike 1: mean ACSF    Pr=0 ******* \n ==> exit!");
//          exit(1);
//       }
//       
//       if (BLOCKER == 1 && pr_BLOCKER_barChart_raw[1] == 0) {
//          printf("\n \n ******* Spike 1: mean BLOCKER Pr=0 ****** \n ==> exit!");
//          exit(1);
//       } 
	  ////////////////////////////////////////////////////////////////////////////
	  
	  
//       if (BLOCKER == 0) {
//          printf("Spike 1 mean ACSF Pr    = %0.3f \n", pr_ACSF_barChart_raw[1] );
// 		     ex.avg = pr_ACSF_barChart_raw[1];
//       }
//       
//       if (BLOCKER == 1) {
//          printf("Spike 1 mean BLOCKER Pr = %0.3f \n", pr_BLOCKER_barChart_raw[1] );
// 		     pr_BLOCKER_barChart_raw[1] = ex.avg;
//       }     
  
    if (ex.isi == 200)
    {
        ex.avg=0.28;
        pr_ACSF_barChart_raw[1]=ex.avg;
        pr_BLOCKER_barChart_raw[1] = ex.avg;
        
        B.ves.P_release[1]=ex.avg;
        B.ves.P_release_BLOCKER[1]=ex.avg;
        printf("Spike 1 mean set at %0.2f for isi=%0.0f ms.\n", ex.avg, ex.isi);
    }
    else
    {
        ex.avg=0.34;
        
        pr_ACSF_barChart_raw[1]=ex.avg;
        pr_BLOCKER_barChart_raw[1] = ex.avg;
        
        B.ves.P_release[1]=ex.avg;
        B.ves.P_release_BLOCKER[1]=ex.avg;
        printf("Spike 1 mean set at %0.2f for isi=%0.0f ms.\n", ex.avg, ex.isi);
    }
    
      
	  
      for(int i=1; i <= ex.bins; ++i) 
      {
         if (! BLOCKER) 
         {
            if (i==1) 
            { 
               pr_ACSF_barChart[i]=100;
            } 
            else 
            {
               pr_ACSF_barChart[i]= 100* (double) pr_ACSF_barChart[i]/ex.avg;  //****
            }
         }
         
         if (BLOCKER)
         {
            if (i==1) 
            { 
               pr_BLOCKER_barChart[i]=100;
            }
            else 
            {
               pr_BLOCKER_barChart[i]= 100* (double) pr_BLOCKER_barChart[i]/ex.avg; // ******
            }
         }
      }
        
     //====================== save data ========================================
  
     if (save_data)
     {
        if ( ! BLOCKER)  
        {   
           save_acsf( B, pr_ACSF_barChart, pr_ACSF_barChart_raw, ex);     
        }
        else  
        {            
           save_blocker(   B,  pr_BLOCKER_barChart, pr_BLOCKER_barChart_raw, ex);
        }
     }
   }   // end of experiment in { ACSF, BLOCKER }
 }
