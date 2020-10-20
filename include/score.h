

double score(double * pr_ACSF_barChart, double * pr_BLOCKER_barChart, EX ex)
{
   // first 10 spikes are estimates based on Fig 10 of McGuinness 2010.
   double pr_ACSF_data_1[11]     = { 10000,  100, 120, 120, 120, 120, 120, 120,  120, 120, 120 }; 
   double pr_BLOCKER_data_1[11]  = { 10000,  100, 120, 120, 120, 120, 120, 120,  120, 120, 120 }; 
   
   double pr_ACSF_data_5[11]     = { 10000,  100, 160, 240, 300, 300, 300, 300,  300, 300, 300  }; 
   double pr_BLOCKER_data_5[11]  = { 10000,  100, 150, 150, 150, 150, 150, 150,  150, 150, 150  }; 
   
   double pr_ACSF_data_20[11]    = { 10000,  100, 150, 190, 200, 250, 250, 250,  250, 250, 250 }; 
   double pr_BLOCKER_data_20[11] = { 10000,  100, 150, 190, 190, 190, 190, 190,  190, 190, 190 }; 
   
   double minMSE=10000000;
   double MSE=0, SE=0;
    
   double start_choice=3, end_choice=3;

   int maxSpikeNumber= ( sizeof(pr_ACSF_data_1)/sizeof(double) )-1;
   
      
   for(int i=1; i < ex.bins && i <= maxSpikeNumber; ++i) {
    
      if (ex.isi == 1000) {
         SE += pow( 100*(pr_ACSF_barChart[i]     - pr_ACSF_data_1[i])/pr_ACSF_data_1[i], 2); 
         SE += pow( 100*(pr_BLOCKER_barChart[i]  - pr_BLOCKER_data_1[i])/pr_BLOCKER_data_1[i],   2);
      }
      else if (ex.isi == 200) {
         SE += pow( 100*(pr_ACSF_barChart[i]     - pr_ACSF_data_5[i])/pr_ACSF_data_5[i],         2); 
         SE += pow( 100*(pr_BLOCKER_barChart[i]  - pr_BLOCKER_data_5[i])/pr_BLOCKER_data_5[i],   2);
      }   
      else if (ex.isi == 50) {
         SE += pow( 100*(pr_ACSF_barChart[i]     - pr_ACSF_data_20[i])/pr_ACSF_data_20[i], 2); 
         SE += pow( 100*(pr_BLOCKER_barChart[i]  - pr_BLOCKER_data_20[i])/pr_BLOCKER_data_20[i], 2);
      }  
      else {
         ; //printf("score.h:  no code for scoring for 75 ms PPF yet \n");    
      } 
   }
 
  int sets=2;
  
  MSE = SE/(10*sets);  // 10 bins per experiment * 2 conditions
   
  printf("\n\n ===> mean MSE=%0.3f  for isi=%0.0f\n\n", MSE, ex.isi);
  
  return MSE;
}

