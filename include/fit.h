


EX fit(double * pr_ACSF_barChart, double * pr_BLOCKER_barChart, EX &ex)
{
   // first 10 spikes are estimates based on Fig 10 of McGuinness 2010.
   //
   double pr_ACSF_data_1[11]     = { 10000,  100, 120, 120, 120, 120, 120, 120,  120, 120, 120 }; 
   double pr_BLOCKER_data_1[11]  = { 10000,  100, 120, 120, 120, 120, 120, 120,  120, 120, 120 }; 
   
   double pr_ACSF_data_5[11]     = { 10000,  100, 160, 240, 300, 300, 300, 300,  300, 300, 300  }; 
   double pr_BLOCKER_data_5[11]  = { 10000,  100, 150, 150, 150, 150, 150, 150,  150, 150, 150  }; 
   
   double pr_ACSF_data_20[11]    = { 10000,  100, 150, 190, 200, 250, 250, 250,  250, 250, 250 }; 
   double pr_BLOCKER_data_20[11] = { 10000,  100, 150, 190, 190, 190, 190, 190,  190, 190, 190 }; 
   
   double minMSE=100000000;
   double MSE=0, SE=0;  // mean squared error and squared error
 
   double start_n1=4, end_n1=5, stepSize=1.0;
   
   double start_choice=3, end_choice=3;

   int maxSpikeNumber= ( sizeof(pr_ACSF_data_1)/sizeof(double) )-1;
    
   printf("max spike num %d \n", maxSpikeNumber);
   
   EX ex1,ex2,ex3, optimal;
   
   buildTrain(ex1);
   buildTrain(ex2);
   buildTrain(ex3);

   double n1=0;
   
   for(n1=start_n1; n1 <= end_n1; n1 = n1 + stepSize)  
   { 
       printf("testing n1=%0.1f \n", n1);
       
       for(int choice=start_choice; choice <= end_choice; ++choice) 
       { 
         if(choice == 1) {  
          
            printf("choice=%d \n", choice);
            ex1.n1=n1;
              
            //sim(pr_ACSF_barChart, pr_BLOCKER_barChart, ex1, 0);
         }
         else if(choice == 2) {
            
            printf("choice=%d \n", choice);
            ex2.n1=n1;
            
            //sim(pr_ACSF_barChart, pr_BLOCKER_barChart, ex2, 0);
         }
         else if(choice == 3) {  
          
            printf("choice=%d \n", choice);
            ex3.n1=n1;   
            
            sim(pr_ACSF_barChart, pr_BLOCKER_barChart, ex3, 0);
         } 
         else {
            printf("\n Error: choice not available \n");
            exit(0);
         }  
            
         // mean squared error for this Hill coefficient n1
         SE=0;
         
         for(int i=1; i <= ex.bins && i <= maxSpikeNumber; ++i) 
         { 
            if (choice==1) {
               SE += pow( 100*(pr_ACSF_barChart[i]     - pr_ACSF_data_1[i])/pr_ACSF_data_1[i], 2); 
               SE += pow( 100*(pr_BLOCKER_barChart[i]  - pr_BLOCKER_data_1[i])/pr_BLOCKER_data_1[i],   2);
            }
            else if (choice==2) {
               SE += pow( 100*(pr_ACSF_barChart[i]     - pr_ACSF_data_5[i])/pr_ACSF_data_5[i],         2); 
               SE += pow( 100*(pr_BLOCKER_barChart[i]  - pr_BLOCKER_data_5[i])/pr_BLOCKER_data_5[i],   2);
            }   
            else if (choice==3) {
               SE += pow( 100*(pr_ACSF_barChart[i]     - pr_ACSF_data_20[i])/pr_ACSF_data_20[i], 2); 
               SE += pow( 100*(pr_BLOCKER_barChart[i]  - pr_BLOCKER_data_20[i])/pr_BLOCKER_data_20[i], 2);
            }      
         } 
      } // end of choices, frequencies 1, 5, 20Hz
      
      int sets=6;
      
      MSE = SE/(ex.bins*sets);  // * 6 because 3 frequencies, same number of bins, 2 conditions each
      
      if (MSE < minMSE) { 
          minMSE=MSE;
          optimal.n1 = n1;  
      }
      printf("minMSE=%0.3f, MSE=%0.3f:  n1=%0.2f \n", minMSE, MSE, ex.n1);
   } // end of n1 loop
      
   printf("Best MSE=%0.3f,  best n1=%0.3f \n", minMSE, optimal.n1); 
   
   return optimal;
}

