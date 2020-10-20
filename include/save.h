



void save(double * data, int tn, const char * filename)
{
    FILE * fp = fopen( filename, "w+" ); // Open file for writing
    
    int i;
    for(i=1;i < tn; ++i )
    {
        fprintf(fp, "%f,", data[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void save_bins(double * data, int bins, const char * filename)
{
    FILE * fp = fopen( filename, "w+" ); // Open file for writing
    
    int i;
    for(i=1;i <= bins; ++i )
    {
        fprintf(fp, "%f,", data[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void save_int(int * data, int tn, const char * filename)
{
    FILE * fp = fopen( filename, "w+" ); // Open file for writing
    
    int i;
    for(i=1;i < tn; ++i )
    {
        fprintf(fp, "%d,", data[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}



void save_acsf(Bouton B, double * pr_ACSF_barChart, double * pr_ACSF_barChart_raw, EX ex) 
{
     save( B.v,         ex.tn, (const char *) "csv/bv.csv" ); 
     save( B.ca_global, ex.tn, (const char *) "csv/bc.csv"); 
     save( B.ves.G_syn, ex.tn, (const char *) "csv/bg.csv");
     
     save( B.ves.VR_event, ex.tn, (const char *)  "csv/b_vr.csv");
     
     save( B.Ivgcc,       ex.tn, (const char *)    "csv/b_ca_Ivgcc.csv"     ); 
     save( B.Inmda,       ex.tn, (const char *)    "csv/b_ca_Inmda.csv"     ); 
                 
     save( B.ca_VGCC,     ex.tn, (const char *)      "csv/b_ca_VGCC.csv"     );
     save( B.ca_PreNMDAR, ex.tn, (const char *)      "csv/b_ca_PreNMDAR.csv" );
     save( B.ca_PreNMDAR_mean, ex.tn, (const char *) "csv/b_ca_PreNMDAR_mean.csv" );
     
     save( B.ves.Ca_MD,  ex.tn, (const char *) "csv/b_ca_MD.csv");             
     
     save( B.ca_RyR, ex.tn, (const char *)     "csv/b_ca_RyR.csv" );
     save( B.er.cer,    ex.tn, (const char *)  "csv/b_cer.csv" ); 
     
     char fn[50];
     char fn_raw[50];
     
     int nnn=((int)1000/ex.isi);
     sprintf(fn,     "csv/pr%dHz.csv", nnn );
     sprintf(fn_raw, "csv/pr%dHz_raw.csv", nnn );
                 
     save_bins( pr_ACSF_barChart,     ex.bins, (const char *)  fn);  
     save_bins( pr_ACSF_barChart_raw, ex.bins, (const char *)  fn_raw); 
     save_bins( pr_ACSF_barChart,     ex.bins, (const char *)  "csv/pr.csv"); 
     save_bins( pr_ACSF_barChart_raw, ex.bins, (const char *)  "csv/pr_raw.csv"); 
                    
     save(B.ca_VGCC_RyR,             ex.tn,   (const char *)  "csv/b_ca_vgcc_ryr.csv"); 

     for(int i=1; i <= ex.tn; ++i) { 
      
       B.ves.RRP[i] = B.ves.RRP[i]/ex.trials;
     }
             
     save(B.ves.P_release, ex.tn, (const char *)  "csv/b_ves_P_release.csv"); 
     
//     save(A.a_ip3, ex.tn, (const char *)  "csv/a_ip3.csv");
//     save(A.ca, ex.tn,    (const char *)  "csv/a_ca.csv");
//     save(A.aG_syn, ex.tn,(const char *)  "csv/a_Gsyn.csv");  
//     save(S.Vm, ex.tn,    (const char *)  "csv/s_Vm.csv");
}



void save_blocker(Bouton B, double * pr_BLOCKER_barChart, double * pr_BLOCKER_barChart_raw, EX ex)
{
   save(B.ves.Ca_MD,             ex.tn, (const char *)  "csv/b_ca_MD_BLOCKER.csv"); 
   save(B.ves.P_release_BLOCKER, ex.tn, (const char *)  "csv/b_ves_P_release_BLOCKER.csv"); 
   char fn[50];
   char fn_raw[50];
   int nnn=((int)1000/ex.isi);
   sprintf(fn,    "csv/prBLOCKER%dHZ.csv", nnn );
   sprintf(fn_raw,"csv/prBLOCKER%dHZ_raw.csv", nnn );
   
   save_bins( pr_BLOCKER_barChart,      ex.bins, (const char *) fn);
   save_bins( pr_BLOCKER_barChart_raw,  ex.bins, (const char *) fn_raw);  
   save_bins( pr_BLOCKER_barChart,      ex.bins, (const char *) "csv/prBLOCKER.csv"); 
   save_bins( pr_BLOCKER_barChart_raw,  ex.bins, (const char *) "csv/prBLOCKER_raw.csv"); 
}





