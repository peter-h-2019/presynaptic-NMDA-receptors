#include <math.h>
#include <stdlib.h>  

struct EX {
               // For Hill equation based calcium sensor.
  double n1;   // Hill coefficient for sensor 1:  activator 
  double n2;   // Hill coefficient for sensor 2:  repressor
  double Kd1;  // Dissociation constant: low affinity
  double Kd2;  // Dissociation constant: very low affinity
 
 
  double avg;  // Use same Pr for spike 1 and normalisation to make valid comparisons.
 
  double isi;     // Inter spike interval.
  double seconds; // Total time in seconds.
  double trials;  // Number of trials of each experiment.
  double deltaT;  // Fixed Time-step; unit: ms,   probably 0.05 ms.
  int astro;      // couple astrocyte {0,1}? 
 
  double base;
  double Tmax;    // Time-duration of simulation; unit: ms

  int tn;         // Total number of time-points
  
  double bins;    // number of bins (spike time points), probably 10, but just 2 in case of PPF
  
  double * Iapp;  // Applied current density.
  double * t;     // because init_double allocates +2
  double tme;
  double beg_pad;  // padding in ms before stimulus
  double end_pad;  // padding in ms after stimulus
  
  
  double lastSpike;         // time point of last spike (most recent spike)
  double * lastSpikeTimes;
  int    * spikes;
  int    spikeCount;
  
  int AP5_exp, RY_exp;
  
  double  rIP3;

  double Ca_ex;  // extracellular [Ca] mM
  double   vca;  // 130.65 if [Ca]ex = 3mM, used Nernst Eq., 125 if 2mM
  
  double ACSF_1000_isi_base_Pr;
  double ACSF_200_isi_base_Pr;
  double ACSF_50_isi_base_Pr;
  
  double base_Pr;
};


//void buildTrain(Ex &ex);
EX buildPoissonTrain(double isi, double seconds, double trials, double deltaT, int calciumTest); // fix me!

double * init_double(int Len);
int    * init_int(int Len);

int      Poisson(double mean);
int      Poisson2(double lambda);

inline double rnd();

int markov(double pp[]);
double heaviside(double d);

// regular frequency spike train
//
void buildTrain(EX &ex) { 
                  // Parameters for Hill function.
    ex.n1=4.0;    // Hill coefficient for vesicle's "calcium sensor" (synaptotagmin)
    ex.n2=2.0;
    
    ex.Kd1=10;    // Dissociation constant for Ca2+ binding to sensor (synaptotagmin)
    ex.Kd2=15;
    
    ex.avg=0.37;  // assumed average Pr   *********************************
	

    ex.beg_pad=5;   // ms, arbitrary amount of time to show state in plots before first spike
    ex.end_pad=100; // padding 

    ex.Tmax = ex.beg_pad + ex.end_pad + (ex.seconds * 1000);  // Time-duration of simulation; unit: ms
  
    ex.tn=(int) ((double) ex.Tmax/ex.deltaT) + 1;             // Total number of time-points, rounded up
    
    //printf("\n utils 99, Tmax=%f deltaT=%f  tn=%d \n", ex.Tmax, ex.deltaT, ex.tn);
  
    ex.bins= 10; 
     
    ex.t   =init_double(ex.tn);  
    ex.Iapp=init_double(ex.tn);  // Applied current density
    
    ex.lastSpike = -1000;
    ex.lastSpikeTimes = init_double(ex.tn);
    ex.spikes = init_int(ex.tn);
    ex.spikeCount = 0;
    
    ex.rIP3 = 0.5;
    ex.Ca_ex= 3;       // Extracellular [Ca] mM,  2 mM is typical.
    ex.vca  = 130.65;  // 130.65 if [Ca]ex = 3 mM as in McGuinness 2010, used Nernst Eq., 125 if 2 mM   

    double tme=0;
    
    for(int i=1; i <= ex.tn; ++i)
    {
       ex.t[i]=tme;         // at t[1] time = 0
       tme+=ex.deltaT;
    } 
    
    for(int i=1; i <= ex.tn; ++i)
    {
        // Applied Current Frequency
        // Regular Spike: 5 Hz with duration of 5 * deltaT, t in 'ms'.  
        
        if (ex.t[i] < ex.beg_pad || ex.t[i] > (ex.Tmax - ex.end_pad)  ) 
        {
            ex.Iapp[i]=0; 
        }
        else if ( fmod(ex.t[i], ex.isi) >= 0 && fmod(ex.t[i],ex.isi) <=4  ) 
        {
            ex.Iapp[i]=10;  // Applied current density
          
            if ( ex.t[i] - ex.lastSpike > 6.34  )  //6.34 restricts max freq to ~157 Hz
            {   
               ex.spikes[i]=1;
               ex.lastSpike=ex.t[i];
               ex.spikeCount++;
            }
        }
        else  
        {
            ex.Iapp[i]=0; 
        }
        ex.lastSpikeTimes[i]=ex.lastSpike;
    }
    return;
}

// Poisson frequency spike train
//
EX buildPoissonTrain(double isi, double seconds, double trials, double deltaT, int calciumTest) 
{
  EX ex;
  ex.seconds=seconds;
  ex.isi=isi;
  ex.trials=trials;
  ex.deltaT=deltaT;
  
  ex.Tmax = ex.seconds * 1000;    // Time-duration of simulation; unit: ms
  ex.tn = (int) ex.Tmax/ex.deltaT;  // Total number of time-points
  // ex.bins=ex.Tmax/ex.isi;
  
  ex.t = init_double(ex.tn);   // because init_double allocates +2
  ex.Iapp=init_double(ex.tn);  // Applied current density
  
  ex.lastSpike = -1000;
  ex.lastSpikeTimes = init_double(ex.tn);
  ex.spikes = init_int(ex.tn);
  ex.spikeCount = 0;
  
  ex.rIP3=0.5;
  ex.Ca_ex=3;     // extracellular [Ca] mM
  ex.vca=130.65;  // 130.65 if [Ca]ex = 3mM, used Nernst Eq., 125 if 2mM   
  
  double tme=0;
  
  for(int i=1; i <= ex.tn; ++i)
  {
     ex.t[i]=tme;         // at t[1] time = 0
     tme+=ex.deltaT;
  } 
  double pSpike=0;
  
  double lambda = 1/(isi/deltaT);
  double pulseWidth = 4/deltaT;
  
  for(int i=1; i <= ex.tn; ++i)
  {    
      pSpike = Poisson(lambda);
    
      // want first spike at time zero, sames as regular spike trains
      if (ex.spikeCount == 0 || pSpike > 0) {          
          for(int j=0; j <= pulseWidth && i+j <= ex.tn; ++j, ++i) 
          { 
             ex.Iapp[i]=10;                            // Applied current density
           
             if ( ex.t[i] - ex.lastSpike > 6.34  )  {  // 6.34 means max ~157 Hz
                
                ex.spikes[i]=1;
                ex.lastSpike=ex.t[i];
                ex.spikeCount++;
                //printf("%f, ", pSpike);
             }
          }
          --i;
      }
      ex.bins=ex.spikeCount;
      ex.lastSpikeTimes[i]=ex.lastSpike;
  }
  //printf("bins = %f \n ", ex.bins);
  return ex;
}



int * init_int(int Len) 
{
  int * temp = new int[Len+2];  
  
  for(int i=0; i <= Len; ++i) 
  {
    temp[i]=0;
  }
  return temp;
}



double * init_double(int Len) 
{
  double * temp = new double[Len+2]; 
  
  for(int i=0; i <= Len; ++i) 
  {
    temp[i]=0;
  }
  return temp;
}


int Poisson(double mean) //Special technique required: Box-Muller method...
{
  double R;
  double sum = 0;
  int i;
  i=-1;
  double z;
  
  while(sum <=mean)
  {
     R = 0;   //  (double) rand()/(double)(RAND_MAX);
     z = -log(R);
     sum+= z;
     i++;
  }
  return i;
}


int Poisson2(double lambda)
{
  int k=0;
  double L=exp(-lambda), p=1;
  
  do {
    ++k;
    p *= (double) rand()/(double)(RAND_MAX);
  } while (p > L);
  
  return --k;
}



// inline small frequently called functions for speed
//
inline double rnd()
{
  return (double)rand() / (double)RAND_MAX ;
}

inline int markov(double pp[])
{
  double u=rnd();    // uniform random number (0 to 1)
  int i=1;           // Initial value of i
  double s=pp[1];    // Probability to stay in the present state
    
  while ((s < u ) && (i < 7))
  {
    i=i+1; 
    s=s+pp[i];
  }
     
  return i;  // present state.
}


inline double heaviside(double d)
{
  if (d < 0) 
  {
     return 0;
  }
  else if (d > 0) 
  {
    return 1;
  }
  else 
  {
   return 0.5;
  }
}
 






