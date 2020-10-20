
class Vesicle_Markov
{
public:
 
// Synaptic vesicle parameters 
static const double nv=2;  // Number of docked vesicles; nv=2 in Nikonenko & Skibo (2006)
static const double gv=60; // Glutamate concentration in single vesicle;  mM;   Montana et al. (2006)

// NOTE:  rates are in nM per ms
                              // Changes: 10 * kno and 10 * ca have same effect
static const double kon=0.3;  // Calcium binding rate to vesicle; 0.3 /uM /ms in Bollman et al (2000)
                               
static const double koff=3.0;  // Calcium dissociation rate from vesicle;  3 per ms in Bollman et al (2000)
                               // Note:  Kd = koff/kon = 10uM or 10,000 nM
static const double Gamma=30;
static const double delta=8;  // Isomerization constants;  per ms;   Bollman et al (2000)

static const double degG=10;  // Extracellular glutamate degradation rate;  per ms;   Destexhe et al (1998)
 
// Vesicle Regulatory time constants ms,  Tsodyks & Markram (1998)
static const double tau_rec=800;  // Vesicle recovery time constant;  ms;   Tsodyks & Markram (1997)
static const double tau_inact=3;  // Vesicle inactivation time constant;  ms;   Tsodyks & Markram (1997)
 
 
int * x1;  // Temporal evolution of synaptic vesicle 1
int * x2;  // Temporal evolution of synaptic vesicle 2

double * G_syn;  // Synaptic glutamate concentration
double * R_syn;  // Fraction of releasable vesicles
double * E_syn;  // Fraction of effective vesicles in synaptic cleft
double * I_syn;  // Fraction of inactivated vesicles 
double * RRP;    // Vesicles released from RRP

double lastRelease;      // Time of most recent vesicle release (ms)

int spikes;
double lastSpike;

int vesiclesReleased;

double * P_release;
double * P_release_BLOCKER;

double * VR_event;
double * AP_event;

double * Ca_MD;



Vesicle_Markov::Vesicle_Markov()
{
    ;
}


Vesicle_Markov::Vesicle_Markov(int tn)
{    
    // allocate memory
    x1=init_int(tn+2);  
    x2=init_int(tn+2);    
    
    Ca_MD =init_double(tn+2); 
    
    G_syn=init_double(tn+2); 
    R_syn=init_double(tn+2);  
    E_syn=init_double(tn+2);  
    I_syn=init_double(tn+2);  
    RRP=init_double(tn+2);  
    
    P_release        =init_double(tn+2);
    P_release_BLOCKER=init_double(tn+2);

    VR_event=init_double(tn+2); 
    AP_event=init_double(tn+2); 
}

// set values for next trial
void Vesicle_Markov::set(int tn) {
 
    for(int i = 1; i < tn+1; ++i)
    {
      R_syn[i]=1;     // Releasable fraction of vesicles
      E_syn[i]=0;     // Effective fraction of vesicles in synaptic cleft
      I_syn[i]=0;     // Inactivated fraction of vesicles
      G_syn[i]=0;     // Glutamate concentration in cleft;  mM

      Ca_MD[i]=0;
      
      P_release[i]=0;
      P_release_BLOCKER[i]=0;

      VR_event[i]=0; 
      AP_event[i]=0;  
    }
    
    R_syn[1]=1;     // Releasable fraction of vesicles
    E_syn[1]=0;     // Effective fraction of vesicles in synaptic cleft
    I_syn[1]=0;     // Inactivated fraction of vesicles
    G_syn[1]=1e-3;  // Glutamate concentration in cleft;  mM
    
    lastRelease=-10;
    spikes=0;
    lastSpike=-10;
    lastSpike=0;
    vesiclesReleased=0;
    
    // initialise values
    double mu[8]={0, 1, 0, 0, 0, 0, 0};  // Inital vector for vesicle '1' & '2.' 
    
    x1[1]=markov(mu);
    x2[1]=markov(mu);  // Initial state of synaptic vesicle
}


// VGCC, preNMDAR and RyR calcium are included in [Ca] at vesicle's calcium sensor.
//
double Vesicle_Markov::release(int i, EX ex, double Vm, double Vrest, double ca, double AP5)
{
double x_factor=2000;

ca = ca + x_factor;

Ca_MD[i]= ca;   // in nM, and saved as nM for plotting  

ca = ca/1000;   // convert from nM to uM because Markov model assumes uM;
 
// Transition probability matrix 'PM', using 1 as first element, so ignoring PM[0][0]
double a21=5* kon * ca * ex.deltaT;
double a12=koff*ex.deltaT;
double a32=4*kon* ca *ex.deltaT;
double a23=2*koff*ex.deltaT;
double a43=3*kon* ca *ex.deltaT;
double a34=3*koff*ex.deltaT;
double a54=2*kon* ca *ex.deltaT;
double a45=4*koff*ex.deltaT;
double a65=kon*   ca *ex.deltaT;
double a56=5*koff*ex.deltaT;
double a76=Gamma*ex.deltaT;
double a67=delta*ex.deltaT;
// Diagonal elements
double D0=1-a21;
double D1=1-a12-a32;
double D2=1-a23-a43;
double D3=1-a34-a54;
double D4=1-a45-a65;
double D5=1-a56-a76;
double D6=1-a67;


double PM[8][8] = {{0,0,0,0,     0,0,0,0},
                  {0,D0,a12,0,  0,0,0,0},
                  {0,a21,D1,a23,0,0,0,0},
                  {0,0,a32,D2,  a34,0,0,0},
                  {0,0,0,a43,   D3,a45,0,0},
                  {0,0,0,0,     a54,D4,a56,0},
                  {0,0,0,0,     0,a65,D5,a67},
                  {0,0,0,0,     0,0,a76,D6 } };

int j=0;

double mm1[8];
double mm2[8];

for(j=1; j < 8; ++j)   // row j,  col x1[i]
{
   mm1[j] = PM[j][x1[i]];
   mm2[j] = PM[j][x2[i]];
   //printf("%f,  %f\n", mm1[j], mm2[j]);
}

x1[i+1]=markov( mm1 );  // State vector for 1st vesicle
x2[i+1]=markov( mm2 );  // State vector for 2nd vesicle
    

// Refractory period of 6.34 ms. 

if ( ex.t[i] - lastSpike > 6.34  &&   Vm > -40 )  
{ 
    ++spikes;
    lastSpike = ex.t[i];
}



double window = 5; // synchronous vesicle release must be within short time window after spike
int synch = 0;

if ( (ex.t[i] - lastSpike) <= window )
{
   synch = 1;
}



// ASSUMPTIONS:  no vesicle depletion, no spontaneous release.

if ( synch > 0 && x1[i]==7 && x2[i] != 7  &&   ex.t[i] - lastRelease  > 6.34 ) { 
    RRP[i]=0.5;  // one vesicle is released
    lastRelease = ex.t[i];    
    vesiclesReleased +=1;
    VR_event[i]=1;
}
else if ( synch > 0 && x2[i]==7 && x1[i] != 7 &&  ex.t[i] - lastRelease > 6.34 ) { 
    RRP[i]=0.5;  // one vesicle is released
    lastRelease = ex.t[i];
    vesiclesReleased +=1;
    VR_event[i]=1;
}
else if ( synch > 0  && x1[i]==7 && x2[i]==7 &&  ex.t[i] - lastRelease > 6.34 ) { 
    RRP[i]=0.5;    // the simulation assumes 0 or 1 of two docked vesicles
    lastRelease = ex.t[i];
    vesiclesReleased +=1;
    VR_event[i]=1;
}
else {
    RRP[i]=0;
    VR_event[i]=0;
}

 
// Fraction of Neuronal Synaptic vesicles in releasable, effective and inactive states respectively.
//
R_syn[i+1]=1;    // R_syn[i]+ex.deltaT*(((I_syn[i])/tau_rec)-((RRP[i])*R_syn[i]));
E_syn[i+1]=E_syn[i]+ex.deltaT*(((RRP[i])*R_syn[i])-(E_syn[i]/tau_inact));  // RRP[i] = {0,0.5, 1.0}
//I_syn[i+1]=1-R_syn[i+1]-E_syn[i+1];

// Glutamate in the synaptic cleft:
G_syn[i+1]=G_syn[i]+ex.deltaT*(nv*gv*E_syn[i]-degG*(G_syn[i]));

return G_syn[i+1];
}
};