
class Vesicle_Markov_6
{
public:
 
// Synaptic vesicle parameters 
static const double nv=2;  // Number of docked vesicles; nv=2 in Nikonenko & Skibo (2006)
static const double gv=60; // Glutamate concentration in single vesicle;  mM;   Montana et al. (2006)

// NOTE:  rates are in nM per ms
                              // Changes: 10 * kon and 10 * ca have same effect
static const double kon=0.3;  // Calcium binding rate to vesicle; 0.3 /uM /ms in Bollman et al (2000)
                               
static const double koff=3.0;  // Calcium dissociation rate from vesicle;  3 per ms in Bollman et al (2000)
                               // Note:  Kd = koff/kon = 10uM or 10,000 nM

static const double degG=10;  // Extracellular glutamate degradation rate;  per ms;   Destexhe et al (1998)
 
// Vesicle Regulatory time constants ms,  Tsodyks & Markram (1998)
static const double tau_rec=800;  // Vesicle recovery time constant;  ms;   Tsodyks & Markram (1997)
static const double tau_inact=3;  // Vesicle inactivation time constant;  ms;   Tsodyks & Markram (1997)
 

double Xn;

double * G_syn;  // Synaptic glutamate conceVesicle_Markov_6ntration
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


Vesicle_Markov_6()
{
    ;
};


Vesicle_Markov_6(int tn)
{    
    Xn=0;
    
    Ca_MD =init_double(tn+2); 
    
    G_syn=init_double(tn+2); 
    R_syn=init_double(tn+2);  
    E_syn=init_double(tn+2);  
    I_syn=init_double(tn+2);  
    RRP  =init_double(tn+2);  
    
    P_release        =init_double(tn+2);
    P_release_BLOCKER=init_double(tn+2);

    VR_event=init_double(tn+2); 
    AP_event=init_double(tn+2); 
};

// set values for next trial
void set(int tn) {
 
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
    
    Xn=0;
};


// VGCC, preNMDAR and RyR calcium are included in [Ca] at vesicle's calcium sensor.
//
//
double release(int i, EX ex, double Vm, double Vrest, double ca, double AP5)
{

double x_factor=0;

ca = ca + x_factor;

Ca_MD[i]= ca;   // in nM, and saved as nM for plotting  

ca = ca/1000;   // convert to uM

double right=0; double left=0;      // Transition probabilities
 
 
if (Xn == 0) {
    right = 5 * kon * ca * ex.deltaT;
    right = pow(right,1)/(0.7 + pow(right,1));
}   
else if (Xn > 0 && Xn < 5) {

    right = (5 - Xn) * kon * ca * ex.deltaT;
    right = pow(right,1)/(0.7 + pow(right,1));
    
    left  = Xn  * koff* ex.deltaT;
}
else {
    left = 5*koff* ex.deltaT;
}




double window = 5; // synchronous vesicle release must be within short time window after spike

double rn=rnd();

if (rn < right) {
   Xn = Xn + 1;
}
else if (rn < (right + left) )  {
   Xn = Xn - 1;
}


if (Xn < 0 || Xn > 5) {
   printf(" ==> ERROR Xn=%0.2f\n ", Xn);
   exit(1);
}

// Refractory period of 6.34 ms. 

if ( ex.t[i] - lastSpike > 6.34  &&   Vm > -40 )  
{ 
    ++spikes;
    lastSpike = ex.t[i];
}


int synch = 0;

if ( (ex.t[i] - lastSpike) <= window )
{
   synch = 1;
}


if (AP5 == 0 ) {
    P_release[i] = right;
}
else {
    P_release_BLOCKER[i] = left;
}


// ASSUMPTIONS:  no vesicle depletion, no spontaneous release.

if ( synch > 0 && Xn == 5  && ex.t[i] - lastRelease  > 6.34 ) { 
    RRP[i]=0.5;  // one vesicle is released
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
R_syn[i+1]=1;    //R_syn[i]+ex.deltaT*(((I_syn[i])/tau_rec)-((RRP[i])*R_syn[i]));
E_syn[i+1]=E_syn[i]+ex.deltaT*(((RRP[i])*R_syn[i])-(E_syn[i]/tau_inact));  // RRP[i] = {0,0.5, 1.0}
//I_syn[i+1]=1-R_syn[i+1]-E_syn[i+1];

// Glutamate in the synaptic cleft:
G_syn[i+1]=G_syn[i]+ex.deltaT*(nv*gv*E_syn[i]-degG*(G_syn[i]));

return G_syn[i+1];
};


};