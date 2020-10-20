//! Allosteric model of the calcium sensor.
/*! 
  This allosteric model of a calcium sensor is based on the paper by Lou et al.,
  "Allosteric modulation of the presynaptic Ca2+ sensor for vesicle fusion", 
  Nature, 2005.
  
  In this model, Ca2+ may bind to 5 sites (6 states).  In their paper, release 
  was assumed to occur from a single homogeneous pool of 2000 vesicles. 
  Each vesicle was defined as being in one of six states (0 to 5).  
  Their model predicted a vesicle fusion rate based on the number of vesicles in 
  each state.  
  
  In this implementation, the model predicts the probability of vesicle release 
  by calculating the sensor's "activation level", i.e. current fusion rate as a 
  fraction of the maximum fusion rate per ms according to the number of vesicles
  in each state.   The fusion rate for these vesicles increases with the number 
  of Ca2+ ions bound per vesicle.
*/

class Vesicle_Allosteric
{
public:
 
//  Ca2+ binding rates:  convert fom per s to per ms, e.g. 4000 --> 4
static const double kon=1e5;   // 1e5 /M /ms;  1e8 /M /s 
static const double koff=4;      // 4      /ms;      4000 /s

static const double b=0.25; // cooperativity factor;  b=0.5 in Lou (2005)

static const double f=31.3;  // multiplier constant to I, representing the 
                             // greater rate of fusion from faster states

static const double I=2e-7; // the instantaneous rate of fusion per ms, 
                            // representing the constant exocytosis of
                            // vesicles, or basal fusion "willingness", 
                            // see Lou 2005, p500,  p501:  2e-4 /s
        
 
// In the calyx of Held, a giant synapse in the auditory system, 
// there are about 600 active zones and thousands of
// vesicles.   Each action potential evokes the release of many vesicles.   
// Most CA3-CA1 synapses have only one
// active zone and release at most one vesicle per action potential. 
// 
// Number of vesicles in each state.
double V0Ca;   
double V1Ca; 
double V2Ca; 
double V3Ca; 
double V4Ca;  
double V5Ca;  
 
double primed; // put vesicle state counts in the steady state for initial [Ca]i 

// Synaptic vesicle parameter 
static const double nv=2;    // Number of vesicles; nv=2 in Nikonenko & Skibo (2006)
static const double gv=60;   // Glutamate concentration in single vesicle;  mM;   Montana et al. (2006)
static const double degG=10; // Extracellular glutamate degradation rate; per ms; Destexhe et al (1998)

// Vesicle Regulatory time constants ms,  Tsodyks & Markram (1998)
static const double tau_rec=800;  // Vesicle recovery time constant;  ms;   Tsodyks & Markram (1997)
static const double tau_inact=3;  // Vesicle inactivation time constant;  ms;   Tsodyks & Markram (1997)

double * G_syn;  // Synaptic glutamate concentration
double * R_syn;  // Fraction of releasable vesicles
double * E_syn;  // Fraction of effective vesicles in synaptic cleft
double * I_syn;  // Fraction of inactivated vesicles 
double * RRP;    // Vesicles released 

double lastRelease;         // Time of most recent vesicle release (ms) 

int spikes;
double lastSpike;
int vesiclesReleased;


double * P_release;         // Pr per ms, control condition
double * P_release_BLOCKER; // Pr per ms, when blocker applied

double * VR_event;
double * AP_event;

double * Ca_MD;


Vesicle_Allosteric::Vesicle_Allosteric() { ; }


Vesicle_Allosteric::Vesicle_Allosteric(int tn)
{    
    // allocate memory 
    
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

// set values at the start of each trial (N trials per experimental condition
void Vesicle_Allosteric::set(int tn) {
 
    for(int i = 1; i < tn+1; ++i)
    {
      G_syn[i]=0; 
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
    
    lastRelease=-100;      // Time of last vesicle release, initialised to 100 ms before start of simulation.
    lastSpike=-100;
    spikes=0;
    
    vesiclesReleased=0;
    
    V0Ca=10;   
    
    V1Ca=2; 
    V2Ca=1; 
    V3Ca=1; 
    V4Ca=1;  
    V5Ca=1;
    
    primed=0;
}



double Vesicle_Allosteric::release(int i, EX ex, double Vm, double Vrest, double Ca, double AP5)
{
 
double x_factor=2000; 
 
Ca_MD[i]= Ca + x_factor;  // nM

// Lou 2005: Figure 4(a) shows a plot where the fusion rate is very [Ca2+] 
// dependent from 2-8 uM  
//
// Ca is in nM, so convert to Moles because the allosteric model in Lou 2005 assumes Moles.
Ca = (Ca + x_factor)* pow(10,-9);      
 

double fr=0;

// ASSUMPTIONS:  no vesicle depletion, no spontaneous release.
    
do
{

double dV0Ca =                      koff*V1Ca             - 5*kon *Ca*V0Ca;

double dV1Ca = 5* kon *Ca*V0Ca -    koff * V1Ca           - 4*kon*Ca*V1Ca + 2* koff * pow(b,1) *V2Ca;

double dV2Ca = 4* kon *Ca*V1Ca - 2* koff * pow(b,1) *V2Ca - 3*kon*Ca*V2Ca + 3* koff * pow(b,2) * V3Ca;

double dV3Ca = 3* kon *Ca*V2Ca - 3* koff * pow(b,2) *V3Ca - 2*kon*Ca*V3Ca + 4* koff * pow(b,3) * V4Ca;

double dV4Ca = 2* kon *Ca*V3Ca - 4* koff * pow(b,3) *V4Ca -   kon*Ca*V4Ca + 5* koff * pow(b,4) * V5Ca;

double dV5Ca =    kon *Ca*V4Ca - 5* koff * pow(b,4) *V5Ca;


V0Ca = V0Ca+  ex.deltaT *dV0Ca;
V1Ca = V1Ca+  ex.deltaT *dV1Ca;    
V2Ca = V2Ca+  ex.deltaT *dV2Ca;
V3Ca = V3Ca+  ex.deltaT *dV3Ca;  
V4Ca = V4Ca+  ex.deltaT *dV4Ca;
V5Ca = V5Ca+  ex.deltaT *dV5Ca;
    
if (V0Ca < 0)  { V0Ca=0; }
if (V1Ca < 0)  { V1Ca=0; }
if (V2Ca < 0)  { V2Ca=0; }
if (V3Ca < 0)  { V3Ca=0; } 
if (V4Ca < 0)  { V4Ca=0; }  
if (V5Ca < 0)  { V5Ca=0; }


// fr is the fusion rate per ms.  Max fr is 6 times the total number of vesicles.
fr = I*V0Ca + I*f*V1Ca  +I*pow(f,2)*V2Ca +I*pow(f,3)*V3Ca + I*pow(f,4)*V4Ca + I*pow(f,5)*V5Ca;

primed+=1;   // prime vesicles for 15 ms to put them in steady state for initial [Ca]i

} while ( primed < 300 ); 
    
    
// convert fr to percent of max
double max_rate = 6*( V0Ca + V1Ca + V2Ca + V3Ca + V4Ca + V5Ca );
double activation = fr / max_rate;
double pr_max = 0.90;


double pr = pr_max * activation;  


if(AP5 == 0) {
    P_release[i] = pr;       // Pr per ms
}
else {
    P_release_BLOCKER[i] = pr; // Pr per ms
}

// Refractory period of 6.34 ms. 

if ( ex.t[i] - lastSpike > 6.34  &&   Vm > -40 )  
{ 
    ++spikes;
    lastSpike = ex.t[i];
}


double window = 2; // synchronous vesicle release must be with time window after spike
int synch = 0;

if ( (ex.t[i] - lastSpike) <= window )
{
   synch = 1;
}

pr = (pr * ex.deltaT);     // per ms

double rv = rnd();


if ( synch > 0 && rv < pr  &&  ex.t[i] - lastRelease  > 6.34 ) {
    RRP[i]=0.5;       // one vesicle is released
    lastRelease = ex.t[i];    
    vesiclesReleased +=1;
    VR_event[i]=1;
}
else {
    RRP[i]=0;         // zero vesicles are released
    VR_event[i]=0;
}

// Fraction of Neuronal Synaptic vesicles in releasable, effective
// and inactive states respectively
R_syn[i+1]=1; // R_syn[i]+ex.deltaT*(((I_syn[i])/tau_rec)-((RRP[i])*R_syn[i]));
E_syn[i+1]=E_syn[i]+ex.deltaT*(((RRP[i])*R_syn[i])-(E_syn[i]/tau_inact));  // RRP[i] = {0,0.5, 1.0}
//I_syn[i+1]=1-R_syn[i+1]-E_syn[i+1];

// Glutamate in the synaptic cleft
G_syn[i+1]=G_syn[i]+ex.deltaT*(nv*gv*E_syn[i]-degG*(G_syn[i]));

return G_syn[i+1];
}
};