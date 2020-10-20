/*
A dissociation constant is a specific type of equilibrium constant that measures 
the propensity of a larger object to separate (dissociate) reversibly into smaller 
components, as when a complex falls apart into its component molecules, or when 
a salt splits up into its component ions. 

The smaller the dissociation constant, the more tightly bound the ligand is, or 
the higher the affinity between ligand and protein. 

Kd:  The dissociation constant is commonly used to describe the affinity between 
a ligand (L) (such as a drug) and a protein (P), i.e. how tightly a ligand binds 
to a particular protein.

Complex C <-->  P + L

Kd=[P][L]/[C]

The dissociation constant has molar units (M), which correspond to the 
concentration of ligand [L] at which the binding site on a particular protein 
is half occupied,  i.e. the concentration of ligand, at which the concentration 
of protein with ligand bound [C], equals the concentration of protein with no 
ligand bound [P].

High affinity (typically 0.1-10 nM rather than uM range) means that the drug 
binds tightly and less drug is required to bind to 50% of the target protein.

A low affinity Ca2+ sensor can be activated by high [Ca2+], e.g. synaptotagmin. 

Ka = 1/Kd    // the association constant is the inverse of the dissociation 
constant.

Chemical equilibrium is also the ratio of the on-rate (kforward) and 
off-rate (kback) constants.

Ka = on-rate/off-rate
*/


class Vesicle_Hill
{
public:
 
// Synaptic vesicle parameter 

double max_docked;
double num_docked;    
double Vpr;
    
static const double nv=1;  // max number of docked vesicles that can be released per AP;  also see Nikonenko & Skibo (2006)
static const double gv=60; // 60; // Glutamate concentration in single vesicle;  mM; Montana et al. (2006)
static const double degG=10; //10 Extracellular glutamate degradation rate; per ms; p4 Destexhe et al (1998)
 
// Vesicle Regulatory time constants ms,  Tsodyks & Markram (1998)

double tau_rec;  // Vesicle recovery time constant; ms; Tsodyks & Markram (1997)

static const double tau_inact=3;  // Vesicle inactivation time constant; ms; Tsodyks & Markram (1997)

double *     R_syn;     // Releasable fraction of vesicles
double *     E_syn;     // Effective fraction of vesicles in synaptic cleft
double *     I_syn;     // Inactivated fraction of vesicles
double *     G_syn;     // Glutamate concentration in cleft;  mM

double * REL;  // vesicles released from the RRP
double * RRP;

double * P_release;
double * P_release_BLOCKER;

double * VR_event;
double * AP_event;

double * Ca_MD;   // [Ca] in the microdomain of the vesicle

double lastRelease;          // Time of most recent vesicle release (ms)

int vesiclesReleased;

int spikeNumber;

double lastSpike;
int spikes;


Vesicle_Hill(){ ; }

Vesicle_Hill(int tn)
{ 
    R_syn=init_double(tn+2); 
    E_syn=init_double(tn+2); 
    I_syn=init_double(tn+2); 
    G_syn=init_double(tn+2);  
    
    REL=init_double(tn+2);
    RRP=init_double(tn+2);
    
    P_release=init_double(tn+2); 
    P_release_BLOCKER=init_double(tn+2); 
     
    VR_event=init_double(tn+2); 
    AP_event=init_double(tn+2);   
    
    Ca_MD=init_double(tn+2); 
};


void set(int tn) {
 
    for(int i = 1; i < tn+1; ++i)
    {
      R_syn[i]=1;     // Releasable fraction of vesicles
      E_syn[i]=0;     // Effective fraction of vesicles in synaptic cleft
      I_syn[i]=0;     // Inactivated fraction of vesicles
      G_syn[i]=0;     // Glutamate concentration in cleft;  mM
 
      REL[i]=0;
      
      P_release[i]=0;
      P_release_BLOCKER[i]=0;

      VR_event[i]=0; 
      AP_event[i]=0;  
      
      Ca_MD[i]=0;
    }
    
    max_docked=10;
    num_docked=5;
    Vpr=0.10;
    
    tau_rec=800;    // Vesicle recovery time constant; ms; Tsodyks & Markram (1997)
    
    G_syn[1]=1e-3;  // Glutamate concentration in cleft;  mM
    
    lastRelease=-100;
    vesiclesReleased=0;
    
    spikeNumber=0;
    
    spikes=0;
    lastSpike=-100;
};


    
double release(int i, EX ex, double Vm, double vr, double Ca, double AP5)
{

double x_factor=000;

Ca_MD[i+1] = Ca + x_factor; // nM

Ca = (Ca + x_factor)/1000.0; // convert from nM to uM, Ca may include Ca2+ from vgcc, preNMDARs, RyRs


// IP3 and RyR:    
// Within the physiological range of cytoplasmic calcium, the IP3-gated channel itself allows 
// positive feedback and then negative feedback for calcium release, whereas the ryanodine 
// receptor/channel behaves solely as a calcium-activated channel.   The existence in the 
// same cell of two channels with different responses to calcium and different ligand 
// sensitivities provides a basis for complex patterns of intracellular calcium regulation. 
// [Bezprozvanny 1991]. 

// Kd=10 uM in Calyx of Held.
// Schneggenburger, R. & Neher, E. Intracellular calcium dependence of transmitter 
// release rates at a fast central synapse. Nature(Lond.) 406, 889-893 (2000).
//
// Pr of vesicle fusion given the [Ca] and Hill function.

double n=ex.n1,  Kd=ex.Kd1;

RRP[i] = 1;  // floor(num_docked);            

double Pr_max = 0.90;  //  1 - pow((1-Vpr), RRP[i]);

// Synaptotagmin 1: low affinity calcium sensor triggers vesicle fusion and release
double Pr1 =  pow(Ca,n) / (pow(Kd,n) + pow(Ca,n));

// Second sensor: a very low affinity calcium sensor (opposes release)
n = ex.n2;   Kd = ex.Kd2; 

double Pr2 = pow(Kd,n) / (pow(Kd,n) + pow(Ca,n));

double Pr = Pr_max * Pr1; // * Pr2;

//printf("Pr=%0.1f   Pr1=%0.1f    Pr2=%0.1f \n", Pr, Pr1, Pr2);

double fusion_rate  = Pr;  //  

if (AP5 == 0 ) {
    P_release[i]=fusion_rate;
}
else {
    P_release_BLOCKER[i]=fusion_rate;
}


// Refractory period of 6.34 ms. 
if ( ex.t[i] - lastSpike > 6.34  &&   Vm > -40 )  
{ 
    ++spikes;
    lastSpike = ex.t[i];   // vesicle release window starts at beginning of most recent spike
}


double window=5;      // p678 Meinrenken 2003:  [Ca]avg_vesicle peaks at 8 uM and decays to 400 nM, 
                      // predicted avg Pr (after 5ms) is 25%
int synch = 0;

if ( (ex.t[i] - lastSpike) <= window )
{
   synch = 1;
}

Pr = (Pr * ex.deltaT);  // Pr per time step

// ASSUMPTIONS:  no spontaneous release.   vesicle depletion versus no vesicle depletion?

if ( synch == 1  && (ex.t[i] - lastRelease) > 6.34 && RRP[i] >= 1 && rnd() < Pr )  
{
   lastRelease = ex.t[i];    
   vesiclesReleased +=1;  
   
   VR_event[i]=1;
   REL[i]=1.0;        // vesicles released from REL: 0 or 1
   num_docked = num_docked - 1;
} 
else
{
   VR_event[i]=0; 
   REL[i]=0;
}


num_docked = num_docked + ex.deltaT * ( (max_docked - num_docked)/tau_rec ) * (Ca - 0.100)/(Ca + 0.100);

// Glutamate in the synaptic cleft
//
// Glutamate peaked at 1.1 millimolar and decayed with a time constant of 
// 1.2 milliseconds at cultured hippocampal synapses (Clements 1992).
//
// In Tewari 2012a, p10, synaptic [Glu] peaked at ~= 0.275 mM
// while extrasynaptic [Glu]= peaked at about 1.75 mM.
// Time course of glutamate in the cleft is 2 ms.
//
//
// Fraction of Neuronal Synaptic vesicles in releasable, effective
// and inactive states respectively
R_syn[i+1]= 1;  // R_syn[i]+ex.deltaT*(((I_syn[i])/tau_rec)-((REL[i])*R_syn[i])); // 1
E_syn[i+1]= E_syn[i]+ex.deltaT*( ( REL[i]*R_syn[i] )-(E_syn[i]/tau_inact) );  // REL[i] = {0,1.0}
I_syn[i+1]= 1-R_syn[i+1]-E_syn[i+1];

// Glutamate in the synaptic cleft
G_syn[i+1]= G_syn[i]+ex.deltaT*(nv*gv*E_syn[i]-degG*(G_syn[i]));


return G_syn[i+1];

};

};