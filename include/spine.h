

class Spine {

public:

//  Spine Parameters
static const double R_in=0.7985e+8;  // Input resistance of dendrite spine; k-ohm; Calculated

static const double tau_mem=50;      // Post-synaptic membr time constant; ms; Tsodyks & Markram (1997)
static const double rad_spine=0.6e-4;  // Radius of spine-head; cm; Dumitriu et al (2010)

double vspine;  // Volume of spine-head; unit: L; Calculated
double sspine;  // Surface area of spine-head; unit: cm^2; Calculated

static const double Vrest=-70;  // Resting membrane potential of spine membrane; unit: mV

                                     
// Post-synaptic spine variables
double * Vm;    // Membrane potential
double * Iampa;  // Current through AMPAR
double * Inmda;
double * Ip2x;

AMPA_spine ampaR;
NMDA_spine nmdaR;
P2X_spine  p2xR;


// Volume of spine-head (sphere); L; Calculated 
// Surface area of spine-head (sphere); cm^2; Calculated
 public:
 
     Spine() {
        ;
     }
 
     Spine(int tn) {
      
         vspine=1e-3*(4.0/3.0)*PI*pow(rad_spine,3);  
         sspine=4*PI*pow(rad_spine,2);      
         
         Vm    =init_double(tn+2);  // Membrane potential
         Iampa =init_double(tn+2);  // Current through AMPAR         
         Inmda =init_double(tn+2);
         Ip2x  =init_double(tn+2);  
         
         ampaR=AMPA_spine();
         nmdaR=NMDA_spine();
         p2xR=P2X_spine(1);
    }
    
    void set(int tn) {
       for(int i = 2; i < tn+2; ++i)
       {
         Vm[i]=0;  // Membrane potential
         Iampa[i]=0;  // Current through AMPAR         
         Inmda[i]=0;
         Ip2x[i]=0;  
       }
       
       // initial conditions
       Vm[1]=-70;  // Resting post-membrane potential; unit: mV

       Iampa[1]=0;  // Current through AMPAR; unit: uA        
       Inmda[1]=0;
       Ip2x[1]=0;
    }
    
                  
    // Euler's method      i,   time t[i],     deltaT,      B.tdr,    B.G_syn[i]);
    void spine_model_1(int i, double t, double deltaT, double tdr, double glu) 
    {
         // Vm: membrane voltage of spine
         Iampa[i]=ampaR.syn(deltaT, glu, Vm[i]);             
         Inmda[i]=nmdaR.syn(deltaT, glu, Vm[i]);
         Ip2x[i] =p2xR.syn(t, tdr); 
               
         Vm[i+1]=Vm[i]+deltaT*((-(Vm[i]- Vrest)-R_in*(Iampa[i]))/tau_mem);  // Membrane potential
    } 
    
    // Heun's method
    void spine_model_2(int i, double t, double deltaT, double tdr, double glu) 
    {
        // Vm: membrane voltage of spine
        Iampa[i]=ampaR.syn(deltaT, glu, Vm[i]);                 
        Inmda[i]=nmdaR.syn(deltaT, glu, Vm[i]);
        Ip2x[i] =p2xR.syn(t, tdr); 
        
        // temp values
        double k1 = deltaT*((-(Vm[i]- Vrest)-R_in*(Iampa[i]))/tau_mem);  // Membrane potential
        
        double Iampa_2 = ampaR.syn(deltaT, glu,  Vm[i]+k1);                 
        //double Inmda_2 = nmdaR.syn(deltaT, glu,  Vm[i]+k1);
        //double Ip2x_2  = p2xR.syn(t, tdr); 
        
        double k2 = deltaT*((-(Vm[i]+k1  - Vrest)-R_in*(Iampa_2))/tau_mem);  // Membrane potential
        
        Vm[i+1] = Vm[i]+ (k1 + k2)/2;
    } 
    
};


