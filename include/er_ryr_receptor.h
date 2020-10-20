//! Ryanodine receptor (RyR)

/*
RyRs reside in the membrane of the ER. 

At nanomolar concentrations, ryanodine locks the receptor in a half-open state. 
At micromolar concentration, it fully closes them.     The effect of the 
nanomolar-level binding is that ryanodine causes release of calcium from 
calcium stores [ Hille, Ionic Channels of Excitable Membranes, 2nd edition].

In studies of central neurons where presynaptic and post-synaptic elements have 
not been studied in tandem, presynaptic CICR has been inferred to contribute 
anywhere from 10 to 80% of the Ca2+ required for neurotransmission or Ca2+ transients.

[Bouchard 2003, Presence and functional significance of presynaptic ryanodine receptors]

In Emptage 2001: 75 ms ISI for paired pulse facilitation (PPF) experiment, inferred ~69%.
*/


class RyR 
{
//  Jcicr: De Schutter and Smolen 1998 (dendrites), RyR mediated CICR in Purkinje Cells
//
//  uptake into stores:
//  Vmax of 4.10^-9 mM^-2  msec^-1 and Kd of 1 mM   Errata?
// 
//  CICR
//  Vmax of   10^-8 mM^-2  msec^-1,  or 3.10^-8 mM^-2 msec^-1, 
//
//  time constant of 1.2 msec
//  and Kd of 0.3 mM.


// Keener and Sneyd 2009, pp302-303, provides a similar model of CICR, 
// which is also based on a sigmoid fx 
//
// k3 = k2*c^n/( Kd^n + c^n)
// J  = k3*(ce - c)
// k3 is rate of release where Kd=500 nM, n=3, k2 = 0.13/s,  
//
// It is necessary only for the ryanodine receptors to be activated by Ca2+ 
// to generate physiological oscillations; inactivation by Ca2+ is not necessary, p303.


// 
// See also Keizer and Levine 1996.
//

public:

double * J_flux;


RyR() { ; }

RyR(int tn) 
{
   J_flux = init_double(tn+1);
}

double Jcicr(int i, double ca, double cer) {
     
    ca = ca/1000;       // convert from nM to uM
    cer = cer/1000;     // Note: [Ca2+]er range is 100 uM to 5 mM

     
    double n=1;         // Hill coefficient.
                        // n=1 makes short stubby Ca spikes, as in De Schutter and Smolen 1998 
                        // n=3 makes sharp tall spikes
                        // low threshold model:    80nM
                        // high threshold model:  200nM,  Vmax=3.8 ...
                        
                        // Parameters from De Schutter 1998, Chapter 6, page 3,
                        // In Methods of Neuronal Modeling by Koch & Segev.
                        //
                        // Note change:  multiply by 3
                        //
						// double x=1e1;   // == 10
                        // double y=10e1;  // == 100
						
    double Vcicr    = 5e-6; // 10^-8 /cm^2 /ms   or 3.8 * 10^-8  /cm^2 /ms     
    double tau_cicr = 1.2;          // 1.2 ms
    
    double Kcicr = 0.3;     // 0.3 uM
    double KT =    0.2;     // 0.2 uM
      
    // inactivation = K / (K + [Ca])     // see Koch 1989, p 105
    
    double J_infinity = 0;
    
    if (ca > KT) {
        J_infinity = Vcicr * ( pow(ca,n)/(pow(ca,n) + pow(Kcicr,n)) ) * (cer - ca);  
    }
    
    J_flux[i + 1] = (J_infinity - J_flux[i])/tau_cicr;
    
    return J_flux[i] * 1000;  // convert back to nM
}
};

// Not used. 
// Based on the model in Gabbiani and Cox, Mathematics for Neuroscientists, p202-203
//
//class RyR_gated
//{
//    // ryanodine receptor flux 
//    //
//    public:
//    static const double nu_r = 10^(-10); // 10^(-6) cm/ms in book
//    
//    static const double K_a = 0.372;     // 0.372 uM in book
//    static const double K_b = 0.636;     // 0.636 uM in book
//    static const double K_c = 0.057;       
//    
//    double * w_r;
//    double w_r_infty;
//    double tau_r_w;
//
//    RyR_gated() 
//    {
//       ;
//    }
//        
//    RyR_gated(int tn) 
//    {
//       w_r = init_double(tn+1);
//    }
//    
//    double Jcicr(int i, double deltaT, double ca, double cer)
//    {
//       ca  =  ca/1000; // convert from nM to uM
//       cer = cer/1000;
//       
//       w_r_infty = (1 + pow(K_a/ca, 4) + pow(ca/K_b, 3))/(1 + (1/K_c) + pow(K_a/ca, 4) + pow(ca/K_b, 3) ); 
//       
//       tau_r_w   =  100;                 // pow(10,4) *  w_r_infty in book 
//       w_r[i+1]  = w_r[i] + (w_r_infty - w_r[i])/tau_r_w;
//   
//       double j_ryr = nu_r * w_r[i] * (1 + pow(ca/K_b,3))/(1 + pow(K_a/ca, 4) + pow(ca/K_b, 3)) * (cer - ca);
//   
//       return -j_ryr*1000;  // convert back to nM 
//    }  




