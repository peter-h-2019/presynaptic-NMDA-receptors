/*  
NMDARs
Ermentrout and Terman, Mathematical Foundations of Neuroscience, 2010, p161

NMDARs produce EPSCs that are slower and last longer that those of AMPARs
In general, Mg++ ions block the channel until the membrane is depolarized 
by open AMPARs.   Both glutamate and glycine or D-serine must bind to 
NMDARs to activate them.

Ermentrout and Terman 2010, p163, Eq (8.8):
The NMDA current is modeled as: I_NMDA = gNMDA * s * B(V)*(V -’ VNMDA) 
where s obeys equations (8.5,8.6) and B(V) represents the magnesium block
(Jahr and Stevens 1990, J. Neuroscience 10, 1830-183  ).

At [Mg++]=2 mM,  V_T= -10 mV,  at 1mM,  -20mV
double V_T = 16.13 * log(Mg/3.37);
*/  


class PreNMDAR
{    
    /* 
    The modulation of transmitter release by presynaptic receptors is 
    an accepted signaling pathway.   The focus of attention was initially 
    on metabotropic G-protein coupled receptors [52,94], but it became 
    clear that numerous populations of presynaptic ionotropic receptors are 
    equally important [53,58,67]. 
    Van Dongen, Biology of the NMDA Receptor, 2009,
    see Chapter 14 Presynaptic NMDA Receptors by Ian C. Duguid, Trevor G. Smart. 
    http://www.ncbi.nlm.nih.gov/books/NBK5275/
    */
    public:
    //  units are mM  and  ms,  Tmax ~ 1mM
    //  Gerstner, p53, max gNMDA = 1.2 nS (granule cell) per synapse.
    static const double a_r=0.072;     // p163, Ermentrout 0.072 /mM/ms, in mM but here using uM !
                                       // a:alpha, larger a_r, a_d means faster rise, decay
    static const double a_d=0.0066;    // p163, Ermentrout 0.0066 /ms
                                       // Destexhe 1998: max gNMDA = 0.01 .. 0.6 nano Siemens per synapse
    static const double Vnmda=0;   
    static const double gNMDA=0.3;    //  .21    0.7 * 10; 
                                      //  max conductance=0.7 for excitatory neuron, p382 of Traub 1994, 
    static const double Mg=2.0;       //  extracellular [Mg++] in mM
    //double Vca;                     //  E_Ca ~ 120 mV
                                      //  See struct Ex in utilities.h 
    double Vca=130.65;                //  130.65 if [Ca]ex = 3 mM as in McGuinness 2010, 
                                      //  used Nernst Eq., 125 if 2 mM 
    double * syn;                    
    
PreNMDAR() 
{
   ;   
}

PreNMDAR(int tn, double v_ca) 
{
   Vca=v_ca;

   syn = init_double(tn);
}

double I_Ca(int i, double deltaT, double glu, double Vm, double ca)
{
   // Mg2+ blocks channel unless membrane is depolarised
   double B =  1/( 1 + exp(-0.062 * Vm) * (Mg/3.57) );     // 3.57 mM     // p163 Ermentrout,2010  
   
   // printf("BR 61 \n");
   // larger a_r and a_d cause fast changing conductance, 
   // syn[i] = fraction of open channels at time t=i 
   syn[i+1] = syn[i] + deltaT * (  a_r * glu * (1 - syn[i]) - a_d * syn[i]);
   
   /* 
     True [Ca] in microdomain near preNMDARs would be much higher than 
     global [Ca], perhaps 10 to 100 times higher.
  
     Below: Hill equation based Ca2+ sensor for Ca2+ dependent inactivation 
     of preNMDARs.  Without inactivation, [Ca] would get too high.
   */
   double n = 4, Kd = 10000;     // nM,  or 10 uM
   double sensor = pow(Kd,n)/(pow(Kd,n) + pow(ca,n));     
                                                                 
   double I_Ca = gNMDA * sensor * syn[i] * B * (Vm - Vca); // Vca=125  if [Ca]ex=2 mM, ~130 if 3 mM
      
   return I_Ca;
};
};


/**
"Presynaptic N-type and P/Q-type Ca2+ channels mediating
synaptic transmission at the calyx of Held of mice", Ishikawa (2005)

*/
// VGCC: inactivated by a Hill function at very high [Ca].

class VGCC_bouton {
public:
    // Calcium current parameters for P-type Ca2+ channels.
    //
    static const double v_half=-17;       // Half-activation voltage for wild type mice; mV; Ishikawa (2005)
    static const double slope_factor=8.4; // Slope factor for wild type mice; mV; Ishikawa (2005)
    
    static const double tau_mc=10;       // Tewari 2012b, the LTP paper;   Ishikawa et al. 2005 

    
    // The fast and slow time constants of inactivation
    // were 28 and 609;  ms; p202 Ishikawa (2005).  Note: 77 and 775 ms for knock out mice.
    
    // Rev potential for Calcium ion determined through Nernst equation. 
    // Assumes extracellular [Ca2+]=2 mM as in Perea & Araque (2007)
    double Vca;    
    
    // Both N-type and P-type channels are present in the nerve terminals at the 
    // end of axons (as well as in cell bodies and in dendrites) and are known 
    // to provide the calcium influx that triggers release of synaptic vesicles 
    // from many nerve terminals (Regehr and Mintz, 1994)
    //
    // p=presynaptic
    // IpCa calcium current: 90% from P/Q-type channels, calyx of Held from immature mouse, Ishikawa (2005)
    //
    //                                  See Nakamura (2015) p152
    static const double g_ca=3.3;    // 3.3 pS;  conductance of single Ca_v_2.1 at calyx,  Sheng (2012) 
                                     // 
                                     // Ca_v_2.2:  2.7 pS; Weber (2010) 
                                     
                                     // if this is changed to 2.3, then must make changes to nmdaR too.
                                     // 2.3e-9 is related to the units used for surface area in bouton.h
                                     
                                     // single channel conductance: 2.2 pS;  Vyleta (2014)
                                     // Vyleta assumed 23 open P/Q-type Ca2+ at the peak 
                                     // (consistent with experimental data from MF bouton).
                                     // See page 7 of Vyleta (2014) supplement.
                                     
                                     // The conductance of the CaV2.2 channel.
                                     // At 2 mM [Ca2+]ex, conductance is 2.76 ± 0.24 pS.
                                     
                                     // single channel current amplitude of 0.35 pA, based on
                                     // 3.3 pS in 2 mM [Ca], measured for Ca_v_2.1 at calyx, 
                                     // p152, Sheng (2012).
                                                     
                                     
    static const double rho_ca=1;  // number of channels in cluster was 32 and 0.03 for surface area in bouton.h
                                   // literature says Density of N-type channels; 3.2 per um^2, or 32e7 cm^2
                                         
                                   // 400 /um^2 ==> 400e8, or 40e7

    double   gc;     // Calcium channel conductance density;  mS / cm^2
    double * mc;     // VGCC gating variable
    
   
    // see book by Liu2012,  p141:  uses HH formalism by Chay and Keizer
   
   
VGCC_bouton()
{    
     ;   
};


VGCC_bouton(int tn, double v_ca) 
{ 
    Vca=v_ca;
    gc= g_ca * rho_ca;    // max single channel conductance * channel density
    
    mc=init_double(tn); // VGCC gating variable for calcium channel
    mc[1]=0;             
};


double I_Ca(int i, double deltaT, double v, double ca_VGCC) 
{
   double mcinf=1/(1+exp((v_half-v)/slope_factor));   // current activation 
    
   mc[i+1]=mc[i]+deltaT*((mcinf - mc[i])/tau_mc);     // VGCC gating variable tau_mc ????

   double n = 2, Kd = 2000;
   
   double sensor = pow(Kd,n)/(pow(Kd,n) + pow(ca_VGCC,n));  // Ca2+ dependent inactivation
   
   // equations due to Erler 2004, except sensor has been added. 
   double I_Ca = gc * sensor * pow(mc[i],2) * (v-Vca);  // VGCC current;  uA per cm**2
   
   return I_Ca;
}

};




