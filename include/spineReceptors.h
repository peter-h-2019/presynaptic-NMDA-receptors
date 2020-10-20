/*   

This module contains a collection of receptors and indirect currents.

Glutamate binds to AMPA, NMDA, GABA A, and GABA B receptors.

Calcium elevation in the astrocyte can be mapped to a repolarizing neuronal 
current: Iastro (Nadkarni and Jung 2004,2007).


Gluatamate:
The neurotransmitter glutamate activates two different kinds of receptors: 
AMPA/kainate which are very fast and NMDA which is implicated in memory and 
long-term potentiation of synapses. Both of these receptors lead to 
excitation of the membrane.

AMPA and NMDA Rs are often co-located:
Weak stimulus (Glu) normally only activates the AMPA Rs, entry of Na+, which 
results in a partial depolarization of the Post Syn Neuron to -35mV.

If stimulus is large enough, depolarization will be large enough to expel 
Mg++ from the NMDA channel.    Then NMDA channels admit not only Na+, 
but Ca2+ as well.

Ca2+ acts as a 2nd messenger, activating several signalling cascades.
This leads to phosphorylation of existing AMPA Rs and the installation of 
    more AMPA Rs.   Result: synaptic enhancement.

AMPA_R: 
mM  ms;        ms:   in typical cortical cells, the rise time is 0.4 to 0.8 ms                                           
a_r=1.1;  a_d=0.19  with [Glu] of 1mM, rise time 1/(1.1+0.19)=0.8ms,decay=5ms

But Zito 2009 says:
  rise, decay times:
    AMPA: 0.2 to 0.4 ms;  ~2 ms
    NMDA: 10 to 50 ms;    ~50 to 500 ms 

AMPA receptors on inhibitory neurons have rise and fall times ~2x faster
than those on excitatory neurons.

AMPA receptors show strong depression, amplitude of Iampa decreases with 
each subsequent spike (ST plasticity).

NMDA R:  
Tglu = Tmax/(1 + exp(-(Vp - V_T)/Kp)

where Tmax is the maximal concentration of transmitter in the synaptic cleft 
and Kp,Vp determine the stiffness and threshold for the release. A good value 
for these constants is Vp=2 mV and Kp=5 mV.  Typically, assume that 1 mM 
of transmitter is the maximum concentration released. 


AMPA, NMDA, and kainate Rs are named after their agonists (AMPA, NMDA, and 
kainic acid).   They always produce excitatory postsynaptic responses 
(Purves, p118).


ATP:
ATP binds to ionotropic P2X receptors and metabotropic P2Y receptors
*/



/*
Ermentrout and Terman, Mathematical Foundations of Neuroscience, 2010, p161 ...

AMPA Rs produce larger EPSCs than other types of ionotropic glu receptors. 
*/
class AMPA_spine
{
 public:
    //   units are mM  and  ms        
    
    // alpha rise and decay, a_r and a_d,  are voltage independent.
    // Ermentrout and Terman 2010.
    static const double a_r=1.1;   // per mM per ms, alpha rise rate of syn conductance
                                   // 1.1e6 per M per second
    static const double a_d=0.67;  // per ms, alpha decay rate of syn conductance
                                   // original 0.67;
    
    static const double gAMPA = 0.35e-6;// max conductance of AMPAR at single synapse, mS
                                        // range is 0.35 to 1.0 nS, p10 Destexhe
    static const double Eampa=0;                                   
    /*   
    In typical cortical cells, the rise time is 0.4 to 0.8 milliseconds.
    Ex.  Tmax/(1.1 + 0.19) ==> 0.8 rise time, 5 ms decay time.
    
    AMPA receptors onto inhibitory interneurons are about twice as fast 
    in rise and fall times as those onto excitatory neurons.
    Real AMPA synapses show quite strong depression. That is, the 
    peak amplitude of the AMPA current decreases with each subsequent spike. 
    */
 

    double s;   // fraction of receptors in open state (gating variable)
    
    AMPA_spine(){
     
        s=0;  
    }
    
    double syn(double deltaT, double glu, double Vm)
    {   
        // gAMPA: maximal conductance
      
        double Iampa=gAMPA * s * (Vm - Eampa);  // AMPAR current (uA)
        s = s + deltaT * (a_r * glu *(1 - s) - a_d*s);  
        
        return Iampa;
    }
};    

/*   Derived from AMPA_R of Ermentrout and Terman, Mathematical Foundations of 
     Neuroscience, 2010, p161
      
     Kainate receptors: in some cases they are found on pre-synaptic terminals
     and serve as feedback mechanisms to regulate glu release.   When found on 
     postsynaptic cells, they generate smaller EPSCs that **rise quickly but decay 
     more slowly** than those mediated by AMPA Rs, (Purves 5th Edition, 2012, p118).
*/ 
    
class Kainate_spine {

    //   units are mM  and  ms        
    static const double Tmax=1,      Kp=5;
    static const double a_r=1.1,     a_d=0.01;   
    static const double Vkainate=0,  gKainate=0.9e-6;  // mS
        
    double s;

    Kainate_spine() {
     s=0;
    }

    double syn(double deltaT, double glu, double Vm) 
    {   
        double Ikainate = gKainate * s * (Vm - Vkainate);
        s = s+ deltaT * (a_r * glu * (1 - s) - a_d * s);
  
        return Ikainate;
    }
};



class NMDA_spine
{    
    public:
    //  units are mM  and  ms,  Tmax ~ 1mM
    //  Gerstner, p53, max gNMDA = 1.2 nS (granule cell)
    static const double a_r=0.072;   // rise and decays times per mM per ms
    static const double a_d=0.0066;      
    static const double Enmda=0;       // Destexhe 1998: max gNMDA = 0.01 .. 0.6 nS per syn
    static const double gNMDA=0.7e-6;  // conductance mS;   0.7 nS from Traub 1994, p382, excitatory neuron
    static const double Mg=1.0;        // conc of extracellular Mg in mM
    
    double s;
    
    NMDA_spine( ) 
    {
       s=0;
    }
    
    // VT half activatio:  VT= 16.13 ln([Mg]/3.57)
    // if [Mg]=2mM, VT = -10mV
    // if [Mg]=1mM, VT = -20mV
    // B = 1/(1 + exp(-(Vm-VT)/16.13) )

    double syn(double deltaT, double glu, double Vm)
    {   
        double B =  1/( 1 + exp(-0.062 * Vm) * Mg/3.57 );  // p163 Ermentrout, 2010
 
        s = s + deltaT * (a_r * glu * (1 - s) - a_d * s);
        
        double Inmda = gNMDA * s * B * (Vm - Enmda);     //    add Inmda to Vd_dt
        
        return Inmda;
    }
    
    
    // assuming interneuron diameter of 15 microns, 6.0 mu A/cm**2
    
    double astroPost(double Ca_astro)
    {
        double Iastro=0;
        double z=(Ca_astro*1000 - 196.69);   //   convert Ca_astro from uM to nM
        if(z > 0) {
            Iastro = (6.0 * heaviside(log(z)) * log(z)); 
        }
        
        return Iastro;                //   time in ms
    }
};


//   GABAergic synapses employ three types of postsynaptic receptors:
//   GABA A,B,and C.  A and C are ionotropic, B receptors are metabotropic.
//  
//   Ermentrout and Terman, Mathematical Foundations of Neuroscience, 2010, p162
//  


class GABA_A_spine {
 
    public:
    //   units are mM  and  ms    
    static const double a_r = 5;  
    static const double a_d = 0.18;  
    static const double V_GABA=-75;   //   reversal potential for Cl- 
    static const double g_GABA=0.5e-6;   // mS,  Destexhe 1998: max g_gaba_a = 0.25 .. 1.2 nS per syn
    
    double s;
    
    GABA_A_spine() {
        s=0;
    }

    double syn(double deltaT, double glu, double Vm)
    {
        double I_gaba = g_GABA * s * (Vm - V_GABA);
        s = s + (a_r * glu * (1 - s) - a_d * s);
        return I_gaba;
    }
};


//   Ermentrout and Terman, Mathematical Foundations of Neuroscience, 2010, p163.
//   There are several models for the activation of GABA_B synapses.  Their book
//   presented the simplest.     Destexhe 1998 presents a more complex one.
//  


class GABA_B_spine
{
    public:
    //   units are mM  and  ms
    static const double Tmax= 1;    
    static const double V_T = 2;          
    static const double Kp  = 5;
    static const double a_r = 0.09;  
    static const double a_d = 0.0012; 
    static const double E_K =-95;    //   -90 ... -105, Ermentrout, FMN p163
    static const double n   = 4;  
    static const double Kd  = 100; 
    static const double K3  = 0.18;  
    static const double K4  = 0.034; 
    static const double g_GABA=1.0e-6;   // mS   Destexhe 1998: max g_gaba_b = 1.0 nS per syn
    
    double r, s;
    
    GABA_B_spine() {
     r=0;
     s=0;
    }    


    double syn(double glu, double Vm)
    {
        //Tglu = Tmax/(1 + exp(-(V_pre - V_T)/Kp) )
         
        //   the current I_GABA is a nonlinear saturating function of s 
        double I_GABA = g_GABA * ( pow(s,n)/( pow(s,n) + Kd )) * (Vm - E_K);
        
        r = r + a_r * glu * (1 - r) - a_d * r;
 
        s = s + K3 * r - K4 * s;
        
        return I_GABA;
    }
};


        

