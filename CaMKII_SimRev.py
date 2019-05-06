# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 13:33:27 2019

@author: yasudar
"""
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from matplotlib.backends.backend_pdf import PdfPages

"""
############################################
CaMKII_Sim(phospho_rate=1, phosphatase=0, autonomous=1, binding_To_PCaMK=0, subtract_baseline=False, SpikeTiming=False, ReleaseProbability=1)
#############################################
phospho_rate = 1 #Relative phosphorylation rate. 1 for WT, 0 for T286A phsophorylation.
phosphatase = 0 #Relative phosphatase activity. 1 for WT, 0 for T286D.

autonomous = 1 #conformation at autonomous states compared to CaM-bound states. In the simplest model, default is 1.
binding_To_PCaMK = 0 #Binding of Ca/CaM to phosphorylated CaMKII compared to non-phosphorylated CaMKII. 
    #In many models, this is considered to be zero. However, we observed
    #binding of CaM to CaMKII-T286D/T305D/T306D

subtract_baseline = False #If you want to plot changes from baseline, put true.

SpikeTiming = True #False for usual uncaging. 
ReleaseProbability = 1 #Only for spike timing.

Return value:
    (t1, activeCaMKII, CaM_Bound, PhosphoCaMKII,CaMKP, CaMKP2)
    t1: time
    activeCaMKII: Active CaMKII
    CaM_Bound: Fraction of CaMKII in CaM-bound state.
    PhosphoCaMKII: Fraction of CaMKII in phosphorylated state (at T286).
    CaMKP: Fraction of CaMKII in P state (see Fig. 6 Chang et al., Nat Commu)
    CaMKP2: Fraction of CaMKII in P2 state (see Fig. 6 Chang et al., Nat Commu)
##############################################
"""
def CaMKII_Sim(phospho_rate=1, phosphatase=1, autonomous=1, binding_To_PCaMK=0, subtract_baseline=False, SpikeTiming=False, ReleaseProbability=1):
    ##########################################
    baseLineTime = 200 #seconds.#Wait for equilibration.
    endTime = 200 #seconds
    dt = 1e-4 #time step in sconds.Needs at least 1e-4.
    
    nSamples = round((baseLineTime + endTime)/dt) #Number of samples.
    t1 = np.array(range(nSamples))*dt #Time. 
    
    #########################################
    #Ca2+ time course.
    #########################################
    if not SpikeTiming:
    ### Usual uncaging protocol.
        print("Calculating Ca2+...(Uncaging)")
        pulseN = 30 #number of uncaging pulses
        Ca_decay = 0.1 #seconds.
        pulseInterval = 2 #seconds
        peakCa = 4e-6 #molar. From Evans et al 2018.
        Ca_decay = 0.1 #seconds. NMDA receptor opening. Evans et al., 2018.
        
        #Ca2+ time course.
        pulseI = round(pulseInterval/dt) #Conert to time step
        Ca = np.zeros(nSamples)
        CaCore = peakCa*np.exp(-t1/Ca_decay) #Ca2+ in response to single uncaging
        baseLine = round(baseLineTime/dt) #Convert to time step
        
        for i in range(pulseN):
            CaCore2 = (0.5*np.exp(-i/5)+0.5)*CaCore
            Ca[baseLine+i*pulseI:] = Ca[baseLine+i*pulseI:] + CaCore2[:-i*pulseI -baseLine]
    else:
    #Spike-timning protocol LTP.
        print("Calculating Ca2+...(STDP)")
        Frequency_STDP = 2; #Hz only for spike timing.
        pulseN = 120 #number of uncaging pulses. pulse interval = 2s.
        peakCaAP = 0.8e-6 # Sabatini et al., 2000, ~1uM
        peakCaPair = 2.4e-6 #When paired, 3x more. Siller et al.
        pulseInterval = 1/Frequency_STDP #seconds
        
        pulseI = round(pulseInterval/dt) 
        releaseEvent = np.random.random_sample(pulseI,) < ReleaseProbability
        peakCa = (peakCaPair - peakCaAP) * releaseEvent + peakCaAP #Only when release occurs, Ca2+ is 3 times higher.
        Ca_decay = 0.02 #seconds. #Sabatini et al., 2002
        
        Ca = np.zeros(nSamples) #Initilize Ca2+. 
        
        baseLine = round(baseLineTime/dt)
        
        for i in range(pulseN): #1--29.
           CaCore = peakCa[i] * np.exp(-t1/Ca_decay)
           Ca[baseLine+i*pulseI: ] = Ca[baseLine+i*pulseI: ] + CaCore[:-i*pulseI - baseLine]
           
    Ca = Ca + 50e-9 #Add resting Ca2+. 50nM is enough to activate CaMKII a little bit.
    
    #########################################
    #Kinetic Parameters
    #########################################
    
    CaMT = 30e-6 #Total calmodulin concentration.
    CaMKII_T = 70e-6 #Total CaMKII concentration.
    
    CaMK = np.zeros(nSamples) #Concentration of K.
    CaMKP = np.zeros(nSamples) #Concentration of P.
    CaMKP2 = np.zeros(nSamples) #Concentration of P2.
    
    NBindingSite = 4
    
    #Array is concentrations of 0 Ca2+, 2 Ca2+ to C-lobe, 2 Ca2+ to N-lobe and 4 Ca.
    CaMCa = np.zeros((NBindingSite, nSamples)) #Concentration of calcium bound CaM
    CaM_CaMK = np.zeros((NBindingSite, nSamples)) #Concentration of calcium bound KCaM.
    CaM_CaMKP = np.zeros((NBindingSite, nSamples)) #Concentraiton of calcium bound PCaM.
    
    #Initial value.
    CaMCa[0] = CaMT
    CaMKP[0] = CaMKII_T
    
    #CaM-CaMKII binding model from:
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2820514/#pcbi.1000675.s001
    #
    
    #Ca2+ binding to CaM.
    #C-lobe
    k_1C_on = 5e6 #1.2-9.6uM-1s-1
    k_1C_off = 50 #10-70 s-1
    k_2C_on = 10e6 #5-25uM-1s-1.
    k_2C_off = 10 #8.5-10s-1.
    
    #N-lobe
    k_1N_on = 100e6 #25-260uM-1s-1
    k_1N_off = 2000 #1000-4000 s-1
    k_2N_on = 200e6 #50-300uM-1s-1.
    k_2N_off = 500 #500->1000.-1
    
    #Ca2+ binding to KCaM
    #C-lobe
    k_K1C_on = 44e6
    k_K1C_off = 33
    k_K2C_on = 44e6
    k_K2C_off = 0.8 #0.49-4.9s-1
    
    #N-lobe
    k_K1N_on = 76e6
    k_K1N_off = 300
    k_K2N_on = 76e6 
    k_K2N_off = 20 #6-60-1
    #
    
    #CaM binding to CaMKII.
    kCaM0_on = 3.8e3 
    kCaM0_off = 5.5
    #kCaM1C_on = 59e3 #Not necessary for this simulation.
    #kCaM1C_off = 6.1 #Not necessary for this simulation.
    kCaM2C_on = 0.92e6
    kCaM2C_off = 6.8
    kCaM2N_on = 0.12e6
    kCaM2N_off = 1.7
    kCaM4_on = 30e6 #14-60
    kCaM4_off = 1.5 #1.1-2.3
    
    #CaM binding to CaMKII array
    kon_KCaM = np.array([kCaM0_on, kCaM2C_on, kCaM2N_on, kCaM4_on]) #CaM -> CaMCaMKII
    koff_KCaM = np.array([kCaM0_off, kCaM2C_off, kCaM2N_off, kCaM4_off])
    #
    #Calmodulin on / off from CaMKII-P. k2
    decay_CaM = 3 #seconds
    kon_PCaM = kon_KCaM * binding_To_PCaMK 
    koff_PCaM = 1/decay_CaM*np.ones(NBindingSite) #Same decay time constant (3 s) for all Ca2+ bound states.
    
    k_phosCaM = 30 * phospho_rate #s-1 #KCaM --> PCaM event. 6.3s-1 Lucicci et al.2008, but a larger value is better. 
    
    #De phosphorylation pathways
    k_dephospho = 1/6 * phosphatase #k3, P --> K Rate of dephosphorylation = 6s. 
    k_P1_P2 = 1/60 #k5, P --> P2, Holding for 60 s. 
    k_P2_P1 = 1/6*0.25 #k4. K --> P2, 0.25 times slower than P --> K.
    
    #########################################
    #Simulation body. Use dA = k * A * dt
    #########################################
    print("Starting simulation...")
    for i in range(1,nSamples):
        #V, VK, VP... are the rate of CaM, K, P increases. We first calculate the rate, and 
        #then apply to the value.
        #X[0]: 0 Ca2+
        #X[1]: 2 Ca2+ in C lobe.
        #X[2]: 2 Ca2+ in N lobe
        #X[3]: 4 Ca2+
        
        #Based on Pepke et al., 2010. Two Ca2+ ions bind to C or N-lobe.
        #2 Ca + CaM0 <-> Ca2CaM-C
        #2 Ca + CaM0 <-> Ca2CaM-N
        #2 Ca + Ca2CaM-C <-> Ca4CaM
        #2 Ca + Ca2CaM-N <-> Ca4CaM
        
        C_on = k_1C_on*k_2C_on/(k_1C_off+k_2C_on*Ca[i-1])
        N_on = k_1N_on*k_2N_on/(k_1N_off+k_2N_on*Ca[i-1])
        C_off = k_1C_off*k_2C_off/(k_1C_off+k_2C_on*Ca[i-1])
        N_off = k_1N_off*k_2N_off/(k_1N_off+k_2N_on*Ca[i-1])
        #
        kon = [C_on,  N_on,  N_on, C_on] 
        koff = [C_off, N_off, N_off, C_off]
    
        V = np.zeros(NBindingSite) #Rate CaM
        V[0] = -kon[0]*Ca[i-1]**2*CaMCa[0,i-1] + koff[0]*CaMCa[1,i-1] + koff[1]*CaMCa[2,i-1] - kon[1]*Ca[i-1]**2*CaMCa[0,i-1] #CaM0
        V[1] = -kon[2]*Ca[i-1]**2*CaMCa[1,i-1] + koff[2]*CaMCa[3,i-1] - koff[0]*CaMCa[1,i-1] + kon[0]*Ca[i-1]**2*CaMCa[0,i-1] #Ca2-C
        V[2] = -kon[3]*Ca[i-1]**2*CaMCa[2,i-1] + koff[3]*CaMCa[3,i-1] - koff[1]*CaMCa[2,i-1] + kon[1]*Ca[i-1]**2*CaMCa[0,i-1]  #Ca2-N
        V[3] = +kon[2]*Ca[i-1]**2*CaMCa[1,i-1] - koff[2]*CaMCa[3,i-1] - koff[3]*CaMCa[3,i-1] + kon[3]*Ca[i-1]**2*CaMCa[2,i-1]  #Ca4
    
        #Pepke et al., 2010. Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII complex.
        #2 Ca + CaM0-CaMKII <-> Ca2CaM-C-CaMKII 
        #2 Ca + CaM0-CaMKII <-> Ca2CaM-N-CaMKII 
        #2 Ca + Ca2CaM-C-CaMKII <-> Ca4CaM-CaMKII
        #2 Ca + Ca2CaM-N-CaMKII <-> Ca4CaM-CaMKII
        
        KC_on = k_K1C_on*k_K2C_on/(k_K1C_off+k_K2C_on*Ca[i-1])
        KN_on = k_K1N_on*k_K2N_on/(k_K1N_off+k_K2N_on*Ca[i-1])
        KC_off = k_K1C_off*k_K2C_off/(k_K1C_off+k_K2C_on*Ca[i-1])
        KN_off = k_K1N_off*k_K2N_off/(k_K1N_off+k_K2N_on*Ca[i-1])
        #
        kon_K = [KC_on,  KN_on,  KN_on, KC_on] 
        koff_K = [KC_off, KN_off, KN_off, KC_off]
    
        VK = np.zeros(NBindingSite) #Rate of CaM-CaMKII
        VK[0] = -kon_K[0]*Ca[i-1]**2*CaM_CaMK[0,i-1] + koff_K[0]*CaM_CaMK[1,i-1] + koff_K[1]*CaM_CaMK[2,i-1] - kon_K[1]*Ca[i-1]**2*CaM_CaMK[0,i-1]
        VK[1] = -kon_K[2]*Ca[i-1]**2*CaM_CaMK[1,i-1] + koff_K[2]*CaM_CaMK[3,i-1] - koff_K[0]*CaM_CaMK[1,i-1] + kon_K[0]*Ca[i-1]**2*CaM_CaMK[0,i-1]    
        VK[2] = -kon_K[3]*Ca[i-1]**2*CaM_CaMK[2,i-1] + koff_K[3]*CaM_CaMK[3,i-1] - koff_K[1]*CaM_CaMK[2,i-1] + kon_K[1]*Ca[i-1]**2*CaM_CaMK[0,i-1]
        VK[3] = +kon_K[2]*Ca[i-1]**2*CaM_CaMK[1,i-1] - koff_K[2]*CaM_CaMK[3,i-1] - koff_K[3]*CaM_CaMK[3,i-1] + kon_K[3]*Ca[i-1]**2*CaM_CaMK[2,i-1]  
        
        #Binding of Ca to CaM-CaMKIIP. The scheme and values are the same as K.
        kon_P = kon_K
        koff_P = koff_K
        #
        VP = np.zeros(NBindingSite) #Rate of CaM-CaMII-P
        VP[0] = -kon_P[0]*Ca[i-1]**2*CaM_CaMKP[0,i-1] + koff_P[0]*CaM_CaMKP[1,i-1] + koff_P[1]*CaM_CaMKP[2,i-1] - kon_P[1]*Ca[i-1]**2*CaM_CaMKP[0,i-1]
        VP[1] = -kon_P[2]*Ca[i-1]**2*CaM_CaMKP[1,i-1] + koff_P[2]*CaM_CaMKP[3,i-1] - koff_P[0]*CaM_CaMKP[1,i-1] + kon_P[0]*Ca[i-1]**2*CaM_CaMKP[0,i-1]    
        VP[2] = -kon_P[3]*Ca[i-1]**2*CaM_CaMKP[2,i-1] + koff_P[3]*CaM_CaMKP[3,i-1] - koff_P[1]*CaM_CaMKP[2,i-1] + kon_P[1]*Ca[i-1]**2*CaM_CaMKP[0,i-1]
        VP[3] = +kon_P[2]*Ca[i-1]**2*CaM_CaMKP[1,i-1] - koff_P[2]*CaM_CaMKP[3,i-1] - koff_P[3]*CaM_CaMKP[3,i-1] + kon_P[3]*Ca[i-1]**2*CaM_CaMKP[2,i-1]  
    
        #Binding of CaM to CaMKII (rate = VK) or CaMII-P (rate = VP)
        #CaM0 <-> CaM0CaMKII
        #2CaCaM-C <-> Ca2CaM-C-CaMKII
        #2CaCaM-N <-> Ca2CaM-N-CaMKII
        #4CaCaM <-> Ca4CaM-CaMKII
        #
        #CaMKP, CaMK = free CaMKII / CaMKII-P
        CaMK_V = 0 #rate of free CaMK increase.
        CaMKP_V = 0 #rate of free CaMKP increase.
        CaMKP2_V = 0 #rate of free CaMKP2 increase.
        for j in range(NBindingSite): #0-3
            V[j] = V[j]   - kon_KCaM[j]*CaMK[i-1]*CaMCa[j,i-1]  - kon_PCaM[j]*CaMKP[i-1]*CaMCa[j,i-1] + koff_KCaM[j]*CaM_CaMK[j,i-1] + koff_PCaM[j]*CaM_CaMKP[j,i-1]
            VK[j] = VK[j] + kon_KCaM[j]*CaMK[i-1]*CaMCa[j,i-1]  - koff_KCaM[j]*CaM_CaMK[j,i-1]
            VP[j] = VP[j] + kon_PCaM[j]*CaMKP[i-1]*CaMCa[j,i-1] - koff_PCaM[j]*CaM_CaMKP[j,i-1]
            
            CaMK_V =  CaMK_V  - kon_KCaM[j]*CaMK[i-1]*CaMCa[j,i-1]  + koff_KCaM[j]*CaM_CaMK[j,i-1]
            CaMKP_V = CaMKP_V - kon_PCaM[j]*CaMKP[i-1]*CaMCa[j,i-1] + koff_PCaM[j]*CaM_CaMKP[j,i-1]
        
        #Phosphorylation CaMXCaMKII -> CaMXCaMKIIP. 
        #Assume no dephosphorylation while binding to CaM. 
        #Fraction active = (CaMKP + sum(CaM_CaMKP,1) + sum(CaM_CaMK,1))/CaMKII_T
        #Rate ~ proportional to the probability that the next subunit is
        #active.
        FA = (CaMKP[i-1] + CaMKP2[i-1] + np.sum(CaM_CaMKP[:,i-1],axis=0) + np.sum(CaM_CaMK[:,i-1],axis=0))/CaMKII_T 
        
        for j in range(1,NBindingSite - 1): 
            VK[j] = VK[j] - k_phosCaM * FA * CaM_CaMK[j,i-1] #K decreases by phosphorylation
            VP[j] = VP[j] + k_phosCaM * FA * CaM_CaMK[j,i-1] #P increases by phosphorylation
        
        #Dephosphorylation.
        # CaMKP -> CaMK
        CaMKP_V = CaMKP_V - k_dephospho * CaMKP[i-1]
        CaMK_V  = CaMK_V + k_dephospho * CaMKP[i-1]
                
        #Second phosphorylation state (P2)
        #CaMKP <-> CaMKP2
        CaMKP2_V = CaMKP2_V + k_P1_P2 * CaMKP[i-1] - k_P2_P1 * CaMKP2[i-1]
        CaMKP_V  = CaMKP_V  - k_P1_P2 * CaMKP[i-1] + k_P2_P1 * CaMKP2[i-1]    
        
        #Finally applying new values 
        #X[i] = X[i-1] + dt*V
        for j in range(NBindingSite):
            CaMCa[j,i] = CaMCa[j,i-1] + dt*V[j]
            CaM_CaMK[j,i] = CaM_CaMK[j,i-1] + dt*VK[j]
            CaM_CaMKP[j,i] = CaM_CaMKP[j,i-1] + dt*VP[j]
        
        CaMK[i] = CaMK[i-1] + dt*CaMK_V
        CaMKP[i] = CaMKP[i-1] + dt*CaMKP_V
        CaMKP2[i] = CaMKP2[i-1] + dt*CaMKP2_V    
        
        PercentDone = i * 100 / nSamples;
        if PercentDone == np.round(PercentDone):
            sys.stdout.write("\rPercent done = %i" % PercentDone)
            sys.stdout.flush()
        
    
    #CaM_CaMKP adn CaM_CaMK are all active.
    #Assuming that autonomnous CaMKII activity is smaller.
    CaMKP = CaMKP * autonomous / CaMKII_T
    CaMKP2 = CaMKP2 * autonomous / CaMKII_T
    PhosphoCaMKII = CaMKP + CaMKP2
    
    CaM_Bound = (np.sum(CaM_CaMKP,axis=0) + np.sum(CaM_CaMK,axis=0))/CaMKII_T
    
    activeCaMKII = CaM_Bound + PhosphoCaMKII
    
    #Subtracting baseline.
    if subtract_baseline:
        activeCaMKII = activeCaMKII - activeCaMKII[baseLine-1]
        CaM_Bound = CaM_Bound - CaM_Bound[baseLine-1]
        PhosphoCaMKII = PhosphoCaMKII - PhosphoCaMKII[baseLine-1]
        CaMKP = CaMKP - CaMKP[baseLine-1]
        CaMKP2 = CaMKP2 - CaMKP2[baseLine-1]
    
    t1 = t1-baseLineTime #Time 0 = first pulse.
    
    #Downsample to save memoty for plotting.
    Downsample = 100;
    nS = int(nSamples/Downsample);
    t1 = np.mean(np.reshape(t1, (nS, Downsample)), axis = 1)
    activeCaMKII = np.mean(np.reshape(activeCaMKII, (nS, Downsample)), axis = 1)
    CaM_Bound = np.mean(np.reshape(CaM_Bound, (nS, Downsample)), axis = 1)
    PhosphoCaMKII = np.mean(np.reshape(PhosphoCaMKII, (nS, Downsample)), axis = 1)
    CaMKP = np.mean(np.reshape(CaMKP, (nS, Downsample)), axis = 1)
    CaMKP2 = np.mean(np.reshape(CaMKP2, (nS, Downsample)), axis = 1)
    
    print("")
    return (t1, activeCaMKII, CaM_Bound, PhosphoCaMKII, CaMKP, CaMKP2)

#####################################################
#Plotting
#####################################################    
if __name__ == "__main__":
    (t1, activeCaMKII, CaM_Bound, PhosphoCaMKII,CaMKP, CaMKP2) = CaMKII_Sim(1, 1, 1, 0, False, False, 1)
    windowResolution = 100
    plotsize = (8,4)
    filename = r'CaMKII_Sim.pdf'
    
    with PdfPages(filename) as export_pdf:
        
        f = plt.figure(figsize = plotsize, dpi=windowResolution) #define the size of the figure.
        
        a = f.add_subplot(1,2,1)
        a.plot(t1, activeCaMKII, '-k')
        a.plot(t1, CaM_Bound, '-r')
        a.plot(t1, PhosphoCaMKII, '-b')
        a.plot(t1, CaMKP, '-c')
        a.plot(t1, CaMKP2, '-m')
        
        a.set_xlim(-2, 20)
        a.set_ylim(0, 0.8)
        
        a = f.add_subplot(1,2,2)
        a.plot(t1, activeCaMKII, '-k')
        a.plot(t1, CaM_Bound, '-r')
        a.plot(t1, PhosphoCaMKII, '-b')
        a.plot(t1, CaMKP, '-c')
        a.plot(t1, CaMKP2, '-m')
        
        a.set_xlim(-10,160)
        a.set_ylim(0, 0.8)
        
        export_pdf.savefig()
        plt.close()
        
    os.startfile(filename)