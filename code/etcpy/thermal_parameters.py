import numpy as np
import pandas as pd
from scipy.optimize import fsolve


def calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt):
    '''
    # Using Trainsition state theory to calculate kcat at temperature T.
    # dHTH, dSTS: entropy and enthalpy at comergence temperatures. Protein
    # unfolding process.
    # dCpu, heat capacity change unpon unfolding.
    # kcatTopt: kcat values at optimal temperature
    # Topt, optimal temperature of the enzyme, in K
    # T, temperature, in K
    #
    '''
    # Constants
    R = 8.314
    TH = 373.5
    TS = 385
    T0 = 30+273.15

    # Use the equation from solvedHT.m and re-organized
    dGuTopt = dHTH +dCpu*(Topt-TH) -Topt*dSTS-Topt*dCpu*np.log(Topt/TS)
    dHt = dHTH+dCpu*(Topt-TH)-dCpt*(Topt-T0)-R*Topt-(dHTH+dCpu*(Topt-TH))/(1+np.exp(-dGuTopt/(R*Topt)))

    # Calculate kcat at reference Temperautre
    kcat0 = kcatTopt/np.exp(np.log(Topt/T0)-(dHt+dCpt*(Topt-T0))/R/Topt+dHt/R/T0+dCpt*np.log(Topt/T0)/R)

    # Calculate kcat at given temperature
    kcatT = kcat0*np.exp(np.log(T/T0)-(dHt+dCpt*(T-T0))/R/T+dHt/R/T0+dCpt*np.log(T/T0)/R)

    return kcatT


def get_dH_dS_dCpu_from_TmT90(Tm,T90):
    '''
    # With knowing Tm and T90, get dHTH, dSTS and dCpu by solving
    # 
    # dHTH = slope*dSTS+intercept
    # deltaG(Tm) = 0
    # deltaG(T90) = -RTln9
    # 
    # to get dSTS, dCpu.
    # 
    # Tm, T90 are in K
    # dHTH, is in J/mol
    # dSTS is in J/mol/K
    # 
    '''
    TH = 373.5
    TS = 385
    R = 8.314
    slope = 299.58
    intercept = 20008
    
    a = np.array([[1,-slope,0],
                 [1,-Tm,Tm-TH-Tm*np.log(Tm/TS)],
                 [1,-T90,T90-TH-T90*np.log(T90/TS)]])
    b = np.array([intercept,0,-R*T90*np.log(9)])
    
    [dHTH,dSTS,dCpu] = list(np.linalg.solve(a,b))

    return dHTH,dSTS,dCpu

def get_dH_dS_dCpu_from_TmLength(Tm,N):
    '''
    # In case of negative obtained from get_dH_dS_dCpu_from_TmT90(Tm,T90), or this is no T90 data available, 
    # using the same euqations from  Sawle and Ghosh, Biophysical Journal, 2011 for delatH* and deltaS*.
    # Then caculate dCpu by solving deltaG(Tm) =0
    # Tm is in K
    # 
    # Tm, T90 are in K
    # dHTH, is in J/mol
    # dSTS is in J/mol/K
    # 
    '''
    dHTH = (4*N+143)*1000
    dSTS = 13.27*N+448
    
    TH = 373.5
    TS = 385
    
    def func(dCp):
        dGTm = dHTH + dCp*(Tm-TH)-Tm*dSTS -Tm*dCp*np.log(Tm/TS)
        return dGTm
    dCpu = fsolve(func,10000)[0]
    return dHTH,dSTS,dCpu

def get_dGu(T: float,dHTH: float,dSTS: float,dCpu: float):
    '''
    # calculate the deltaG of unfolding process at temperature T
    # dHTH, dSTS are enthalpy and entropy at TH and TS, repspectively
    # dCpu is the heat capacity change unpoin unfolding
    '''
    TH = 373.5
    TS = 385

    dGu: float = dHTH +dCpu*(T-TH) -T*dSTS-T*dCpu*np.log(T/TS);
    return dGu


def get_fNT(T: float,dHTH: float,dSTS: float,dCpu: float):
    '''
    # Calculate the fraction of enzyme in native state
    # dHTH, dSTS are enthalpy and entropy at TH and TS, repspectively
    # dCpu is the heat capacity change unpoin unfolding
    # T, temeperature, a single value or a numpy array, in K
    '''
    R = 8.314;
    dGu = get_dGu(T,dHTH,dSTS,dCpu);
    f: float = 1/(1+np.exp(-dGu/R/T));
    return f

def get_dGu(T: float,dHTH: float,dSTS: float,dCpu: float):
    '''
    # calculate the deltaG of unfolding process at temperature T
    # dHTH, dSTS are enthalpy and entropy at TH and TS, repspectively
    # dCpu is the heat capacity change unpoin unfolding
    '''
    TH = 373.5;
    TS = 385;

    dGu: float = dHTH +dCpu*(T-TH) -T*dSTS-T*dCpu*np.log(T/TS);
    return dGu


def get_fNT(T: float,dHTH: float,dSTS: float,dCpu: float):
    '''
    # Calculate the fraction of enzyme in native state
    # dHTH, dSTS are enthalpy and entropy at TH and TS, repspectively
    # dCpu is the heat capacity change unpoin unfolding
    # T, temeperature, a single value or a numpy array, in K
    '''
    R = 8.314;
    dGu = get_dGu(T,dHTH,dSTS,dCpu);
    f: float = 1/(1+np.exp(-dGu/R/T));
    return f

def getNGAMT(T):
    # T is in K, a single value
    def NGAM_function(T):
        return 0.740 + 5.893/(1+np.exp(31.920-(T-273.15))) + 6.12e-6*(T-273.15-16.72)**4

    lb = 5+273.15
    ub = 40+273.15
    if T < lb: NGAM_T = NGAM_function(lb)
    elif T > ub: NGAM_T =NGAM_function(ub)
    else: NGAM_T = NGAM_function(T)

    return NGAM_T

def calculate_thermal_params(params):
    '''
    # params, a dataframe with at least following columns: Tm,T90,Length,dCpt,Topt,dCpt. All are in standard units.
    # 
    # The script will return a dataframe with following columns: dHTH,dSTS,dCpu,Topt,dCpt
    # 
    '''
    thermalparams = pd.DataFrame()
    
    # step 1: calculate dHTH,dSTS,dCpu from tm, t90/length
    for ind in params.index:
        tm,t90 = params.loc[ind,'Tm'],params.loc[ind,'T90']
        if t90 is not None and t90>tm and ~np.isnan(t90):
            dHTH,dSTS,dCpu = get_dH_dS_dCpu_from_TmT90(tm,t90)
            if dCpu <0: dHTH,dSTS,dCpu = get_dH_dS_dCpu_from_TmLength(tm,params.loc[ind,'Length'])
        else: dHTH,dSTS,dCpu = get_dH_dS_dCpu_from_TmLength(tm,params.loc[ind,'Length'])
         
        thermalparams.loc[ind,'dHTH'] = dHTH
        thermalparams.loc[ind,'dSTS'] = dSTS
        thermalparams.loc[ind,'dCpu'] = dCpu
        
        
    # step 2. copy columns Topt and dCpt
    thermalparams['Topt'] = params['Topt']
    thermalparams['dCpt'] = params['dCpt']
    
    return thermalparams
