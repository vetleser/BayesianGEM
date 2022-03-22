import cobra
from cobra import Model
import numpy as np
import pandas as pd

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

def change_rxn_coeff(rxn: cobra.Reaction,met: cobra.Metabolite,new_coeff: float):
    '''
    # This is based on the rxn.add_metabolites function. If there the metabolite is already in the reaction,
    # new and old coefficients will be added. For example, if the old coeff of metA is 1, use
    # rxn.add_metabolites({metA:2}), After adding, the coeff of metA is 1+2 = 3
    #
    '''

    diff_coeff = new_coeff-rxn.metabolites[met]
    rxn.add_metabolites({met:diff_coeff})



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


def map_fNT(model: Model,T: float,df: pd.DataFrame,Tadj: float=0, reference_model: Model=None):
    '''
    # apply the fraction of enzymes in native state to each protein.
    # model, cobra model
    # T, temperature, in K, float
    # Tadj, This is to adjust the orginal denaturation curve by moving to left by
    # Tadj degrees.
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu
    # reference_model In the case of warm-starting FBA, the parameters to modify are fetched from the
    # reference model cobra model.
    #
    #
    # Gang Li, 2019-05-03


    # in the enzyme constrained model, protein usage was describe as reaction
    # draw_prot_P0CW41 68.69778 prot_pool --> prot_P0CW41. 68.69778 is the molecular weight of protein P0CW41
    # the total enzyme pool was constrained by the upper bound of reaction prot_pool_exchange:
    #          --> prot_pool
    # To map fNT, change the coefficient of prot_pool in reactions like MW*prot_pool --> prot_P0CW41,
    #                                MW/fNT prot_pool --> prot_P0CW41

    '''

    met: cobra.Metabolite = model.metabolites.prot_pool
    if reference_model is not None:
        reference_met: cobra.Metabolite = reference_model.metabolites.prot_pool

    rxn: cobra.Reaction
    for rxn in met.reactions:
        # this is to ignore reaction 'prot_pool_exchange': --> prot_pool
        if len(rxn.metabolites)<2: continue

        uniprot_id: str = rxn.id.split('_')[-1]
        cols = ['dHTH', 'dSTS','dCpu','Topt']
        [dHTH, dSTS,dCpu,topt] = df.loc[uniprot_id,cols]
        fNT = get_fNT(T+Tadj,dHTH,dSTS,dCpu)
        if fNT < 1e-32: fNT = 1e-32
        if reference_model is not None:
            rxn_id = rxn.id
            reference_rxn: float = reference_model.reactions.get_by_id(rxn_id)
            new_coeff: float = reference_rxn.metabolites[reference_met]/fNT
        else: 
            new_coeff: float = rxn.metabolites[met]/fNT
        change_rxn_coeff(rxn,met,new_coeff)



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
    R = 8.314;
    TH = 373.5;
    TS = 385;
    T0 = 30+273.15;

    # Use the equation from solvedHT.m and re-organized
    dGuTopt = dHTH +dCpu*(Topt-TH) -Topt*dSTS-Topt*dCpu*np.log(Topt/TS);
    dHt = dHTH+dCpu*(Topt-TH)-dCpt*(Topt-T0)-R*Topt-(dHTH+dCpu*(Topt-TH))/(1+np.exp(-dGuTopt/(R*Topt)));

    # Calculate kcat at reference Temperautre
    kcat0 = kcatTopt/np.exp(np.log(Topt/T0)-(dHt+dCpt*(Topt-T0))/R/Topt+dHt/R/T0+dCpt*np.log(Topt/T0)/R);

    # Calculate kcat at given temperature
    kcatT = kcat0*np.exp(np.log(T/T0)-(dHt+dCpt*(T-T0))/R/T+dHt/R/T0+dCpt*np.log(T/T0)/R);

    return kcatT


def map_kcatT(model: Model,T: float,df: pd.DataFrame, reference_model: Model=None):
    '''
    # Apply temperature effect on enzyme kcat.
    # based on trainsition state theory
    # model, cobra model
    # T, temperature, in K
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # reference_model In the case of warm-starting FBA, the parameters to modify are fetched from the
    # reference model cobra model.
    # Ensure that Topt is in K. Other parameters are in standard units.
    #
    # Gang Li, 2019-05-03
    #
    '''
    warm_start = reference_model is not None
    met: cobra.Metabolite
    for met in model.metabolites:
        

        # look for those metabolites: prot_uniprotid
        if not met.id.startswith('prot_'): continue

        # ingore metabolite: prot_pool
        if met.id == 'prot_pool': continue
        uniprot_id = met.id.split('_')[1]

        if reference_model is not None:
            met_id = met.id
            reference_met = reference_model.metabolites.get_by_id(met_id)

        # Change kcat value.
        # pmet_r_0001 + 1.8518518518518518e-07 prot_P00044 + 1.8518518518518518e-07 prot_P32891 -->
        # 2.0 s_0710 + s_1399
        #
        # 1.8518518518518518e-07 is correponding to 1/kcat
        # change the kcat to kcat(T)
        # In some casese, this coefficient could be 2/kcat or some other values. This doesn't matter.
        #
        # a protein could be involved in several reactions
        cols = ['dHTH', 'dSTS','dCpu','Topt','dCpt']
        [dHTH, dSTS,dCpu,Topt,dCpt]=df.loc[uniprot_id,cols]

        rxn: cobra.Reaction
        for rxn in met.reactions:
            if rxn.id.startswith('draw_prot'): continue

            if warm_start:
                # assume that Topt in the original model is measured at Topt
                rxn_id = rxn.id
                reference_rxn = reference_model.reactions.get_by_id(rxn_id)
                kcatTopt = -1/reference_rxn.metabolites[reference_met]
            else: 
                # assume that Topt in the original model is measured at Topt
                kcatTopt = -1/rxn.metabolites[met]
            kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
            if kcatT < 1e-32: kcatT = 1e-32
            new_coeff = -1/kcatT
            if warm_start:
                rxn._metabolites[met] = new_coeff
                {
                    rxn.forward_variable: new_coeff,
                    rxn.reverse_variable: -new_coeff,
                    }
                model.constraints
            else:
                change_rxn_coeff(rxn,met,new_coeff)


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

def set_NGAMT(model,T):
    # T is in K
    NGAM_T = getNGAMT(T)
    rxn = model.reactions.NGAM
    #ori_lb,ori_ub = rxn.lower_bound,rxn.upper_bound
    rxn.lower_bound = NGAM_T
    rxn.upper_bound = NGAM_T


def set_sigma(model,sigma):
    rxn = model.reactions.prot_pool_exchange
    #ori_ub_sigma = rxn.upper_bound
    rxn.upper_bound = 0.17866*sigma

