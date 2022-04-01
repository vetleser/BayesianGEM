import cobra
from cobra import Model
import numpy as np
import pandas as pd
from .thermal_parameters import *


def solve_unboundedness(model: Model, extreme_bound: float = 1000):
    """This function ensures no bound has an infinite magnitude, but are set to a threshold
    The input model is changed in-place

    Args:
        model (Model): A Cobra model to the modified
        extreme_bound (float): The extreme absolute value for bounds. Any bound exceeding this will be set
        to this value respecting the its sign
    """
    for reaction in model.reactions:
        if reaction.lower_bound < extreme_bound:
            reaction.lower_bound = - extreme_bound
        if reaction.upper_bound > extreme_bound:
            reaction.upper_bound = extreme_bound  


def change_rxn_coeff(rxn: cobra.Reaction,met: cobra.Metabolite,new_coeff: float):
    '''
    # This is based on the rxn.add_metabolites function. If there the metabolite is already in the reaction,
    # new and old coefficients will be added. For example, if the old coeff of metA is 1, use
    # rxn.add_metabolites({metA:2}), After adding, the coeff of metA is 1+2 = 3
    #
    '''

    diff_coeff = new_coeff-rxn.metabolites[met]
    rxn.add_metabolites({met:diff_coeff})






def map_fNT(model: Model,T: float,df: pd.DataFrame,Tadj: float=0):
    '''
    # apply the fraction of enzymes in native state to each protein.
    # model, cobra model
    # T, temperature, in K, float
    # Tadj, This is to adjust the orginal denaturation curve by moving to left by
    # Tadj degrees.
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu
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

    rxn: cobra.Reaction
    for rxn in met.reactions:
        # this is to ignore reaction 'prot_pool_exchange': --> prot_pool
        if len(rxn.metabolites)<2: continue

        uniprot_id: str = rxn.id.split('_')[-1]
        cols = ['dHTH', 'dSTS','dCpu','Topt']
        [dHTH, dSTS,dCpu,topt] = df.loc[uniprot_id,cols]
        fNT = get_fNT(T+Tadj,dHTH,dSTS,dCpu)
        if fNT < 1e-32: fNT = 1e-32
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


def map_kcatT(model: Model,T: float,df: pd.DataFrame):
    '''
    # Apply temperature effect on enzyme kcat.
    # based on trainsition state theory
    # model, cobra model
    # T, temperature, in K
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ensure that Topt is in K. Other parameters are in standard units.
    #
    # Gang Li, 2019-05-03
    #
    '''
    met: cobra.Metabolite
    for met in model.metabolites:
        

        # look for those metabolites: prot_uniprotid
        if not met.id.startswith('prot_'): continue

        # ingore metabolite: prot_pool
        if met.id == 'prot_pool': continue
        uniprot_id = met.id.split('_')[1]

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

            # assume that Topt in the original model is measured at Topt
            kcatTopt = -1/rxn.metabolites[met]
            kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
            if kcatT < 1e-32: kcatT = 1e-32
            new_coeff = -1/kcatT
            change_rxn_coeff(rxn,met,new_coeff)


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

