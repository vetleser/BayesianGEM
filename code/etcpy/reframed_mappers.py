import reframed
from reframed import CBModel
import numpy as np
import pandas as pd
from .thermal_parameters import *


def map_fNT(model: CBModel,T: float,df: pd.DataFrame,Tadj: float=0):
    '''
    # apply the fraction of enzymes in native state to each protein.
    # model, reframed model
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

    met: str = model.metabolites.prot_pool.id

    rxn: reframed.CBReaction
    for rxn in model.reactions.values():
        # We would preferingly only iterate over reactions where
        # the protein pool is present, but unlike cobrapy, we do not
        # have any sound way to do this exect looping over all reactions
        if met not in rxn.stoichiometry: continue
        # this is to ignore reaction 'prot_pool_exchange': --> prot_pool
        if len(rxn.stoichiometry)<2: continue
        uniprot_id: str = rxn.id.split('_')[-1]
        cols = ['dHTH', 'dSTS','dCpu','Topt']
        [dHTH, dSTS,dCpu,topt] = df.loc[uniprot_id,cols]
        fNT = get_fNT(T+Tadj,dHTH,dSTS,dCpu)
        if fNT < 1e-32: fNT = 1e-32
        new_coeff: float = rxn.stoichiometry[met]/fNT
        rxn.stoichiometry[met] = new_coeff
    return



def map_kcatT(model: CBModel,T: float,df: pd.DataFrame):
    '''
    # Apply temperature effect on enzyme kcat.
    # based on trainsition state theory
    # model, reframed model
    # T, temperature, in K
    # df, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ensure that Topt is in K. Other parameters are in standard units.
    #
    # Gang Li, 2019-05-03
    #
    '''
    cols = ['dHTH', 'dSTS','dCpu','Topt','dCpt']
    for rxn in model.reactions.values():
        rxn: reframed.CBReaction
        if rxn.id.startswith('draw_prot'): continue
        for met in rxn.stoichiometry:
            if not met.startswith('prot_'): continue
            # ingore metabolite: prot_pool
            if met == 'prot_pool': continue
            uniprot_id = met.split('_')[1]
            [dHTH, dSTS,dCpu,Topt,dCpt]=df.loc[uniprot_id,cols]
            # Change kcat value.
            # pmet_r_0001 + 1.8518518518518518e-07 prot_P00044 + 1.8518518518518518e-07 prot_P32891 -->
            # 2.0 s_0710 + s_1399
            #
            # 1.8518518518518518e-07 is correponding to 1/kcat
            # change the kcat to kcat(T)
            # In some casese, this coefficient could be 2/kcat or some other values. This doesn't matter.
            #
            # a protein could be involved in several reactions
            # assume that Topt in the original model is measured at Topt
            kcatTopt = -1/rxn.stoichiometry[met]
            kcatT = calculate_kcatT(T,dHTH,dSTS,dCpu,kcatTopt,dCpt,Topt)
            if kcatT < 1e-32: kcatT = 1e-32
            new_coeff = -1/kcatT
            rxn.stoichiometry[met] = new_coeff
    return
        


def set_NGAMT(model,T):
    # T is in K
    NGAM_T = getNGAMT(T)
    rxn = model.reactions.NGAM
    #ori_lb,ori_ub = rxn.lower_bound,rxn.upper_bound
    rxn.lb = NGAM_T
    rxn.ub = NGAM_T


def set_sigma(model,sigma):
    rxn = model.reactions.prot_pool_exchange
    #ori_ub_sigma = rxn.upper_bound
    rxn.ub = 0.17866*sigma

