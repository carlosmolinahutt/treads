# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The University of British Columbia, Vancouver, BC, Canada.
#
# This file is part of downtime assessment framework.
#
# This code is developed to support the framework to estimate earthquake induced downtime 
# and recovery trajectory of residential buildings 
# proposed by Molina Hutt et al. (2021).
#
# The proposed framework can be found in the manuscript entitled 
# Molina Hutt, C., Vahanvaty, T., and Kourehpaz, P. (2021) 
# "an analytical framework to assess earthquake induced downtime and model recovery of buildings", Earthquake Spectra.
# 
# This code will be publicly available on GitHub; however, it should not be redistributed or modified without permission.
#
# Contributor(s):
# Pouria Kourehpaz

"""
This module calculates component repair classes

"""

def RC_calc(DMG_input, RCtable_input, DL_summary_input):
    import sys
    import numpy as np
    import pandas as pd
    from more_itertools import locate
    
    #determine if DMG file contains "nan" entries or not
    DMG = pd.read_csv(DMG_input, low_memory=False, header=None)
    if DMG[4:].isnull().values.any()==True:
        sys.exit("DMG.csv file contains blank entries")
   
    RC_table = pd.read_csv(RCtable_input)
    FG_all = RC_table["Fragility"]
    RC_all = RC_table[["RC-DS0","RC-DS1","RC-DS2","RC-DS3", "RC-DS4", "RC-DS5"]]
    RS_all = RC_table["Repair Sequence"]
    sip_DMG_bld = RC_table["SiP bld level dmg threshold"]
    sip_DMG_str = RC_table["SiP story level dmg threshold"]
    stab_DMG_bld = RC_table["Stab bld level dmg threshold"]
    stab_DMG_str = RC_table["Stab story level dmg threshold"]
    PG = DMG.loc[1, 1:] #performance groups indicating the location of fragility
    
    FG = DMG.loc[0, 1:] #fragility groups
    n_PG = len(FG)
    
    #check if all fragilities are defined in the repair class table
    RS_FG =[]
    for i in range(n_PG):
        indx_FG = FG_all[FG_all==FG[i+1]].index
        if indx_FG.values.size==0:
            sys.exit("Fragility not found. Check Repair_Class_Table.csv file")
        RS_FG.append(RS_all[indx_FG[0]])
    
    RS = list(set(RS_FG))
    RS.sort()
    
    story_FG = []
    for i in range(n_PG):
        story_FG.append((PG[i+1][len(PG[i+1])-4:len(PG[i+1])-1]))    
    story = list(set(story_FG))
    story.sort()
    
    #determine indices for collpase, irreparable, and repairable realizations    
    DL_summary=pd.read_csv(DL_summary_input)
    cases_collapse = DL_summary["collapses/collapsed"]
    indx_collapse = cases_collapse[cases_collapse==1].index
    cases_irreparable = DL_summary["reconstruction/irreparable"]
    indx_irreparable = cases_irreparable[cases_irreparable==1].index
    indx_all = DL_summary["#Num"]
    indx_repairable =  indx_all.drop(indx_collapse.union(indx_irreparable))
    
    
    FG = DMG.loc[0, 1:] #fragility groups
    n_PG = len(FG) #number of perfomance groups
    
    DS=[]
    DSG_DS = DMG.loc[2, 1:] #damage state groups
    for i in range(len(DSG_DS)):
        DS.append(float(DSG_DS[i+1][0]))
    
    RC_matrix = np.zeros((len(indx_repairable),n_PG))
    
    #create DMG file for repairable realizations only
    DMG_all = DMG.loc[4:, 1:]
    for i in range(len(DMG_all)):
        if int(DMG.loc[i+4, 0]) in indx_irreparable:            
            DMG_all.iloc[i,:]=np.nan
    DMG_repairable = DMG_all.dropna()
    DMG_repairable = DMG_repairable.astype('float')
    DMG_repairable = np.array(DMG_repairable)
    
    #create repair class matrix
    DS=pd.Series(DS)
    for j in range(n_PG):
        indx_FG = FG_all[FG_all==FG[j+1]].index
        DS_num = DS[j]
        for i in range(len(indx_repairable)):
            if DMG_repairable[i,j] != 0:
                RC_matrix[i,j] = RC_all.iloc[indx_FG[0], int(DS_num)]
    
    #modify the repair class matrix per the building and story level thresholds for stability recovery state                
    FG_uniq = FG.unique()
    for j in range(len(FG_uniq)): 
        indx_FG = list(locate(FG, lambda x: x == FG_uniq[j]))
        indx_FG_table = list(locate(FG_all, lambda x: x == FG_uniq[j]))
        for i in range(len(indx_repairable)):
            RC_FG = RC_matrix[i][indx_FG]
            indx_RC5_1 = list(locate(RC_matrix[i][indx_FG], lambda x: x == 5))
            indx_RC5 = np.array(indx_FG)[indx_RC5_1]
            loc_FG = np.array(story_FG)[indx_RC5]
            for k in range(len(loc_FG)):
                if (sum(DMG_repairable[i, (indx_RC5)]) <= float(stab_DMG_bld[indx_FG_table]) * np.max(DMG_repairable[:,(indx_FG)],initial=0) * len(set(np.array(PG)[indx_FG])) and DMG_repairable[i, (indx_RC5)][k] <= float(stab_DMG_str[indx_FG_table])*np.max(DMG_repairable[:,(indx_FG)],initial=0)):
                    RC_matrix[i][indx_RC5[k]] = 4
    
    #adjust repair class matrix for components where RC=4 do not exist
    RC_s1 = RC_all[(RC_all == 5).any(axis=1)]
    RC_s2 = RC_s1[(RC_s1 == 4).any(axis=1)]
    RC_s = RC_s1[~RC_s1.isin(RC_s2)].dropna() #find any component that has RC=5 but not RC=4
    comp = list(FG_all[RC_s.index])
    indx_comp=[]
    for i in range(len(comp)):
        indx_comp.append(list(locate(FG, lambda x: x == comp[i])))  
    for i in range(len(indx_repairable)):
        for j in range(len(np.concatenate(indx_comp))):
            if (RC_matrix[i,int(np.concatenate(indx_comp)[j])])==4:
                RC_matrix[i, int(np.concatenate(indx_comp)[j])] = 3
    
    #modify the repair class matrix per the building and story level thresholds for shelter-in-place recovery state
    for j in range(len(FG_uniq)): 
        indx_FG = pd.Series(locate(FG, lambda x: x == FG_uniq[j]))
        indx_FG_table = pd.Series(locate(FG_all, lambda x: x == FG_uniq[j]))
        for i in range(len(indx_repairable)):
            RC_FG = RC_matrix[i][indx_FG]
            indx_RC4_1 = pd.Series(locate(RC_FG, lambda x: x == 4))
            indx_RC4 = np.array(indx_FG)[indx_RC4_1]
            loc_FG = np.array(story_FG)[indx_RC4]
            for k in range(len(loc_FG)):
                if (sum(DMG_repairable[i, (indx_RC4)]) <= float(sip_DMG_bld[indx_FG_table]) * np.max(DMG_repairable[:,(indx_FG)],initial=0) * len(set(np.array(PG)[indx_FG]))) and (DMG_repairable[i, (indx_RC4)][k] <= float(sip_DMG_str[indx_FG_table])*np.max(DMG_repairable[:,(indx_FG)],initial=0)):
                    RC_matrix[i, indx_RC4[k]] = 3
    
    RC_matrix4 = np.c_[np.arange(0,len(indx_repairable),1), RC_matrix]
    RC = pd.concat([DMG.loc[:3, :], pd.DataFrame(RC_matrix4)]) #repair class matrix consistant with DMG file
    RC.to_csv('RC_component.csv', header=False, index=False)
    print('Repair class matrix is generated')
    
    #create a matrix that contains the maximum repair class per repair sequence per realization across the building 
    RCmax_RS=[]
    for i in range(len(indx_repairable)): 
        RC_realization = RC_matrix[i,:]
        RCmax_RS.append([])
        for j in range(len(story)):
            for k in range(len(RS)):
                indx_ST = list(locate(story_FG, lambda x: x == story[j]))
                indx_RS = list(locate(RS_FG, lambda x: x == RS[k]))
                indx_intersect = list(set.intersection(set(indx_ST),set(indx_RS)))
                RCmax_RS[i].append(np.max((RC_realization[indx_intersect]),initial=0))

    print('Repair class calculation is completed')
    return RCmax_RS, RC, indx_repairable, indx_irreparable, indx_collapse, story
