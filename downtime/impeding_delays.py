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
This module calculates impeding factor delays

"""

def IF_calc(total_cost, coeff_financing, RCmax_RS, N_DMG_RC5, N_DMG_RC3_RS, DMG_input, DL_summary_input, n_structural, t_structural, n_facade, t_facade, IF_delays_input):
    import numpy as np
    import pandas as pd
    
    DMG = pd.read_csv(DMG_input, low_memory=False, header=None)
    PG = DMG.loc[1, 1:] #performance groups indicating the location of fragility
    
    FG = DMG.loc[0, 1:] #fragility groups
    n_PG = len(FG)
    DL_summary=pd.read_csv(DL_summary_input)
    IF_delays=pd.read_csv(IF_delays_input)
    med_RC1=IF_delays['Median_RC=1']
    med_RC2=IF_delays['Median_RC>1']
    disp_RC1=IF_delays['Dispersion_RC=1']
    disp_RC2=IF_delays['Dispersion_RC>1']


    cases_collapse = DL_summary["collapses/collapsed"]
    indx_collapse = cases_collapse[cases_collapse==1].index
    cases_irreparable = DL_summary["reconstruction/irreparable"]
    indx_irreparable = cases_irreparable[cases_irreparable==1].index
    indx_all = DL_summary["#Num"]
    indx_repairable =  indx_all.drop(indx_collapse.union(indx_irreparable))
    
    story_FG = []
    for i in range(n_PG):
        story_FG.append((PG[i+1][len(PG[i+1])-4:len(PG[i+1])-1]))
    
    story = list(set(story_FG))
    story.sort()
    
    #%% repairable buildings -- impeding factors
    #inspection
    IF_inspection = np.zeros(len(indx_repairable))
    
    for i in range (len(RCmax_RS)):
        if (max(RCmax_RS[i,:])>=1):
           IF_inspection[i]= np.random.lognormal(np.log(med_RC2[0]),disp_RC2[0])
    IF_inspection_mat2 = np.append(IF_inspection, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #engineering
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(indx_repairable)):
        x=0
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_eng = (n1*np.random.lognormal(np.log(med_RC1[1]),disp_RC1[1],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[1]),disp_RC2[1],len(indx_repairable)))/len(story)
    IF_eng_mat2 = np.append(IF_eng, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #permitting
    IF_permit = (n1*np.random.lognormal(np.log(med_RC1[2]),disp_RC1[2],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[2]),disp_RC2[2],len(indx_repairable)))/len(story)
    IF_permit_mat2 = np.append(IF_permit, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #financing
    repair_cost_all = DL_summary[DL_summary.iloc[:,3]==0]
    repair_cost = repair_cost_all.iloc[:,4]
    ELR = (repair_cost/total_cost).to_numpy()
    deductible = 0.1
    
    IF_finance = np.zeros(len(indx_repairable))
    rand_num = np.random.rand(len(indx_repairable))
    for i in range(len(indx_repairable)):
         #financing is suppoorted by insurance
        if rand_num[i] <=coeff_financing[0]:
            if (ELR[i]>0.05) and (ELR[i]<=deductible):
                IF_finance[i] = 0.5*np.random.lognormal(np.log(med_RC2[4]),disp_RC2[4])
            elif (ELR[i]>deductible):
                IF_finance[i] = max(0.5*np.random.lognormal(np.log(med_RC2[4]),disp_RC2[4]) , np.random.lognormal(np.log(med_RC2[3]),disp_RC2[3]))
         #financing is suppoorted by private loans
        elif coeff_financing[0]< rand_num[i] <=coeff_financing[1]+coeff_financing[0]:
            if (ELR[i]>0.05) and (ELR[i]<=0.1):
                IF_finance[i] = 0.5*np.random.lognormal(np.log(med_RC2[4]),disp_RC2[4])
            elif (ELR[i]>0.1):
                IF_finance[i] = np.random.lognormal(np.log(med_RC2[4]),disp_RC2[4])
        #financing is suppoorted by public loans  
        else: 
            if (ELR[i]>0.1):
                IF_finance[i] = np.random.lognormal(np.log(med_RC2[5]),disp_RC2[5])
            else:
                IF_finance[i] = 0.5*np.random.lognormal(np.log(med_RC2[5]),disp_RC2[5])
    
    IF_finance_mat2 = np.append(IF_finance, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #contractor mobilization for repair sequence 1
    IF_cm_rs1 = (n1*np.random.lognormal(np.log(med_RC1[6]),disp_RC1[6],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[6]),disp_RC2[6],len(indx_repairable)))/len(story)
    IF_cm_rs1_mat2 = np.append(IF_cm_rs1, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #contractor mobilization for repair sequence 2
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(RCmax_RS)):
        x=1
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_cm_rs2 = (n1*np.random.lognormal(np.log(med_RC1[7]),disp_RC1[7],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[7]),disp_RC2[7],len(indx_repairable)))/len(story)
    IF_cm_rs2_mat2 = np.append(IF_cm_rs2, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #contractor mobilization for repair sequence 3
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(RCmax_RS)):
        x=2
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_cm_rs3 = (n1*np.random.lognormal(np.log(med_RC1[8]),disp_RC1[8],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[8]),disp_RC2[8],len(indx_repairable)))/len(story)  
    IF_cm_rs3_mat2 = np.append(IF_cm_rs3, np.zeros(len(indx_irreparable)+len(indx_collapse)))
        
    #contractor mobilization for repair sequence 4
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(RCmax_RS)):
        x=3
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_cm_rs4 = (n1*np.random.lognormal(np.log(med_RC1[9]),disp_RC1[9],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[9]),disp_RC2[9],len(indx_repairable)))/len(story) 
    IF_cm_rs4_mat2 = np.append(IF_cm_rs4, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #contractor mobilization for repair sequence 5
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(RCmax_RS)):
        x=4
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_cm_rs5 = (n1*np.random.lognormal(np.log(med_RC1[10]),disp_RC1[10],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[10]),disp_RC2[10],len(indx_repairable)))/len(story)
    IF_cm_rs5_mat2 = np.append(IF_cm_rs5, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #contractor mobilization for repair sequence 6
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(RCmax_RS)):
        x=5
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_cm_rs6 = (n1*np.random.lognormal(np.log(med_RC1[11]),disp_RC1[11],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[11]),disp_RC2[11],len(indx_repairable)))
    IF_cm_rs6_mat2 = np.append(IF_cm_rs6, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    #contractor mobilization for repair sequence 7
    n1=np.zeros(len(indx_repairable))
    n2=np.zeros(len(indx_repairable))
    for i in range (len(RCmax_RS)):
        x=6
        for j in range(len(story)):
            if ((RCmax_RS[i,x])==1):
                n1[i]=n1[i]+1
            elif ((RCmax_RS[i,x])>1):
                n2[i]=n2[i]+1
            x = x+7
    IF_cm_rs7 = (n1*np.random.lognormal(np.log(med_RC1[12]),disp_RC1[12],len(indx_repairable)) + n2*np.random.lognormal(np.log(med_RC2[12]),disp_RC2[12],len(indx_repairable)))/len(story)
    IF_cm_rs7_mat2 = np.append(IF_cm_rs7, np.zeros(len(indx_irreparable)+len(indx_collapse)))
    
    # Stabilization
    IF_stab_RC5 = np.zeros((len(RCmax_RS))) #structrual repairs
    IF_stab_mob = 0
    N_DMG_RC5_tot=sum(N_DMG_RC5.transpose())
    k1 = (t_structural[0]-t_structural[1])/(n_structural[1]-n_structural[0])
    k2 = (t_facade[0]-t_facade[1])/(n_facade[1]-n_facade[0])
    t1 = k1*n_structural[0]+t_structural[0]
    t2 = k2*n_facade[0]+t_facade[0]
    for i in range(len(RCmax_RS)):
        if N_DMG_RC5_tot[i]>0 and N_DMG_RC5_tot[i] <= n_structural[0]:
            IF_stab_RC5[i] = IF_stab_mob + N_DMG_RC5_tot[i]*np.random.lognormal(np.log(t_structural[0]),0.4)
        elif (N_DMG_RC5_tot[i] >n_structural[0]) and (N_DMG_RC5_tot[i] <n_structural[1]):
           IF_stab_RC5[i] = IF_stab_mob + N_DMG_RC5_tot[i]*np.random.lognormal(np.log(-k1*N_DMG_RC5_tot[i]+t1),0.4)
        elif N_DMG_RC5_tot[i] >=n_structural[1]:
           IF_stab_RC5[i] = IF_stab_mob + (N_DMG_RC5_tot[i]*np.random.lognormal(np.log(t_structural[1]),0.4))/2
           
    IF_stab_RC3 = np.zeros((len(RCmax_RS))) #curtain walls
    N_DMG_RC3_RS_mat = np.squeeze(N_DMG_RC3_RS)[:,np.arange(2,len(np.squeeze(N_DMG_RC3_RS).transpose()),7)]
    N_DMG_RC3_RS3 = sum(N_DMG_RC3_RS_mat.transpose())
    
    for i in range(len(RCmax_RS)):
        if N_DMG_RC3_RS3[i]>0 and N_DMG_RC3_RS3[i] <=n_facade[0]:
            IF_stab_RC3[i] = IF_stab_mob + N_DMG_RC3_RS3[i]*np.random.lognormal(np.log(t_facade[0]),0.4)
        elif N_DMG_RC3_RS3[i]>n_facade[0] and N_DMG_RC3_RS3[i] <n_facade[1]:
            IF_stab_RC3[i] = IF_stab_mob + N_DMG_RC3_RS3[i]*np.random.lognormal(np.log(-k2*N_DMG_RC3_RS3[i]+t2),0.4)
        elif N_DMG_RC3_RS3[i]>= n_facade[1]:
            IF_stab_RC3[i] = IF_stab_mob + N_DMG_RC3_RS3[i]*np.random.lognormal(np.log(t_facade[1]),0.4)/2
    
    IF_stab = np.maximum(IF_stab_RC3, IF_stab_RC5) #total stabilization time
    IF_stab_mat2 = np.append(IF_stab, np.zeros(len(indx_irreparable)+len(indx_collapse)))
            
        
    #%% impeding factor delays for irreparable buildings
    IF_reconst_eng = np.random.lognormal(np.log(med_RC2[13]),disp_RC2[13],len(indx_irreparable)+len(indx_collapse))
    IF_reconst_insur = np.random.lognormal(np.log(med_RC2[14]),disp_RC2[14],len(indx_irreparable)+len(indx_collapse))
    IF_reconst_demol = np.random.lognormal(np.log(med_RC2[15]),disp_RC2[15],len(indx_irreparable)+len(indx_collapse))
    
    IF_reconst = np.maximum(IF_reconst_eng, IF_reconst_insur, IF_reconst_demol)
    IF_reconst_mat2 = np.append(np.zeros(len(indx_repairable)),IF_reconst)

    IF_matrix = [IF_inspection_mat2.astype(object), IF_eng_mat2.astype(object), IF_permit_mat2.astype(object), IF_finance_mat2.astype(object), IF_cm_rs1_mat2.astype(object), 
                 IF_cm_rs2_mat2.astype(object), IF_cm_rs3_mat2.astype(object), IF_cm_rs4_mat2.astype(object), IF_cm_rs5_mat2.astype(object), 
                 IF_cm_rs6_mat2.astype(object), IF_cm_rs7_mat2.astype(object), IF_stab_mat2.astype(object), IF_reconst_mat2.astype(object)]
    
    (pd.DataFrame(IF_matrix)).T.to_csv('IF_delays.csv', header=['IF_inspection','IF_eng','IF_permit','IF_finance',
                                                                'IF_cm_RS1','IF_cm_RS2','IF_cm_RS3','IF_cm_RS4',
                                                                'IF_cm_RS5','IF_cm_RS6','IF_cm_RS7','IF_stab','IF_reconstruct'], index_label='#Num')

    print('Impeding factor calculation is completed')
    output = [IF_inspection, IF_eng, IF_permit, IF_finance, IF_cm_rs1, IF_cm_rs2, IF_cm_rs3, IF_cm_rs4, IF_cm_rs5, IF_cm_rs6, IF_cm_rs7, IF_stab, IF_reconst]
    return output