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
This module calculates component repair time and performs worker allocation adjusments

"""

def RT_calc(n_elev, rep_phases, fl_area, RCtable_input, RC_input, DMG_input, DV_rec_time_input, DL_summary_input, max_num_workers):
    import sys
    import numpy as np
    import pandas as pd
    from more_itertools import locate

    RC_table = pd.read_csv(RCtable_input) 
    DMG = pd.read_csv(DMG_input, low_memory=False, header=None)

    RC = RC_input
    DV_rec_time = pd.read_csv(DV_rec_time_input, low_memory=False, header=None)
    rec_time = DV_rec_time.loc[4:,1:]
    DMG_all = DMG.loc[4:, 1:]
    RC_data = RC.iloc[4:,1:]
    FG_all = RC_table["Fragility"]
    Qty_FEMA = RC_table['Qty unit FEMA P-58']
    Wrk_DMG = RC_table['Workers per damaged Qty']
    
    
    FG = DMG.loc[0, 1:] #fragility groups
    n_PG = len(FG)  #number of performance groups
    DL_summary=pd.read_csv(DL_summary_input)
    
    #determine indices for collpase, irreparable, and repairable realizations
    cases_collapse = DL_summary["collapses/collapsed"]
    indx_collapse = cases_collapse[cases_collapse==1].index
    cases_irreparable = DL_summary["reconstruction/irreparable"]
    indx_irreparable = cases_irreparable[cases_irreparable==1].index
    indx_all = DL_summary["#Num"]
    indx_repairable =  indx_all.drop(indx_collapse.union(indx_irreparable))
    
    for i in range(len(DMG_all)):
        if int(DMG.loc[i+4, 0]) in indx_irreparable:            
            DMG_all.iloc[i,:]=np.nan
    DMG_repairable = DMG_all.dropna()
    DMG_repairable = DMG_repairable.astype('float')
    DMG_repairable = np.array(DMG_repairable)
    
    RT_RC2 = np.zeros((len(indx_repairable),n_PG))  #repair times for functional recovery 
    DMG_RC2 = np.zeros((len(indx_repairable),n_PG))
    RT_RC3 = np.zeros((len(indx_repairable),n_PG))  #repair times for reoccupancy
    DMG_RC3 = np.zeros((len(indx_repairable),n_PG))
    RT_RC4 = np.zeros((len(indx_repairable),n_PG))  #repair times for shelter-in-place
    DMG_RC4 = np.zeros((len(indx_repairable),n_PG))
    
    DMG_RC5 = np.zeros((len(indx_repairable),n_PG))
    
    #determine the number of damaged components and repair times for each recovery state
    for i in range(len(indx_repairable)):
        for j in range(n_PG):
            if DMG_repairable[i,j] != 0 and float(RC_data.iloc[i,j]) >= 2:
               RT_RC2[i,j] = float(rec_time.iloc[i,j])
               DMG_RC2[i,j] = DMG_repairable[i,j]
               
               if float(RC_data.iloc[i,j]) >= 3:
                   RT_RC3[i,j] = float(rec_time.iloc[i,j])
                   DMG_RC3[i,j] = DMG_repairable[i,j]
                   
                   if float(RC_data.iloc[i,j]) >= 4:
                       RT_RC4[i,j] = float(rec_time.iloc[i,j])
                       DMG_RC4[i,j] = DMG_repairable[i,j]
                       
                       if float(RC_data.iloc[i,j]) >= 5:
                           DMG_RC5[i,j] = DMG_repairable[i,j]
    
    #create matrix of FEMA P-58 coefficients for fragilities
    Qty_norm=[]
    Wrk_norm=[]
    for i in range(n_PG):
        indx_FG = FG_all[FG_all==FG[i+1]].index
        Qty_norm.append(Qty_FEMA[indx_FG])
        Wrk_norm.append(Wrk_DMG[indx_FG])
    Qty_norm = np.stack( Qty_norm, axis=0)
    Wrk_norm = np.stack( Wrk_norm, axis=0)
    Qty_norm_mat = np.repeat(Qty_norm.T, len(indx_repairable), axis=0)
    Wrk_norm_mat = np.repeat(Wrk_norm.T, len(indx_repairable), axis=0)
    
    #number of damaged components matrix for each recovery state
    N_DMG_RC2 = np.ceil(np.divide(DMG_RC2,Qty_norm_mat)) 
    N_DMG_RC3 = np.ceil(np.divide(DMG_RC3,Qty_norm_mat)) 
    N_DMG_RC4 = np.ceil(np.divide(DMG_RC4,Qty_norm_mat)) 
    N_DMG_RC5 = np.ceil(np.divide(DMG_RC5,Qty_norm_mat))
    
    #edit the number of damaged FEMA P-58 elevator components (simultaneous damage states)
    PG = DMG.loc[1, 1:]
    story_FG = []
    for i in range(n_PG):
        story_FG.append((PG[i+1][len(PG[i+1])-4:len(PG[i+1])-1]))
    
    story = list(set(story_FG))
    story.sort()
    
    indx = (FG[FG=='D1014.011'].index)-1
    for k in range(len(indx_repairable)):
        maxx = np.max(N_DMG_RC2[k][indx], initial=0.0001)
        a = N_DMG_RC2[k][indx]
        a[a<maxx]=0
        N_DMG_RC2[k][indx] = (a*n_elev)/np.max([n_elev,np.sum(a)])
        
    indx = (FG[FG=='D1014.012'].index)-1
    for k in range(len(indx_repairable)):
        maxx = np.max(N_DMG_RC2[k][indx], initial=0.0001)
        a = N_DMG_RC2[k][indx]
        a[a<maxx]=0
        N_DMG_RC2[k][indx] = (a*n_elev)/np.max([n_elev,np.sum(a)])
        
    indx = (FG[FG=='D1014.021'].index)-1
    for k in range(len(indx_repairable)):
        maxx = np.max(N_DMG_RC2[k][indx], initial=0.0001)
        a = N_DMG_RC2[k][indx]
        a[a<maxx]=0
        N_DMG_RC2[k][indx] = (a*n_elev)/np.max([n_elev,np.sum(a)])
    
    indx = (FG[FG=='D1014.022'].index)-1
    for k in range(len(indx_repairable)):
        maxx = np.max(N_DMG_RC2[k][indx], initial=0.0001)
        a = N_DMG_RC2[k][indx]
        a[a<maxx]=0
        N_DMG_RC2[k][indx] = (a*n_elev)/np.max([n_elev,np.sum(a)])     
    
    #an original matrix for number of workers is generated using the values defined in the repair class table
    N_Wrk_RC2 = (np.multiply(N_DMG_RC2,Wrk_norm_mat)) #number of workers matrix
    N_Wrk_RC3 = (np.multiply(N_DMG_RC3,Wrk_norm_mat)) #number of workers matrix
    N_Wrk_RC4 = (np.multiply(N_DMG_RC4,Wrk_norm_mat)) #number of workers matrix
    
    #%% repair time calculation per sequence per story
    
    FG_all = RC_table["Fragility"]
    RS_all = RC_table["Repair Sequence"]
    
    RS_FG =[]
    
    for i in range(n_PG):    
        indx_FG = FG_all[FG_all==FG[i+1]].index
        RS_FG.append(RS_all[indx_FG[0]])
    
    RS = list(set(RS_FG))
    RS.sort()
    
    header=[]
    for i in range(len(story)):
        for j in range(len(RS)):
            header.append(str(story[i])+'_'+str(RS[j]))
    header_ = np.vstack((header, [""]*len(header)))
    col = np.append('St_RSeq', np.append('#Num', np.arange(0,len(indx_repairable),1)))
    
    
    RT_RC2_RS=[]
    N_Wrk_RC2_RS=[]
    RT_RC3_RS=[]
    N_Wrk_RC3_RS=[]
    N_DMG_RC3_RS=[]
    RT_RC4_RS=[]
    N_Wrk_RC4_RS=[]
    N_DMG_RC5_RS=[]
    for i in range(len(indx_repairable)):
        RT_RC2_realization = RT_RC2[i,:]
        N_Wrk_RC2_realization = N_Wrk_RC2[i,:]
        RT_RC2_RS.append([])
        N_Wrk_RC2_RS.append([])
        RT_RC3_realization = RT_RC3[i,:]
        N_Wrk_RC3_realization = N_Wrk_RC3[i,:]
        N_DMG_RC3_realization = N_DMG_RC3[i,:]
        RT_RC3_RS.append([])
        N_Wrk_RC3_RS.append([])
        N_DMG_RC3_RS.append([])
        RT_RC4_realization = RT_RC4[i,:]
        N_Wrk_RC4_realization = N_Wrk_RC4[i,:]
        RT_RC4_RS.append([])
        N_Wrk_RC4_RS.append([])
        N_DMG_RC5_realization = N_DMG_RC5[i,:]
        N_DMG_RC5_RS.append([])
        for j in range(len(story)):
            for k in range(len(RS)):
                indx_ST = list(locate(story_FG, lambda x: x == story[j]))
                indx_RS = list(locate(RS_FG, lambda x: x == RS[k]))
                indx_intersect = list(set.intersection(set(indx_ST),set(indx_RS)))
                RT_RC2_RS[i].append(np.sum(RT_RC2_realization[indx_intersect]))
                N_Wrk_RC2_RS[i].append(np.sum(N_Wrk_RC2_realization[indx_intersect]))
                
                RT_RC3_RS[i].append(np.sum(RT_RC3_realization[indx_intersect]))
                N_Wrk_RC3_RS[i].append(np.sum(N_Wrk_RC3_realization[indx_intersect]))
                N_DMG_RC3_RS[i].append(np.sum(N_DMG_RC3_realization[indx_intersect]))
                
                RT_RC4_RS[i].append(np.sum(RT_RC4_realization[indx_intersect]))
                N_Wrk_RC4_RS[i].append(np.sum(N_Wrk_RC4_realization[indx_intersect]))
                
                N_DMG_RC5_RS[i].append(np.sum(N_DMG_RC5_realization[indx_intersect]))
    
    
    #%% worker allocation limits
  
    N_Wrk_RC2_RS_org = np.squeeze(N_Wrk_RC2_RS)
    N_Wrk_RC3_RS_org = np.squeeze(N_Wrk_RC3_RS) 
    N_Wrk_RC4_RS_org = np.squeeze(N_Wrk_RC4_RS) 
    N_Wrk_RC2_RS_adj = np.squeeze(N_Wrk_RC2_RS)
    N_Wrk_RC3_RS_adj = np.squeeze(N_Wrk_RC3_RS)
    N_Wrk_RC4_RS_adj = np.squeeze(N_Wrk_RC4_RS)
    
    
    max_worker_fl = np.zeros((len(indx_repairable),len(story)))
    for i in range(len(indx_repairable)):
        for j in range(len(story)):
            if (float(np.max(RC.iloc[i+4,:])) >= 4):
                max_worker_fl[i,j] = (1/500)*fl_area[j]
            else: max_worker_fl[i,j] = (1/1000)*fl_area[j]
    
    if (len(story) <= 5):
        max_worker_seq = max_num_workers
    elif (len(story) > 5) and (len(story) <= 20):
        max_worker_seq = [2 * k for k in max_num_workers]
    else:
        max_worker_seq = [3 * k for k in max_num_workers]
    
    #worker limit per repair sequence based on the floor area
    for i in range(len(indx_repairable)):
        for j in range(len(story)):
            N_Wrk_RC2_RS_adj[i,j*7] = min(N_Wrk_RC2_RS_org[i,j*7], max_worker_fl[i,j])
            N_Wrk_RC3_RS_adj[i,j*7] = min(N_Wrk_RC3_RS_org[i,j*7], max_worker_fl[i,j])
            N_Wrk_RC4_RS_adj[i,j*7] = min(N_Wrk_RC4_RS_org[i,j*7], max_worker_fl[i,j])
            
            N_Wrk_RC2_RS_adj[i,j*7+2] = min(N_Wrk_RC2_RS_org[i,j*7+2], max_worker_fl[i,j])
            N_Wrk_RC3_RS_adj[i,j*7+2] = min(N_Wrk_RC3_RS_org[i,j*7+2], max_worker_fl[i,j])
            N_Wrk_RC4_RS_adj[i,j*7+2] = min(N_Wrk_RC4_RS_org[i,j*7+2], max_worker_fl[i,j])
            
            N_Wrk_RC2_RS_adj[i,j*7+5] = min(N_Wrk_RC2_RS_org[i,j*7+5], max_worker_fl[i,j])
            N_Wrk_RC3_RS_adj[i,j*7+5] = min(N_Wrk_RC3_RS_org[i,j*7+5], max_worker_fl[i,j])
            N_Wrk_RC4_RS_adj[i,j*7+5] = min(N_Wrk_RC4_RS_org[i,j*7+5], max_worker_fl[i,j])
            
            N_Wrk_RC2_RS_adj[i,j*7+6] = min(N_Wrk_RC2_RS_org[i,j*7+6], max_worker_fl[i,j])
            N_Wrk_RC3_RS_adj[i,j*7+6] = min(N_Wrk_RC3_RS_org[i,j*7+6], max_worker_fl[i,j])
            N_Wrk_RC4_RS_adj[i,j*7+6] = min(N_Wrk_RC4_RS_org[i,j*7+6], max_worker_fl[i,j])
            
            summ2 =  N_Wrk_RC2_RS_org[i,j*7+1]+N_Wrk_RC2_RS_org[i,j*7+3]+N_Wrk_RC2_RS_org[i,j*7+4]
            summ3 =  N_Wrk_RC3_RS_org[i,j*7+1]+N_Wrk_RC3_RS_org[i,j*7+3]+N_Wrk_RC3_RS_org[i,j*7+4]
            summ4 =  N_Wrk_RC4_RS_org[i,j*7+1]+N_Wrk_RC4_RS_org[i,j*7+3]+N_Wrk_RC4_RS_org[i,j*7+4]
            
            if (summ2 > max_worker_fl[i,j]):
                N_Wrk_RC2_RS_adj[i,j*7+1] = np.round((N_Wrk_RC2_RS_org[i,j*7+1])*(max_worker_fl[i,j]/summ2))
                N_Wrk_RC2_RS_adj[i,j*7+3] = np.round((N_Wrk_RC2_RS_org[i,j*7+3])*(max_worker_fl[i,j]/summ2))
                N_Wrk_RC2_RS_adj[i,j*7+4] = np.round((N_Wrk_RC2_RS_org[i,j*7+4])*(max_worker_fl[i,j]/summ2))
            if (summ3 > max_worker_fl[i,j]):
                N_Wrk_RC3_RS_adj[i,j*7+1] = np.round((N_Wrk_RC3_RS_org[i,j*7+1])*(max_worker_fl[i,j]/summ3))
                N_Wrk_RC3_RS_adj[i,j*7+3] = np.round((N_Wrk_RC3_RS_org[i,j*7+3])*(max_worker_fl[i,j]/summ3))
                N_Wrk_RC3_RS_adj[i,j*7+4] = np.round((N_Wrk_RC3_RS_org[i,j*7+4])*(max_worker_fl[i,j]/summ3))
            if (summ4 > max_worker_fl[i,j]):
                N_Wrk_RC4_RS_adj[i,j*7+1] = np.round((N_Wrk_RC4_RS_org[i,j*7+1])*(max_worker_fl[i,j]/summ4))
                N_Wrk_RC4_RS_adj[i,j*7+3] = np.round((N_Wrk_RC4_RS_org[i,j*7+3])*(max_worker_fl[i,j]/summ4))
                N_Wrk_RC4_RS_adj[i,j*7+4] = np.round((N_Wrk_RC4_RS_org[i,j*7+4])*(max_worker_fl[i,j]/summ4))
    
    
    #worker limit per repair sequence based on the numbers defined in Paul et al. (2018) table 1
    N_Wrk_RC2_RS_adj2 = N_Wrk_RC2_RS_adj
    N_Wrk_RC3_RS_adj2 = N_Wrk_RC3_RS_adj
    N_Wrk_RC4_RS_adj2 = N_Wrk_RC4_RS_adj
    for i in range(len(indx_repairable)):
        x=0
        for j in range(len(rep_phases)):
            for k in range(len(RS)):
                summ = sum(N_Wrk_RC2_RS_adj[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)]) 
                if (summ > max_worker_seq[k]):
                    N_Wrk_RC2_RS_adj2[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)] = np.round(N_Wrk_RC2_RS_adj[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)] * (max_worker_seq[k]/summ))
            x = rep_phases[j]*7+x
            
    for i in range(len(indx_repairable)):
        x=0
        for j in range(len(rep_phases)):
            for k in range(len(RS)):
                summ = sum(N_Wrk_RC3_RS_adj[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)]) 
                if (summ > max_worker_seq[k]):
                    N_Wrk_RC3_RS_adj2[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)] = np.round(N_Wrk_RC3_RS_adj[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)] * (max_worker_seq[k]/summ))
            x = rep_phases[j]*7+x
            
    for i in range(len(indx_repairable)):
        x=0
        for j in range(len(rep_phases)):
            for k in range(len(RS)):
                summ = sum(N_Wrk_RC4_RS_adj[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)]) 
                if (summ > max_worker_seq[k]):
                    N_Wrk_RC4_RS_adj2[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)] = np.round(N_Wrk_RC4_RS_adj[i,k+x:len(RS)*rep_phases[j]+k+x:len(RS)] * (max_worker_seq[k]/summ))
            x = rep_phases[j]*7+x
    
    #worker limit based on the total floor area defined in the REDi guidelines
    bld_area_ft2 = sum(fl_area)
    max_worker_tot1 = max(2.5*0.0001*bld_area_ft2 + 10, 20)
    max_worker_tot = min(max_worker_tot1, 260)
    
    N_Wrk_RC2_RS_adj3 = N_Wrk_RC2_RS_adj2
    N_Wrk_RC3_RS_adj3 = N_Wrk_RC3_RS_adj2
    N_Wrk_RC4_RS_adj3 = N_Wrk_RC4_RS_adj2
    for i in range(len(indx_repairable)):
        y=0
        for j in range(len(rep_phases)):
            summ = sum(N_Wrk_RC2_RS_adj2[i,y:y+rep_phases[j]*7])
            if (summ > max_worker_tot):
                N_Wrk_RC2_RS_adj3[i,y:y+rep_phases[j]*7] = np.round(N_Wrk_RC2_RS_adj2[i,y:y+rep_phases[j]*7] * max_worker_tot/summ)
            y = y+rep_phases[j]*7
            
    for i in range(len(indx_repairable)):
        y=0
        for j in range(len(rep_phases)):
            summ = sum(N_Wrk_RC3_RS_adj2[i,y:y+rep_phases[j]*7])
            if (summ > max_worker_tot):
                N_Wrk_RC3_RS_adj3[i,y:y+rep_phases[j]*7] = np.round(N_Wrk_RC3_RS_adj2[i,y:y+rep_phases[j]*7] * max_worker_tot/summ)
            y = y+rep_phases[j]*7
    
    for i in range(len(indx_repairable)):
        y=0
        for j in range(len(rep_phases)):
            summ = sum(N_Wrk_RC4_RS_adj2[i,y:y+rep_phases[j]*7])
            if (summ > max_worker_tot):
                N_Wrk_RC4_RS_adj3[i,y:y+rep_phases[j]*7] = np.round(N_Wrk_RC4_RS_adj2[i,y:y+rep_phases[j]*7] * max_worker_tot/summ)
            y = y+rep_phases[j]*7
    
    
    #elevator components adjustment to ensure the repair time is divided equally across stories 
    elev_array = np.maximum.reduce(N_Wrk_RC2_RS_adj3[:,5:-1:7].T)
    n_wrk_elev = max(RC_table[RC_table["Repair Sequence"]==6]["Workers per damaged Qty"])
    elev_array = np.minimum(elev_array, np.ones((1,len(indx_repairable)))*n_elev*n_wrk_elev)
    indx_elev_array = np.argmax(np.sum((N_Wrk_RC2_RS_adj3[:,5:-1:7]), axis=0))
    N_Wrk_RC2_RS_adj3[:,5+7*indx_elev_array]=elev_array
    if len(np.nonzero(np.sum((N_Wrk_RC2_RS_adj3[:,5:-1:7]), axis=0)))>1:
        sys.exit("Elevator fragility should be assigned to one floor only")
    RT_RC2_RS_adj4 = np.divide(np.squeeze(RT_RC2_RS), N_Wrk_RC2_RS_adj3, out=np.zeros_like(np.squeeze(RT_RC2_RS)), where=N_Wrk_RC2_RS_adj3!=0) 
    RT_RC2_RS_days = RT_RC2_RS_adj4
    RT_elev_adj = np.maximum.reduce(RT_RC2_RS_days[:,5:-1:7].T)/len(story) #find the max elevator repair time across stories
    RT_RC2_RS_days[:,np.arange(5,len(RT_RC2_RS_days.T),7)] = np.tile(RT_elev_adj,(len(story),1)).T
    
    RT_RC3_RS_days = np.divide(np.squeeze(RT_RC3_RS), N_Wrk_RC3_RS_adj3, out=np.zeros_like(np.squeeze(RT_RC3_RS)), where=N_Wrk_RC3_RS_adj3!=0) 
    RT_RC4_RS_days = np.divide(np.squeeze(RT_RC4_RS), N_Wrk_RC4_RS_adj3, out=np.zeros_like(np.squeeze(RT_RC4_RS)), where=N_Wrk_RC4_RS_adj3!=0) 
    
    
    RT_RC2_RS_days_mat=np.vstack((header_, (RT_RC2_RS_days)))
    pd.DataFrame(np.c_[col.T, RT_RC2_RS_days_mat]).to_csv('RT_RSeq_FR.csv', header=False, index=False)
    RT_RC3_RS_days_mat=np.vstack((header_, (RT_RC3_RS_days)))
    pd.DataFrame(np.c_[col.T, RT_RC3_RS_days_mat]).to_csv('RT_RSeq_RO.csv', header=False, index=False)
    RT_RC4_RS_days_mat=np.vstack((header_, (RT_RC4_RS_days)))
    pd.DataFrame(np.c_[col.T, RT_RC4_RS_days_mat]).to_csv('RT_RSeq_SiP.csv', header=False, index=False)
       
    print('Repair time calculation is completed')
    return RT_RC2_RS_days, RT_RC3_RS_days, RT_RC4_RS_days, N_DMG_RC5, N_DMG_RC3_RS       
                