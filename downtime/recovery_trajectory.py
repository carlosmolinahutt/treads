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
This module calculates recovery trajectories

"""

def RecTr_calc(story, rep_phases, Qt_facade, reconst_time, indx_repairable, indx_irreparable, indx_collapse, IF_output, RT_RC2_RS_days, RT_RC3_RS_days, RT_RC4_RS_days, RCmax_RS, N_DMG_RC3_RS): 
    import numpy as np
    import pandas as pd
    np.seterr(divide='ignore', invalid='ignore')
    
    IF_inspection = IF_output[0]
    IF_eng = IF_output[1]
    IF_permit = IF_output[2]
    IF_finance = IF_output[3]
    IF_cm_rs1 = IF_output[4]
    IF_cm_rs2 = IF_output[5]
    IF_cm_rs3 = IF_output[6]
    IF_cm_rs4 = IF_output[7]
    IF_cm_rs5 = IF_output[8]
    IF_cm_rs6 = IF_output[9]
    IF_cm_rs7 = IF_output[10]
    IF_stab = IF_output[11]
    IF_reconst = IF_output[12]
    
    story_bm = rep_phases[len(rep_phases)-1] #number of basement stories
    story_gr = sum(rep_phases) - story_bm #number of above grade stories
    
    usability_repairable = np.append([1,0],np.linspace(0, 1, story_gr+1))
    usability_irreparable_1 = np.append([1], np.zeros(len(usability_repairable)-2))
    usability_irreparable = np.append(usability_irreparable_1, [1])
    usability = np.vstack((usability_repairable.T,usability_irreparable.T))
    
    #downtime for irreparable and collapse scenarios   
    DT_final_irreparable = np.zeros((len(indx_irreparable)+len(indx_collapse), len(usability_repairable)))
    DT_irr_tot = IF_reconst + np.ones(len(indx_irreparable)+len(indx_collapse))*reconst_time*sum(rep_phases)
    DT_final_irreparable[:,-1] = DT_irr_tot
    DT_final_irreparable[:,-2] = DT_irr_tot
    
    
##downtime to functional recovery
    RT_RC_RS_days = RT_RC2_RS_days
    
    #downtime calculation for repair path
    DT_A1 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A2 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A4 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A5 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_B = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_C = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_D = np.zeros((len(indx_repairable), len(usability_repairable)))
    
    max_RTbm_2_4_5 = np.maximum.reduce([RT_RC_RS_days[:,-6-7*(story_bm-1):-5:7], RT_RC_RS_days[:,-4-7*(story_bm-1):-3:7], RT_RC_RS_days[:,-3-7*(story_bm-1):-2:7]])
    DT_A1[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs1, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-7-7*(story_bm-1):-6:7]+max_RTbm_2_4_5, axis=1)
    DT_A2[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs2, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-6-7*(story_bm-1):-5:7], axis=1)
    DT_A4[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs4, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-4-7*(story_bm-1):-3:7], axis=1)
    DT_A5[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs5, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-3-7*(story_bm-1):-2:7], axis=1)
    DT_A = np.maximum.reduce([DT_A1, DT_A2, DT_A4, DT_A5])
    DT_B[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs3, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-5-7*(story_bm-1):-4:7], axis=1)
    DT_C[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs6, IF_eng+IF_permit]) + sum(RT_RC_RS_days[:,-2-7*(story_bm-1):-1:7].T) #2 workers per elevator for the entire bld
    DT_D[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs7, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-1-7*(story_bm-1):len(RT_RC_RS_days.T):7], axis=1)
    
    #dowtime calculation for each rapair phase assuming rapair is peformed every 1, 2, or 3 stories
    RT_RS1 = np.zeros((len(indx_repairable), story_gr))
    RT_A1 = np.zeros((len(indx_repairable), story_gr))
    RT_A2 = np.zeros((len(indx_repairable), story_gr))
    RT_A4 = np.zeros((len(indx_repairable), story_gr))
    RT_A5 = np.zeros((len(indx_repairable), story_gr))
    RT_B = np.zeros((len(indx_repairable), story_gr))
    RT_C = np.zeros((len(indx_repairable), story_gr))
    RT_D = np.zeros((len(indx_repairable), story_gr))
    for i in range(len(indx_repairable)):
        n=0
        m=0
        for j in range(len(rep_phases)-1):
            if rep_phases[j]==1:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])

                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                                    
                m=m+1
                n=n+rep_phases[j]*7

            elif rep_phases[j]==2:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])

                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                
                RT_RS1[i,m+1] = min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n], max_RT_RS1-RT_RC_RS_days[i,0+n]),max_RT_RS1)
                RT_A1[i,m+1] = max(min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n]+max_RTgr_2_4_5[1], max_RT_A1-RT_RC_RS_days[i,0+n]-max_RTgr_2_4_5[0]),max_RT_A1), RT_A1[i,m])
                RT_A2[i,m+1] = min(RT_A2[i,m] + min(RT_RC_RS_days[i,8+n], max_RT_A2-RT_RC_RS_days[i,1+n]),max_RT_A2)
                RT_A4[i,m+1] = min(RT_A4[i,m] + min(RT_RC_RS_days[i,10+n], max_RT_A4-RT_RC_RS_days[i,3+n]),max_RT_A4)
                RT_A5[i,m+1] = min(RT_A5[i,m] + min(RT_RC_RS_days[i,11+n], max_RT_A5-RT_RC_RS_days[i,4+n]),max_RT_A5)
                RT_B[i,m+1] = min(RT_B[i,m] + min(RT_RC_RS_days[i,9+n], max_RT_B-RT_RC_RS_days[i,2+n]),max_RT_B)
                RT_C[i,m+1] = RT_C[i,m] + RT_RC_RS_days[i,12+n]
                RT_D[i,m+1] = min(RT_D[i,m] + min(RT_RC_RS_days[i,13+n], max_RT_D-RT_RC_RS_days[i,6+n]),max_RT_D)
                     
                m=m+2
                n=n+rep_phases[j]*7
  
            elif rep_phases[j]==3:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])
                
                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                
                RT_RS1[i,m+1] = min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n], max_RT_RS1-RT_RC_RS_days[i,0+n]),max_RT_RS1)
                RT_A1[i,m+1] = max(min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n]+max_RTgr_2_4_5[1], max_RT_A1-RT_RC_RS_days[i,0+n]-max_RTgr_2_4_5[0]),max_RT_A1), RT_A1[i,m])
                RT_A2[i,m+1] = min(RT_A2[i,m] + min(RT_RC_RS_days[i,8+n], max_RT_A2-RT_RC_RS_days[i,1+n]),max_RT_A2)
                RT_A4[i,m+1] = min(RT_A4[i,m] + min(RT_RC_RS_days[i,10+n], max_RT_A4-RT_RC_RS_days[i,3+n]),max_RT_A4)
                RT_A5[i,m+1] = min(RT_A5[i,m] + min(RT_RC_RS_days[i,11+n], max_RT_A5-RT_RC_RS_days[i,4+n]),max_RT_A5)
                RT_B[i,m+1] = min(RT_B[i,m] + min(RT_RC_RS_days[i,9+n], max_RT_B-RT_RC_RS_days[i,2+n]),max_RT_B)
                RT_C[i,m+1] = RT_C[i,m] + RT_RC_RS_days[i,12+n]
                RT_D[i,m+1] = min(RT_D[i,m] + min(RT_RC_RS_days[i,13+n], max_RT_D-RT_RC_RS_days[i,6+n]),max_RT_D)
                
                RT_RS1[i,m+2] = min(RT_RS1[i,m+1] + min(RT_RC_RS_days[i,14+n], max_RT_RS1-RT_RC_RS_days[i,7+n]),max_RT_RS1)
                RT_A1[i,m+2] = max(min(RT_RS1[i,m+1] + min(RT_RC_RS_days[i,14+n]+max_RTgr_2_4_5[2], max_RT_A1-RT_RC_RS_days[i,7+n]-max_RTgr_2_4_5[1]),max_RT_A1), RT_A1[i,m+1])
                RT_A2[i,m+2] = min(RT_A2[i,m+1] + min(RT_RC_RS_days[i,15+n], max_RT_A2-RT_RC_RS_days[i,8+n]),max_RT_A2)
                RT_A4[i,m+2] = min(RT_A4[i,m+1] + min(RT_RC_RS_days[i,17+n], max_RT_A4-RT_RC_RS_days[i,10+n]),max_RT_A4)
                RT_A5[i,m+2] = min(RT_A5[i,m+1] + min(RT_RC_RS_days[i,18+n], max_RT_A5-RT_RC_RS_days[i,11+n]),max_RT_A5)
                RT_B[i,m+2] = min(RT_B[i,m+1] + min(RT_RC_RS_days[i,16+n], max_RT_B-RT_RC_RS_days[i,9+n]),max_RT_B)
                RT_C[i,m+2] = RT_C[i,m+1] + RT_RC_RS_days[i,19+n]
                RT_D[i,m+2] = min(RT_D[i,m+1] + min(RT_RC_RS_days[i,20+n], max_RT_D-RT_RC_RS_days[i,13+n]),max_RT_D)
                
                m=m+3
                n=n+rep_phases[j]*7

    for i in range(len(indx_repairable)):
        for j in range(len(rep_phases)-1):
            
            RT_RS1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_RS1[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_RS1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A1[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A2[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A2[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A2[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A4[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A4[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A4[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A5[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A5[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A5[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_B[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_B[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_B[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_C[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_C[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_C[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_D[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_D[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_D[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
    RT_A = np.maximum.reduce([RT_A1, RT_A2, RT_A4, RT_A5])
    RT_RS2 = RT_A2
    RT_RS4 = RT_A4
    RT_RS5 = RT_A5
    RT_RS3 = RT_B
    RT_RS6 = RT_C
    RT_RS7 = RT_D
    
    #generate repair time stepping fuctions for functional recovery
    with pd.ExcelWriter('RT_stepfunc_FR.xlsx') as writer:  
        pd.DataFrame(RT_RS1).to_excel(writer, sheet_name='RSeq1', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS2).to_excel(writer, sheet_name='RSeq2', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS3).to_excel(writer, sheet_name='RSeq3', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS4).to_excel(writer, sheet_name='RSeq4', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS5).to_excel(writer, sheet_name='RSeq5', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS6).to_excel(writer, sheet_name='RSeq6', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS7).to_excel(writer, sheet_name='RSeq7', header=story[0:story_gr], index_label='#Num')
    
    a=np.zeros(len(indx_repairable))
    b=np.zeros(len(indx_repairable))
    c=np.zeros(len(indx_repairable))
    d=np.zeros(len(indx_repairable))
    for i in range(len(indx_repairable)):
        if max(RT_A[i,:]) != 0:
            a[i]=1
        if max(RT_B[i,:]) != 0:
            b[i]=1
        if max(RT_C[i,:]) != 0:
            c[i]=1
        if max(RT_D[i,:]) != 0:
            d[i]=1
    aa=np.tile(a,(len(usability_repairable),1)).T
    bb=np.tile(b,(len(usability_repairable),1)).T
    cc=np.tile(c,(len(usability_repairable),1)).T
    dd=np.tile(d,(len(usability_repairable),1)).T
    
    RT_A = RT_A + np.tile(DT_A[:,2],(story_gr,1)).T
    RT_B = RT_B + np.tile(DT_B[:,2],(story_gr,1)).T
    RT_C = RT_C + np.tile(DT_C[:,2],(story_gr,1)).T
    RT_D = RT_D + np.tile(DT_D[:,2],(story_gr,1)).T
    
    DT_A[:,3:]=RT_A
    DT_B[:,3:]=RT_B
    DT_C[:,3:]=RT_C
    DT_D[:,3:]=RT_D
    DT_final_repairable = np.maximum.reduce([DT_A*aa, DT_B*bb, DT_C*cc, DT_D*dd])

    #utility time consideration for downtime to functional recovery calculation
    DT_utility = np.zeros((len(DT_final_repairable),len(DT_final_repairable.T)))
    k = np.maximum(np.random.lognormal(np.log(10),1,len(DT_final_repairable)),np.random.lognormal(np.log(4),.55,len(DT_final_repairable)),np.random.lognormal(np.log(3),1.2,len(DT_final_repairable)))
    for i in range(len(DT_final_repairable)):
        DT_utility[i,2:]=k[i]
    
    DT_A = DT_A*aa
    DT_B = DT_B*bb
    DT_C = DT_C*cc
    DT_D = DT_D*dd
    
    mat_adj_A = np.where(np.divide(DT_A,np.transpose(np.repeat([DT_A[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_B = np.where(np.divide(DT_B,np.transpose(np.repeat([DT_B[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_C = np.where(np.divide(DT_C,np.transpose(np.repeat([DT_C[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_D = np.where(np.divide(DT_D,np.transpose(np.repeat([DT_D[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj = np.where(np.divide(DT_final_repairable,np.transpose(np.repeat([DT_final_repairable[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_A2 = np.concatenate((np.ones((len(DT_A),3)), mat_adj_A), axis=1)
    mat_adj_B2 = np.concatenate((np.ones((len(DT_B),3)), mat_adj_B), axis=1)
    mat_adj_C2 = np.concatenate((np.ones((len(DT_C),3)), mat_adj_C), axis=1)
    mat_adj_D2 = np.concatenate((np.ones((len(DT_D),3)), mat_adj_D), axis=1)
    mat_adj2 = np.concatenate((np.ones((len(DT_final_repairable),3)), mat_adj), axis=1)
    DT_A = DT_A * mat_adj_A2
    DT_B = DT_B * mat_adj_B2
    DT_C = DT_C * mat_adj_C2
    DT_D = DT_D * mat_adj_D2
    DT_final_repairable = DT_final_repairable * mat_adj2
    
    # adjustment for downtime stepping functions to ensure the usability can be restored if no repair is required in lower stories
    for i in range(len(DT_final_repairable)):
        if DT_final_repairable[i,3]==0:
            indx = np.asarray(np.where(DT_final_repairable[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_final_repairable[i,indx_max]=DT_final_repairable[i,2]
            DT_final_repairable[i,2]=0
    #for i in range(len(DT_A)):
        if DT_A[i,3]==0:
            indx = np.asarray(np.where(DT_A[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_A[i,indx_max]=DT_A[i,2]
            DT_A[i,2]=0
    #for i in range(len(DT_B)):
        if DT_B[i,3]==0:
            indx = np.asarray(np.where(DT_B[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_B[i,indx_max]=DT_B[i,2]
            DT_B[i,2]=0
    #for i in range(len(DT_C)):
        if DT_C[i,3]==0:
            indx = np.asarray(np.where(DT_C[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_C[i,indx_max]=DT_C[i,2]
            DT_C[i,2]=0
    #for i in range(len(DT_D)):
        if DT_D[i,3]==0:
            indx = np.asarray(np.where(DT_D[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_D[i,indx_max]=DT_D[i,2]
            DT_D[i,2]=0   
    
    #ensure that the downtime is not less than the inspection time in each repair phase
    for i in range(len(DT_final_repairable)):
        if DT_final_repairable[i,2]==0:
            DT_final_repairable[i,2:][DT_final_repairable[i,2:]==0]=IF_inspection[i]      
    #for i in range(len(DT_A)):
        if DT_A[i,2]==0:
            DT_A[i,2:][DT_A[i,2:]==0]=IF_inspection[i] 
    #for i in range(len(DT_B)):
        if DT_B[i,2]==0:
            DT_B[i,2:][DT_B[i,2:]==0]=IF_inspection[i] 
    #for i in range(len(DT_C)):
        if DT_C[i,2]==0:
            DT_C[i,2:][DT_C[i,2:]==0]=IF_inspection[i] 
    #for i in range(len(DT_D)):
        if DT_D[i,2]==0:
            DT_D[i,2:][DT_D[i,2:]==0]=IF_inspection[i] 
    
    #compare the utility repair time vs the total downtime only if FR is triggered
    for i in range(len(DT_final_repairable)):
        if (DT_final_repairable[i,-1]<DT_utility[i,-1]) and (DT_final_repairable[i,-1]!=DT_final_repairable[i,2]):
            DT_final_repairable[i,:]=DT_utility[i,:]
    
    row_id_rep=[]
    for i in range(len(indx_repairable)):
        row_id_rep.append('real_'+str(i)+'_repairable')
    row_id_irr=[]
    for i in range(len(indx_irreparable) + len(indx_collapse)):
        row_id_irr.append('real_'+str(i)+'_irreparable')
    row_id_all = row_id_rep + row_id_irr
    row_id_all=['usability_repairable','usability_irreparable']+row_id_all
    
    #generate downtime to functional recovery stepping functions
    DT_final_RC2 = np.concatenate((DT_final_repairable, DT_final_irreparable), axis=0)
    DT_final_RC2_use = np.concatenate((usability, DT_final_RC2), axis=0)
    DT_final_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_final_RC2_use]
    pd.DataFrame(DT_final_RC2_use).to_csv('DT_stepfunc_FR.csv', header=None, index=None)
    
    DT_A_RC2 = np.concatenate((DT_A, DT_final_irreparable), axis=0)
    DT_B_RC2 = np.concatenate((DT_B, DT_final_irreparable), axis=0)
    DT_C_RC2 = np.concatenate((DT_C, DT_final_irreparable), axis=0)
    DT_D_RC2 = np.concatenate((DT_D, DT_final_irreparable), axis=0)
    DT_utility = np.concatenate((DT_utility, DT_final_irreparable), axis=0)
    
    DT_A_RC2_use = np.concatenate((usability, DT_A_RC2), axis=0)
    DT_B_RC2_use = np.concatenate((usability, DT_B_RC2), axis=0)
    DT_C_RC2_use = np.concatenate((usability, DT_C_RC2), axis=0)
    DT_D_RC2_use = np.concatenate((usability, DT_D_RC2), axis=0)
    DT_utility_use = np.concatenate((usability, DT_utility), axis=0)
    
    DT_A_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_A_RC2_use]
    DT_B_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_B_RC2_use]
    DT_C_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_C_RC2_use]
    DT_D_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_D_RC2_use]
    DT_utility_use = np.c_[np.squeeze(row_id_all).T, DT_utility_use]
    
    with pd.ExcelWriter('DT_path_FR.xlsx', options= {'strings_to_numbers': True}) as writer:  
        pd.DataFrame(DT_A_RC2_use).to_excel(writer, sheet_name='A', header=None, index_label=None, index=False)
        pd.DataFrame(DT_B_RC2_use).to_excel(writer, sheet_name='B', header=None, index_label=None, index=False)
        pd.DataFrame(DT_C_RC2_use).to_excel(writer, sheet_name='C', header=None, index_label=None, index=False)
        pd.DataFrame(DT_D_RC2_use).to_excel(writer, sheet_name='D', header=None, index_label=None, index=False)
        pd.DataFrame(DT_utility_use).to_excel(writer, sheet_name='utility', header=None, index_label=None, index=False)
        
    
##downtime to reoccupancy
    RT_RC_RS_days = RT_RC3_RS_days
    
    #downtime calculation for repair path
    DT_A1 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A2 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A4 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A5 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_B = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_C = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_D = np.zeros((len(indx_repairable), len(usability_repairable)))
    
    max_RTbm_2_4_5 = np.maximum.reduce([RT_RC_RS_days[:,-6-7*(story_bm-1):-5:7], RT_RC_RS_days[:,-4-7*(story_bm-1):-3:7], RT_RC_RS_days[:,-3-7*(story_bm-1):-2:7]])
    DT_A1[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs1, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-7-7*(story_bm-1):-6:7]+max_RTbm_2_4_5, axis=1)
    DT_A2[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs2, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-6-7*(story_bm-1):-5:7], axis=1)
    DT_A4[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs4, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-4-7*(story_bm-1):-3:7], axis=1)
    DT_A5[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs5, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-3-7*(story_bm-1):-2:7], axis=1)
    DT_A = np.maximum.reduce([DT_A1, DT_A2, DT_A4, DT_A5])
    DT_B[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs3, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-5-7*(story_bm-1):-4:7], axis=1)
    DT_C[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs6, IF_eng+IF_permit]) + sum(RT_RC_RS_days[:,-2-7*(story_bm-1):-1:7].T) #2 workers per elevator for the entire bld
    DT_D[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs7, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-1-7*(story_bm-1):len(RT_RC_RS_days.T):7], axis=1)
    
    RT_RS1 = np.zeros((len(indx_repairable), story_gr))
    RT_A1 = np.zeros((len(indx_repairable), story_gr))
    RT_A2 = np.zeros((len(indx_repairable), story_gr))
    RT_A4 = np.zeros((len(indx_repairable), story_gr))
    RT_A5 = np.zeros((len(indx_repairable), story_gr))
    RT_B = np.zeros((len(indx_repairable), story_gr))
    RT_C = np.zeros((len(indx_repairable), story_gr))
    RT_D = np.zeros((len(indx_repairable), story_gr))
    for i in range(len(indx_repairable)):
        n=0
        m=0
        for j in range(len(rep_phases)-1):
            if rep_phases[j]==1:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])

                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                                    
                m=m+1
                n=n+rep_phases[j]*7

            elif rep_phases[j]==2:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])

                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                
                RT_RS1[i,m+1] = min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n], max_RT_RS1-RT_RC_RS_days[i,0+n]),max_RT_RS1)
                RT_A1[i,m+1] = max(min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n]+max_RTgr_2_4_5[1], max_RT_A1-RT_RC_RS_days[i,0+n]-max_RTgr_2_4_5[0]),max_RT_A1), RT_A1[i,m])
                RT_A2[i,m+1] = min(RT_A2[i,m] + min(RT_RC_RS_days[i,8+n], max_RT_A2-RT_RC_RS_days[i,1+n]),max_RT_A2)
                RT_A4[i,m+1] = min(RT_A4[i,m] + min(RT_RC_RS_days[i,10+n], max_RT_A4-RT_RC_RS_days[i,3+n]),max_RT_A4)
                RT_A5[i,m+1] = min(RT_A5[i,m] + min(RT_RC_RS_days[i,11+n], max_RT_A5-RT_RC_RS_days[i,4+n]),max_RT_A5)
                RT_B[i,m+1] = min(RT_B[i,m] + min(RT_RC_RS_days[i,9+n], max_RT_B-RT_RC_RS_days[i,2+n]),max_RT_B)
                RT_C[i,m+1] = RT_C[i,m] + RT_RC_RS_days[i,12+n]
                RT_D[i,m+1] = min(RT_D[i,m] + min(RT_RC_RS_days[i,13+n], max_RT_D-RT_RC_RS_days[i,6+n]),max_RT_D)
                     
                m=m+2
                n=n+rep_phases[j]*7
  
            elif rep_phases[j]==3:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])
                
                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                
                RT_RS1[i,m+1] = min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n], max_RT_A1-RT_RC_RS_days[i,0+n]),max_RT_RS1)
                RT_A1[i,m+1] = max(min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n]+max_RTgr_2_4_5[1], max_RT_A1-RT_RC_RS_days[i,0+n]-max_RTgr_2_4_5[0]),max_RT_A1), RT_A1[i,m])
                RT_A2[i,m+1] = min(RT_A2[i,m] + min(RT_RC_RS_days[i,8+n], max_RT_A2-RT_RC_RS_days[i,1+n]),max_RT_A2)
                RT_A4[i,m+1] = min(RT_A4[i,m] + min(RT_RC_RS_days[i,10+n], max_RT_A4-RT_RC_RS_days[i,3+n]),max_RT_A4)
                RT_A5[i,m+1] = min(RT_A5[i,m] + min(RT_RC_RS_days[i,11+n], max_RT_A5-RT_RC_RS_days[i,4+n]),max_RT_A5)
                RT_B[i,m+1] = min(RT_B[i,m] + min(RT_RC_RS_days[i,9+n], max_RT_B-RT_RC_RS_days[i,2+n]),max_RT_B)
                RT_C[i,m+1] = RT_C[i,m] + RT_RC_RS_days[i,12+n]
                RT_D[i,m+1] = min(RT_D[i,m] + min(RT_RC_RS_days[i,13+n], max_RT_D-RT_RC_RS_days[i,6+n]),max_RT_D)
                
                RT_RS1[i,m+2] = min(RT_A1[i,m+1] + min(RT_RC_RS_days[i,14+n], max_RT_RS1-RT_RC_RS_days[i,7+n]),max_RT_RS1)
                RT_A1[i,m+2] = max(min(RT_RS1[i,m+1] + min(RT_RC_RS_days[i,14+n]+max_RTgr_2_4_5[2], max_RT_A1-RT_RC_RS_days[i,7+n]-max_RTgr_2_4_5[1]),max_RT_A1), RT_A1[i,m+1])
                RT_A2[i,m+2] = min(RT_A2[i,m+1] + min(RT_RC_RS_days[i,15+n], max_RT_A2-RT_RC_RS_days[i,8+n]),max_RT_A2)
                RT_A4[i,m+2] = min(RT_A4[i,m+1] + min(RT_RC_RS_days[i,17+n], max_RT_A4-RT_RC_RS_days[i,10+n]),max_RT_A4)
                RT_A5[i,m+2] = min(RT_A5[i,m+1] + min(RT_RC_RS_days[i,18+n], max_RT_A5-RT_RC_RS_days[i,11+n]),max_RT_A5)
                RT_B[i,m+2] = min(RT_B[i,m+1] + min(RT_RC_RS_days[i,16+n], max_RT_B-RT_RC_RS_days[i,9+n]),max_RT_B)
                RT_C[i,m+2] = RT_C[i,m+1] + RT_RC_RS_days[i,19+n]
                RT_D[i,m+2] = min(RT_D[i,m+1] + min(RT_RC_RS_days[i,20+n], max_RT_D-RT_RC_RS_days[i,13+n]),max_RT_D)
                
                m=m+3
                n=n+rep_phases[j]*7

    for i in range(len(indx_repairable)):
        for j in range(len(rep_phases)-1):
            
            RT_RS1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_RS1[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_RS1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A1[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A2[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A2[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A2[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A4[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A4[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A4[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A5[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A5[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A5[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_B[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_B[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_B[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_C[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_C[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_C[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_D[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_D[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_D[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
    RT_A = np.maximum.reduce([RT_A1, RT_A2, RT_A4, RT_A5])
    RT_RS2 = RT_A2
    RT_RS4 = RT_A4
    RT_RS5 = RT_A5
    RT_RS3 = RT_B
    RT_RS6 = RT_C
    RT_RS7 = RT_D
    
    with pd.ExcelWriter('RT_stepfunc_RO.xlsx') as writer:  
        pd.DataFrame(RT_RS1).to_excel(writer, sheet_name='RSeq1', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS2).to_excel(writer, sheet_name='RSeq2', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS3).to_excel(writer, sheet_name='RSeq3', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS4).to_excel(writer, sheet_name='RSeq4', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS5).to_excel(writer, sheet_name='RSeq5', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS6).to_excel(writer, sheet_name='RSeq6', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS7).to_excel(writer, sheet_name='RSeq7', header=story[0:story_gr], index_label='#Num')
    
    a=np.zeros(len(indx_repairable))
    b=np.zeros(len(indx_repairable))
    c=np.zeros(len(indx_repairable))
    d=np.zeros(len(indx_repairable))
    for i in range(len(indx_repairable)):
        if max(RT_A[i,:]) != 0:
            a[i]=1
        if max(RT_B[i,:]) != 0:
            b[i]=1
        if max(RT_C[i,:]) != 0:
            c[i]=1
        if max(RT_D[i,:]) != 0:
            d[i]=1
    aa=np.tile(a,(len(usability_repairable),1)).T
    bb=np.tile(b,(len(usability_repairable),1)).T
    cc=np.tile(c,(len(usability_repairable),1)).T
    dd=np.tile(d,(len(usability_repairable),1)).T
    
    RT_A = RT_A + np.tile(DT_A[:,2],(story_gr,1)).T
    RT_B = RT_B + np.tile(DT_B[:,2],(story_gr,1)).T
    RT_C = RT_C + np.tile(DT_C[:,2],(story_gr,1)).T
    RT_D = RT_D + np.tile(DT_D[:,2],(story_gr,1)).T
    
    DT_A[:,3:]=RT_A
    DT_B[:,3:]=RT_B
    DT_C[:,3:]=RT_C
    DT_D[:,3:]=RT_D
    
    DT_final_repairable = np.maximum.reduce([DT_A*aa, DT_B*bb, DT_D*dd])

    DT_A = DT_A*aa
    DT_B = DT_B*bb
    DT_D = DT_D*dd
    
    mat_adj_A = np.where(np.divide(DT_A,np.transpose(np.repeat([DT_A[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_B = np.where(np.divide(DT_B,np.transpose(np.repeat([DT_B[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_D = np.where(np.divide(DT_D,np.transpose(np.repeat([DT_D[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj = np.where(np.divide(DT_final_repairable,np.transpose(np.repeat([DT_final_repairable[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_A2 = np.concatenate((np.ones((len(DT_A),3)), mat_adj_A), axis=1)
    mat_adj_B2 = np.concatenate((np.ones((len(DT_B),3)), mat_adj_B), axis=1)
    mat_adj_D2 = np.concatenate((np.ones((len(DT_D),3)), mat_adj_D), axis=1)
    mat_adj2 = np.concatenate((np.ones((len(DT_final_repairable),3)), mat_adj), axis=1)
    DT_A = DT_A * mat_adj_A2
    DT_B = DT_B * mat_adj_B2
    DT_D = DT_D * mat_adj_D2
    DT_final_repairable = DT_final_repairable * mat_adj2
                
    for i in range(len(DT_final_repairable)):
        if DT_final_repairable[i,3]==0:
            indx = np.asarray(np.where(DT_final_repairable[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_final_repairable[i,indx_max]=DT_final_repairable[i,2]
            DT_final_repairable[i,2]=0    
    #for i in range(len(DT_A)):
        if DT_A[i,3]==0:
            indx = np.asarray(np.where(DT_A[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_A[i,indx_max]=DT_A[i,2]
            DT_A[i,2]=0
    #for i in range(len(DT_B)):
        if DT_B[i,3]==0:
            indx = np.asarray(np.where(DT_B[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_B[i,indx_max]=DT_B[i,2]
            DT_B[i,2]=0
    #for i in range(len(DT_D)):
        if DT_D[i,3]==0:
            indx = np.asarray(np.where(DT_D[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_D[i,indx_max]=DT_D[i,2]
            DT_D[i,2]=0
    
    for i in range(len(DT_final_repairable)):
        if DT_final_repairable[i,2]==0:
            DT_final_repairable[i,2:][DT_final_repairable[i,2:]==0]=IF_inspection[i]        
    #for i in range(len(DT_A)):
        if DT_A[i,2]==0:
            DT_A[i,2:][DT_A[i,2:]==0]=IF_inspection[i] 
    #for i in range(len(DT_B)):
        if DT_B[i,2]==0:
            DT_B[i,2:][DT_B[i,2:]==0]=IF_inspection[i] 
    #for i in range(len(DT_D)):
        if DT_D[i,2]==0:
            DT_D[i,2:][DT_D[i,2:]==0]=IF_inspection[i]     
    
    DT_final_RC3 = np.concatenate((DT_final_repairable, DT_final_irreparable), axis=0)
    DT_final_RC3_use = np.concatenate((usability, DT_final_RC3), axis=0)
    DT_final_RC3_use = np.c_[np.squeeze(row_id_all).T, DT_final_RC3_use]
    pd.DataFrame(DT_final_RC3_use).to_csv('DT_stepfunc_RO.csv', header=None, index=None)
    
    DT_A_RC2 = np.concatenate((DT_A, DT_final_irreparable), axis=0)
    DT_B_RC2 = np.concatenate((DT_B, DT_final_irreparable), axis=0)
    DT_D_RC2 = np.concatenate((DT_D, DT_final_irreparable), axis=0)
    
    DT_A_RC2_use = np.concatenate((usability, DT_A_RC2), axis=0)
    DT_B_RC2_use = np.concatenate((usability, DT_B_RC2), axis=0)
    DT_D_RC2_use = np.concatenate((usability, DT_D_RC2), axis=0)
    
    #elevator (path C) is not required to be repaired to achieve reoccupancy
    DT_A_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_A_RC2_use]
    DT_B_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_B_RC2_use]
    DT_D_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_D_RC2_use]
    
    with pd.ExcelWriter('DT_path_RO.xlsx', options= {'strings_to_numbers': True}) as writer:  
        pd.DataFrame(DT_A_RC2_use).to_excel(writer, sheet_name='A', header=None, index_label=None, index=False)
        pd.DataFrame(DT_B_RC2_use).to_excel(writer, sheet_name='B', header=None, index_label=None, index=False)
        pd.DataFrame(DT_D_RC2_use).to_excel(writer, sheet_name='D', header=None, index_label=None, index=False)
    
##downtime to shelter-in-place
    RT_RC_RS_days = RT_RC4_RS_days
    
    #repair paths:
    DT_A1 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A2 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A4 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_A5 = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_B = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_C = np.zeros((len(indx_repairable), len(usability_repairable)))
    DT_D = np.zeros((len(indx_repairable), len(usability_repairable)))
    
    max_RTbm_2_4_5 = np.maximum(RT_RC_RS_days[:,-6-7*(story_bm-1):-5:7], RT_RC_RS_days[:,-4-7*(story_bm-1):-3:7], RT_RC_RS_days[:,-3-7*(story_bm-1):-2:7])
    DT_A1[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs1, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-7-7*(story_bm-1):-6:7]+max_RTbm_2_4_5, axis=1)
    DT_A = DT_A1
    DT_B[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs3, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-5-7*(story_bm-1):-4:7], axis=1)
    DT_C[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs6, IF_eng+IF_permit]) + sum(RT_RC_RS_days[:,-2-7*(story_bm-1):-1:7].T) #2 workers per elevator for the entire bld
    DT_D[:,2] = IF_inspection + np.maximum.reduce([IF_stab, IF_finance, IF_cm_rs7, IF_eng+IF_permit]) + np.amax(RT_RC_RS_days[:,-1-7*(story_bm-1):len(RT_RC_RS_days.T):7], axis=1)
    
    RT_RS1 = np.zeros((len(indx_repairable), story_gr))
    RT_A1 = np.zeros((len(indx_repairable), story_gr))
    RT_A2 = np.zeros((len(indx_repairable), story_gr))
    RT_A4 = np.zeros((len(indx_repairable), story_gr))
    RT_A5 = np.zeros((len(indx_repairable), story_gr))
    RT_B = np.zeros((len(indx_repairable), story_gr))
    RT_C = np.zeros((len(indx_repairable), story_gr))
    RT_D = np.zeros((len(indx_repairable), story_gr))
    for i in range(len(indx_repairable)):
        n=0
        m=0
        for j in range(len(rep_phases)-1):
            if rep_phases[j]==1:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])

                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                                    
                m=m+1
                n=n+rep_phases[j]*7

            elif rep_phases[j]==2:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])

                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                
                RT_RS1[i,m+1] = min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n], max_RT_RS1-RT_RC_RS_days[i,0+n]),max_RT_RS1)
                RT_A1[i,m+1] = max(min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n]+max_RTgr_2_4_5[1], max_RT_A1-RT_RC_RS_days[i,0+n]-max_RTgr_2_4_5[0]),max_RT_A1), RT_A1[i,m])
                RT_A2[i,m+1] = min(RT_A2[i,m] + min(RT_RC_RS_days[i,8+n], max_RT_A2-RT_RC_RS_days[i,1+n]),max_RT_A2)
                RT_A4[i,m+1] = min(RT_A4[i,m] + min(RT_RC_RS_days[i,10+n], max_RT_A4-RT_RC_RS_days[i,3+n]),max_RT_A4)
                RT_A5[i,m+1] = min(RT_A5[i,m] + min(RT_RC_RS_days[i,11+n], max_RT_A5-RT_RC_RS_days[i,4+n]),max_RT_A5)
                RT_B[i,m+1] = min(RT_B[i,m] + min(RT_RC_RS_days[i,9+n], max_RT_B-RT_RC_RS_days[i,2+n]),max_RT_B)
                RT_C[i,m+1] = RT_C[i,m] + RT_RC_RS_days[i,12+n]
                RT_D[i,m+1] = min(RT_D[i,m] + min(RT_RC_RS_days[i,13+n], max_RT_D-RT_RC_RS_days[i,6+n]),max_RT_D)
                     
                m=m+2
                n=n+rep_phases[j]*7
  
            elif rep_phases[j]==3:
                max_RTgr_2_4_5 = np.maximum.reduce([RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7],RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7],RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7]])
                max_RT_A1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7]+max_RTgr_2_4_5)
                max_RT_RS1 = np.amax(RT_RC_RS_days[i,0+n:0+7*rep_phases[j]+n:7])
                max_RT_A2 = np.amax(RT_RC_RS_days[i,1+n:1+7*rep_phases[j]+n:7])
                max_RT_A4 = np.amax(RT_RC_RS_days[i,3+n:3+7*rep_phases[j]+n:7])
                max_RT_A5 = np.amax(RT_RC_RS_days[i,4+n:4+7*rep_phases[j]+n:7])
                max_RT_B = np.amax(RT_RC_RS_days[i,2+n:2+7*rep_phases[j]+n:7])
                max_RT_D = np.amax(RT_RC_RS_days[i,6+n:6+7*rep_phases[j]+n:7])
                
                RT_RS1[i,m] = min(RT_RC_RS_days[i,0+n],max_RT_RS1)
                RT_A1[i,m] = min(RT_RC_RS_days[i,0+n]+max_RTgr_2_4_5[0],max_RT_A1)
                RT_A2[i,m] = min(RT_RC_RS_days[i,1+n],max_RT_A2)
                RT_A4[i,m] = min(RT_RC_RS_days[i,3+n],max_RT_A4)
                RT_A5[i,m] = min(RT_RC_RS_days[i,4+n],max_RT_A5)
                RT_B[i,m] = min(RT_RC_RS_days[i,2+n],max_RT_B)
                RT_C[i,m] = RT_RC_RS_days[i,5+n]
                RT_D[i,m] = min(RT_RC_RS_days[i,6+n],max_RT_D)
                
                RT_RS1[i,m+1] = min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n], max_RT_A1-RT_RC_RS_days[i,0+n]),max_RT_RS1)
                RT_A1[i,m+1] = max(min(RT_RS1[i,m] + min(RT_RC_RS_days[i,7+n]+max_RTgr_2_4_5[1], max_RT_A1-RT_RC_RS_days[i,0+n]-max_RTgr_2_4_5[0]),max_RT_A1), RT_A1[i,m])
                RT_A2[i,m+1] = min(RT_A2[i,m] + min(RT_RC_RS_days[i,8+n], max_RT_A2-RT_RC_RS_days[i,1+n]),max_RT_A2)
                RT_A4[i,m+1] = min(RT_A4[i,m] + min(RT_RC_RS_days[i,10+n], max_RT_A4-RT_RC_RS_days[i,3+n]),max_RT_A4)
                RT_A5[i,m+1] = min(RT_A5[i,m] + min(RT_RC_RS_days[i,11+n], max_RT_A5-RT_RC_RS_days[i,4+n]),max_RT_A5)
                RT_B[i,m+1] = min(RT_B[i,m] + min(RT_RC_RS_days[i,9+n], max_RT_B-RT_RC_RS_days[i,2+n]),max_RT_B)
                RT_C[i,m+1] = RT_C[i,m] + RT_RC_RS_days[i,12+n]
                RT_D[i,m+1] = min(RT_D[i,m] + min(RT_RC_RS_days[i,13+n], max_RT_D-RT_RC_RS_days[i,6+n]),max_RT_D)
                
                RT_RS1[i,m+2] = min(RT_A1[i,m+1] + min(RT_RC_RS_days[i,14+n], max_RT_RS1-RT_RC_RS_days[i,7+n]),max_RT_RS1)
                RT_A1[i,m+2] = max(min(RT_RS1[i,m+1] + min(RT_RC_RS_days[i,14+n]+max_RTgr_2_4_5[2], max_RT_A1-RT_RC_RS_days[i,7+n]-max_RTgr_2_4_5[1]),max_RT_A1), RT_A1[i,m+1])
                RT_A2[i,m+2] = min(RT_A2[i,m+1] + min(RT_RC_RS_days[i,15+n], max_RT_A2-RT_RC_RS_days[i,8+n]),max_RT_A2)
                RT_A4[i,m+2] = min(RT_A4[i,m+1] + min(RT_RC_RS_days[i,17+n], max_RT_A4-RT_RC_RS_days[i,10+n]),max_RT_A4)
                RT_A5[i,m+2] = min(RT_A5[i,m+1] + min(RT_RC_RS_days[i,18+n], max_RT_A5-RT_RC_RS_days[i,11+n]),max_RT_A5)
                RT_B[i,m+2] = min(RT_B[i,m+1] + min(RT_RC_RS_days[i,16+n], max_RT_B-RT_RC_RS_days[i,9+n]),max_RT_B)
                RT_C[i,m+2] = RT_C[i,m+1] + RT_RC_RS_days[i,19+n]
                RT_D[i,m+2] = min(RT_D[i,m+1] + min(RT_RC_RS_days[i,20+n], max_RT_D-RT_RC_RS_days[i,13+n]),max_RT_D)
                
                m=m+3
                n=n+rep_phases[j]*7

    for i in range(len(indx_repairable)):
        for j in range(len(rep_phases)-1):
            
            RT_RS1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_RS1[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_RS1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A1[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A1[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A2[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A2[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A2[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A4[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A4[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A4[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_A5[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_A5[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_A5[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_B[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_B[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_B[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_C[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_C[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_C[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
            RT_D[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])] = np.amax(RT_D[i,sum(rep_phases[:j]):sum(rep_phases[:j+1])]) + RT_D[i,sum(rep_phases[:j+1]):sum(rep_phases[:j+2])]
    RT_A = RT_A1
    RT_RS2 = RT_A2
    RT_RS4 = RT_A4
    RT_RS5 = RT_A5
    RT_RS3 = RT_B
    RT_RS6 = RT_C
    RT_RS7 = RT_D
    
    with pd.ExcelWriter('RT_stepfunc_SiP.xlsx', options= {'strings_to_numbers': True}) as writer:  
        pd.DataFrame(RT_RS1).to_excel(writer, sheet_name='RSeq1', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS2).to_excel(writer, sheet_name='RSeq2', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS3).to_excel(writer, sheet_name='RSeq3', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS4).to_excel(writer, sheet_name='RSeq4', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS5).to_excel(writer, sheet_name='RSeq5', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS6).to_excel(writer, sheet_name='RSeq6', header=story[0:story_gr], index_label='#Num')
        pd.DataFrame(RT_RS7).to_excel(writer, sheet_name='RSeq7', header=story[0:story_gr], index_label='#Num')
    
    a=np.zeros(len(indx_repairable))
    b=np.zeros(len(indx_repairable))
    c=np.zeros(len(indx_repairable))
    d=np.zeros(len(indx_repairable))
    for i in range(len(indx_repairable)):
        if max(RT_A[i,:]) != 0:
            a[i]=1
        if max(RT_B[i,:]) != 0:
            b[i]=1
        if max(RT_C[i,:]) != 0:
            c[i]=1
        if max(RT_D[i,:]) != 0:
            d[i]=1
    aa=np.tile(a,(len(usability_repairable),1)).T
    bb=np.tile(b,(len(usability_repairable),1)).T
    cc=np.tile(c,(len(usability_repairable),1)).T
    dd=np.tile(d,(len(usability_repairable),1)).T
    
    RT_A = RT_A + np.tile(DT_A[:,2],(story_gr,1)).T
    RT_B = RT_B + np.tile(DT_B[:,2],(story_gr,1)).T
    RT_C = RT_C + np.tile(DT_C[:,2],(story_gr,1)).T
    RT_D = RT_D + np.tile(DT_D[:,2],(story_gr,1)).T
    
    DT_A[:,3:]=RT_A
    DT_B[:,3:]=RT_B
    DT_C[:,3:]=RT_C
    DT_D[:,3:]=RT_D
    
    DT_final_repairable = np.maximum.reduce([DT_A*aa, DT_D*dd])

    DT_A = DT_A*aa
    DT_D = DT_D*dd
    
    mat_adj_A = np.where(np.divide(DT_A,np.transpose(np.repeat([DT_A[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_D = np.where(np.divide(DT_D,np.transpose(np.repeat([DT_D[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj = np.where(np.divide(DT_final_repairable,np.transpose(np.repeat([DT_final_repairable[:,2]],story_gr+3,axis=0)))[:,3:]==1,0,1)
    mat_adj_A2 = np.concatenate((np.ones((len(DT_A),3)), mat_adj_A), axis=1)
    mat_adj_D2 = np.concatenate((np.ones((len(DT_D),3)), mat_adj_D), axis=1)
    mat_adj2 = np.concatenate((np.ones((len(DT_final_repairable),3)), mat_adj), axis=1)
    DT_A = DT_A * mat_adj_A2
    DT_D = DT_D * mat_adj_D2
    DT_final_repairable = DT_final_repairable * mat_adj2
                
    for i in range(len(DT_final_repairable)):
        if DT_final_repairable[i,3]==0:
            indx = np.asarray(np.where(DT_final_repairable[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_final_repairable[i,indx_max]=DT_final_repairable[i,2]
            DT_final_repairable[i,2]=0    
    #for i in range(len(DT_A)):
        if DT_A[i,3]==0:
            indx = np.asarray(np.where(DT_A[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_A[i,indx_max]=DT_A[i,2]
            DT_A[i,2]=0
    #for i in range(len(DT_D)):
        if DT_D[i,3]==0:
            indx = np.asarray(np.where(DT_D[i,:]==0))
            indx_max = max(np.squeeze(indx))
            DT_D[i,indx_max]=DT_D[i,2]
            DT_D[i,2]=0            

    for i in range(len(DT_final_repairable)):
        if DT_final_repairable[i,2]==0:
            DT_final_repairable[i,2:][DT_final_repairable[i,2:]==0]=IF_inspection[i]        
    #for i in range(len(DT_A)):
        if DT_A[i,2]==0:
            DT_A[i,2:][DT_A[i,2:]==0]=IF_inspection[i] 
    #for i in range(len(DT_D)):
        if DT_D[i,2]==0:
            DT_D[i,2:][DT_D[i,2:]==0]=IF_inspection[i]
            
    DT_final_RC4 = np.concatenate((DT_final_repairable, DT_final_irreparable), axis=0)
    DT_final_RC4_use = np.concatenate((usability, DT_final_RC4), axis=0)
    DT_final_RC4_use = np.c_[np.squeeze(row_id_all).T, DT_final_RC4_use]
    pd.DataFrame(DT_final_RC4_use).to_csv('DT_stepfunc_SiP.csv', header=None, index=None)
    
    DT_A_RC2 = np.concatenate((DT_A, DT_final_irreparable), axis=0)
    DT_D_RC2 = np.concatenate((DT_D, DT_final_irreparable), axis=0)
    
    DT_A_RC2_use = np.concatenate((usability, DT_A_RC2), axis=0)
    DT_D_RC2_use = np.concatenate((usability, DT_D_RC2), axis=0)
    
    #structural (path A) and staricase (path D) repairs are only required to achieve sheltering capacity
    DT_A_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_A_RC2_use]
    DT_D_RC2_use = np.c_[np.squeeze(row_id_all).T, DT_D_RC2_use]
    
    with pd.ExcelWriter('DT_path_SiP.xlsx', options= {'strings_to_numbers': True}) as writer:  
        pd.DataFrame(DT_A_RC2_use).to_excel(writer, sheet_name='A', header=None, index_label=None, index=False)
        pd.DataFrame(DT_D_RC2_use).to_excel(writer, sheet_name='D', header=None, index_label=None, index=False)
    
    #%% summary stats
    zero_DT_RC2 = np.percentile(DT_final_RC2[:,-1],0)
    tenth_DT_RC2 = np.percentile(DT_final_RC2[:,-1],10)
    med_DT_RC2 = np.median(DT_final_RC2[:,-1])
    mean_DT_RC2 = np.mean(DT_final_RC2[:,-1])
    ninety_DT_RC2 = np.percentile(DT_final_RC2[:,-1],90)
    hundred_DT_RC2 = np.percentile(DT_final_RC2[:,-1],100)
    
    zero_DT_RC3 = np.percentile(DT_final_RC3[:,-1],0)
    tenth_DT_RC3 = np.percentile(DT_final_RC3[:,-1],10)
    med_DT_RC3 = np.median(DT_final_RC3[:,-1])
    mean_DT_RC3 = np.mean(DT_final_RC3[:,-1])
    ninety_DT_RC3 = np.percentile(DT_final_RC3[:,-1],90)
    hundred_DT_RC3 = np.percentile(DT_final_RC3[:,-1],100)
    
    zero_DT_RC4 = np.percentile(DT_final_RC4[:,-1],0)
    tenth_DT_RC4 = np.percentile(DT_final_RC4[:,-1],10)
    med_DT_RC4 = np.median(DT_final_RC4[:,-1])
    mean_DT_RC4 = np.mean(DT_final_RC4[:,-1])
    ninety_DT_RC4 = np.percentile(DT_final_RC4[:,-1],90)
    hundred_DT_RC4 = np.percentile(DT_final_RC4[:,-1],100)
    
    DT_summary_RC2 =  [zero_DT_RC2, tenth_DT_RC2, med_DT_RC2, mean_DT_RC2, ninety_DT_RC2, hundred_DT_RC2]
    DT_summary_RC3 =  [zero_DT_RC3, tenth_DT_RC3, med_DT_RC3, mean_DT_RC3, ninety_DT_RC3, hundred_DT_RC3]
    DT_summary_RC4 =  [zero_DT_RC4, tenth_DT_RC4, med_DT_RC4, mean_DT_RC4, ninety_DT_RC4, hundred_DT_RC4]
    
    row_id = ["Minimum", "10th Percentile", "Median", "Mean", "90th Percentile", "Maximum"]
    col_id = ["Downtime","Functional Recovery", "Re-Occupancy", "Shelter-in-Place"]
    DT_summ_1 = np.c_[(row_id), np.array(DT_summary_RC2), np.array(DT_summary_RC3), np.array(DT_summary_RC4)]
    DT_summ = np.vstack((col_id,DT_summ_1))
    pd.DataFrame(DT_summ).to_csv('DT_summary.csv', header=False, index=False)
    
    #determine the proability of hindering a recovery state based on the max repair class & damaged facade components for the stability recovery state
    RCmax_repairable = np.max(np.squeeze(RCmax_RS), axis=1)
    N_DMG_RC3_RS_mat = np.squeeze(N_DMG_RC3_RS)[:,np.arange(2,len(np.squeeze(N_DMG_RC3_RS).transpose()),7)]
    N_DMG_RC3_RS3 = sum(N_DMG_RC3_RS_mat.transpose())
    #stability is hindered if damage facade components exceed 50% of the total
    for i in range(len(indx_repairable)):
        if N_DMG_RC3_RS3[i] > 0.5*Qt_facade:
            RCmax_repairable[i]=5
    RCmax = np.append(RCmax_repairable, 5*np.ones(len(indx_irreparable)+len(indx_collapse)))
    RS_stats = ['prob (RS not achieved)', len(RCmax[RCmax>=2])/len(RCmax), len(RCmax[RCmax>=3])/len(RCmax), len(RCmax[RCmax>=4])/len(RCmax)]
    index_label = ['Recovery State','Functional Recovery','Reoccupancy','Shelter-in-Place']
    pd.DataFrame(np.c_[index_label, RS_stats]).to_csv('RS_stats.csv', header=False, index=False)
    
    print('Downtime calculations for "Functional Recovery", "Re-Occupancy", and "Shelter-in-Place" recovery states are completed')
    return DT_final_RC2, DT_final_RC3, DT_final_RC4, DT_summary_RC2, DT_summary_RC3, DT_summary_RC4, RCmax