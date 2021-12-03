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
This module generates downtime assessment outputs

"""
def run_analysis(parameters, files):
    import numpy as np
    import sys
    import repair_class
    import main
    import impeding_delays
    import recovery_trajectory
    
    
    RC_results = repair_class.RC_calc(files[1], files[0], files[2])
    RCmax_RS = np.squeeze(RC_results[0])
    RC = RC_results[1]
    indx_repairable = RC_results[2]
    indx_irreparable = RC_results[3]
    indx_collapse = RC_results[4]
    
    if sum(parameters[1])!=len(RCmax_RS.T)/7:
        sys.exit("Sum of repair phases is not equal to the number of stories")
    
    #repair time calculation
    RT_results = main.RT_calc(parameters[0], parameters[1], parameters[2], files[0], RC, files[1], files[3], files[2], parameters[11])
    
    N_DMG_RC5 = RT_results[3]
    N_DMG_RC3_RS = RT_results[4]
    
    #impeding factor calculation
    IF_results = impeding_delays.IF_calc(parameters[3], parameters[4], RCmax_RS, N_DMG_RC5, N_DMG_RC3_RS, files[1], files[2], parameters[7], parameters[8], parameters[9], parameters[10], files[4])
    
    IF_output = IF_results
    RT_RC2_RS_days = RT_results[0]
    RT_RC3_RS_days = RT_results[1] 
    RT_RC4_RS_days = RT_results[2]
    story = RC_results[5]
    
    #recovery trajectory calculation
    DT_results = recovery_trajectory.RecTr_calc(story, parameters[1], parameters[5], parameters[6], indx_repairable, indx_irreparable, indx_collapse, IF_output, RT_RC2_RS_days, RT_RC3_RS_days, RT_RC4_RS_days, RCmax_RS, N_DMG_RC3_RS)
    return
