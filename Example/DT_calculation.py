# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Pouria Kourehpaz.
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
#
# Contributor(s):
# Pouria Kourehpaz

"""
This module generates downtime assessment outputs

"""

def run_treads(input_params_json, RCtable_input, IF_delays_input, DMG_input, DL_summary_input, DV_rec_time_input, output_path):
    
    import numpy as np
    import sys
    import treads.repair_class as repair_class
    import treads.main as main
    import treads.impeding_delays as impeding_delays
    import treads.recovery_trajectory as recovery_trajectory
    import json
    with open(input_params_json) as f:
        input_param = json.load(f)

    RC_results = repair_class.RC_calc(DMG_input, RCtable_input, DL_summary_input, output_path)
    RCmax_RS = np.squeeze(RC_results[0])
    RC = RC_results[1]
    indx_repairable = RC_results[2]
    indx_irreparable = RC_results[3]
    indx_collapse = RC_results[4]
    
    if sum(input_param['repair_phases'])!=len(RCmax_RS.T)/7:
        sys.exit("Sum of repair phases is not equal to the number of stories")
    
    #repair time calculation
    RT_results = main.RT_calc(input_param['elevator_quantity'], input_param['repair_phases'], input_param['floor_area'], RCtable_input, RC, DMG_input, DV_rec_time_input, DL_summary_input, input_param['max_number_workers'], output_path)
    
    N_DMG_RC5 = RT_results[3]
    N_DMG_RC3_RS = RT_results[4]
    
    #impeding factor calculation
    IF_results = impeding_delays.IF_calc(input_param['total_cost'], input_param['financing_coeff'], RCmax_RS, N_DMG_RC5, N_DMG_RC3_RS, DMG_input, DL_summary_input, input_param['stabilization']['limit_structural'], input_param['stabilization']['time_structural'], 
                                         input_param['stabilization']['limit_facade'], input_param['stabilization']['time_facade'], IF_delays_input, output_path)
    
    IF_output = IF_results
    RT_RC2_RS_days = RT_results[0]
    RT_RC3_RS_days = RT_results[1] 
    RT_RC4_RS_days = RT_results[2]
    story = RC_results[5]
    
    #recovery trajectory calculation
    DT_results = recovery_trajectory.RecTr_calc(story, input_param['repair_phases'], input_param['facade_quantity_eq'], input_param['time_story_reconstruction'], indx_repairable, indx_irreparable, indx_collapse, IF_output, RT_RC2_RS_days, RT_RC3_RS_days, RT_RC4_RS_days, RCmax_RS, N_DMG_RC3_RS, output_path)
    return DT_results, IF_results
