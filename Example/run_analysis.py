# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Pouria Kourehpaz.
#
#
# This code is developed to support the framework to estimate earthquake induced downtime 
# and recovery trajectory of residential buildings 
# proposed by Molina Hutt et al. (2021).
#
# The proposed framework can be found in the manuscript entitled 
# Molina Hutt, C., Vahanvaty, T., and Kourehpaz, P. (2021) 
# "An analytical framework to assess earthquake induced downtime and model recovery of buildings", Earthquake Spectra.
# 
#
# Contributor(s):
# Pouria Kourehpaz

"""
1. populate input parameters json file
2. define input files
3. define output path
4. execute this script to perform DT assessments

"""

import DT_calculation

input_parameters = 'input_parameters.json'
RCtable_input = 'Repair_Class_Table.csv'
IF_delays_input = 'IF_delays_input.csv'
DMG_input = 'DMG.csv' #pelicun output
DL_summary_input = 'DL_summary.csv' #pelicun output
DV_rec_time_input = 'DV_rec_time.csv' #pelicun output
output_path = '**insert output directory here**'

DT_calculation.run_treads(input_parameters, RCtable_input, IF_delays_input, DMG_input, DL_summary_input, DV_rec_time_input, output_path)
