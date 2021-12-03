# -*- coding: utf-8 -*-
#
# This file is part of DT assessment framework.
#
# This code is developed to support the framework to estimate earthquake induced downtime 
# and recovery trajectory of buildings 
# proposed by C. Molina Hutt, T. Vahanvaty, P. Kourehpaz. 
#
# The proposed framework can be found in the manuscript entitled 
# "An analytical framework to assess earthquake induceddowntime and model recovery of buildings"
# The manuscript will be submitted to Earthquake Spectra for publication.
# 
# This code will be publicly available on GitHub after the manuscript is published, 
# and it should not be distributed until then.
#
# Contributor(s):
# Pouria Kourehpaz

"""
Define input parameters and files and execute this script to perform DT assessment

"""
import timeit
start = timeit.default_timer()
import DT_calculation

rep_phases = [3,3,3,3,4]    #last element should indicate # of basements, if none, insert 0. Note: max # of stories in each phase <=3 (except # of basements)
fl_area = 16*[10000] #per floor area in ft2 ; note: basements come after the above grade floors
total_cost = 47564000 #building replacement cost
coeff_financing = [0.16,0.72,0.12]    #entry 1,2,3: contribution from (1) insurance, (2) private loan, (3) public loan
Qt_facade = 12*2010/30 #number of facade units normalized by the FEMA P-58 component unit
reconst_time = 14 #reconstruction time per story in days
n_structural = [3,7] #lower and upper limits for the number of structural components required to peform stabilization times
t_structural = [6,4] #stabilization times for a structural component; lower and upper limit values in days
n_facade = [20,100] #lower and upper limits for the number of facade components required to peform stabilization times
t_facade = [0.1429,0.076] #stabilization times for a facade component; lower and upper limit values in days
n_elev = 2 #number of elevators
max_num_workers = [15, 15, 15, 9, 9, 6, 6] #max number of workers per repair sequence in each phase for low-rise buildings; note: these numbers will be doubled for mid-rise and tripled for high-rise buildings


RCtable_input = 'Repair_Class_Table.csv'
IF_delays_input = 'IF_delays_input.csv'
DMG_input = 'DMG.csv'
DL_summary_input = 'DL_summary.csv'
DV_rec_time_input = 'DV_rec_time.csv'

parameters=[n_elev, rep_phases, fl_area, total_cost, coeff_financing, Qt_facade, reconst_time, n_structural, t_structural, n_facade, t_facade, max_num_workers]
files = [RCtable_input, DMG_input, DL_summary_input, DV_rec_time_input, IF_delays_input]
DT_calculation.run_analysis(parameters, files)

stop = timeit.default_timer()
print('Time: ', stop - start)


