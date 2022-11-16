<p align="center">
<img src="https://github.com/carlosmolinahutt/ESR-Lab-Repositories/blob/master/dt_framework.png" 
	height="350"/>
</p>

<h1 align = "center"> TREADS: Tool for Recovery Estimation And Downtime Simulation of buildings</h1>
	    
</p>

## Introduction

`treads` is a Python package to evaluate earthquake induced downtime and model recovery of buildings. This tool implements the framework presented in: 

Molina Hutt, C., Vahanvaty, T. and Kourehpaz, P. (2022). “An analytical framework to assess earthquake induced downtime and model recovery of buildings.” *Earthquake Spectra*, 38(2): 1283-1320. https://doi.org/10.1177%2F87552930211060856 

This tool is fully compatible with SimCenter’s tool for loss assessment, i.e., pelicun (https://github.com/NHERI-SimCenter/pelicun) 

## Requirements

`treads` runs under Python 3.6+. The following packages are required for it to work properly:

`numpy`  `pandas` `os` `sys` `more_itertools` `json` 

You can install these using `pip`.

## Installation

`treads` is available at the Python Package Index (PyPI). You can simply install it using `pip` as follows:

```
pip install treads
```
## Basic Demo
```python
import DT_calculation 	# refer to "Example" folder

input_parameters = 'input_parameters.json'
RCtable_input = 'Repair_Class_Table.csv'
IF_delays_input = 'IF_delays_input.csv'

DMG_input = 'DMG.csv' 	# pelicun output
DL_summary_input = 'DL_summary.csv' 	# pelicun output
DV_rec_time_input = 'DV_rec_time.csv' 	# pelicun output

output_path = '**insert output directory here**'

DT_calculation.run_treads(input_parameters, RCtable_input, IF_delays_input, DMG_input, DL_summary_input, DV_rec_time_input, output_path)
```
## Outputs

`treads` estimates earthquake-induced downtime to achieve Functional Recovery (FR), Re-Occupancy (RO), and Shelter-in-Place (SiP) post-earthquake recovery states for residential buildings. The following output files will be generated once you run `treads`: 

- **RC_component.csv:**  Component repair class matrix. 
- **DT_summary.csv:**  10th percentile, 90th percentile, median, and mean downtime estimates.
- **RS_stats.csv:**  Probability of a building not achieving different recovery states immediately after an earthquake. 
- **DT_stepfunc_xx.csv:**  Governing recovery trajectories to each recovery state (xx= FR, RO, SiP).
- **DT_path_xx.xlsx:**  Recovery trajectories to each recovery state for each repair path (xx= FR, RO, SiP).
- **RT_stepfunc_xx.xlsx:**  Repair time stepping functions for each repair sequence when each recovery state is achieved (xx= FR, RO, SiP).
- **RT_RSeq_xx.csv:**  Repair time per story for each repair sequence when each recovery state is achieved (xx= FR, RO, SiP).
- **IF_delays.csv:**  Impeding factor delays. 

## Contact

Pouria Kourehpaz, University of British Columbia, Vancouver, BC, Canada. email: pouria.kourehpaz@ubc.ca
