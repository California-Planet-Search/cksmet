# CKS Metallicity

## Calibrate the LAMOST metallicities to the CKS scale

run_cksmet.py calibrate-lamo

## Calc Completeness

run_cksmet.py calc-comp

## Produce occurrence bins

run_cksmet.py calc-occur

## Fit various regions of the per-prad-smet plane with

run_cksmet.py fit-occur

## Create plots

run_cksmet.py create-plot all

## Output fit statistics

run_cksmet.py create-val all

## Output sample statistics

run_cksmet.py create-table all

## 

run_cksmet.py calibrate-lamo
run_cksmet.py calc-comp # issue with big endian, workaround run twice
run_cksmet.py calc-occur
run_cksmet.py fit-occur
run_cksmet.py create-plot all
run_cksmet.py create-table all
run_cksmet.py create-val all
run_cksmet.py update-paper

