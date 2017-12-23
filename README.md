# CKS Metallicity

## Calibrate the LAMOST metallicities to the CKS scale

run_cksmet.py calibrate-lamo
run_cksmet.py create-plot lamo-on-cks

## Apply cuts and create samples

run_cksmet.py create-samples
run_cksmet.py create-plot prad-smet-cut

## Calc Completeness

run_cksmet.py calc-comp

## Produce occurrence bins

run_cksmet.py calc-occur

## Occurrence surface

run_cksmet.py calc-occur-surface

## Fit various regions of the per-prad-smet plane with

run_cksmet.py fit-occur

## Create plots

run_cksmet.py create-plot all

## Sample the population

run_cksmet.py calc-population

## Output fit statistics

run_cksmet.py create-val all

## Output sample statistics

run_cksmet.py create-table all


## Cookbook to do a fresh build of paper 

Takes about 5min (with bootstrap, about 20 min)

```bash
rm data/*.pkl
rm load_table_cache.hdf
run_cksmet.py calibrate-lamo
#run_cksmet.py calibrate-lamo-bootstrap
run_cksmet.py create-samples
run_cksmet.py calc-comp
run_cksmet.py calc-occur
run_cksmet.py calc-occur-surface
run_cksmet.py fit-occur
run_cksmet.py calc-population
run_cksmet.py create-table all
run_cksmet.py create-plot all
run_cksmet.py create-val all
run_cksmet.py update-paper 
```

## Cookbook to perturb metallicities

We tested the extent to which zero-point offsets between the LAMOST
and CKS metallicities could affect the final results. Go into
cksmet/io.py and add or subtract a given metallicity from the lamost-cal
