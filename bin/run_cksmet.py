#!/usr/bin/env python
from argparse import ArgumentParser
import glob
import os

import pandas as pd
from matplotlib.pylab import *

import cksphys.iso
import cksphys.io
import cksphys.calc
import cksphys._isoclassify
import cksphys._isochrones
import cksphys.plotting.compare
import cksphys.plotting.hr_diagram

from cksphys.config import ISO_CSVFN
import cksphys.tables

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr_stats = subpsr.add_parser('stats', parents=[psr_parent])
    psr_stats.set_defaults(func=stats)

    psr_plots = subpsr.add_parser('create-plots', parents=[psr_parent], )
    psr_plots.set_defaults(func=create_plots)

    psr_paper = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr_paper.set_defaults(func=update_paper)

    psr_cal = subpsr.add_parser('calibrate-lamo', parents=[psr_parent])
    psr_cal.set_defaults(func=calibrate_lamo)

    args = psr.parse_args()
    args.func(args)

def stats(args):
    import cksmet.stats
    cksmet.stats.stats()

def calibrate_lamo(args):
    import cksmet.calibrate
    import cksmet.plotting.calibrate
    cksmet.calibrate.calibrate_lamo()
    cksmet.plotting.calibrate.validation_lamo()
    gcf().savefig('fig_lamo-on-cks.pdf')
    #cksspec.plotting.compare.comparison_three('cks-lamo')
    #gcf().savefig('fig_sm-sx.pdf')

def create_plots(args):
    import cksmet.plotting.smet
    '''
    cksmet.plotting.smet.cuts()
    gcf().savefig('fig_prad-smet-cuts.pdf')

    cksmet.plotting.smet.prad_fe()
    gcf().savefig('fig_prad-fe.pdf')

    cksmet.plotting.smet.prad_fe_percentiles()
    gcf().savefig('fig_prad-fe-percentiles.pdf')

    cksmet.plotting.samples.samples()
    gcf().savefig('fig_stellar-samples.pdf')

    cksmet.plotting.smet.period_prad_slices(mode='four-equal-stars')
    gcf().savefig('fig_per-prad-slices-equal-stars.pdf')

    cksmet.plotting.smet.period_prad_slices(mode='four-equal-smet')
    gcf().savefig('fig_per-prad-slices-equal-smet.pdf')

    fig1, fig2 = cksmet.plotting.samples.lamo_detectability() 
    fig1.savefig('fig_lamo-smet-hr.pdf')
    fig2.savefig('fig_lamo-smet-kepmag-cdpp.pdf')
    '''
    fig = cksmet.plotting.samples.smet_snr() 
    fig.savefig('fig_smet-snr.pdf')

def update_paper(args):
    files = [
        'fig_prad-smet-cuts.pdf',
        'fig_prad-fe.pdf',
        'fig_prad-fe-percentiles.pdf',
        'fig_stellar-samples.pdf',
        'fig_per-prad-slices-equal-stars.pdf',
        'fig_period_prad_slices',
        'fig_lamo-on-cks.pdf',
        'fig_lamo-smet-hr.pdf',
        'fig_lamo-smet-kepmag-cdpp.pdf',
        'fig_smet-snr.pdf'
    ]

    for _file in files:
        cmd = 'cp {} paper/'.format(_file)
        print cmd
        os.system(cmd)
        
if __name__=="__main__":
    main()

