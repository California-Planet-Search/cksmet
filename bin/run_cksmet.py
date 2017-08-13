#!/usr/bin/env python
from argparse import ArgumentParser
import glob
import os
from matplotlib.pylab import *
import cPickle as pickle

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('calibrate-lamo', parents=[psr_parent])
    psr2.set_defaults(func=calibrate_lamo)

    psr2 = subpsr.add_parser('calc-comp', parents=[psr_parent])
    psr2.set_defaults(func=calc_comp)

    psr2 = subpsr.add_parser('calc-occur', parents=[psr_parent])
    psr2.set_defaults(func=calc_occur)

    psr2 = subpsr.add_parser('stats', parents=[psr_parent])
    psr2.set_defaults(func=stats)

    psr2 = subpsr.add_parser('create-plots', parents=[psr_parent], )
    psr2.set_defaults(func=create_plots)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

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

def calc_comp(args):
    import cksmet.analysis
    import cksmet.io
    comp = cksmet.analysis.calc_completeness()
    fn = cksmet.io.COMPLETENESS_FILE
    print "computing completeness"
    with open(fn,'w') as f:
        pickle.dump(comp,f)
        print "saved {}".format(fn) 

def calc_occur(args):
    import cksmet.io
    print "calc occur"
    df = cksmet.io.load_table('occur-nper=2-nsmet=5',cache=2)
    df = cksmet.io.load_table('occur-nsmet=5',cache=2)
    df = cksmet.io.load_table('occur-nsmet=2',cache=2)

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

