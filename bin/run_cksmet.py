#!/usr/bin/env python
from argparse import ArgumentParser
import glob
import os
from matplotlib.pylab import *
import cPickle as pickle
from collections import OrderedDict
import cksmet.plotting.smet
import cksmet.plotting.occur
import cksmet.plotting.comp
import cksmet.tables
import cksmet.calibrate
import cksmet.plotting.calibrate
import cksmet.analysis
import cksmet.io
import cksmet.kstest

def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('calibrate-lamo', parents=[psr_parent])
    psr2.set_defaults(func=calibrate_lamo)

    psr2 = subpsr.add_parser('create-samples', parents=[psr_parent])
    psr2.set_defaults(func=create_samples)

    psr2 = subpsr.add_parser('calc-comp', parents=[psr_parent])
    psr2.set_defaults(func=calc_comp)

    psr2 = subpsr.add_parser('calc-occur', parents=[psr_parent])
    psr2.set_defaults(func=calc_occur)

    psr2 = subpsr.add_parser('fit-occur', parents=[psr_parent])
    psr2.set_defaults(func=fit_occur)

    # Products for paper
    psr2 = subpsr.add_parser('create-val', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_val)

    psr2 = subpsr.add_parser('create-plot', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_plot)

    psr2 = subpsr.add_parser('create-table', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_table)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)


def calibrate_lamo(args):
    cksmet.calibrate.calibrate_lamo()

def create_samples(args):
    df = cksmet.io.load_table('cks-cuts',cache=2)
    df = cksmet.io.load_table('field-cuts',cache=2)
    df = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=2)

def calc_comp(args):
    cksmet.io.load_object('comp',cache=2)

def calc_occur(args):
    print "calc occur"
    df = cksmet.io.load_object('occur-test',cache=2) # just to make sure things are working
    df = cksmet.io.load_object('occur-hot-jup',cache=2)
    df = cksmet.io.load_object('occur-nper=2-nsmet=5',cache=2)
    df = cksmet.io.load_object('occur-nsmet=5',cache=2)
    df = cksmet.io.load_object('occur-nsmet=2',cache=2)
    df = cksmet.io.load_object('occur-nsmet=1',cache=2)

def fit_occur(args):
    fits = [
        'fit_per-sub-se',
        'fit_per-sup-se',
        'fit_per-sub-sn',
        'fit_per-sup-sn',
#        'fit_per-sup-ss',

        'fit_smet-hot-se',
        'fit_smet-warm-se',

        'fit_smet-hot-sn',
        'fit_smet-warm-sn',

        'fit_smet-warm-ss',
        'fit_smet-hot-jup',
    ]

    for fit in fits:
        cksmet.io.load_object(fit,cache=2)

def create_table(args):
    w = Workflow()
    w.create_file('table', args.name ) 

def create_plot(args):
    w = Workflow()
    w.create_file('plot', args.name ) 

def create_val(args):
    w = Workflow()
    w.create_file('val',args.name) 

def update_paper(args):
    w = Workflow()
    w.update_paper() 


class Workflow(object):
    def __init__(self):
        d = OrderedDict()
        d['lamo-on-cks'] = cksmet.plotting.calibrate.validation_lamo
        d['prad-smet-cuts'] = cksmet.plotting.smet.cuts
        d['stellar-samples'] = cksmet.plotting.samples.samples
        d['prad-fe'] = cksmet.plotting.smet.prad_fe
        d['prad-fe-percentiles'] = cksmet.plotting.smet.prad_fe_percentiles
        d['per-prad-slices-equal-stars'] = lambda : cksmet.plotting.smet.period_prad_slices(mode='four-equal-stars')
        d['smet-snr'] = cksmet.plotting.samples.smet_snr
        d['checkerboard'] =  cksmet.plotting.occur.fig_checkerboard
        d['prob-detect-transit'] =  cksmet.plotting.comp.fig_prob_detect_transit
        d['per-small2'] = cksmet.plotting.occur.fig_per_small2
        d['smet-large4'] = cksmet.plotting.occur.fig_smet_large4
        d['smet-small4'] = cksmet.plotting.occur.fig_smet_small4
        self.plot_dict = d

        d = OrderedDict()
        d['cuts-lamost'] = cksmet.tables.cuts_lamost
        d['cuts-planets'] = cksmet.tables.cuts_planets
        d['cuts-field'] = cksmet.tables.cuts_field
        d['smet-stats'] = cksmet.tables.smet_stats
        d['kstest'] = cksmet.kstest.kstest_region
        self.table_dict = d

        d = OrderedDict()
        d['samp'] = cksmet.tables.val_samp
        d['fit'] = cksmet.tables.val_fit
        self.val_dict = d

        d = OrderedDict()
        d['table'] = self.table_dict
        d['plot'] = self.plot_dict
        d['val'] = self.val_dict
        self.all_dict = d

    def key2fn(self, key, kind):
        if kind=='plot':
            return 'fig_'+key+'.pdf'
        if kind=='table':
            return 'tab_'+key+'.tex'
        if kind=='val':
            return 'val_'+key+'.tex'
            
    def create_file(self, kind, name):
        i = 0

        for key, func in self.all_dict[kind].iteritems():
            if kind=='plot':
                if name=='all':
                    func()
                elif name==key:
                    func()
                else:
                    continue

                    
                fn = self.key2fn(key, 'plot')
                gcf().savefig(fn)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                fn = self.key2fn(key, 'table')
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='val':
                fn = self.key2fn(key, 'val')
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                lines1 = [
                    "\\newcommand{\%s}[1]{%%" % key,
                    "\IfEqCase{#1}{%",
                ]

                lines2 = [
                    "}[\PackageError{tree}{Undefined option to tree: #1}{}]%",
                    "}%"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            i+=1

        if i==0:
            assert False, key + " not a valid key"

    def update_paper(self):
        for kind, d in self.all_dict.iteritems():
            for key, val in d.iteritems():
                fn = self.key2fn(key, kind)
                cmd = 'cp {} paper/'.format(fn)
                print cmd
                os.system(cmd)

if __name__=="__main__":
    main()

