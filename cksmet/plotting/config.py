from matplotlib.pylab import * 

def label():
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    xt = [1,3,10,30,100,300]
    yt = [0.5,1,2,4,8,16]
    xticks(xt,xt)
    yticks(yt,yt)
    fig = gcf()
    rename()

def rename():
    """
    Plotting routines will give short machine readable names, here we
    change to human-readable names

    """
    fig = gcf()
    namemap = {
        'prob_det':'Prob Det.',
        'prob_tr':'Prob Tr.',
        'prob_trdet':'Prob Tr. and Det.',
        'plnt_occur':'Planet Occurrence per Bin',
        'prad':'Planet Size (Earth-radii)',
        'fstars':'Completeness',
    }
    for o in fig.findobj(matplotlib.text.Text):
        for oldname, newname in namemap.iteritems():
            if o.get_text()==oldname:
                o.set_text(newname)
