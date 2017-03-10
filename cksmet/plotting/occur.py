from matplotlib.pylab import *
import seaborn as sns
sns.set_style('whitegrid')

def mesfac(st):
    fig,axL = subplots(ncols=2,figsize=(8,4))
    sca(axL[0])
    loglog()
    plot(st.plnt_mes.mes,st.plnt_mes.mes_formula,'.')
    xlabel('MES (TPS)')
    ylabel('MES (Formula)')
    s = '{:.2f} $\pm$ {:.2f} (dex)'.format(
        st.logmes_diff_std, st.logmes_diff_med
    )
    text(0.6,0.1,s,transform=axL[0].transAxes)
    sca(axL[1])
    loglog()
    plot(st.plnt_mes.mes,st.plnt_mes.mes_formula_scaled,'.')
    xlabel('MES (TPS)')
    ylabel('MES (Formula-Scaled)')
    x = [10,1e4]
    for ax in axL:   
        sca(ax)
        ax.plot(x,x)
        xlim(*x)
        ylim(*x)
        
    fig.set_tight_layout(True)

