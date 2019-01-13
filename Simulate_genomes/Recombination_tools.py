

import numpy as np
import plotly.graph_objs as go


def Check_cop(Cop_func,Cop_kwargs,range_diff= [-10,10],steps= 100):
    
    Cop_windows= {}
    
    for angle in np.linspace(range_diff[0],range_diff[1],steps):
        
        Cop, ID= Cop_func(angle,range_diff,**Cop_kwargs)
        
        Cop_windows[angle]= Cop
    
    
    probs= [x for x in Cop_windows.keys()]
    winds= [Cop_windows[x] for x in Cop_windows.keys()]
    
    fig_data= [go.Scatter(
    x= probs,
    y= winds,
    mode= 'markers'
    )]
    
    layout = go.Layout(
        title= 'Reccombination map. {}'.format(ID),
        yaxis=dict(
            title='Reccombination probability',
            range= [min(winds)-1e-4,max(winds) * 2]
            ),
        xaxis=dict(
            title='Proxy genome position')
    )

    fig= go.Figure(data=fig_data, layout=layout)
    
    return fig



#######
#######

def Rec_prob_uniform(angle,range_windows,Cop):
    ID= 'uniform'
    return Cop, ID


def Rec_prob_sinusoid(angle,range_windows,Cop_range):
    ID= 'sinusoid'
    
    progress= (angle - range_windows[0]) #/ (range_windows[1] - range_windows[0])
    #print(progress)
    Cop= (sin(progress) + 1) * Cop_range
    
    return Cop, ID


def Rec_prob_modal(angle,range_windows,modes,multiplier,N= 100,bandwidth= .2):
    ID= 'modal. modes: {}'.format(modes)
    N_modes= len(modes)
    progress= (angle - range_windows[0]) / (range_windows[1] - range_windows[0])
    
    sam= np.repeat(modes,N)
    sam= np.array(sam).reshape(-1,1)
    
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(sam)
    
    log_dens = kde.score_samples(progress)
    
    Cop= np.exp(log_dens)[0] * multiplier
    #print(Cop)
    return Cop, ID


def Rec_prob_rdist(angle,range_windows,multiplier,c,loc,scale):
    ID= 'rdist'
    #Size= (range_windows[1] - range_windows[0])
    
    progress= (angle - range_windows[0]) / ((range_windows[1] - range_windows[0]) / 2) - 1 
    #print(progress)
    
    Cop= multiplier * scipy.stats.rdist.pdf(progress,c,loc,scale)
    
    return Cop, ID
