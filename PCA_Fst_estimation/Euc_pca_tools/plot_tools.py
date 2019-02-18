import plotly
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import itertools as it

def plot_fst(fst_x,y_true):

    trace1= go.Scatter(
        x= [np.exp(x) for x in fst_x],
        y= [np.exp(x) for x in y_true],
        mode= 'markers'
        )

    trace2= go.Scatter(
        x= fst_x,
        y= y_true,
        mode= 'markers'
        )

    fig = tools.make_subplots(rows=1, cols=2,
                             subplot_titles=('Euc to Fst','log of relationship'))

    fig.append_trace(trace1, 1, 1)
    fig.append_trace(trace2, 1, 2)

    fig['layout']['xaxis1'].update(title='Fst')
    fig['layout']['xaxis2'].update(title='log Fst')

    fig['layout']['yaxis1'].update(title='Euc')
    fig['layout']['yaxis2'].update(title='log Euc')


    layout = go.Layout(
        title= 'Euclidian distances to fst',
        yaxis=dict(
            title='fst'),
        xaxis=dict(
            title='fst')
    )

    fig= go.Figure(data=fig, layout=layout)
    iplot(fig)



def plot_3D(feats,color_indexes,colors= {},var_comps= [],Names= [],ref_names= []):
    
    if len(colors) == 0:
        colors= {x:x for x in color_indexes}
    
    if len(var_comps) == 0:
        var_comps= [0 for x in range(feats.shape[1])]
    
    if len(Names) == 0:
        Names= ['ind' + str(x).zfill(int()) for x in range(feats.shape[0])]
    
    if len(ref_names) == 0:
        ref_names= {x:str(x) for x in color_indexes}
    
    fig_data= [go.Scatter3d(
            x = feats[color_indexes[i],0],
            y = feats[color_indexes[i],1],
            z = feats[color_indexes[i],2],
            mode= "markers",
            name= ref_names[i],
            text= ['ID: {}, gp: {}'.format(Names[x], ref_names[i]) for x in color_indexes[i]],
            marker= {
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
            'color': colors[i],
            "opacity": 1
          }
        ) for i in color_indexes.keys()]

    layout = go.Layout(
        scene= go.Scene(
        yaxis=dict(
            title='{}'.format(round(var_comps[1],3))),
        xaxis=dict(
        title= '{}'.format(round(var_comps[0],3))),
        zaxis=dict(
        title= '{}'.format(round(var_comps[2],3))))
    )


    fig = go.Figure(data=fig_data,layout= layout)
    iplot(fig)



def plot_vertix(feats,label_select,Vertices,m_coeff= 1,b= 0, color= '#1f77b4'):

    fig_vertix= []
    d= 0

    for pair in it.combinations(Vertices,2):
        coordinates= [np.mean(feats[label_select[z],:],axis= 0) for z in pair]
        coordinates= np.array(coordinates)

        MS_pair_dist= pairwise_distances(coordinates,metric= 'euclidean')
        iu_select= np.triu_indices(2,1)
        MS_pair_dist= MS_pair_dist[iu_select][0]
        MS_pair_fst= np.exp(m_coeff * np.log(MS_pair_dist) + b) 

        fig_vertix.append(go.Scatter3d(
        x= coordinates[:,0],
        y= coordinates[:,1], 
        z= coordinates[:,2],
        text= ['group: {}'.format(x) for x in pair], 
        name= 'pred Fst {}: {}'.format(pair, round(MS_pair_fst,2)),
        marker=dict(
            size=5,
            color= '#1f77b4'
        ),
        line=dict(
            color= d,
            width= 2,
        )))

        d += 1


    fig_vertix.extend([go.Scatter3d(
    x= feats[label_select[i],0],
    y= feats[label_select[i],1], 
    z= feats[label_select[i],2],
    mode= 'markers',
    name= 'Cl: {}'.format(i),
    marker=dict(
        size=5,
        color= '#1f77b4'
    )
    ) for i in Vertices])


    layout = go.Layout(
        scene= go.Scene(
        yaxis=dict(
            showgrid= False,
            title='PC2'),
        xaxis=dict(
            showgrid= False,
            title= 'PC1'),
        zaxis=dict(
            showgrid= False,
            title= 'PC3')
    ))


    fig = go.Figure(data=fig_vertix,layout= layout)
    iplot(fig)
