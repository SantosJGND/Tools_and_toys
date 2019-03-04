from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth

import numpy as np
import itertools as it
import plotly.graph_objs as go
from plotly import tools
from plotly.offline import iplot

from structure_tools.vcf_geno_tools import check_densities


def Struct_subplots_lib(Struct_lib,vector_lib,Fsts_test,Npops= 4,Ncols= 2,range_diff= [0,10],func= 'prior_func',kwargs= 'prior_kwargs',xaxis= '',yaxis= '',title= ''):
    
    from structure_tools.Generate_samples import Check_Path_Rfreqs
    
    Ncols= 2
    titles= list(Struct_lib.keys())
    print(titles)

    fig_subplots = tools.make_subplots(rows= int(len(Struct_lib) / float(Ncols)) + (len(Struct_lib) % Ncols > 0), cols=Ncols,
                             subplot_titles=tuple(titles))

    #####
    for gp in range(len(titles)):

        pos1= int(float(gp) / Ncols) + 1
        pos2= gp - (pos1-1)*Ncols + 1

        title= titles[gp]

        prior_func= Struct_lib[title]['prior_func']
        prior_kwargs= Struct_lib[title]['prior_kwargs']

        fig, Pops, prior= Check_Path_Rfreqs(Npops,vector_lib,Fsts_test,prior_func,prior_kwargs,Pops= [],random= True,n_comp= vector_lib.shape[1],range_diff= range_diff,steps= range_diff[1])

        for trace1 in fig['data']:
            
            fig_subplots.append_trace(trace1, pos1, pos2)
        
        if yaxis:
            fig_subplots['layout']['yaxis' + str(gp + 1)].update(title=yaxis)
        if xaxis:
            fig_subplots['layout']['xaxis' + str(gp + 1)].update(title=xaxis)


    layout = go.Layout(
        title= title
    )

    fig= go.Figure(data=fig_subplots, layout=layout)
    iplot(fig_subplots)


def Rec_subplots_lib(Cop_lib,vector_lib,Npops= 4,Ncols= 2,range_diff= [0,10],steps= 100,func= 'cop_func',kwargs= 'cop_kwargs',xaxis= '',yaxis= '',title= ''):
    
    from plotly import tools
    from structure_tools.Recombination_tools import Check_cop

    titles= list(Cop_lib.keys())
    print(titles)

    fig_subplots = tools.make_subplots(rows= int(len(Cop_lib) / float(Ncols)) + (len(Cop_lib) % Ncols > 0), cols=Ncols,
                             subplot_titles=tuple(titles))


    #####
    for gp in range(len(titles)):

        pos1= int(float(gp) / Ncols) + 1
        pos2= gp - (pos1-1)*Ncols + 1

        title= titles[gp]

        Cop_func= Cop_lib[title]['cop_func']
        Cop_kwargs= Cop_lib[title]['cop_kwargs']

        fig= Check_cop(Cop_func,Cop_kwargs,steps= steps,range_diff= range_diff)

        trace1= fig[0]

        fig_subplots.append_trace(trace1, pos1, pos2)
        
        if yaxis:
            fig_subplots['layout']['yaxis' + str(gp + 1)].update(title=yaxis)
        if xaxis:
            fig_subplots['layout']['xaxis' + str(gp + 1)].update(title=xaxis)


    layout = go.Layout(
        title= title
    )

    fig= go.Figure(data=fig_subplots, layout=layout)
    iplot(fig_subplots)

def Admixture_subplots_lib(Geneo,Ncols= 2,xaxis= '',yaxis= '',title= ''):

    snmf_gps= sorted(Geneo.keys())

    Ncols= 2

    fig_box_subplots = tools.make_subplots(rows=int(len(snmf_gps) / float(Ncols)) + (len(snmf_gps) % Ncols > 0), cols=2,
                             subplot_titles=tuple(['Gp: {}'.format(x) for x in snmf_gps]))

    for gp in range(len(snmf_gps)):

        pos1= int(float(gp) / Ncols) + 1
        pos2= gp - (pos1-1)*Ncols + 1

        for prop in range(len(snmf_gps) - 1):
            trace= go.Box(
                y= [x[prop] for x in Geneo[gp]],
                name= 'gp: {}, Anc: {}'.format(gp,prop)
            )

            fig_box_subplots.append_trace(trace, pos1, pos2)

        if yaxis:
            fig_box_subplots['layout']['yaxis' + str(gp + 1)].update(title=yaxis)
        if xaxis:
            fig_box_subplots['layout']['xaxis' + str(gp + 1)].update(title=xaxis)


    layout = go.Layout(
        title= title,
    )

    fig= go.Figure(data=fig_box_subplots, layout=layout)
    iplot(fig)


def plot_global_pca(feats,label_select,PCA_color_ref,title= '',height= 500,width= 950):
    ##

    from plotly import tools

    fig_pca_subplots = tools.make_subplots(rows=1, cols=2,subplot_titles=tuple([title]*2))

    for subp in range(2):

        n_plot= subp

        coords= label_select

        for i in coords.keys():
            trace= go.Scatter(
            x = feats[coords[i],0],
            y = feats[coords[i],subp + 1],
            mode= "markers",
            name= str(i),
            marker= {
            'color': PCA_color_ref[i],
            'line': {'width': 0},
            'size': 6,
            'symbol': 'circle',
            "opacity": .8})

            fig_pca_subplots.append_trace(trace, int(n_plot/float(2)) + 1, subp + 1)

        fig_pca_subplots['layout']['yaxis' + str(n_plot + 1)].update(title='PC{}'.format(subp + 2))
        fig_pca_subplots['layout']['xaxis' + str(n_plot + 1)].update(title='PC1')

    fig_pca_subplots['layout'].update(height= height,width= width)

    iplot(fig_pca_subplots)


def window_sample_plot(Windows,label_select,PCA_color_ref,plot_who= [],shade= [],Chr= 1,windows_pick= 4,height= 1500,width= 1000):

    Ncols= 2
    opac= [.8] * len(label_select)
    
    if shade:
        for mit in shade:
            opac[mit]= .2

    windows_pick= np.random.choice(list(Windows[1].keys()),windows_pick,replace= False)
    titles= ['window: ' + str(x) for x in windows_pick]
    titles= np.repeat(titles,2)

    fig_pca_subplots = tools.make_subplots(rows= int(len(titles) / float(Ncols)) + (len(titles) % Ncols > 0), cols=Ncols,
                             subplot_titles=tuple(titles))

    n_plot= 1
    for row in range(len(windows_pick)):
        pca_local= PCA(n_components=3, whiten=False,svd_solver='randomized')

        feats_local= pca_local.fit_transform(Windows[1][windows_pick[row]])

        for subp in range(2):

            coords= label_select
            paint= list(coords.keys())
            
            if plot_who:
                paint= plot_who

            for i in paint:
                trace= go.Scatter(
                x = feats_local[coords[i],0],
                y = feats_local[coords[i],subp + 1],
                mode= "markers",
                name= str(i),
                marker= {
                'color': PCA_color_ref[i],
                'line': {'width': 0},
                'size': 6,
                'symbol': 'circle',
                "opacity": opac[i]})

                fig_pca_subplots.append_trace(trace, row + 1, subp + 1)

            fig_pca_subplots['layout']['yaxis' + str(n_plot)].update(title='PC{}'.format(subp + 2))
            fig_pca_subplots['layout']['xaxis' + str(n_plot)].update(title='PC1')
            n_plot += 1

    fig_pca_subplots['layout'].update(height= height,width= width)

    iplot(fig_pca_subplots)


def PC_analysis_plot(pc_density,pc_coords,kde_class_labels,PCA_color_ref,range_windows= [],plot_choice= 'coords'):

    x_range= range(len(pc_density))
    
    if range_windows:
        x_range= range(range_windows[0],range_windows[1])
        
    if plot_choice == 'density':
        x_coords= [z for z in it.chain(*[[x] * 100 for x in x_range])]
        y_coords= [z for z in it.chain(*[list(np.linspace(-8,8,100)) for x in range(len(pc_density))])]
        z_coords= [z for z in it.chain(*pc_density)]


        fig_data= [go.Scatter(
        x= x_coords,
        y= y_coords,
        mode= 'markers',
        marker= {
            'color': z_coords,
            'colorscale':'Viridis',
            'line': {'width': 0},
            'size': 8,
            'symbol': 'circle',
            "opacity": .6
        })
        ]

        layout = go.Layout(
            title= 'PC1 density',
            yaxis=dict(
                title='PC1 density of projections across data sets'),
            xaxis=dict(
                title='Ideogram')
        )


        fig= go.Figure(data=fig_data, layout=layout)

    if plot_choice== 'coords':
        x_coords= [z for z in it.chain(*[[x] * len(kde_class_labels) for x in x_range])]
        z_coords= [z for z in it.chain(*pc_coords)]

        class_colors= [z for z in it.chain(*[kde_class_labels for x in range(len(pc_coords))])]
        class_colors= [PCA_color_ref[x-1] for x in class_colors]

        fig_data= [go.Scatter(
        x= x_coords,
        y= z_coords,
        mode= 'markers',
        marker= dict(
            color= class_colors,
            size= 4,
            opacity= .8
        )
        )
        ]

        layout = go.Layout(
            title= 'PC1 coordinates',
            yaxis=dict(
                title='Individual positions along PC1 across data sets'),
            xaxis=dict(
                title='Ideogram')
        )

        fig= go.Figure(data=fig_data, layout=layout)

    iplot(fig)



def fst_window_plot(freq_matrix,ref_labels,sort= True,window_range= [],y_range= [0,.3]):
    
    tuples= list(it.combinations(ref_labels,2))
    x_range= list(range(len(freq_matrix)))
    
    if sort:
        freq_matrix= sorted(freq_matrix)
        
    if window_range:
        x_range= list(range(window_range[0],window_range[1]))
    
    freq_matrix= np.array(freq_matrix)
        
    fig_fst= [go.Scatter(
        x= x_range,
        y= freq_matrix[:,i],
        name= str(tuples[i])
    ) for i in range(freq_matrix.shape[1])]
    
    layout = go.Layout(
        title= 'ref Fst,sorted= {}'.format(sort),
        yaxis=dict(
            title='Fst',
            range= [0,.3]),
        xaxis=dict(
            title= ''.format(['extraction order','sorted'][int(sort)]))
    )

    fig= go.Figure(data=fig_fst, layout=layout)

    iplot(fig)


def freq_dist_plt(freqs_matrix,n_chose= 100,height= 500,width= 900):

    title= ['across','individual']

    fig_freq_subplots = tools.make_subplots(rows=2, cols=1,subplot_titles=tuple(title))


    Chose= np.random.choice(range(freqs_matrix.shape[0]),n_chose)

    Across= list(it.chain(*[freqs_matrix[x] for x in Chose]))

    X_plot = np.linspace(0, 1, 1000)

    freq_kde = KernelDensity(kernel='gaussian', bandwidth=0.02).fit(np.array(Across).reshape(-1,1))

    log_dens = freq_kde.score_samples(X_plot.reshape(-1,1))

    trace= go.Scatter(x=X_plot, y=np.exp(log_dens), 
                                mode='lines', fill='tozeroy', name= 'freq',
                                line=dict(color='blue', width=2))
    ##

    fig_freq_subplots.append_trace(trace, 1, 1)
    fig_freq_subplots['layout']['yaxis' + str(1)].update(title='density')
    fig_freq_subplots['layout']['xaxis' + str(1)].update(title='frequency')


    dist_freq= check_densities(freqs_matrix,n_chose)

    trace1= go.Scatter(
        x= np.linspace(0, 1, dist_freq.shape[1]),
        y= np.mean(dist_freq,axis= 0),
        mode= 'markers+lines',
        name= 'mean'
    )

    trace2= go.Scatter(
        x= np.linspace(0, 1, dist_freq.shape[1]),
        y= np.std(dist_freq,axis= 0),
        mode= 'markers+lines',
        name= 'sd'
    )

    fig_freq_subplots.append_trace(trace1, 2, 1)
    fig_freq_subplots.append_trace(trace2, 2, 1)
    fig_freq_subplots['layout']['yaxis2'].update(title='density')
    fig_freq_subplots['layout']['xaxis2'].update(title='frequency')


    fig_freq_subplots['layout'].update(height= height,width= width)

    iplot(fig_freq_subplots)

