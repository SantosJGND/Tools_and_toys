import plotly.graph_objs as go
from plotly import tools
from plotly.offline import iplot


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


def plot_global_pca(feats,PCA_color_ref,title= '',height= 500,width= 950):
    ## perform MeanShift clustering.
    bandwidth = estimate_bandwidth(feats, quantile=0.15)

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=True, min_bin_freq=15)
    ms.fit(feats)
    labels1 = ms.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}

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

def window_sample_plot(Windows,PCA_color_ref,Chr= 1,windows_pick= 4,height= 1500,width= 1000):

    windows_pick= 4
    Chr= 1
    Ncols= 2
    height= 1500
    width= 1000

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

            for i in coords.keys():
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
                "opacity": .8})

                fig_pca_subplots.append_trace(trace, row + 1, subp + 1)

            fig_pca_subplots['layout']['yaxis' + str(n_plot)].update(title='PC{}'.format(subp + 2))
            fig_pca_subplots['layout']['xaxis' + str(n_plot)].update(title='PC1')
            n_plot += 1

    fig_pca_subplots['layout'].update(height= height,width= width)

    iplot(fig_pca_subplots)
