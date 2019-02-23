import plotly.graph_objs as go
from plotly import tools
from plotly.offline import iplot


def Struct_subplots_lib(Struct_lib,vector_lib,Npops= 4,Ncols= 2,range_diff= [0,10],func= 'prior_func',kwargs= 'prior_kwargs',xaxis= '',yaxis= '',title= ''):
    
    from synth_tools.Generate_samples import Check_Path
    
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

        fig, Pops, prior= Check_Path(Npops,vector_lib,prior_func,prior_kwargs,Pops= [],random= True,n_comp= vector_lib.shape[1],range_diff= range_diff,steps= range_diff[1])

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
    from synth_tools.Recombination_tools import Check_cop

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
