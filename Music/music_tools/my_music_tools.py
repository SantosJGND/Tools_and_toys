import pandas as pd
import numpy as np
import itertools as it
import collections
import re
import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

import plotly.graph_objs as go
from plotly import tools
from plotly.offline import iplot

from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale

import music_tools.thinkplot as thinkplot
import music_tools.thinkdsp as thinkdsp

from IPython.display import clear_output

from structure_tools.Modules_tools import return_fsts2

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)



def get_sound_coords(pc_coords_s,sampleRate = 44100,
                                frequency = 1500,
                                length = 5,
                                group= 2,
                                PC_select= 1,
                                qtl= 0.95
                                ):
    inters= []
    for windl in pc_coords_s:

        ci = scipy.stats.norm.interval(qtl, loc=np.mean(windl), scale=np.std(windl))
        inters.append(ci)

    inters= np.array(inters)

    from scipy.interpolate import interp1d
    q= inters[:,1] - inters[:,0]

    max_q= max(q)
    q= q / max_q

    ###
    ###
    ###
    
    t = np.linspace(0, length, inters.shape[0])
    f2= interp1d(t, q, kind='cubic')

    xfloor = np.linspace(0, length, sampleRate * length)
    roof= f2(xfloor)
    y = np.sin(xfloor * frequency * 2) * roof
    y= y * max_q / 2
    
    print_some= np.linspace(0,len(xfloor)-1,1000)
    print_some= [int(x) for x in print_some]

    fig_test= [go.Scatter(
        x= [xfloor[e] for e in print_some],
        y= [y[e] for e in print_some]
    )]
    
    layout= go.Layout(
        title= 'gp{}_pc{}_sR{}_Hz{}_l{}'.format(group,PC_select,sampleRate,frequency,length),
        xaxis= dict(
            title= 'seconds'
        ),
        yaxis= dict(
            title= 'amp'
        )
    )
    
    fig_freqs= go.Figure(data= fig_test,layout= layout)


    from scipy.io import wavfile
    

    wavfile.write('gp{}_pc{}_sR{}_Hz{}_l{}.wav'.format(group,PC_select,sampleRate,frequency,length), sampleRate, y)
    
    wave= thinkdsp.Wave(y,framerate=sampleRate)

    wave.make_audio()

    return fig_freqs, wave


#### some functions to get the spectrogram from thinkdsp.
def find_index(x, xs):
    """Find the index corresponding to a given value in an array."""
    n = len(xs)
    start = xs[0]
    end = xs[-1]
    i = round((n-1) * (x - start) / (end - start))
    return int(i)

def return_spec_array(spect,high= 5000):
        """
        high: highest frequency component to plot
        """
        fs = spect.frequencies()
        i = None if high is None else find_index(high, fs)
        fs = fs[:i]
        ts = spect.times()

        # make the array
        size = len(fs), len(ts)
        array = np.zeros(size, dtype=np.float)

        # copy amplitude from each spectrum into a column of the array
        for j, t in enumerate(ts):
            spectrum = spect.spec_map[t]
            array[:, j] = spectrum.amps[:i]

        return array, ts, fs

def plot_frame(fr_array,frame_plot,spec_fs,height= 0,width= 0):
    frame_plot= 57

    probs= fr_array[:,frame_plot]
    probs[probs < 0] = 0

    fig= [go.Scatter(
        x= spec_fs,
        y= probs,
        mode= 'lines'
    )]
    

    layout= go.Layout(
        xaxis= dict(
            title= 'frequency'
        ),
        yaxis= dict(
            title= 'amp'
        )
    )
    
    if height:
        layout.update(height= height)
    if width:
        layout.update(width= width)

    iplot(fig)




def cluster_threshold(center_list,t,scor_func):
    trim_cl= {}

    for cml in center_list:
        if not trim_cl:
            trim_cl[0]= np.array([cml])
            continue

        d= 0

        for clamp in trim_cl.keys(): 
            dists= abs(trim_cl[clamp] - cml)

            if min(dists) < 200:
                trim_cl[clamp]= np.array([*list(trim_cl[clamp]),cml])
                d += 1

        if d == 0:
            trim_cl[len(trim_cl)]= np.array([cml])
    
    pval_trim= [[scor_func(x) for x in trim_cl[z]] for z in trim_cl.keys()]
    #pval_trim= [[np.exp(sklearn_scor.score_samples(x.reshape(-1,1))) for x in trim_cl[z]] for z in trim_cl.keys()]
    trim_cl= [trim_cl[x][np.argmax(pval_trim[x])] for x in trim_cl.keys()]

    return trim_cl



def frame_peaks(array_spec,
                spec_fs,
                spec_ts,
                frame= 0,
                Sample_N= 500,
                p_threshold= 0.0004,
                amp_cap= 4,
                peak_cap= .7,
                peak_iso= 200,
                band_qtl= 0.02,
                frame_plot= False):
    
    ## get probs from
    probs= list(array_spec[:,frame]) 
    probs= np.array(probs)
    probs[probs > amp_cap]= amp_cap

    prob_sum= np.sum(probs)
    probs= probs / prob_sum
        
    
    # #############################################################################
    # Compute clustering with MeanShift

    # The following bandwidth can be automatically detected using
    new_freqs= np.random.choice(list(spec_fs),Sample_N,replace= True,p= probs)

    new_freqs= new_freqs.reshape(-1,1)

    bandwidth = estimate_bandwidth(new_freqs, quantile=band_qtl, n_samples=Sample_N)

    if bandwidth== 0:
        bandwidth= peak_iso

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True,cluster_all= False).fit(new_freqs)

    labels = ms.labels_
    cluster_centers = ms.cluster_centers_

    cluster_centers= list(it.chain(*cluster_centers))

    ## trim_clusters:
    ### interpolation makes it easier to chose between neibhour centroids
    ### that are unlikely to exist as obs. frequency values.
    from scipy.interpolate import interp1d
    f2= interp1d(spec_fs, array_spec[:,frame], kind='cubic')
    
    cluster_centers= cluster_threshold(cluster_centers,peak_iso,f2)
    cluster_centers= cluster_threshold(cluster_centers,peak_iso,f2)

    ####
    #### get amplitudes of peaks and store them
    
    peak_cent= []
    amps_centres= []

    shapes= []
    for cent in cluster_centers:

        closest= abs(spec_fs - cent)

        closet= np.argmin(closest)

        amp_sel= array_spec[closet,frame]

        if amp_sel >= peak_cap:
            peak_cent.append(cent)
            amps_centres.append(amp_sel)
    
    ## get time stamps for each of the peaks.
    time_spec= [spec_ts[frame]]* len(amps_centres)

    if frame_plot:
        
        kde= KernelDensity(kernel='gaussian', bandwidth= bandwidth).fit(new_freqs)
        X_plot = np.linspace(0, max(spec_fs) + 100, 1000)[:, np.newaxis]
        log_dens = kde.score_samples(X_plot)
        
        fig= [go.Scatter(
            x= spec_fs,
            y= array_spec[:,frame],
            mode= 'lines'
        )]
        
        #fig= [go.Scatter(x=X_plot[:, 0], y=np.exp(log_dens), mode='lines', fill='tozeroy', line=dict(color='#AAAAFF', width=2))]

        shapes= []

        for center in peak_cent:

            shapes.append({
                'type': 'line',
                'x0': center,
                'y0': 0,
                'x1': center,
                'y1': max(array_spec[:,frame]),
                'line': {
                    'color': 'red',
                    'width': 4,
                    'dash': 'solid'
                },
            })
        
        layout= go.Layout(
            title= 'frame inx: {}'.format(frame),
            shapes= shapes,
            xaxis= dict(title= 'frequency'),
            yaxis= dict(title= 'amplitude')
        )
        
        figure_frame= go.Figure(data= fig,layout= layout)
        
        return peak_cent, time_spec, amps_centres, figure_frame

    else:
        return peak_cent, time_spec, amps_centres



def filter_output(samps_tracks):

    fig= [go.Scatter(
        x= samps_tracks[:,0],
        y= samps_tracks[:,1],
        mode= 'markers',
        marker= dict(
            color= samps_tracks[:,2],
            colorscale= 'Viridis'
        )
    )]

    layout= go.Layout(
        title= 'Spect. peaks',
        xaxis= dict(
            title= 'time (s)'),
        yaxis= dict(
            title= 'frequency')
    )
    fig= go.Figure(data= fig,layout= layout)

    return fig



def break_down_spec(actual_tracks,N_neigh= 0,
                                    ms_layer2= True,
                                    scale_spec= False,
                                       qtl_I= 0.05,
                                       qtl_II= 0.1,
                                       clst_all_I= True,
                                       clst_all_II= True):
    ###

    if scale_spec:
        samps_tracks= scale(actual_tracks,axis= 0)
    else:
        samps_tracks= actual_tracks

    # #############################################################################

    bandwidth = estimate_bandwidth(samps_tracks, quantile= qtl_I, n_samples=samps_tracks.shape[0])
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True,cluster_all= clst_all_I).fit(samps_tracks)

    labels = ms.labels_
    coords= {
        z: [x for x in range(len(labels)) if labels[x] == z] for z in list(set(labels)) if z != -1
    }
    
    names_plots= ['MS1']
    
    ####
    
    fig= [go.Scatter(
        x= [actual_tracks[x,0] for x in coords[i]],
        y= [actual_tracks[x,1] for x in coords[i]],
        mode= 'markers',
        name= str(i),
        marker= dict(
            color= i
        )
    ) for i in coords.keys()]

    layout= go.Layout(
        title= 'MS clust. I',
        xaxis= dict(
            title= 'time (s)'),
        yaxis= dict(
            title= 'frequency')
    )

    figures= [go.Figure(data= fig, layout= layout)]
    ####

    ## an extra step to clean this up.
    if ms_layer2:
        extra_cls= {}

        for clust in coords.keys():
            subset= samps_tracks[coords[clust],:]
            subset= scale(subset,axis= 0)

            if subset.shape[0] > 10:

                bandwidth = estimate_bandwidth(subset, quantile= qtl_II, n_samples=subset.shape[0])
                if bandwidth > 0:
                    ms = MeanShift(bin_seeding=True,cluster_all= clst_all_II,bandwidth= bandwidth).fit(subset)

                    labels_local = ms.labels_
                    coords_local= {
                        z: [coords[clust][x] for x in range(len(labels_local)) if labels_local[x] == z] for z in list(set(labels_local)) if z != -1
                    }
                    coords_local= {z:coords_local[z] for z in coords_local.keys() if len(coords_local[z]) > 3}
                    coords_keys= sorted(coords_local.keys())

                    if len(coords_keys) > 1:
                        coords[clust]= coords_local[coords_keys[0]]
                        for cl in coords_keys[1:]:
                            extra_cls[len(extra_cls) + len(coords)]= coords_local[cl]


        coords.update(extra_cls)
        names_plots.append('MSII')

    ########################
    ## get just neighbours
    ##
    if N_neigh:
        extra_cls= {}


        for clust in coords.keys():
            subset= samps_tracks[coords[clust],:]

            if subset.shape[0] >= 2*N_neigh:

                t= list(np.arange(0,len(coords[clust]),N_neigh))
                coords_local= list(sorted(coords[clust]))
                if len(t) > 1:
                    if len(coords_local) - t[-1] < N_neigh:
                        t[-1] = len(coords_local)
                    else:
                        t.append(len(coords_local))
                    coords[clust]= coords_local[t[0]:t[1]]
                    for cl in range(1,len(t)-1):
                        extra_cls[len(extra_cls) + len(coords)]= coords_local[t[cl]:t[cl +1]]



        coords.update(extra_cls)
        names_plots[-1]= names_plots[-1] + '_neighs{}'.format(N_neigh)

    fig= [go.Scatter(
        x= [actual_tracks[x,0] for x in coords[i]],
        y= [actual_tracks[x,1] for x in coords[i]],
        mode= 'markers',
        name= str(i),
        marker= dict(
            color= i
        )
    ) for i in coords.keys()]

    layout= go.Layout(
        title= 'MS clust. II',
        xaxis= dict(
            title= 'time (s)'),
        yaxis= dict(
            title= 'frequency')
    )

    figures.append(go.Figure(data= fig, layout= layout))
    
    return coords, figures, names_plots




def meanfreq_plot(samps_tracks,coords,ts_list,peaks,height= 0,width= 0):
    fig= [go.Scatter(
        x= [samps_tracks[x,0] for x in coords[i]],
        y= [np.mean([samps_tracks[x,1] for x in coords[i]])] * len(coords[i]),
        mode= 'markers',
        name= str(i),
        marker= dict(
            color= i
        )
    ) for i in coords.keys() if (max(coords[i]) - min(coords[i])) / len(coords[i]) < 100]

    layout= go.Layout(
        title= 'MS clust. I',
        xaxis= dict(
            title= 'time (s)'),
        yaxis= dict(
            title= 'frequency'),
        height= height,
        width= width
    )

    fig_mean= go.Figure(data= fig, layout= layout)
    return fig_mean





