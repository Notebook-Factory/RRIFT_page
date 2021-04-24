# Reference region and input function tail (RRIFT) method

<img src="images/NotebookFactory.png">

# Introduction <a name="introduction"></a> 

This is a Jupyter notebook that introduces interactive plots created in plotly from some of figures in the paper [Pharmacokinetic modeling of dynamic contrast‐enhanced MRI using a reference region and input function tail](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27913). 

Link to the paper's GitHub repository: https://github.com/MPUmri/RRIFT

GitHub link to this notebook: https://github.com/Notebook-Factory/RRIFT_notebooks


# Abstract <a name="abstract"></a> 

##  Purpose
Quantitative analysis of dynamic contrast‐enhanced MRI (DCE‐MRI) requires an arterial input function (AIF) which is difficult to measure. We propose the reference region and input function tail (RRIFT) approach which uses a reference tissue and the washout portion of the AIF.

## Methods
RRIFT was evaluated in simulations with 100 parameter combinations at various temporal resolutions (5‐30 s) and noise levels (σ = 0.01‐0.05 mM). RRIFT was compared against the extended Tofts model (ETM) in 8 studies from patients with glioblastoma multiforme. Two versions of RRIFT were evaluated: one using measured patient‐specific AIF tails, and another assuming a literature‐based AIF tail.

## Results
RRIFT estimated the transfer constant $K^{trans}$ and interstitial volume $V_{e}$ with median errors within 20% across all simulations. RRIFT was more accurate and precise than the ETM at temporal resolutions slower than 10 s. The percentage error of $K^{trans}$ had a median and interquartile range of −9 ± 45% with the ETM and −2 ± 17% with RRIFT at a temporal resolution of 30 s under noiseless conditions. RRIFT was in excellent agreement with the ETM in vivo, with concordance correlation coefficients (CCC) of 0.95 for $K^{trans}$, 0.96 for $V_{e}$, and 0.73 for the plasma volume $V_{p}$ using a measured AIF tail. With the literature‐based AIF tail, the CCC was 0.89 for $K^{trans}$, 0.93 for $V_{e}$ and 0.78 for $V_{p}$.

## Conclusions 
Quantitative DCE‐MRI analysis using the input function tail and a reference tissue yields absolute kinetic parameters with the RRIFT method. This approach was viable in simulation and in vivo for temporal resolutions as low as 30 s.

import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from numpy import array
import sys
from scipy import stats
from scipy.io import loadmat
import statsmodels.api as sm
import pandas as pd
from itertools import cycle
import itertools
from pandas.core.common import flatten
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from IPython.core.display import display, HTML
init_notebook_mode(connected=True)
config={'showLink': False, 'displayModeBar': False}



## Figure 1 <a name="figure-1"></a>


* (A) Reference Tissue <a name="reference-tissue"></a>

file = loadmat('fig1vars.mat')

Crr = np.squeeze(np.asarray((np.reshape(file['Crr'], 120))))
Ct = file['Ct']
Cp = np.squeeze(np.asarray((np.reshape(file['Cp'], 120))))
t = np.squeeze(np.asarray((np.reshape(file['t'], 120))))
idx = np.squeeze(np.asarray((np.reshape(file['idx'], 4))))
estKepRR = np.squeeze(np.asarray((np.reshape(file['estKepRR'], 1))))
pkCE = file['pkCE']
denum = np.squeeze(np.asarray((np.reshape(file['denum'], 84))))
num = np.squeeze(np.asarray((np.reshape(file['num'], 84))))

# ======================================= First subplot ===================================================================================

sub1 = go.Figure()

# sub1.add_trace(go.Scatter(x=t,y=Crr,mode='lines',line_color='#41B3A8')

sub1.add_trace(go.Scatter(x=t, y=Crr,
                    mode='lines', line_color='#41B3A8'))

sub1.update_layout(title='(A) Reference Tissue',
                   xaxis_title='Time [min]',
                   yaxis_title=r'Concentration [mM]',
                   plot_bgcolor="#fff",
                   xaxis=dict(range=[0, 10],
                              tickvals=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
                   yaxis=dict(range=[-0.01, 0.25], #the range from zero cuts out some values off
                              tickvals=[0, 0.05, 0.10, 0.15, 0.20, 0.25]),
                   
                 )

sub1.add_annotation(text='K<sub>RR</sub><sup>trans</sup> 0.007 [min<sup>-1</sup>] K<sub>ep, RR</sub> 0.50 [min<sup>-1</sup>]',
                  xref="x", yref="y",
                  x=5, y=0.02, showarrow=False, font = dict(size = 26))

sub1.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
sub1.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
               
# sub1.show()

plot(sub1, filename = 'figures/fig1-1.html', config = config)

display(HTML('figures/fig1-1.html'))

* (B) Tissue of Interest <a name="tissue-of-interest"></a>

file = loadmat('fig1vars.mat')

Ct = file['Ct']
t = np.squeeze(np.asarray((np.reshape(file['t'], 120))))
idx = np.squeeze(np.asarray((np.reshape(file['idx'], 4))))
estKepRR = np.squeeze(np.asarray((np.reshape(file['estKepRR'], 1))))
pkCE = file['pkCE']


# Matlab's indexing starts at 1, whereas Python's starts from 0, so we need to
# lower all the indexes in idx by one to match the approprate data from the array

for i in range(idx.size):
    idx[i] = idx[i]-1

    
fig = go.Figure()

palette_lines = cycle(px.colors.qualitative.Plotly)
palette_labels = cycle(px.colors.qualitative.Plotly)
line_colours = cycle(['#292348','#585E9A','#88A8D5','#83C9C4'])
pkCE_index = cycle([4, 9, 14, 19])

pkCE_reshaped = np.ravel(pkCE)


for el in idx:
    reshaped_Ct = np.reshape(Ct[:,el], 120)
    final_Ct = np.squeeze(np.asarray(reshaped_Ct))
    trace_label = round(pkCE_reshaped[next(pkCE_index)], 3)
    fig.add_trace(go.Scatter(name=str(trace_label),
                             x=t, 
                             y=final_Ct, 
                             mode='lines',
                             line_color=next(line_colours)))

fig.update_layout(title='(B) Tissue of Interest',
                   xaxis=dict(title='Time [min]',
                              range=[0, 10],
                              tickvals=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
                   yaxis=dict(title='Concentration [mM]',
                              range=[-0.1, 1.5],
                              tickvals=[0, 0.5, 1, 1.5]),
                   plot_bgcolor="#fff",
)

fig.add_annotation(text='Estimate k<sub>ep, RR</sub> (min<sup>-1</sup>)',
                  xref="x", yref="y",
                  x=6, y=1.19, showarrow=False, font = dict(size = 21))

fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig1-2.html', config = config)

display(HTML('figures/fig1-2.html'))

* (C) Input Function <a name="input-function"></a>

file = loadmat('fig1vars.mat')

Cp = np.squeeze(np.asarray((np.reshape(file['Cp'], 120))))
t = np.squeeze(np.asarray((np.reshape(file['t'], 120))))

fig = go.Figure()

fig.add_trace(go.Scatter(x=t, 
                         y=Cp, 
                         mode='lines',
                         line_color='red'))

fig.update_layout(title='(C) Input Function',
                   xaxis_title='Time [min]',
                   yaxis_title='Concentration [mM]',
                   xaxis_range=[0,10],
                   yaxis_range=[-0.1,11],
                   plot_bgcolor="#fff",
)

#Adding a shape region
fig.add_vrect(
    x0=3, x1=10,
    fillcolor="LightGray", opacity=0.5,
    layer="below", line_width=0,
    annotation=dict(text="Input Function Tail",
                  xref="paper", yref="paper",
                    x=7.6, y=0.7,
#                   xanchor="center", yanchor="center", 
#                   xshift=-140, yshift=-100,
                  showarrow=False, font=dict(size=20, color='black'))
)

# T_{tail}
# An arrow for T_{tail} that points to the right side of the highlighted region
fig.add_annotation(
    x=10,  # arrows' head
    y=5,  # arrows' head
    ax=3, # arrows' tail
    ay=5,  # arrows' tail
    xref='x',
    yref='y',
    axref='x',
    ayref='y',
    text='',  # if you want only the arrow
    showarrow=True,
    arrowhead=3,
    arrowsize=2,
    arrowwidth=2,
    arrowcolor='black'
)

# An arrow for T_{tail} that points to the left side of the highlighted region
fig.add_annotation(
    x=3,  # arrows' head
    y=5,  # arrows' head
    ax=10, # arrows' tail
    ay=5,  # arrows' tail
    xref='x',
    yref='y',
    axref='x',
    ayref='y',
    text='',  # if you want only the arrow
    showarrow=True,
    arrowhead=3,
    arrowsize=2,
    arrowwidth=2,
    arrowcolor='black'
)

# Text annotation for T_{tail} below the arrows
fig.add_annotation(  
    x=6.5,
    y=4.6,
    text='T<sub>tail</sub>',
    xanchor="center", yanchor="top", 
    showarrow=False, font=dict(size=20, color='black'))


# T_{start}
# An arrow for T_{start}
fig.add_annotation(
    x=3,  # arrows' head
    y=0,  # arrows' head
    ax=3, # arrows' tail
    ay=2,  # arrows' tail
    xref='x',
    yref='y',
    axref='x',
    ayref='y',
    text='T<sub>start</sub>',
    font=dict(size=16, color='black'),
    showarrow=True,
    arrowhead=3,
    arrowsize=1,
    arrowwidth=2,
    arrowcolor='black'
)

# T_{end}
# An arrow for T_{end}
fig.add_annotation(
    x=10,  # arrows' head
    y=0,  # arrows' head
    ax=10, # arrows' tail
    ay=2,  # arrows' tail
    xref='x',
    yref='y',
    axref='x',
    ayref='y',
    text='T<sub>end</sub>',
    font=dict(size=16, color='black'),
    showarrow=True,
    arrowhead=3,
    arrowsize=1,
    arrowwidth=2,
    arrowcolor='black'
)

fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black', tick0=0, dtick=1)
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
# fig.show()

plot(fig, filename = 'figures/fig1-3.html', config = config)

display(HTML('figures/fig1-3.html'))


* (D) RRIFT Fit <a name="rrift-fit"></a> 

file = loadmat('fig1vars.mat')

denum = np.squeeze(np.asarray((np.reshape(file['denum'], 84))))
num = np.squeeze(np.asarray((np.reshape(file['num'], 84))))

line_fit = sm.OLS(num ,sm.add_constant(denum)).fit().fittedvalues

fig=go.Figure()
fig.add_trace(go.Scatter(name="Markers",
                         x=denum, 
                         y=num, 
                         mode='markers',
                         marker_symbol="circle-open", 
                         marker_size=14, 
                         marker_color='red'))
fig.add_trace(go.Scatter(name="Linear Fit",
                         x=denum, 
                         y=line_fit, 
                         mode='lines',
                         line_color="black"))

fig.add_annotation(text='Slope | K<sub>RR</sub><sup>trans</sup>: 0.071 min<sup>-1</sup><br>R<sub>2</sub>: 0.9996',
                  xref="x", yref="y",
                  x=4.7, y=0.1, showarrow=False, font = dict(size = 22))

fig.update_layout(title='(D) RRIFT Fit',
                   xaxis_title='Denominator [mM * min]',
                   yaxis_title='Numerator [mM]',
                   xaxis=dict(showline=True,
                              range=[-0.1,7],
                              tickvals=[0, 1, 2, 3, 4, 5, 6, 7]),
                   yaxis=dict(showline=True,
                              range=[-0.015,0.5],
                              tickvals=[0, 0.1, 0.2, 0.3, 0.4, 0.5]),
                   plot_bgcolor="#fff",
)

fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')

# fig.show()

plot(fig, filename = 'figures/fig1-4.html', config = config)

display(HTML('figures/fig1-4.html'))

## Figure 2 <a name="figure-2"></a> 

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))

curErr3 = file['curErr3']
errQt3 = file['errQt3']
errMd3 = file['errMd3']

curErr4 = file['curErr4']
errQt4 = file['errQt4']
errMd4 = file['errMd4']

curErr5 = file['curErr5']
errQt5 = file['errQt5']
errMd5 = file['errMd5']

fig = make_subplots(rows=1, cols=3, 
                    subplot_titles=('K<sup>trans</sup>', 'V<sup>e</sup>', 'V<sup>p</sup>'))

legend_values=cycle(['5','10','15','30'])

# First subplot 
line_colours = cycle(['#253494','#2C7fB8','#41B6C4','#A1DAB4'])

for i in range(0,4):
    finalErrMd3=np.squeeze(np.asarray((np.reshape(errMd3[i,:], 6))))
    # take another look here - is it errorQt[i,:,1] because MATLAB indexing starts from 1,
    # or is it just like in MATLAB?
    finalErrorDiff3=np.squeeze(np.asarray((np.reshape(abs(errQt3[i,:,0]-errMd3[i,:]), 6))))
    finalErrorDiffMinus3=np.squeeze(np.asarray((np.reshape(abs(errQt3[i,:,1]-errMd3[i,:]), 6))))
    fig.append_trace(go.Scatter(
        x=listSigmaC,
        y=finalErrMd3,
        showlegend=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff3,
            arrayminus=finalErrorDiffMinus3,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4
    ),row=1, col=1)
    


# Second subplot

line_colours = cycle(['#253494','#2C7fB8','#41B6C4','#A1DAB4'])

for i in range(0,4):
    finalErrMd4=np.squeeze(np.asarray((np.reshape(errMd4[i,:], 6))))
    
    finalErrorDiff4=np.squeeze(np.asarray((np.reshape(abs(errQt4[i,:,0]-errMd4[i,:]), 6))))
    finalErrorDiffMinus4=np.squeeze(np.asarray((np.reshape(abs(errQt4[i,:,1]-errMd4[i,:]), 6))))
    fig.append_trace(go.Scatter(
        x=listSigmaC,
        y=finalErrMd4,
        showlegend=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff4,
            arrayminus=finalErrorDiffMinus4,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4       
    ),row=1, col=2)


# Third subplot
    


line_colours = cycle(['#253494','#2C7fB8','#41B6C4','#A1DAB4'])

for i in range(0,4):
    finalErrMd5=np.squeeze(np.asarray((np.reshape(errMd5[i,:], 6))))
    finalErrorDiff5=np.squeeze(np.asarray((np.reshape(abs(errQt5[i,:,0]-errMd5[i,:]), 6))))
    finalErrorDiffMinus5=np.squeeze(np.asarray((np.reshape(abs(errQt5[i,:,1]-errMd5[i,:]), 6))))
    fig.add_trace(go.Scatter(
        x=listSigmaC,
        y=finalErrMd5,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff5,
            arrayminus=finalErrorDiffMinus5,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        showlegend=True,
        line_color=next(line_colours),
        line_width=4
    ),row=1, col=3)



fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]')
fig.update_yaxes(range=[-50,50], tickvals=[-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50], row=1, col=1)
fig.update_yaxes(range=[-50,50], tickvals=[-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50], row=1, col=2)
fig.update_yaxes(range=[-50,50], tickvals=[-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50], row=1, col=3)

fig.update_layout(yaxis=dict(title='Percent Error',
                             showline=True,
                             ), 
                  plot_bgcolor='#ffffff'
                  #legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30) -- change legend position
)


# names = set()
# fig.for_each_trace(
#     lambda trace:
#         trace.update(showlegend=False)
#         if (trace.name in names) else names.add(trace.name))

fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig2.html', config = config)

display(HTML('figures/fig2.html'))

### Subfigures with sliders 

* $K^{trans}$: <a name="fig-2-1"></a> 

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))

curErr3 = file['curErr3']
errQt3 = file['errQt3']
errMd3 = file['errMd3']

# First subplot
line_colours = cycle(['#253494','#2C7fB8','#41B6C4','#A1DAB4'])
legend_values=cycle(['5','10','15','30'])

fig = go.Figure()
for i in range(0,4):
    finalErrMd3=np.squeeze(np.asarray((np.reshape(errMd3[i,:], 6))))
    finalErrorDiff3=np.squeeze(np.asarray((np.reshape(abs(errQt3[i,:,0]-errMd3[i,:]), 6))))
    finalErrorDiffMinus3=np.squeeze(np.asarray((np.reshape(abs(errQt3[i,:,1]-errMd3[i,:]), 6))))
    fig.add_trace(go.Scatter(
        visible=False,
        x=listSigmaC,
        y=finalErrMd3,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff3,
            arrayminus=finalErrorDiffMinus3,
            symmetric=False,
            visible=True,
            thickness=4),
        line_color=next(line_colours),
        line_width=4, showlegend=False, name=next(legend_values)
    ))
    

# Default value on the slider    
fig.data[0].visible = True
legend_values = [5, 10, 15, 30]

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],  # layout attribute
        value=legend_values[i],
        label=legend_values[i]        
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Temporal resolution [seconds]: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='K<sub>trans</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-35,5],
                             tickvals=[-35,-30,-25,-20,-15,-10,-5,0,5],
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders)

fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig2-1.html', config = config)

display(HTML('figures/fig2-1.html'))

* $V_{e}$: <a name="fig-2-2"></a> 

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))
curErr4 = file['curErr4']
errQt4 = file['errQt4']
errMd4 = file['errMd4']

# Second subplot

line_colours = cycle(['#253494','#2C7fB8','#41B6C4','#A1DAB4'])
legend_values=cycle(['5','10','15','30'])
fig = go.Figure()

for i in range(0,4):
    finalErrMd4=np.squeeze(np.asarray((np.reshape(errMd4[i,:], 6))))
    
    finalErrorDiff4=np.squeeze(np.asarray((np.reshape(abs(errQt4[i,:,0]-errMd4[i,:]), 6))))
    finalErrorDiffMinus4=np.squeeze(np.asarray((np.reshape(abs(errQt4[i,:,1]-errMd4[i,:]), 6))))
    fig.add_trace(go.Scatter(
        visible=False,
        x=listSigmaC,
        y=finalErrMd4,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff4,
            arrayminus=finalErrorDiffMinus4,
            symmetric=False,
            visible=True,
            thickness=4),
        line_color=next(line_colours),
        line_width=4, showlegend=False, name=next(legend_values)       
    ))
    
    
# Default value on the slider    
fig.data[0].visible = True

legend_values = [5, 10, 15, 30]

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
#               {"title": "Current value selected on slider: " + str(legend_values[i])}
             ],  # layout attribute
        value=legend_values[i],
        label=legend_values[i]        
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Temporal resolution [seconds]: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='V<sub>e</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-10,10],
                             tickvals=[-10,-5,0,5, 10],
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()
plot(fig, filename = 'figures/fig2-2.html', config = config)

display(HTML('figures/fig2-2.html'))

* $V_{p}$: <a name="fig-2-3"></a>

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))
curErr5 = file['curErr5']
errQt5 = file['errQt5']
errMd5 = file['errMd5']


# Third subplot

fig = go.Figure()
line_colours = cycle(['#253494','#2C7fB8','#41B6C4','#A1DAB4'])
legend_values=cycle(['5','10','15','30'])
for i in range(0,4):
    finalErrMd5=np.squeeze(np.asarray((np.reshape(errMd5[i,:], 6))))
    finalErrorDiff5=np.squeeze(np.asarray((np.reshape(abs(errQt5[i,:,0]-errMd5[i,:]), 6))))
    finalErrorDiffMinus5=np.squeeze(np.asarray((np.reshape(abs(errQt5[i,:,1]-errMd5[i,:]), 6))))
    fig.add_trace(go.Scatter(
        visible=False,
        x=listSigmaC,
        y=finalErrMd5,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff5,
            arrayminus=finalErrorDiffMinus5,
            symmetric=False,
            visible=True,
            thickness=4),
        line_color=next(line_colours),
        line_width=4, showlegend=False, name=next(legend_values)) 
    )
    
    
# Default value on the slider    
fig.data[0].visible = True

legend_values = [5, 10, 15, 30]

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=legend_values[i],
        label=legend_values[i]        
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Temporal resolution [seconds]: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='V<sub>p</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-170,50],
                             tickvals=[-160,-130,-100,-70,-40,-10,20,50],
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')



# fig.show()

plot(fig, filename = 'figures/fig2-3.html', config = config)

display(HTML('figures/fig2-3.html'))

## Figure 3 <a name="figure-3"></a> 

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))

curErr = file['curErr']
errQt = file['errQt']
errMd = file['errMd']

curErr1 = file['curErr1']
errQt1 = file['errQt1']
errMd1 = file['errMd1']

curErr2 = file['curErr2']
errQt2 = file['errQt2']
errMd2 = file['errMd2']

fig = make_subplots(rows=1, cols=3, 
                    subplot_titles=('k&#770;<sub>ep,RR</sub>', 'K<sup>trans</sup><sub>RR</sub>', 'V<sub>e,RR</sub>'),
                    )

# First subplot 
line_colours = cycle(['#006837','#31A354','#78C679','#C2E699'])

legend_values = cycle([5, 10, 15, 30])


for i in range(0,4):
    finalErrMd=np.squeeze(np.asarray((np.reshape(errMd[i,:], 6))))
    finalErrorDiff=np.squeeze(np.asarray((np.reshape(abs(errQt[i,:,0]-errMd[i,:]), 6))))
    finalErrorDiffMinus=np.squeeze(np.asarray((np.reshape(abs(errQt[i,:,1]-errMd[i,:]), 6))))
    fig.add_trace(go.Scatter(
        x=listSigmaC,
        y=finalErrMd,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff,
            arrayminus=finalErrorDiffMinus,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4, showlegend=False
    ),row=1, col=1)
    



# Second subplot

line_colours = cycle(['#006837','#31A354','#78C679','#C2E699'])

for i in range(0,4):
    finalErrMd1=np.squeeze(np.asarray((np.reshape(errMd1[i,:], 6))))
    
    finalErrorDiff1=np.squeeze(np.asarray((np.reshape(abs(errQt1[i,:,0]-errMd1[i,:]), 6))))
    finalErrorDiffMinus1=np.squeeze(np.asarray((np.reshape(abs(errQt1[i,:,1]-errMd1[i,:]), 6))))
    fig.add_trace(go.Scatter(
        x=listSigmaC,
        y=finalErrMd1,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff1,
            arrayminus=finalErrorDiffMinus1,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4, showlegend=False       
    ),row=1, col=2)




# Third subplot

legend_values=cycle(['5','10','15','30'])

    
line_colours = cycle(['#006837','#31A354','#78C679','#C2E699'])

for i in range(0,4):
    finalErrMd2=np.squeeze(np.asarray((np.reshape(errMd2[i,:], 6))))
    
    finalErrorDiff2=np.squeeze(np.asarray((np.reshape(abs(errQt2[i,:,0]-errMd2[i,:]), 6))))
    finalErrorDiffMinus2=np.squeeze(np.asarray((np.reshape(abs(errQt2[i,:,1]-errMd2[i,:]), 6))))
    fig.add_trace(go.Scatter(
        x=listSigmaC,
        y=finalErrMd2,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff2,
            arrayminus=finalErrorDiffMinus2,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4, showlegend=True      
    ),row=1, col=3)

# General


fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


fig.update_yaxes(range=[-35,5], tickvals=[-35,-30,-25,-20,-15,-10,-5,0,5], row=1, col=1)
fig.update_yaxes(range=[-35,5], tickvals=[-35,-30,-25,-20,-15,-10,-5,0,5], row=1, col=2)
fig.update_yaxes(range=[-35,5], tickvals=[-35,-30,-25,-20,-15,-10,-5,0,5], row=1, col=3)

fig.update_xaxes(range=[0,0.05], tickvals=[0,0.01,0.02,0.03,0.04,0.05], row=1, col=1)
fig.update_xaxes(range=[0,0.05], tickvals=[0,0.01,0.02,0.03,0.04,0.05], row=1, col=2)
fig.update_xaxes(range=[0,0.05], tickvals=[0,0.01,0.02,0.03,0.04,0.05], row=1, col=3)

fig.update_layout(yaxis=dict(title='Percent Error',
                             showline=True), 
                  plot_bgcolor='#ffffff',
)

# fig.show()

plot(fig, filename = 'figures/fig3.html', config = config)

display(HTML('figures/fig3.html'))

### Subfigures with sliders 

* $\widehat{k_{ep,RR}}$ <a name="fig-3-1"></a> 

# THIS CAN PROBABLY BE LEFT OUT BUT IT ONLY WORKS HALF OF THE TIME IF I LEAVE IT OUT SO IT'S HERE FOR NOW

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))

curErr = file['curErr']
errQt = file['errQt']
errMd = file['errMd']

fig = go.Figure()

legend_values = cycle([5, 10, 15, 30])

line_colours = cycle(['#006837','#31A354','#78C679','#C2E699'])

for i in range(0,4):
    finalErrMd=np.squeeze(np.asarray((np.reshape(errMd[i,:], 6))))
    finalErrorDiff=np.squeeze(np.asarray((np.reshape(abs(errQt[i,:,0]-errMd[i,:]), 6))))
    finalErrorDiffMinus=np.squeeze(np.asarray((np.reshape(abs(errQt[i,:,1]-errMd[i,:]), 6))))
    fig.add_trace(go.Scatter(
        name=next(legend_values),
        visible=False,
        x=listSigmaC,
        y=finalErrMd,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff,
            arrayminus=finalErrorDiffMinus,
            symmetric=False,
            visible=True,
            thickness=4),
        line_color=next(line_colours),
        line_width=4, showlegend=False
    ))
    
# Default value on the slider    
fig.data[0].visible = True


# Create and add slider
legend_values = [5, 10, 15, 30]
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=legend_values[i],
        label=legend_values[i]        
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Temporal resolution [seconds]: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]


fig.update_layout(title=dict(text='k&#770;<sub>ep,RR</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-35,5],
                             tickvals=[-35,-30,-25,-20,-15,-10,-5,0,5],
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig3-1.html', config = config)

display(HTML('figures/fig3-1.html'))

* $K^{trans}_{RR}$ <a name="fig-3-2"></a> 

# THIS CAN PROBABLY BE LEFT OUT BUT IT ONLY WORKS HALF OF THE TIME IF I LEAVE IT OUT SO IT'S HERE FOR NOW

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))

curErr1 = file['curErr1']
errQt1 = file['errQt1']
errMd1 = file['errMd1']
fig = go.Figure()

legend_values = cycle([5, 10, 15, 30])

line_colours = cycle(['#006837','#31A354','#78C679','#C2E699'])
    
for i in range(0,4):
    finalErrMd1=np.squeeze(np.asarray((np.reshape(errMd1[i,:], 6))))
    
    finalErrorDiff1=np.squeeze(np.asarray((np.reshape(abs(errQt1[i,:,0]-errMd1[i,:]), 6))))
    finalErrorDiffMinus1=np.squeeze(np.asarray((np.reshape(abs(errQt1[i,:,1]-errMd1[i,:]), 6))))
    fig.add_trace(go.Scatter(
        name=next(legend_values),
        visible=False,
        x=listSigmaC,
        y=finalErrMd1,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff1,
            arrayminus=finalErrorDiffMinus1,
            symmetric=False,
            visible=True,
            thickness=4),
        line_color=next(line_colours),
        line_width=4, showlegend=False       
    ))



    
# Default value on the slider    
fig.data[0].visible = True


# Create and add slider
legend_values = [5, 10, 15, 30]
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=legend_values[i],
        label=legend_values[i]        
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Temporal resolution [seconds]: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]


fig.update_layout(title=dict(text='K<sup>trans</sup><sub>RR</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-25,5],
                             tickvals=[-25,-20,-15,-10,-5,0,5],
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig3-2.html', config = config)
display(HTML('figures/fig3-2.html'))

* $V_{e,RR}$ <a name="fig-3-3"></a> 

# THIS CAN PROBABLY BE LEFT OUT BUT IT ONLY WORKS HALF OF THE TIME IF I LEAVE IT OUT SO IT'S HERE FOR NOW

file = loadmat('fig2andfig3vars.mat')

listSigmaC = np.squeeze(np.asarray((np.reshape(file['listSigmaC'], 6))))

curErr2 = file['curErr2']
errQt2 = file['errQt2']
errMd2 = file['errMd2']

fig = go.Figure()

legend_values = cycle([5, 10, 15, 30])

line_colours = cycle(['#006837','#31A354','#78C679','#C2E699'])

for i in range(0,4):
    finalErrMd2=np.squeeze(np.asarray((np.reshape(errMd2[i,:], 6))))
    
    finalErrorDiff2=np.squeeze(np.asarray((np.reshape(abs(errQt2[i,:,0]-errMd2[i,:]), 6))))
    finalErrorDiffMinus2=np.squeeze(np.asarray((np.reshape(abs(errQt2[i,:,1]-errMd2[i,:]), 6))))
    fig.add_trace(go.Scatter(
        name=next(legend_values),
        visible=False,
        x=listSigmaC,
        y=finalErrMd2,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=finalErrorDiff2,
            arrayminus=finalErrorDiffMinus2,
            symmetric=False,
            visible=True,
            thickness=4),
        line_color=next(line_colours),
        line_width=4, showlegend=False      
    ))



    
# Default value on the slider    
fig.data[0].visible = True


# Create and add slider
legend_values = [5, 10, 15, 30]
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=legend_values[i],
        label=legend_values[i]        
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Temporal resolution [seconds]: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='V<sub>e,RR</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-10,5],
                             tickvals=[-10,-5,0,5],
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='σ<sub>noise</sub>[Mm]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig3-3.html', config = config)
display(HTML('figures/fig3-3.html'))

## Figure 4 <a name="figure-4"></a> 

file = loadmat('fig4vars.mat')

trueKtRR = np.squeeze(np.asarray(np.reshape(file['trueKtRR'], 9)))
trueVeRR = np.squeeze(np.asarray(np.reshape(file['trueVeRR'], 9)))

Ktrans_ci1 = file['Ktrans_ci1']
Ktrans_ci2 = file['Ktrans_ci2']
Ktrans_ci3 = file['Ktrans_ci3']
Ktrans_avg1 = np.squeeze(np.asarray(np.reshape(file['Ktrans_avg1'], 9)))
Ktrans_avg2 = np.squeeze(np.asarray(np.reshape(file['Ktrans_avg2'], 9)))
Ktrans_avg3 = np.squeeze(np.asarray(np.reshape(file['Ktrans_avg3'], 9)))

Ve_ci1 = file['Ve_ci1']
Ve_ci2 = file['Ve_ci2']
Ve_ci3 = file['Ve_ci3']
Ve_avg1 = np.squeeze(np.asarray(np.reshape(file['Ve_avg1'], 9)))
Ve_avg2 = np.squeeze(np.asarray(np.reshape(file['Ve_avg2'], 9)))
Ve_avg3 = np.squeeze(np.asarray(np.reshape(file['Ve_avg3'], 9)))

Vp_ci1 = file['Vp_ci1']
Vp_ci2 = file['Vp_ci2']
Vp_ci3 = file['Vp_ci3']
Vp_avg1 = np.squeeze(np.asarray(np.reshape(file['Vp_avg1'], 9)))
Vp_avg2 = np.squeeze(np.asarray(np.reshape(file['Vp_avg2'], 9)))
Vp_avg3 = np.squeeze(np.asarray(np.reshape(file['Vp_avg3'], 9)))


# The subplot titles are different in the paper and in the code (check which ones to use)
fig = make_subplots(rows=1, cols=3, 
                    subplot_titles=('K<sup>trans</sup>', 'v<sub>e</sub>', r'v<sub>p</sub>'))

legend_values=cycle(['ETM','RRIFT','RRM w/ fixed RR params'])
line_colours = cycle(['#bfbfbf', '#292447', '#e45947'])


# First subplot - Ktrans

Ktrans_ci = [Ktrans_ci3, Ktrans_ci2, Ktrans_ci1]
Ktrans_avg = [Ktrans_avg3, Ktrans_avg2, Ktrans_avg1]

for i in range(0,3):
    ci_values_negative = np.squeeze(np.asarray(np.reshape(Ktrans_ci[i][:,0], 9)))
    ci_values_positive = np.squeeze(np.asarray(np.reshape(Ktrans_ci[i][:,1], 9)))
    negative = Ktrans_avg[i]-ci_values_negative
    positive = ci_values_positive-Ktrans_avg[i]
    fig.append_trace(go.Scatter(
        x=trueKtRR,
        y=Ktrans_avg[i],
        showlegend=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=positive,
            arrayminus=negative,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4
    ),row=1, col=1)


# Second subplot - Ve

Ve_ci = [Ve_ci3, Ve_ci2, Ve_ci1]
Ve_avg = [Ve_avg3, Ve_avg2, Ve_avg1]

for i in range(0,3):
    ci_values_negative = np.squeeze(np.asarray(np.reshape(Ve_ci[i][:,0], 9)))
    ci_values_positive = np.squeeze(np.asarray(np.reshape(Ve_ci[i][:,1], 9)))
    negative = Ve_avg[i]-ci_values_negative
    positive = ci_values_positive-Ve_avg[i]
    fig.append_trace(go.Scatter(
        x=trueVeRR,
        y=Ve_avg[i],
        showlegend=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=positive,
            arrayminus=negative,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4
    ),row=1, col=2)


# Third subplot - Vp

Vp_ci = [Vp_ci3, Vp_ci2, Vp_ci1]
Vp_avg = [Vp_avg3, Vp_avg2, Vp_avg1]

for i in range(0,3):
    ci_values_negative = np.squeeze(np.asarray(np.reshape(Vp_ci[i][:,0], 9)))
    ci_values_positive = np.squeeze(np.asarray(np.reshape(Vp_ci[i][:,1], 9)))
    negative = Vp_avg[i]-ci_values_negative
    positive = ci_values_positive-Vp_avg[i]
    fig.append_trace(go.Scatter(
        x=trueKtRR,
        y=Vp_avg[i],
        showlegend=True,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=positive,
            arrayminus=negative,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4
    ),row=1, col=3)
    
    
fig.update_xaxes(title_text='Reference K<sup>trans</sup>[min<sup>-1</sup>]', row=1, col=1)
fig.update_xaxes(title_text='Reference v<sub>e</sub>', row=1, col=2)
fig.update_xaxes(title_text='Reference K<sup>trans</sup>[min<sup>-1</sup>]', row=1, col=3)

fig.update_yaxes(range=[-100,100], tickvals=[-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100], row=1, col=1)
fig.update_yaxes(range=[-100,100], tickvals=[-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100], row=1, col=2)
fig.update_yaxes(range=[-100,100], tickvals=[-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100], row=1, col=3)

    
fig.update_layout(yaxis=dict(title='Percent Error',showline=True), 
                  plot_bgcolor='#ffffff'
                  #legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30) -- change legend position
)

fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')

    
# fig.show()

plot(fig, filename = 'figures/fig4.html', config = config)
display(HTML('figures/fig4.html'))

* $K^{trans}$ <a name="fig-4-1"></a>

# THESE IMPORTS ARE PROBABLY NOT NEEDED BUT IT SEEMS TO NOT ALWAYS WORK FOR ME IF I LEAVE THEM OUT

file = loadmat('fig4vars.mat')

trueKtRR = np.squeeze(np.asarray(np.reshape(file['trueKtRR'], 9)))
trueVeRR = np.squeeze(np.asarray(np.reshape(file['trueVeRR'], 9)))

Ktrans_ci1 = file['Ktrans_ci1']
Ktrans_ci2 = file['Ktrans_ci2']
Ktrans_ci3 = file['Ktrans_ci3']
Ktrans_avg1 = np.squeeze(np.asarray(np.reshape(file['Ktrans_avg1'], 9)))
Ktrans_avg2 = np.squeeze(np.asarray(np.reshape(file['Ktrans_avg2'], 9)))
Ktrans_avg3 = np.squeeze(np.asarray(np.reshape(file['Ktrans_avg3'], 9)))

fig = go.Figure()

legend_values=cycle(['ETM','RRIFT','RRM w/ fixed RR params'])
line_colours = cycle(['#bfbfbf', '#292447', '#e45947'])

Ktrans_ci = [Ktrans_ci3, Ktrans_ci2, Ktrans_ci1]
Ktrans_avg = [Ktrans_avg3, Ktrans_avg2, Ktrans_avg1]


for i in range(0,3):
    ci_values_negative = np.squeeze(np.asarray(np.reshape(Ktrans_ci[i][:,0], 9)))
    ci_values_positive = np.squeeze(np.asarray(np.reshape(Ktrans_ci[i][:,1], 9)))
    negative = Ktrans_avg[i]-ci_values_negative
    positive = ci_values_positive-Ktrans_avg[i]
    fig.add_trace(go.Scatter(
        x=trueKtRR,
        y=Ktrans_avg[i],
        visible=False,
        showlegend=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=positive,
            arrayminus=negative,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        line_color=next(line_colours),
        line_width=4
    ))    
     
        
# Default value on the slider    
        
fig.data[0].visible = True


# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=next(legend_values),
        label=next(legend_values)
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=4,
    currentvalue={"prefix": "Value: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='K<sup>trans</sup>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-40,160],
                             tickvals=[-40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160]
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='Reference K<sup>trans</sup>[min<sup>-1</sup>]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig4-1.html', config = config)
display(HTML('figures/fig4-1.html'))

* $V_e$ <a name="fig-4-2"></a> 

# THESE IMPORTS ARE PROBABLY NOT NEEDED BUT IT SEEMS TO NOT ALWAYS WORK FOR ME IF I LEAVE THEM OUT

file = loadmat('fig4vars.mat')

trueVeRR = np.squeeze(np.asarray(np.reshape(file['trueVeRR'], 9)))

Ve_ci1 = file['Ve_ci1']
Ve_ci2 = file['Ve_ci2']
Ve_ci3 = file['Ve_ci3']

Ve_avg1 = np.squeeze(np.asarray(np.reshape(file['Ve_avg1'], 9)))
Ve_avg2 = np.squeeze(np.asarray(np.reshape(file['Ve_avg2'], 9)))
Ve_avg3 = np.squeeze(np.asarray(np.reshape(file['Ve_avg3'], 9)))


fig = go.Figure()

legend_values=cycle(['ETM','RRIFT','RRM w/ fixed RR params'])
line_colours = cycle(['#bfbfbf', '#292447', '#e45947'])

Ve_ci = [Ve_ci3, Ve_ci2, Ve_ci1]
Ve_avg = [Ve_avg3, Ve_avg2, Ve_avg1]

for i in range(0,3):
    ci_values_negative = np.squeeze(np.asarray(np.reshape(Ve_ci[i][:,0], 9)))
    ci_values_positive = np.squeeze(np.asarray(np.reshape(Ve_ci[i][:,1], 9)))
    negative = Ve_avg[i]-ci_values_negative
    positive = ci_values_positive-Ve_avg[i]
    fig.add_trace(go.Scatter(
        x=trueVeRR,
        y=Ve_avg[i],
        visible=False,
        showlegend=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=positive,
            arrayminus=negative,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4
    ))


# Default value on the slider    
        
fig.data[0].visible = True


# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=next(legend_values),
        label=next(legend_values)
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=3,
    currentvalue={"prefix": "Value: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='v<sub>e</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-30, 50],
                             tickvals=[-30, -20, -10, 0, 10, 20, 30, 40, 50]
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='Reference v<sub>e</sub>', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig4-2.html', config = config)
display(HTML('figures/fig4-2.html'))

* $V_p$ <a name="fig-4-3"></a> 

# THESE IMPORTS ARE PROBABLY NOT NEEDED BUT IT SEEMS TO NOT ALWAYS WORK FOR ME IF I LEAVE THEM OUT

file = loadmat('fig4vars.mat')

trueKtRR = np.squeeze(np.asarray(np.reshape(file['trueKtRR'], 9)))

Vp_ci1 = file['Vp_ci1']
Vp_ci2 = file['Vp_ci2']
Vp_ci3 = file['Vp_ci3']

Vp_avg1 = np.squeeze(np.asarray(np.reshape(file['Vp_avg1'], 9)))
Vp_avg2 = np.squeeze(np.asarray(np.reshape(file['Vp_avg2'], 9)))
Vp_avg3 = np.squeeze(np.asarray(np.reshape(file['Vp_avg3'], 9)))

fig = go.Figure()

legend_values=cycle(['ETM','RRIFT','RRM w/ fixed RR params'])
line_colours = cycle(['#bfbfbf', '#292447', '#e45947'])

Vp_ci = [Vp_ci3, Vp_ci2, Vp_ci1]
Vp_avg = [Vp_avg3, Vp_avg2, Vp_avg1]

for i in range(0,3):
    ci_values_negative = np.squeeze(np.asarray(np.reshape(Vp_ci[i][:,0], 9)))
    ci_values_positive = np.squeeze(np.asarray(np.reshape(Vp_ci[i][:,1], 9)))
    negative = Vp_avg[i]-ci_values_negative
    positive = ci_values_positive-Vp_avg[i]
    fig.add_trace(go.Scatter(
        x=trueKtRR,
        y=Vp_avg[i],
        visible=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=positive,
            arrayminus=negative,
            symmetric=False,
            visible=True,
            thickness=4),
        name=next(legend_values),
        legendgroup=i,
        line_color=next(line_colours),
        line_width=4
    ))
    

# Default value on the slider    
        
fig.data[0].visible = True


# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)}],
        value=next(legend_values),
        label=next(legend_values)
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=3,
    currentvalue={"prefix": "Value: "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]



fig.update_layout(title=dict(text='v<sub>p</sub>',
                             x=0.5,
                             xanchor='center'),
                  yaxis=dict(title='Percent Error',
                             showline=True,
                             range=[-80, 240],
                             tickvals=[-80, -40, 0, 40, 80, 120, 160, 200, 240]
                             ), 
                  plot_bgcolor='#ffffff',
                  legend=dict(yanchor='bottom', y=0.20, xanchor='right', x=0.30),
                  sliders=sliders
)

fig.update_xaxes(title_text='Reference K<sup>trans</sup>[min<sup>-1</sup>]', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


# fig.show()

plot(fig, filename = 'figures/fig4-3.html', config = config)
display(HTML('figures/fig4-3.html'))

## Figure 5 <a name="figure-5"></a> 

file = loadmat('fig5vars.mat')

RRIFT_measured_AIF_tail = file['RRIFT_measured_AIF_tail']
overlay_measured_AIF_tail = np.squeeze(np.asarray((np.reshape(file['overlay_measured_AIF_tail'], 2))))
RRIFT_assumed_AIF_tail = file['RRIFT_assumed_AIF_tail']
overlay_assumed_AIF_tail = np.squeeze(np.asarray((np.reshape(file['overlay_assumed_AIF_tail'], 2))))
RRM_assumed_Ktrans_RR = file['RRM_assumed_Ktrans_RR']
overlay_assumed_Ktrans_RR = np.squeeze(np.asarray((np.reshape(file['overlay_assumed_Ktrans_RR'], 2))))

# RRIFT assumed Ktrans RR
fig = go.Figure()
fig.add_trace(go.Heatmap(z=RRM_assumed_Ktrans_RR,
                          colorscale='jet',
                          zmin=0, zmax=3,
                          hovertemplate="<br>".join(["z: %{z}<extra></extra>",])
                         ))
fig.add_shape(type="line",
    x0=0, y0=0, x1=100, y1=100,
    line=dict(
        color="white",
        width=2,
        dash="dot",
    )          
)

fig.add_annotation(text="CCC: {:.4f} \n ρ: {:.4f}".format(overlay_assumed_Ktrans_RR[0], overlay_assumed_Ktrans_RR[1]),
                    xref="paper", yref="paper",
                    x=0.05, y=0.9, showarrow=False, 
                    font = dict(size = 14, color="#fff"))
                 

# RRIFT with assumed AIF tail
fig.add_trace(go.Heatmap(z=RRIFT_assumed_AIF_tail, colorscale='jet', zmin=0, zmax=3,
                         hovertemplate="<br>".join(["z: %{z}<extra></extra>"]), visible=False))


# fig.add_annotation(text="CCC: {:.4f} \n ρ: {:.4f}".format(overlay_assumed_AIF_tail[0], 
#                                                            overlay_assumed_AIF_tail[1]),
#                     xref="paper", yref="paper",
#                     x=0.05, y=0.9, showarrow=False, 
#                     font = dict(size = 14, color="#fff"), visible=False)



#RRIFT measured AIF tail
fig.add_trace(go.Heatmap(z=RRIFT_measured_AIF_tail, zmin=0, zmax=3, colorscale='jet', 
                         hovertemplate="<br>".join(["z: %{z}<extra></extra>"]), visible=False))
fig.add_shape(type="line",
    x0=0, y0=0, x1=100, y1=100,
    line=dict(
        color="white",
        width=2,
        dash="dot",
    ), visible=False
)


fig.update_layout(title='RRM with assumed K<sub>RR</sub><sup>trans</sup>',
                   xaxis_title='ETM - K<sup>trans</sup>[min<sup>-1</sup>]',
                   yaxis_title='K<sup>trans</sup>[min<sup>-1</sup>]',
                   xaxis_tickmode = "array", 
                   xaxis_tickvals = [0, 25, 50, 75, 100], 
                   xaxis_ticktext = ['0', '0.05', '0.10', '0.15', '0.20'],
                   yaxis_tickmode = "array", 
                   yaxis_tickvals = [0, 25, 50, 75, 100], 
                   yaxis_ticktext = ['0', '0.05', '0.10', '0.15', '0.20'],
                   plot_bgcolor="#fff",
                   width=700, height=700,
)

fig.update_traces(
    colorbar_tickmode="array",
    colorbar_tickvals = [0, 1, 2, 3],
    colorbar_ticktext = ['<b>0</b>', '<b>10<sup>1</sup></b>', '<b>10<sup>2</sup></b>', '<b>10<sup>3</sup></b>']
)


annotation1 = [dict(
    text = "CCC: 0.856 \n ρ: 0.909",
    xref="paper", yref="paper",
    x=0.05, y=0.9, showarrow=False,
    font = dict(size = 14, color="#fff")
)]

annotation2 = [dict(
    text = "CCC: 0.886 \n ρ: 0.949",
    xref="paper", yref="paper",
    x=0.05, y=0.9, showarrow=False,
    font = dict(size = 14, color="#fff")
)]

annotation3 = [dict(
    text = "CCC: 0.947 \n ρ: 0.964",
    xref="paper", yref="paper",
    x=0.05, y=0.9, showarrow=False,
    font = dict(size = 14, color="#fff")
)]

fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0,
        y=1.076,
        xanchor="left",
        yanchor="top",
        direction="right",
        type='buttons',
        buttons=list(
            [
            dict(label = 'RRM with assumed K<sub>RR</sub><sup>trans</sup>',
                  method = 'update',
                  args = [{'visible': [True, False, False]},
                          {'title': 'RRM with assumed K<sub>RR</sub><sup>trans</sup>',
                           'showlegend':True, 'visible': True, 'annotations': annotation1}]),
             dict(label = 'RRIFT with assumed AIF tail',
                  method = 'update',
                  args = [{'visible': [False, True, False]},
                          {'title': 'RRIFT with assumed AIF tail',
                           'showlegend':True, 'visible': True, 'annotations': annotation2}]),
             dict(label = 'RRIFT with measured AIF tail',
                  method = 'update',
                  args = [{'visible': [False, False, True]},
                          {'title': 'RRIFT with measured AIF tail',
                           'showlegend':True, 'visible': True, 'annotations': annotation3}])
            ]
        )
    )
   ]
)


# fig.show()

plot(fig, filename = 'figures/fig5.html', config = config)
display(HTML('figures/fig5.html'))

## Figure 6 <a name="figure-6"></a> 

* $K_{trans}$ <a name="fig-6-1"></a> 

files = ["fig6pat1.mat", 
         "fig6pat2.mat", 
         "fig6pat3.mat", 
         "fig6pat4.mat",
         "fig6pat5.mat", 
         "fig6pat6.mat", 
         "fig6pat7.mat", 
         "fig6pat8.mat"]

# Ktrans
Ktrans_lims = [0, 0.16]

patient = loadmat(files[2])
img1 = patient['Ktrans_image_12']
img2 = patient['Ktrans_image_13']
img3 = patient['Ktrans_image_22']
img4 = patient['Ktrans_image_23']

fig = make_subplots(rows=2, cols=2, column_widths=[0.5, 0.5], subplot_titles=('ETM', 'ETM', 'RRIFT', 'RRIFT'))

fig1 = go.Heatmap(z=img1, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  colorscale='jet', colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
fig.append_trace(fig1, 1, 1)

fig2 = go.Heatmap(z=img2, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  colorscale='jet', colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
fig.append_trace(fig2, 1, 2)

fig3 = go.Heatmap(z=img3, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  colorscale='jet', colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
fig.append_trace(fig3, 2, 1)

fig4 = go.Heatmap(z=img4, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  colorscale='jet', colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
fig.append_trace(fig4, 2, 2)
    
fig.update_layout(width = 700, height = 700, title_text='Patient 3')
fig.update_yaxes(autorange="reversed")
    
# Figures for each patient
# STARTING FROM FOURTH FILE BECAUSE PATIENT 1 AND PATIENT 2 HAVE EMPTY VALUES
for i in range (3, 8):
    patient = loadmat(files[i])
    img1 = patient['Ktrans_image_12']
    img2 = patient['Ktrans_image_13']
    img3 = patient['Ktrans_image_22']
    img4 = patient['Ktrans_image_23']
    
    fig1 = go.Heatmap(z=img1, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
    fig.append_trace(fig1, 1, 1)
    
    fig2 = go.Heatmap(z=img2, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
    fig.append_trace(fig2, 1, 2)
    
    fig3 = go.Heatmap(z=img1, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
    fig.append_trace(fig3, 2, 1)
    
    fig4 = go.Heatmap(z=img2, zmin=Ktrans_lims[0], zmax=Ktrans_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'})
    fig.append_trace(fig4, 2, 2)

fig.update_yaxes(autorange="reversed")
# fig.update_xaxes(rangemode="normal")

            
fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.4,
        y=1.2,
        xanchor="left",
        yanchor="top",
        direction="down",
        type='dropdown',
        buttons=list(
            [
            dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 3',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 4',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 4',
                           'showlegend':False}]),
             dict(label = 'Patient 5',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 5',
                           'showlegend':False}]),
             dict(label = 'Patient 6',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 6',
                           'showlegend':False}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':False}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':False}]),
#              dict(label = 'Patient 7',
#                   method = 'update',
#                   args = [{'visible': [False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       True, True, True, True,
#                                       False, False, False, False]}, # the index of True aligns with the indices of plot traces
#                           {'title': 'Patient 7',
#                            'showlegend':False}]),
#              dict(label = 'Patient 8',
#                   method = 'update',
#                   args = [{'visible': [False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       True, True, True, True]}, # the index of True aligns with the indices of plot traces
#                           {'title': 'Patient 8',
#                            'showlegend':False}]),
            ]
            )
        )
    ])

fig.update_yaxes(showticklabels = False)
fig.update_xaxes(showticklabels = False)

# fig.show()

plot(fig, filename = 'figures/fig6-1.html', config = config)
display(HTML('figures/fig6-1.html'))

* $V_{e}$ <a name="fig-6-2"></a> 

files = ["fig6pat1.mat", 
         "fig6pat2.mat", 
         "fig6pat3.mat", 
         "fig6pat4.mat",
         "fig6pat5.mat", 
         "fig6pat6.mat", 
         "fig6pat7.mat", 
         "fig6pat8.mat"]
# Ve
Ve_lims = [0, 0.6]
patient = loadmat(files[2])

img1 = patient['Ve_image_12']
img2 = patient['Ve_image_12']
img3 = patient['Ve_image_22']
img4 = patient['Ve_image_23']

fig = make_subplots(rows=2, cols=2, column_widths=[0.5, 0.5], subplot_titles=('ETM', 'ETM', 'RRIFT', 'RRIFT'))

fig1 = go.Heatmap(z=img1, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                 colorbar={"title": 'V<sub>e</sub>'})
fig.append_trace(fig1, 1, 1)

fig2 = go.Heatmap(z=img2, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                 colorbar={"title": 'V<sub>e</sub>'})
fig.append_trace(fig2, 1, 2)

fig3 = go.Heatmap(z=img3, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                 colorbar={"title": 'V<sub>e</sub>'})
fig.append_trace(fig3, 2, 1)

fig4 = go.Heatmap(z=img4, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                 colorbar={"title": 'V<sub>e</sub>'})
fig.append_trace(fig4, 2, 2)

fig.update_layout(width = 700, height = 700, title = "Patient 3")
fig.update_yaxes(autorange="reversed")

# STARTING FROM THE FOURTH PATIENT BECAUSE 1 AND 2 HAVE EMPTY VALUES
for i in range (3, 8):
    patient = loadmat(files[i])

    img1 = patient['Ve_image_12']
    img2 = patient['Ve_image_12']
    img3 = patient['Ve_image_22']
    img4 = patient['Ve_image_23']
    
    fig1 = go.Heatmap(z=img1, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', visible=False,
                     colorbar={"title": 'V<sub>e</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
    fig.append_trace(fig1, 1, 1)

    fig2 = go.Heatmap(z=img2, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', visible=False,
                     colorbar={"title": 'V<sub>e</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
    fig.append_trace(fig2, 1, 2)

    fig3 = go.Heatmap(z=img3, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', visible=False,
                     colorbar={"title": 'V<sub>e</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
    fig.append_trace(fig3, 2, 1)

    fig4 = go.Heatmap(z=img4, zmin=Ve_lims[0], zmax=Ve_lims[1], colorscale='jet', visible=False,
                     colorbar={"title": 'V<sub>e</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
    fig.append_trace(fig4, 2, 2)
    

    
fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.4,
        y=1.2,
        xanchor="left",
        yanchor="top",
        direction="down",
        type='dropdown',
        buttons=list(
            [dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 3',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 4',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 4',
                           'showlegend':False}]),
             dict(label = 'Patient 5',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 5',
                           'showlegend':False}]),
             dict(label = 'Patient 6',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 6',
                           'showlegend':False}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':False}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':False}]),
#              dict(label = 'Patient 7',
#                   method = 'update',
#                   args = [{'visible': [False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       True, True, True, True,
#                                       False, False, False, False]}, # the index of True aligns with the indices of plot traces
#                           {'title': 'Patient 7',
#                            'showlegend':False}]),
#              dict(label = 'Patient 8',
#                   method = 'update',
#                   args = [{'visible': [False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       True, True, True, True]}, # the index of True aligns with the indices of plot traces
#                           {'title': 'Patient 8',
#                            'showlegend':False}]),
            ]
            )
        )
    ])

fig.update_layout(width = 700, height = 700)
fig.update_yaxes(autorange="reversed", showticklabels = False)
fig.update_xaxes(showticklabels = False)

# fig.show()
plot(fig, filename = 'figures/fig6-2.html', config = config)
display(HTML('figures/fig6-2.html'))

* $V_{p}$ <a name="fig-6-3"></a> 

files = ["fig6pat1.mat", 
        "fig6pat2.mat", 
        "fig6pat3.mat", 
        "fig6pat4.mat",
        "fig6pat5.mat", 
        "fig6pat6.mat", 
        "fig6pat7.mat", 
        "fig6pat8.mat"]

#Vp
Vp_lims = [0, 0.1]

patient = loadmat(files[2])

img1 = patient['Vp_image_12']
img2 = patient['Vp_image_12']
img3 = patient['Vp_image_22']
img4 = patient['Vp_image_23']

# Filter out NaN values
img1[np.isnan(img1)] = 0
img2[np.isnan(img2)] = 0
img3[np.isnan(img3)] = 0
img4[np.isnan(img4)] = 0

fig = make_subplots(rows=2, cols=2, column_widths=[0.5, 0.5], subplot_titles=('ETM', 'ETM', 'RRIFT', 'RRIFT'))

fig1 = go.Heatmap(z=img1, zmin=Vp_lims[0], zmax=Vp_lims[1], colorscale='jet', 
                  colorbar={"title": 'V<sub>p</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
fig.append_trace(fig1, 1, 1)

fig2 = go.Heatmap(z=img2, zmin=Vp_lims[0], zmax=Vp_lims[1], colorscale='jet', 
                  colorbar={"title": 'V<sub>p</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
fig.append_trace(fig2, 1, 2)

fig3 = go.Heatmap(z=img3, zmin=Vp_lims[0], zmax=Vp_lims[1], colorscale='jet', 
                  colorbar={"title": 'V<sub>p</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
fig.append_trace(fig3, 2, 1)

fig4 = go.Heatmap(z=img4, zmin=Vp_lims[0], zmax=Vp_lims[1], colorscale='jet', 
                  colorbar={"title": 'V<sub>p</sub>'}, hovertemplate="<br>".join(["z: %{z}<extra></extra>"]))
fig.append_trace(fig4, 2, 2)

# STARTING FROM PATIENT 4 BECAUSE PATIENT 1 AND 2 ARE EMPTY 
for i in range (3, 8):
    patient = loadmat(files[i])

    img1 = patient['Vp_image_12']
    img2 = patient['Vp_image_12']
    img3 = patient['Vp_image_22']
    img4 = patient['Vp_image_23']
    
    # Filter out NaN values
    img1[np.isnan(img1)] = 0
    img2[np.isnan(img2)] = 0
    img3[np.isnan(img3)] = 0
    img4[np.isnan(img4)] = 0

    fig1 = go.Heatmap(z=img1, zmin=Ve_lims[0], zmax=Vp_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'V<sub>p</sub>'})
    fig.append_trace(fig1, 1, 1)

    fig2 = go.Heatmap(z=img2, zmin=Ve_lims[0], zmax=Vp_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'V<sub>p</sub>'})
    fig.append_trace(fig2, 1, 2)

    fig3 = go.Heatmap(z=img3, zmin=Ve_lims[0], zmax=Vp_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'V<sub>p</sub>'})
    fig.append_trace(fig3, 2, 1)

    fig4 = go.Heatmap(z=img4, zmin=Ve_lims[0], zmax=Vp_lims[1], colorscale='jet', hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False, colorbar={"title": 'V<sub>p</sub>'})
    fig.append_trace(fig4, 2, 2)

    
fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.4,
        y=1.2,
        xanchor="left",
        yanchor="top",
        direction="down",
        type='dropdown',
        buttons=list(
            [dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 3',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 4',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 4',
                           'showlegend':False}]),
             dict(label = 'Patient 5',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 5',
                           'showlegend':False}]),
             dict(label = 'Patient 6',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 6',
                           'showlegend':False}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':False}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':False}]),
#              dict(label = 'Patient 7',
#                   method = 'update',
#                   args = [{'visible': [False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       True, True, True, True,
#                                       False, False, False, False]}, # the index of True aligns with the indices of plot traces
#                           {'title': 'Patient 7',
#                            'showlegend':False}]),
#              dict(label = 'Patient 8',
#                   method = 'update',
#                   args = [{'visible': [False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       False, False, False, False,
#                                       True, True, True, True]}, # the index of True aligns with the indices of plot traces
#                           {'title': 'Patient 8',
#                            'showlegend':False}]),
            ]
            )
        )
    ])

fig.update_layout(width = 700, height = 700)
fig.update_yaxes(autorange="reversed", showticklabels = False)
fig.update_xaxes(showticklabels = False)

plot(fig, filename = 'figures/fig6-3.html', config = config)
display(HTML('figures/fig6-3.html'))

## Figure 7 <a name="figure-7"></a>

fig = make_subplots(rows=1, cols=3)

files = ["./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-0185-1.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-0185-2.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-0185-3.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-0881-1.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-0881-2.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-1802-1.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-2570-1.mat",
         "./RRIFT/data/TCGA-GBM-Results/c02_postprocessed/TCGA-06-5417-1.mat"]

patient = loadmat(files[0])

t = np.squeeze(np.asarray((np.reshape(patient['t'], 70))))
Cp = np.squeeze(np.asarray((np.reshape(patient['Cp'], 70))))
Crr = np.squeeze(np.asarray((np.reshape(patient['Crr'], 70))))
num = np.squeeze(np.asarray((np.reshape(patient['num'], 38))))
denum = np.squeeze(np.asarray((np.reshape(patient['denum'], 38))))

labels = ['Blood plasma', 'Muscle', 'RRIFT fit', 'RRIFT']

fig1 = go.Scatter(
    name = labels[0],
    x=t,
    y=Cp,
    showlegend=False,
    visible=True,
    line_color='#e3988b',
)

fig.append_trace(fig1, 1, 1)
#     fig.append_trace(fig1, 1, 2)
#     fig.append_trace(fig1, 1, 3)

fig2 = go.Scatter(
    x=t,
    y=Crr,
    name = labels[1],
    showlegend=False,
    visible=True,
    line_color='#97d6d0',
)

#     fig.append_trace(fig2, 1, 1)
fig.append_trace(fig2, 1, 2)
#     fig.append_trace(fig2, 1, 3)

line_fit = sm.OLS(num ,sm.add_constant(denum)).fit().fittedvalues


fig3_1 = go.Scatter(x=denum, 
                         y=num, 
                         mode='markers',
                         marker_symbol="circle-open",
                         name=labels[2],
                         marker_size=14, 
                         marker_color='red', showlegend=False, visible=True)

#     fig.append_trace(fig3_1, 1, 1)
#     fig.append_trace(fig3_1, 1, 2)
fig.append_trace(fig3_1, 1, 3)

fig3_2 = go.Scatter(x=denum, 
                         y=line_fit,
                         name=labels[3],
                         mode='lines',
                         line_color="black", showlegend=False, visible=True)

fig.append_trace(fig3_2, 1, 3)

fig.update_layout(title_text='Patient 1') 
fig.update_xaxes(title_text='Time [min]',
                 range=[-1, 6], 
                 tickvals=[0, 2, 4, 6])

for i in range (1, 8):
    patient = loadmat(files[i])
    t = np.squeeze(np.asarray((np.reshape(patient['t'], 70))))
    Cp = np.squeeze(np.asarray((np.reshape(patient['Cp'], 70))))
    Crr = np.squeeze(np.asarray((np.reshape(patient['Crr'], 70))))
    num = np.squeeze(np.asarray((np.reshape(patient['num'], 38))))
    denum = np.squeeze(np.asarray((np.reshape(patient['denum'], 38))))
    
    fig1 = go.Scatter(
        x=t,
        y=Cp,
        showlegend=False,
        name=labels[0],
        visible=False,
        line_color='#e3988b',
    )
    
    fig.append_trace(fig1, 1, 1)
#     fig.append_trace(fig1, 1, 2)
#     fig.append_trace(fig1, 1, 3)

    fig2 = go.Scatter(
        x=t,
        y=Crr,
        name=labels[1],
        showlegend=False,
        visible=False,
        line_color='#97d6d0',
    )
    
#     fig.append_trace(fig2, 1, 1)
    fig.append_trace(fig2, 1, 2)
#     fig.append_trace(fig2, 1, 3)
    
    line_fit = sm.OLS(num ,sm.add_constant(denum)).fit().fittedvalues


    fig3_1 = go.Scatter(x=denum, 
                             y=num, 
                             mode='markers',
                             name=labels[2],
                             marker_symbol="circle-open", 
                             marker_size=14, 
                             marker_color='red', showlegend=False, visible=False)
    
#     fig.append_trace(fig3_1, 1, 1)
#     fig.append_trace(fig3_1, 1, 2)
    fig.append_trace(fig3_1, 1, 3)
    
    fig3_2 = go.Scatter(x=denum, 
                             y=line_fit,
                             name=labels[3],
                             mode='lines',
                             line_color="black", showlegend=False, visible=False)
    
    fig.append_trace(fig3_2, 1, 3)

fig.update_yaxes(title_text='Concentration [mM]',
                 row=1, col=1)
fig.update_yaxes(title_text='Concentration [mM]',
                 row=1, col=2)
fig.update_yaxes(title_text='Numerator',
                 row=1, col=3)
    
fig.update_layout(
    plot_bgcolor="#fff",
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.1,
        y=1.2,
        xanchor="left",
        yanchor="top",
        direction="right",
        type='buttons',
        buttons=list(
            [dict(label = 'Patient 1',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 1',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 2',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 2',
                           'showlegend':True}]),
             dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 3',
                           'showlegend':True}]),
             dict(label = 'Patient 4',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 4',
                           'showlegend':True}]),
             dict(label = 'Patient 5',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 5',
                           'showlegend':True}]),
             dict(label = 'Patient 6',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 6',
                           'showlegend':True}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':True}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':True}]),
            ]
            )
        )
    ])



fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')



# fig.show()

plot(fig, filename = 'figures/fig7.html', config = config)
display(HTML('figures/fig7.html'))

## Figure 8 <a name="figure-8"></a> 

fig = make_subplots(rows=1, cols=3, 
                    subplot_titles=('k<sub>ep,RR</sub> [min<sup>-1</sup>]', 'K<sub>RR</sub><sup>trans</sup>[min<sup>-1</sup>]', 'V<sub>e,RR</sub>'))

line_colour = '#010101'

file = loadmat('fig8vars.mat')
xv1 = np.squeeze(np.asarray(np.reshape(file['xv1'], 8)))
yv1 = np.squeeze(np.asarray(np.reshape(file['yv1'], 8)))
xr1 = file['xr1']
yr1 = file['yr1']

xv2 = np.squeeze(np.asarray(np.reshape(file['xv2'], 8)))
yv2 = np.squeeze(np.asarray(np.reshape(file['yv2'], 8)))
xr2 = file['xr2']
yr2 = file['yr2']

xv3 = np.squeeze(np.asarray(np.reshape(file['xv3'], 8)))
yv3 = np.squeeze(np.asarray(np.reshape(file['yv3'], 8)))
xr3 = file['xr3']
yr3 = file['yr3']

x = cycle([xv1, xv2, xv3])
y = cycle([yv1, yv2, yv3])



# First subplot vars

y_negative1 = np.squeeze(np.asarray(np.reshape(yr1[:,0], 8)))
y_positive1 = np.squeeze(np.asarray(np.reshape(yr1[:,1], 8)))

x_negative1 = np.squeeze(np.asarray(np.reshape(xr1[:,0], 8)))
x_positive1 = np.squeeze(np.asarray(np.reshape(xr1[:,1], 8)))

# Second subplot vars

y_negative2 = np.squeeze(np.asarray(np.reshape(yr2[:,0], 8)))
y_positive2 = np.squeeze(np.asarray(np.reshape(yr2[:,1], 8)))

x_negative2 = np.squeeze(np.asarray(np.reshape(xr2[:,0], 8)))
x_positive2 = np.squeeze(np.asarray(np.reshape(xr2[:,1], 8)))

# Third subplot vars

y_negative3 = np.squeeze(np.asarray(np.reshape(yr3[:,0], 8)))
y_positive3 = np.squeeze(np.asarray(np.reshape(yr3[:,1], 8)))

x_negative3 = np.squeeze(np.asarray(np.reshape(xr3[:,0], 8)))
x_positive3 = np.squeeze(np.asarray(np.reshape(xr3[:,1], 8)))

x_positive_values = cycle([x_positive1, x_positive2, x_positive3])
x_negative_values = cycle([x_negative1, x_negative2, x_negative3])
y_positive_values = cycle([y_positive1, y_positive2, y_positive3])
y_negative_values = cycle([y_negative1, y_negative2, y_negative3])


for i in range (1, 4):
    fig.add_trace(go.Scatter(
        x=next(x),
        y=next(y),
        showlegend=False,
        hoverinfo='x+y',
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=next(y_positive_values),
            arrayminus=next(y_negative_values),
            symmetric=False,
            visible=False,
            thickness=1),
        error_x=dict(
            type='data', # value of error bar given in data coordinates
            array=next(x_positive_values),
            arrayminus=next(x_negative_values),
            symmetric=False,
            visible=True,
            thickness=2),
        line_color=line_colour,
        mode='markers',
        marker_symbol='circle-open', marker_size=12, marker_line_width=2
    ),row=1, col=i)



# Adding the x=y dashed line on each of the three subplots

fig.append_trace(go.Scatter(x=[0,1], y=[0,1], showlegend=False, mode='lines', line={'dash': 'dash', 'color': '#000400', 'width': 0.5}), row=1, col=1)
fig.append_trace(go.Scatter(x=[0,1], y=[0,1], showlegend=False, mode='lines', line={'dash': 'dash', 'color': '#000400', 'width': 0.5}), row=1, col=2)
fig.append_trace(go.Scatter(x=[0,1], y=[0,1], showlegend=False, mode='lines', line={'dash': 'dash', 'color': '#000400', 'width': 0.5}), row=1, col=3)

# Adding the CCC text on each of the three subplots

fig.add_annotation(  
    x=0.5,
    y=0.14,
    text='CCC: 0.917', 
    showarrow=False, font=dict(size=40, color='black'), row=1, col=1)
fig.add_annotation(  
    x=0.12,
    y=0.013,
    text='CCC: 0.926',
    xref="paper", yref="paper", 
    showarrow=False, font=dict(size=40, color='black'), row=1, col=2)
fig.add_annotation(  
    x=0.22,
    y=0.113,
    text='CCC: 0.877',
    xref="paper", yref="paper", 
    showarrow=False, font=dict(size=40, color='black'), row=1, col=3)

fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_layout(plot_bgcolor="#fff")

fig.update_xaxes(title_text='ETM - k<sub>ep,RR</sub> [min<sup>-1</sup>]',
                 range=[0.1,0.8], 
                 tickvals=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
                 row=1, col=1)
fig.update_xaxes(title_text='ETM - K<sub>RR</sub><sup>trans</sup>[min<sup>-1</sup>]',
                 range=[0, 0.2], 
                 tickvals=[0, 0.1, 0.2],
                 row=1, col=2)
fig.update_xaxes(title_text='ETM - V<sub>e,RR</sub>',
                 range=[0.1, 0.3], 
                 tickvals=[0.1, 0.2, 0.3],
                 row=1, col=3)


fig.update_yaxes(title_text='RRIFT',
                 range=[0.1,0.8], 
                 tickvals=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
                 row=1, col=1)
fig.update_yaxes(title_text='RRIFT',
                 range=[0, 0.2], 
                 tickvals=[0, 0.1, 0.2],
                 row=1, col=2)
fig.update_yaxes(title_text='RRIFT',
                 range=[0.1, 0.3], 
                 tickvals=[0.1, 0.2, 0.3],
                 row=1, col=3)


# fig.show()

plot(fig, filename = 'figures/fig8.html', config = config)
display(HTML('figures/fig8.html'))

### Figure 8 - choose which subplot to view <a name="fig-8-1"></a> 


fig = go.Figure()

line_colour = '#010101'
titles = cycle(['k<sub>ep,RR</sub> [min<sup>-1</sup>]', 'K<sub>RR</sub><sup>trans</sup>[min<sup>-1</sup>]', 'V<sub>e,RR</sub>'])

file = loadmat('fig8vars.mat')
xv1 = np.squeeze(np.asarray(np.reshape(file['xv1'], 8)))
yv1 = np.squeeze(np.asarray(np.reshape(file['yv1'], 8)))
xr1 = file['xr1']
yr1 = file['yr1']

xv2 = np.squeeze(np.asarray(np.reshape(file['xv2'], 8)))
yv2 = np.squeeze(np.asarray(np.reshape(file['yv2'], 8)))
xr2 = file['xr2']
yr2 = file['yr2']

xv3 = np.squeeze(np.asarray(np.reshape(file['xv3'], 8)))
yv3 = np.squeeze(np.asarray(np.reshape(file['yv3'], 8)))
xr3 = file['xr3']
yr3 = file['yr3']

# First subplot vars

y_negative1 = np.squeeze(np.asarray(np.reshape(yr1[:,0], 8)))
y_positive1 = np.squeeze(np.asarray(np.reshape(yr1[:,1], 8)))

x_negative1 = np.squeeze(np.asarray(np.reshape(xr1[:,0], 8)))
x_positive1 = np.squeeze(np.asarray(np.reshape(xr1[:,1], 8)))

# Second subplot vars

y_negative2 = np.squeeze(np.asarray(np.reshape(yr2[:,0], 8)))
y_positive2 = np.squeeze(np.asarray(np.reshape(yr2[:,1], 8)))

x_negative2 = np.squeeze(np.asarray(np.reshape(xr2[:,0], 8)))
x_positive2 = np.squeeze(np.asarray(np.reshape(xr2[:,1], 8)))

# Third subplot vars

y_negative3 = np.squeeze(np.asarray(np.reshape(yr3[:,0], 8)))
y_positive3 = np.squeeze(np.asarray(np.reshape(yr3[:,1], 8)))

x_negative3 = np.squeeze(np.asarray(np.reshape(xr3[:,0], 8)))
x_positive3 = np.squeeze(np.asarray(np.reshape(xr3[:,1], 8)))

x_positive_values = cycle([x_positive1, x_positive2, x_positive3])
x_negative_values = cycle([x_negative1, x_negative2, x_negative3])
y_positive_values = cycle([y_positive1, y_positive2, y_positive3])
y_negative_values = cycle([y_negative1, y_negative2, y_negative3])

x = cycle([xv1, xv2, xv3])
y = cycle([yv1, yv2, yv3])


for i in range (1, 4):
    fig.add_trace(go.Scatter(
        name=next(titles),
        x=next(x),
        y=next(y),
        showlegend=False,
        hoverinfo='x+y',
        visible=False,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=next(y_positive_values),
            arrayminus=next(y_negative_values),
            symmetric=False,
            visible=False,
            thickness=1),
        error_x=dict(
            type='data', # value of error bar given in data coordinates
            array=next(x_positive_values),
            arrayminus=next(x_negative_values),
            symmetric=False,
            visible=True,
            thickness=2),
        line_color=line_colour,
        mode='markers',
        marker_symbol='circle-open', marker_size=12, marker_line_width=2
    ))

    
# Default value on the slider    
fig.data[0].visible = True
legend_values = cycle(['ETM - k<sub>ep,RR</sub> [min<sup>-1</sup>]',
                       'ETM - K<sup>trans</sup><sub>RR</sub>[min<sub>-1</sub>]', 
                       'ETM - V<sub>e,RR</sub>'])

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
             {'title.text': next(titles)}],
        label=next(legend_values)       
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=0,
    currentvalue={"prefix": "ETM - "},
    pad={"t": 70}, #upper padding 
    steps=steps
)]

fig.update_xaxes(title_text='ETM', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(title_text='RRIFT', ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_layout(plot_bgcolor="#fff", sliders=sliders, title={'text': next(titles)})




# fig.show()
plot(fig, filename = 'figures/fig8-1.html', config = config)
display(HTML('figures/fig8-1.html'))

## Figure 9 <a name="figure-9"></a>

* $K_{trans}[min^{-1}]$ <a name="fig-9-1"></a> 

files = ["fig9patient1.mat",
         "fig9patient2.mat",
         "fig9patient3.mat",
         "fig9patient7.mat",
         "fig9patient8.mat"]

patient = loadmat(files[0])

mapKtE1 = patient['mapKtE1']
mapKtR1 = patient['mapKtR1']
mapKtE2 = patient['mapKtE2']
mapKtR2 = patient['mapKtR2']

myS = patient['myS']

fig = make_subplots(rows=2, cols=2, column_widths=[0.5, 0.5], subplot_titles=('ETM', 'RRIFT', 'ETM', 'RRIFT'))

sub1 = go.Heatmap(z=mapKtE1[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15, colorscale='jet',
                  colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub1,1,1)

sub2 = go.Heatmap(z=mapKtR1[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub2,1,2)

sub3 = go.Heatmap(z=mapKtE2[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub3,2,1)

sub4 = go.Heatmap(z=mapKtR2[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub4,2,2)

for i in range(1,5):
    patient = loadmat(files[i])

    mapKtE1 = patient['mapKtE1']
    mapKtR1 = patient['mapKtR1']
    mapKtE2 = patient['mapKtE2']
    mapKtR2 = patient['mapKtR2']

    myS = patient['myS']

    sub1 = go.Heatmap(z=mapKtE1[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub1,1,1)
    
    sub2 = go.Heatmap(z=mapKtR1[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub2,1,2)
    
    sub3 = go.Heatmap(z=mapKtE2[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub3,2,1)
    
    sub4 = go.Heatmap(z=mapKtR2[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'K<sub>trans</sub>[min<sup>-1</sup>]'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub4,2,2)


fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.09,
        y=1.13,
        xanchor="left",
        yanchor="top",
        direction="right",
        type='dropdown',
        buttons=list(
            [dict(label = 'Patient 1',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 1',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 2',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 2',
                           'showlegend':True}]),
             dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 3',
                           'showlegend':True}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':True}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':True}]),
            ]
            )
        )
    ])

# fig.add_annotation(text=r'$\Delta t = 5.5 s$',
#                   xref="x", yref="y",
#                   x=37, y=4, showarrow=False, font = dict(size = 26, color="white"), row=1, col=1)
# fig.add_annotation(text=r'$\Delta t = 5.5 s$',
#                   xref="x", yref="y",
#                   x=37, y=4, showarrow=False, font = dict(size = 26, color="white"), row=1, col=2)
# fig.add_annotation(text=r'$\Delta t = 33 s$',
#                   xref="x", yref="y",
#                   x=37, y=4, showarrow=False, font = dict(size = 26, color="white"), row=2, col=1)
# fig.add_annotation(text=r'$\Delta t = 33 s$',
#                   xref="x", yref="y",
#                   x=37, y=4, showarrow=False, font = dict(size = 26, color="white"), row=2, col=2)

fig.update_layout(yaxis_autorange="reversed", width = 700, height = 700)
fig.update_layout(xaxis = dict(scaleanchor = 'x'))
fig.update_layout(yaxis = dict(scaleanchor = 'y'))

fig.update_layout(title_text="Patient 1", plot_bgcolor='rgba(0,0,0,0)')
fig.update_yaxes(autorange="reversed", showticklabels = False)
fig.update_xaxes(showticklabels = False)
# fig.show()

plot(fig, filename = 'figures/fig9-1.html', config = config)
display(HTML('figures/fig9-1.html'))

* $V_e$ <a name="fig-9-2"></a>

files = ["fig9patient1.mat",
         "fig9patient2.mat",
         "fig9patient3.mat",
         "fig9patient7.mat",
         "fig9patient8.mat"]

patient = loadmat(files[0])

mapVeE1 = patient['mapVeE1']
mapVeR1 = patient['mapVeR1']
mapVeE2 = patient['mapVeE2']
mapVeR2 = patient['mapVeR2']

myS = patient['myS']


fig = make_subplots(rows=2, cols=2, column_widths=[0.5, 0.5], subplot_titles=('ETM', 'RRIFT', 'ETM', 'RRIFT'))

sub1 = go.Heatmap(z=mapVeE1[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.5,
                  colorscale='jet',
                  colorbar={"title": 'Ve'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub1,1,1)

sub2 = go.Heatmap(z=mapVeR1[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.5,
                  colorscale='jet',
                  colorbar={"title": 'Ve'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub2,1,2)

sub3 = go.Heatmap(z=mapVeE2[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.5,
                  colorscale='jet',
                  colorbar={"title": 'Ve'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub3,2,1)

sub4 = go.Heatmap(z=mapVeR2[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.5,
                  colorscale='jet',
                  colorbar={"title": 'Ve'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub4,2,2)


for i in range(1,5):
    patient = loadmat(files[i])

    mapVeE1 = patient['mapVeE1']
    mapVeR1 = patient['mapVeR1']
    mapVeE2 = patient['mapVeE2']
    mapVeR2 = patient['mapVeR2']

    myS = patient['myS']

    sub1 = go.Heatmap(z=mapVeE1[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.5,
                      colorscale='jet',
                      colorbar={"title": 'Ve'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub1,1,1)
    
    sub2 = go.Heatmap(z=mapVeR1[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.5,
                      colorscale='jet',
                      colorbar={"title": 'Ve'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub2,1,2)
    
    sub3 = go.Heatmap(z=mapVeE2[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.5,
                      colorscale='jet',
                      colorbar={"title": 'Ve'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub3,2,1)
    
    sub4 = go.Heatmap(z=mapVeR2[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.5,
                      colorscale='jet',
                      colorbar={"title": 'Ve'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub4,2,2)

# fig.add_annotation(text=r'$\Delta t = 5.5 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=1, col=1)
# fig.add_annotation(text=r'$\Delta t = 5.5 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=1, col=2)
# fig.add_annotation(text=r'$\Delta t = 33 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=2, col=1)
# fig.add_annotation(text=r'$\Delta t = 33 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=2, col=2)
    
fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.09,
        y=1.13,
        xanchor="left",
        yanchor="top",
        direction="right",
        type='dropdown',
        buttons=list(
            [dict(label = 'Patient 1',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 1',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 2',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 2',
                           'showlegend':True}]),
             dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 3',
                           'showlegend':True}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':True}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':True}]),
            ]
            )
        )
    ])


fig.update_layout(title_text="Patient 1", yaxis_autorange="reversed", width = 700, height = 700)
fig.update_yaxes(autorange="reversed", showticklabels = False)
fig.update_xaxes(showticklabels = False)
fig.update_layout(xaxis = dict(scaleanchor = 'x'))
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')

# fig.show()
plot(fig, filename = 'figures/fig9-2.html', config = config)
display(HTML('figures/fig9-2.html'))

* $V_p$ <a name="fig-9-3"></a> 

files = ["fig9patient1.mat",
         "fig9patient2.mat",
         "fig9patient3.mat",
         "fig9patient7.mat",
         "fig9patient8.mat"]

patient = loadmat(files[0])

mapVpE1 = patient['mapVpE1']
mapVpR1 = patient['mapVpR1']
mapVpE2 = patient['mapVpE2']
mapVpR2 = patient['mapVpR2']

myS = patient['myS']


fig = make_subplots(rows=2, cols=2, column_widths=[0.5, 0.5], subplot_titles=('ETM', 'RRIFT', 'ETM', 'RRIFT'))

sub1 = go.Heatmap(z=mapVpE1[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'Vp'},
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  hoverinfo='x+y+z',
                  
                  visible=True)
fig.append_trace(sub1,1,1)

sub2 = go.Heatmap(z=mapVpR1[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'Vp'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub2,1,2)

sub3 = go.Heatmap(z=mapVpE2[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'Vp'},
                  hoverinfo='x+y+z',
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  visible=True)
fig.append_trace(sub3,2,1)

sub4 = go.Heatmap(z=mapVpR2[:,:,myS-1].squeeze(),
                  zmin=0, zmax=0.15,
                  colorscale='jet',
                  colorbar={"title": 'Vp'},
                  hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                  hoverinfo='x+y+z',
                  visible=True)
fig.append_trace(sub4,2,2)

for i in range(1,5):
    patient = loadmat(files[i])

    mapVpE1 = patient['mapVpE1']
    mapVpR1 = patient['mapVpR1']
    mapVpE2 = patient['mapVpE2']
    mapVpR2 = patient['mapVpR2']

    myS = patient['myS']

    sub1 = go.Heatmap(z=mapVpE1[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'Vp'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub1,1,1)
    
    sub2 = go.Heatmap(z=mapVpR1[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'Vp'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub2,1,2)
    
    sub3 = go.Heatmap(z=mapVpE2[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'Vp'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub3,2,1)
    
    sub4 = go.Heatmap(z=mapVpR2[:,:,myS-1].squeeze(),
                      zmin=0, zmax=0.15,
                      colorscale='jet',
                      colorbar={"title": 'Vp'},
                      hoverinfo='x+y+z',
                      hovertemplate="<br>".join(["z: %{z}<extra></extra>"]),
                      visible=False)
    fig.append_trace(sub4,2,2)

# fig.add_annotation(text=r'$\Delta t = 5.5 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=1, col=1)
# fig.add_annotation(text=r'$\Delta t = 5.5 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=1, col=2)
# fig.add_annotation(text=r'$\Delta t = 33 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=2, col=1)
# fig.add_annotation(text=r'$\Delta t = 33 s$',
#                   xref="x", yref="y",
#                   x=32, y=56, showarrow=False, font = dict(size = 26, color="white"), row=2, col=2)
    
fig.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        x=0.09,
        y=1.13,
        xanchor="left",
        yanchor="top",
        direction="right",
        type='dropdown',
        buttons=list(
            [dict(label = 'Patient 1',
                  method = 'update',
                  args = [{'visible': [True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]},
                          {'title': 'Patient 1',
                           'showlegend':True, 'visible': True}]),
             dict(label = 'Patient 2',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 2',
                           'showlegend':True}]),
             dict(label = 'Patient 3',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 3',
                           'showlegend':True}]),
             dict(label = 'Patient 7',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True,
                                      False, False, False, False]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 7',
                           'showlegend':True}]),
             dict(label = 'Patient 8',
                  method = 'update',
                  args = [{'visible': [False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      False, False, False, False,
                                      True, True, True, True]}, # the index of True aligns with the indices of plot traces
                          {'title': 'Patient 8',
                           'showlegend':True}]),
            ]
            )
        )
    ])


fig.update_layout(yaxis_autorange="reversed", title_text="Patient 1", width = 700, height = 700)
fig.update_yaxes(autorange="reversed", showticklabels = False)
fig.update_layout(xaxis = dict(scaleanchor = 'x'))
fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
fig.update_xaxes(showticklabels = False)

# fig.show()

plot(fig, filename = 'figures/fig9-3.html', config = config)
display(HTML('figures/fig9-3.html'))

## Figure 10 <a name="figure-10"></a> 

fig = make_subplots(rows=1, cols=3, 
                    subplot_titles=('K<sup>trans</sup>', 'V<sub>e</sub>', 'V<sub>p</sub>'))


line_colours = ['#585e9a', '#f36a61']
line_names = ['RRIFT', 'ETM']

file = loadmat('fig10vars.mat')
TRes = np.squeeze(np.asarray(np.reshape(file['TRes'], 10)))

errQt1 = file['errQt1']
errMd1 = np.squeeze(np.asarray(np.reshape(file['errMd1'], 10)))

errQt2 = file['errQt2']
errMd2 = np.squeeze(np.asarray(np.reshape(file['errMd2'], 10)))

errQt3 = file['errQt3']
errMd3 = np.squeeze(np.asarray(np.reshape(file['errMd3'], 10)))

errQt4 = file['errQt4']
errMd4 = np.squeeze(np.asarray(np.reshape(file['errMd4'], 10)))

errQt5 = file['errQt5']
errMd5 = np.squeeze(np.asarray(np.reshape(file['errMd5'], 10)))

errQt6 = file['errQt6']
errMd6 = np.squeeze(np.asarray(np.reshape(file['errMd6'], 10)))


# First subplot

err_negative = np.squeeze(np.asarray(np.reshape(abs(errQt1[0,:]-errMd1), 10)))
err_positive = np.squeeze(np.asarray(np.reshape(abs(errQt1[1,:]-errMd1), 10)))

fig.append_trace(go.Scatter(
        name=line_names[0],
        x=TRes,
        y=errMd1,
        showlegend=False,
        legendgroup=1,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=err_positive,
            arrayminus=err_negative,
            symmetric=False,
            thickness=4, width=9),
        line_color=line_colours[0],
    ),row=1, col=1)

err_negative = np.squeeze(np.asarray(np.reshape(abs(errQt2[0,:]-errMd2), 10)))
err_positive = np.squeeze(np.asarray(np.reshape(abs(errQt2[1,:]-errMd2), 10)))

fig.append_trace(go.Scatter(
        name=line_names[1],
        x=TRes,
        y=errMd2,
        showlegend=False,
        legendgroup=2,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=err_positive,
            arrayminus=err_negative,
            symmetric=False,
            thickness=4, width=9),
        line_color=line_colours[1],
    ),row=1, col=1)



# Second subplot 

err_negative = np.squeeze(np.asarray(np.reshape(abs(errQt3[0,:]-errMd3), 10)))
err_positive = np.squeeze(np.asarray(np.reshape(abs(errQt3[1,:]-errMd3), 10)))

fig.append_trace(go.Scatter(
        name=line_names[0],
        x=TRes,
        y=errMd3,
        showlegend=True,
        legendgroup=1,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=err_positive,
            arrayminus=err_negative,
            symmetric=False,
            thickness=4, width=9),
        line_color=line_colours[0],
    ),row=1, col=2)

err_negative = np.squeeze(np.asarray(np.reshape(abs(errQt4[0,:]-errMd4), 10)))
err_positive = np.squeeze(np.asarray(np.reshape(abs(errQt4[1,:]-errMd4), 10)))

fig.append_trace(go.Scatter(
        name=line_names[1],
        x=TRes,
        y=errMd4,
        showlegend=True,
        legendgroup=2,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=err_positive,
            arrayminus=err_negative,
            symmetric=False,
            thickness=4, width=9),
        line_color=line_colours[1],
    ),row=1, col=2)

# Third subplot

err_negative = np.squeeze(np.asarray(np.reshape(abs(errQt5[0,:]-errMd5), 10)))
err_positive = np.squeeze(np.asarray(np.reshape(abs(errQt5[1,:]-errMd5), 10)))

fig.append_trace(go.Scatter(
        name=line_names[0],
        x=TRes,
        y=errMd5,
        showlegend=False,
        legendgroup=1,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=err_positive,
            arrayminus=err_negative,
            symmetric=False,
            thickness=4, width=9),
        line_color=line_colours[0],
    ),row=1, col=3)

err_negative = np.squeeze(np.asarray(np.reshape(abs(errQt6[0,:]-errMd6), 10)))
err_positive = np.squeeze(np.asarray(np.reshape(abs(errQt6[1,:]-errMd6), 10)))

fig.append_trace(go.Scatter(
        name=line_names[1],
        x=TRes,
        y=errMd6,
        showlegend=False,
        legendgroup=2,
        error_y=dict(
            type='data', # value of error bar given in data coordinates
            array=err_positive,
            arrayminus=err_negative,
            symmetric=False,
            thickness=4, width=9),
        line_color=line_colours[1]
    ),row=1, col=3)


fig.add_hline(y=0, line_dash="dot", line_color="black", line_width=1.69, col="all")

fig.update_xaxes(title_text='Temporal Resolutions [s]',
                 range=[0,50], 
                 tickvals=[0, 10, 20, 30, 40, 50],
                 row=1, col=1)
fig.update_xaxes(title_text='Temporal Resolutions [s]',
                 range=[0, 50], 
                 tickvals=[0, 10, 20, 30, 40, 50],
                 row=1, col=2)
fig.update_xaxes(title_text='Temporal Resolutions [s]',
                 range=[0, 50], 
                 tickvals=[0, 10, 20, 30, 40, 50],
                 row=1, col=3)


fig.update_yaxes(title_text='Percent Change',
                 range=[-100, 100], 
                 tickvals=[-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100],
                 row=1, col=1)
fig.update_yaxes(title_text='Percent Change',
                 range=[-100, 100], 
                 tickvals=[-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100],
                 row=1, col=2)
fig.update_yaxes(title_text='Percent Change',
                 range=[-100, 100], 
                 tickvals=[-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100],
                 row=1, col=3)


fig.update_layout(plot_bgcolor='#ffffff')
fig.update_xaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')
fig.update_yaxes(ticks="outside", showline=True, linewidth=2, linecolor='black')


    
# fig.show()
plot(fig, filename = 'figures/fig10.html', config = config)
display(HTML('figures/fig10.html'))