## NOTE: to run the supplementary material, use Jupyter notebook. See README

## Topopt GUI for jupyter notebook
##
## Part of supplementary material of "Overhang control in topology optimization: 
## a comparison of continuous front propagation-based and discrete layer-by-layer
## overhang control", E. van de Ven, R. Maas, C. Ayas, M. Langelaar, F. van Keulen,
## 2019
##
## Disclaimer:                                                              
## The author reserves all rights but does not guarantee that the code is   
## free from errors. Furthermore, the author shall not be liable in any     
## event caused by the use of the program.                                  
##
## Code by Emiel van de Ven, 2020
## emiel@emielvandeven.nl

import topopt
import threading
import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display

# inputs for the optimization: [nelx, nely, filter_rad, volfrac, penal, movelimit]
opt_controls = [widgets.BoundedIntText(value=100,min=1,max=10000,step=1,description='nelx:'),
widgets.BoundedIntText(value=50,min=2,max=10000,step=1,description='nely:'),
widgets.FloatSlider(value=1.2,min=1,max=5,step=0.1,description='Filter radius:'),
widgets.FloatSlider(value=0.5,min=0.1,max=1,step=0.1,description='Volume Frac:'),
widgets.FloatSlider(value=3,min=1,max=8,step=0.1,description='Penalization:'),
widgets.FloatSlider(value=0.05,min=0.01,max=0.3,step=0.001,description='OC movelim.:')]

# inputs for the AM Filter: [Filter type, continuation, k, v_void]
am_controls = [widgets.Dropdown(options=['none','layer-by-layer', 'front-propagation', 'front-propagation-improved'],value='front-propagation-improved',description='AMFilter:'),
              widgets.Checkbox(value=False,description='Continuation on first 10 iters'),
              widgets.FloatSlider(value=2,min=0,max=10,step=0.01,description='k:'),
              widgets.FloatSlider(value=0.25,min=0,max=0.8,step=0.01,description='v_void:')]
am_controls[3].layout.visibility = 'hidden'

# button to start/stop optimization
w_run = widgets.ToggleButton(value=False,description='Start Optimization',
                             layout=Layout(margin='20px 0px 20px 50px',width='50%',height='50px'))

# output label to display optimization stats
w_out = widgets.Label(value="",layout=Layout(fontSize='50px'))

# Group controls together in VBox/HBox and define layout
w_fc = widgets.VBox([widgets.HTML(value="<center><b>Fixed controls:</b></center>"),
                     opt_controls[0],opt_controls[1],opt_controls[2]],layout=Layout(margin='0px 50px'))
w_dc = widgets.VBox([widgets.HTML(value="<center><b>Dynamic controls (can be changed during optimizaton):</b></center>"),
                     opt_controls[3],opt_controls[4],opt_controls[5]])
w_am1 = widgets.VBox([am_controls[0],am_controls[1]],layout=Layout(margin='0px 50px'))       
w_am2 = widgets.VBox([widgets.HTML(value="<center><b>Front propagation specific controls:</b></center>"),
                       am_controls[2],am_controls[3]])
controls = widgets.VBox([widgets.HTML(value="<h4><b>Optimization controls:</b></h4>"),
                         widgets.HBox([w_fc, w_dc]),
                         widgets.HTML(value="<h4><b>AM Filter controls:</b></h4>"),
                         widgets.HBox([w_am1, w_am2]),
                         w_run,
                         w_out])

display(controls)

# callback function when run/stop button is toggled
def run_opt(b):
    if w_run.value==True:
        # optimization is started in seperate thread to allow updating of sliders during optimization
        w_run.description = 'Stop Optimization'
        thread = threading.Thread(target=topopt.main, args=(opt_controls,am_controls,w_run,w_out))
        thread.start()
    else: 
        w_run.description = 'Restart Optimization'

# callback function for AMFilter type dropdown. This function hides/shows the relevant options for the selected filter
def fp_controls(b):
    if am_controls[0].value=='none':
        am_controls[1].layout.visibility = 'hidden'
    else:
        am_controls[1].layout.visibility = 'visible'

    if am_controls[0].value=='front-propagation':
        am_controls[3].layout.visibility = 'visible'
    else:
        am_controls[3].layout.visibility = 'hidden'

    if am_controls[0].value=='front-propagation' or am_controls[0].value=='front-propagation-improved':
        w_am2.layout.visibility = 'visible'
    else:
        w_am2.layout.visibility = 'hidden'

# attach callback functions
w_run.observe(run_opt, 'value')
am_controls[0].observe(fp_controls, 'value')
