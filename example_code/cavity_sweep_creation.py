import phidl
import phidl.geometry as pg
from phidl import quickplot as qp
import cavity_functions as cf
import numpy as np
from phidl import set_quickplot_options 

set_quickplot_options(show_ports=True, show_subports=True,new_window=True,blocking=True)

# Goal of this is to make cleaner code for making gds's, our old code is a bit nasty.

global z
z=0
def cavity(lat_const,lat_def,hx):
    global z
    z += 1
    hx_def = lat_def
    pcc_params = {
        'a': lat_const,
        'cX': hx,
        'cY': 1.4,
        'cW': 2,
        'cH': .5,
        'dA': lat_def,
        'dY': 0,
        'dX': hx_def,
        
        'nleft': 20,
        'nright': 7,
        'ndef': 6,
        'hole_type': 'ellipse',
        'defect_type': 'cubic',
        'min_dim': 1e-7
    }
    cav = cf.make_cavity(pcc_params)
    # text = pg.text(str(z),size= 4).move(destination=(0,5))
    # return pg.boolean(cav,text,'or')
    return cav

lat_const_sweep = [.275]
lat_def_sweep = np.linspace(.16,.19,3)
hx_sweep = np.linspace(.65,.7,3)
E = pg.gridsweep(
    function = cavity,
    param_x= {'lat_const': lat_const_sweep},
    param_y = {'lat_def' : lat_def_sweep,
               'hx': hx_sweep},
    spacing = (5,5))

cross = pg.rectangle((3,1)).move(destination=(-2,1))
cross1 = pg.rectangle((3,1)).rotate(angle=90)
cross = pg.boolean(cross,cross1,'or')
E.add_ref(cross).move(destination = (-10,34.5))
# E.add_ref(cross).move(destination = (141,34.5))
# pg.packer(E,spacing=5)
E.write_gds('thinfilmPhC_final.gds')

#This is really cool!
#We essentially now have made 'building blocks' that we can use to make more complex structures
#What would be cool to build up in phidl.
#1. Have one layer that is 'designed' parameters and one layer is the 'compensated' parameters
# that we actually make the beamer file with.
    # From here, once we fabricate and can roughly approximate actual fabbed parameters with 
    # Katie's fitting, we can compare the 3 different stages of beam profiles, nominally we want
    # the fabricated layer to overlap with the designed layer, so make changes to compensation factors
    # to achieve this

#2. Way to incorporate this with wvgsolver. I.e. if I make a sweep of designed parameters, I want to study
# the designed sensitivtivty + expected wavelength range. I can do this by feeding the designed parameters for
# each device into wvgsolver and extract lambda, Q, etc. This would be really cool

#3. We've made simple cavities, but take advantage of the nice phidl feature for doing electrode routing for 
# strain tuning + other things. It should be way simpler