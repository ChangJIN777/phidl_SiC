import phidl.geometry as pg
from phidl import Device, CrossSection
import phidl.path as pp
import numpy as np
from phidl import quickplot as qp # Rename "quickplot()" to the easier "qp()"
from phidl import set_quickplot_options 
import re

set_quickplot_options(show_ports=True, show_subports=True,new_window=True,blocking=True)

def make_hole(hole_params,layer):
    a = hole_params['a']
    hx = hole_params['hx']
    hy = hole_params['hy']
    w = hole_params['w']
    xpos = hole_params['xpos']
    hole_type = hole_params['hole_type']
    if hole_type == 'ellipse':
        hole = pg.ellipse((hx/2,hy/2)).move(origin=(0,0),destination=(xpos,0))
        return hole
    elif hole_type == 'anvil':
        r1 = hole_params['r2']
        hole = Device('rect')
        xpts = (-hx/2,0,hx/2, 0)
        ypts = (0, hy/2, 0, -hy/2)
        poly1 = hole.add_polygon([xpts, ypts], layer = 0)
        poly1.fillet(r1)
        hole.move(destination=(xpos,0))
        return hole    


def cubic_defect(i, nleft, nright, ndef, maxdef):
    if nleft < i < nleft + 2 * ndef - 1:
        x = np.abs(i - nleft - ndef)
        return 1 - maxdef * (2 * (x / ndef) ** 3 - 3 * (x / ndef) ** 2 + 1)
    else:
        return 1


def quadratic_defect(i, nleft, nright, ndef, maxdef):
    if nleft < i < nleft + 2 * ndef - 1:
        x = np.abs(i - nleft - ndef)
        return 1 - maxdef * (2 * (x / ndef) ** 3 - 3 * (x / ndef) ** 2 + 1)
    else:
        return 1


def generate_cavity_and_check(pcc_params):
    ndef = pcc_params['ndef']
    nleft = pcc_params['nleft']
    nright = pcc_params['nright']
    ntot = nleft + 2 * ndef + nright
    hole_index = np.arange(ntot)

    if pcc_params['defect_type'] == 'cubic':
        defect = cubic_defect
    elif pcc_params['defect_type'] == 'quadratic':
        defect = quadratic_defect

    a = pcc_params['a']
    lat_def = pcc_params['dA']
    lat_def_list = np.array([defect(i, nleft, nright, ndef, lat_def) for i in hole_index])
    lat_list = a * lat_def_list

    hx_scale = pcc_params['cX']
    hx_def = pcc_params['dX']
    hx_def_list = np.array([defect(i, nleft, nright, ndef, hx_def) for i in hole_index])
    hx_list = hx_scale * a * hx_def_list

    hy_scale = pcc_params['cY']
    hy_def = pcc_params['dY']
    hy_def_list = np.array([defect(i, nleft, nright, ndef, hy_def) for i in hole_index])
    hy_list = hy_scale * a * hy_def_list

    # Check fabrication feasibility
    # Check mirror cell
    w0 = pcc_params['cW']*a
    a_mirr = a
    hx_mirr = hx_list[0]
    hy_mirr = hy_list[0]

    a_def = lat_list[nleft + ndef]
    hx_def = hx_list[nleft + ndef]
    hy_def = hy_list[nleft + ndef]
    min_dim = pcc_params['min_dim']
    flag = False
    if a_mirr-hx_mirr < min_dim or hx_mirr < min_dim or (w0 - hy_mirr) < min_dim or hy_mirr < min_dim \
            or (a_def-hx_def) < min_dim or hx_def < min_dim or (w0 - hy_def) < min_dim or hy_def < min_dim:
        # If any of the cavity hole dimensions drop below a minimum fabrication tolerance, throw it out
        flag = True

    xc = np.sum(lat_list[:(nleft+ndef)])
    xpos = [x - xc for x in np.cumsum(lat_list)]
    return [lat_list, hx_list, hy_list, xpos, flag]


def build_holes(pcc_params):
    #Put this for now, but this might not be the way we make it

    a_list, hx_list, hy_list, xpos_list, flag = generate_cavity_and_check(pcc_params)
    a = pcc_params['a']
    w0 = pcc_params['cW'] * a
    hole_list = Device('rect')

    #This is for 300 pA EBL, a nice thing to do would be to make a dictionary of 
    # fabrication compensation numbers and calibrate.
   
    hy_comp = -.03 
    hx_comp = -.03

    if flag:
        return None, flag

    for (lat_const, hx, hy,xpos) in zip(a_list,hx_list,hy_list,xpos_list):
        # Take each unit cell and add it to our model. First update the unit cell parameters according to our cell list
        hole_params = {
            'a': lat_const,
            'hx': hx+hx_comp,
            'hy': hy+hy_comp,
            'w': w0,
            'xpos': xpos,
            'ref_index': 2.4028,
            'hole_type': pcc_params['hole_type']
        }

        hole = make_hole(hole_params,layer=0)
        hole_list.add(hole)

    # This DOES NOT put the cavity center at the (0,0) like we want.
    xL = xpos_list[0]-a_list[0]/2
    xR = xpos_list[-1]+a_list[-1]/2

    return hole_list,xL,xR

def make_tether_1(trenchW, tethL, tethW,w):
    tether1 = pg.taper(tethL/2,width1 = w, width2 = tethW)
    tether2 = pg.taper(tethL/2,width1=tethW,width2=w).move(destination=(tethL/2,0))
    tether = pg.boolean(tether1,tether2,'or')

    vert_teth = Device()
    top_block = pg.taper(.2*trenchW,width1=2,width2 = .3)
    top_block.add_port('top_block',(.2*trenchW,0), orientation= 0)

    main_block = pg.rectangle((.6*trenchW,.3))
    main_block.add_port('main_block_left',(0,.3/2), orientation= 180)
    main_block.add_port('main_block_right',(.6*trenchW,.3/2), orientation= 0)

    tri_ref1 = vert_teth.add_ref(top_block)
    tri_ref2 = vert_teth.add_ref(top_block)
    main_ref = vert_teth.add_ref(main_block)

    main_ref.connect('main_block_left',destination=tri_ref1.ports['top_block'])
    tri_ref2.connect('top_block',destination=main_ref.ports['main_block_right'])
    vert_teth.rotate(90).move(destination=(tethL/2,-trenchW/2))
    tether = pg.boolean(vert_teth,tether,'or')
    tether.add_port('tether_left',(0,0), orientation= 180)
    tether.add_port('tether_right',(tethL,0), orientation= 0)
    return tether

def make_tether_2(trenchW, tethL, tethW,w):
    def my_custom_width_fun(t):
        # Note: Custom width/offset functions MUST be vectorizable--you must be able
        # to call them with an array input like my_custom_width_fun([0, 0.1, 0.2, 0.3, 0.
        #˓→4])
        sigma = 3
        gauss_width = np.exp(-sigma*(t-1/2)**2)
        return gauss_width
    # Create the Path
    P = pp.straight(length = tethL)
    # Create two cross-sections: one fixed width, one modulated by my_custom_offset_fun
    X = CrossSection()
    X.add(width = w,                   offset = 0, layer = 0)
    X.add(width = my_custom_width_fun, offset = 0,  layer = 0)
    # Extrude the Path to create the Device
    tether = P.extrude(X)


    vert_teth = Device()
    top_block = pg.taper(.2*trenchW,width1=2,width2 = .3)
    top_block.add_port('top_block',(.2*trenchW,0), orientation= 0)

    main_block = pg.rectangle((.6*trenchW,.3))
    main_block.add_port('main_block_left',(0,.3/2), orientation= 180)
    main_block.add_port('main_block_right',(.6*trenchW,.3/2), orientation= 0)

    tri_ref1 = vert_teth.add_ref(top_block)
    tri_ref2 = vert_teth.add_ref(top_block)
    main_ref = vert_teth.add_ref(main_block)

    main_ref.connect('main_block_left',destination=tri_ref1.ports['top_block'])
    tri_ref2.connect('top_block',destination=main_ref.ports['main_block_right'])
    vert_teth.rotate(90).move(destination=(tethL/2,-trenchW/2))
    tether = pg.boolean(vert_teth,tether,'or')
    tether.add_port('tether_left',(0,0), orientation= 180)
    tether.add_port('tether_right',(tethL,0), orientation= 0)
    return tether

def make_cavity(pcc_params):

    w = pcc_params['cW']*pcc_params['a']
    trenchW = 6
    wvg_length = 5
    teth_length = 10
    teth_width = 1.3*w

    cav_holes,xL,xR = build_holes(pcc_params)
    
    init_beam = pg.rectangle((xR+abs(xL)+1,w),layer = 0).move(destination=(xL-1,-w/2))
    cav = pg.xor_diff(init_beam,cav_holes)
    cav.add_port('cav_wvg_port',midpoint=(xR,0),width = w, orientation= 0)
    cav.write_gds('test_hole.gds')

    wvg_section = pg.rectangle((wvg_length,w))
    wvg_section.add_port('wvg_left',(0,w/2), orientation= 180)
    wvg_section.add_port('wvg_right',(wvg_length,w/2), orientation= 0)

    tether = make_tether_2(trenchW,teth_length,teth_width,w)

    tap_wvg = pg.taper(30,width1=w,width2=.06)
    tap_wvg.add_port('tap_wvg_left',(0,0), orientation= 180)

    E = Device('full_device')
    cav_ref = E.add_ref(cav)
    wvg_ref1 = E.add_ref(wvg_section)
    wvg_ref1.connect('wvg_left', destination=cav_ref.ports['cav_wvg_port'])

    teth_ref1 = E.add_ref(tether)
    teth_ref1.connect('tether_left', destination=wvg_ref1.ports['wvg_right'])

    # wvg_ref2 = E.add_ref(wvg_section)
    # wvg_ref2.connect('wvg_left', destination=teth_ref1.ports['tether_right'])

    # teth_ref2 = E.add_ref(tether)
    # teth_ref2.connect('tether_left', destination=wvg_ref2.ports['wvg_right'])

    tap_ref = E.add_ref(tap_wvg)
    tap_ref.connect('tap_wvg_left', destination=teth_ref1.ports['tether_right'])
    full_dev_width = 5 + abs(xL)+xR + teth_length+wvg_length + 30 + 10
    trench_rect = pg.rectangle((full_dev_width,trenchW)).move((xL-1,-trenchW/2))

    dev = pg.xor_diff(trench_rect,E)
    dev_text = pg.text(str(int(pcc_params['dA']*100)) + '_' + str(int(pcc_params['a']*1e3)) + '_'+str(int(pcc_params['cX']*100) ),size = 4).movey(5)
    # dev = pg.boolean(dev,dev_text,'or')
    return dev

def import_GDS(fullfilepath):
    """this function import phidl geometry from existing GDS files 

    Args:
        fullfilepath (string): the full file path of the GDS file 

    Returns:
        phidl device: the device imported from the gds file 
    """
    # importing the geometry from existing GDS files 
    c = pg.import_gds(fullfilepath, cellname = None, flatten = False)
    return c 

def make_cavity_v2(fullfilepath):
    """_summary_

    Args:
        gds_file_name (string): the file name of the cavity gds generated by waveguide solver

    Returns:
        _type_: _description_
    """
    
    # import the device from the desired gds file 
    cav = import_GDS(fullfilepath)
    # get the height and the width of the imported device 
    cavity_width = cav.ymax-cav.ymin 
    cavity_length = cav.xmax-cav.xmin
    xR = cav.xmax
    xL = cav.xmin
    cav.add_port('cav_wvg_port',midpoint=(xR,0),width = cavity_width, orientation= 0)
    
    trenchW = 6
    wvg_length = 5
    teth_length = 10
    teth_width = 1.3*cavity_width

    wvg_section = pg.rectangle((wvg_length,cavity_width))
    wvg_section.add_port('wvg_left',(0,cavity_width/2), orientation= 180)
    wvg_section.add_port('wvg_right',(wvg_length,cavity_width/2), orientation= 0)

    tether = make_tether_2(trenchW,teth_length,teth_width,cavity_width)

    # make the taper waveguide region
    tap_wvg = pg.taper(30,width1=cavity_width,width2=.06)
    tap_wvg.add_port('tap_wvg_left',(0,0), orientation= 180)

    E = Device('full_device')
    cav_ref = E.add_ref(cav)
    wvg_ref1 = E.add_ref(wvg_section)
    wvg_ref1.connect('wvg_left', destination=cav_ref.ports['cav_wvg_port'])

    teth_ref1 = E.add_ref(tether)
    teth_ref1.connect('tether_left', destination=wvg_ref1.ports['wvg_right'])

    tap_ref = E.add_ref(tap_wvg)
    tap_ref.connect('tap_wvg_left', destination=teth_ref1.ports['tether_right'])
    full_dev_width = 5 + abs(xL)+xR + teth_length+wvg_length + 30 + 10
    trench_rect = pg.rectangle((full_dev_width,trenchW)).move((xL,-trenchW/2))

    dev = pg.xor_diff(trench_rect,E)
    # dev = pg.boolean(dev,dev_text,'or')
    return dev

def add_alignment_markers(pad):
    
    xmin = pad.xmin
    xmax = pad.xmax
    ymin = pad.ymin
    ymax = pad.ymax
    
    xmin = xmin - 2
    xmax = xmax + 2
    ymin = ymin - 2
    ymax = ymax + 2
    
    xpts1 = (xmin, xmin+1, xmin+1, xmin+1, xmin+2)
    ypts1 = (ymin, ymin+1, ymin, ymin-1, ymin)
    poly1 = pad.add_polygon( [xpts1, ypts1], layer = 0)

    xpts2 = (xmax, 1+xmax, 1+xmax, 1+xmax, 2+xmax)
    ypts2 = (ymin, ymin+1, ymin, ymin-1, ymin)
    poly2 = pad.add_polygon( [xpts2, ypts2], layer = 0)


    xpts3 = (xmax, 1+xmax, 1+xmax, 1+xmax, 2+xmax)
    ypts3 = (ymax, 1+ymax, ymax, -1+ymax, ymax)
    poly3 = pad.add_polygon( [xpts3, ypts3], layer = 0)

    xpts4 = (xmin, xmin+1, xmin+1, xmin+1, xmin+2)
    ypts4 = (ymax, 1+ymax, 0+ymax, -1+ymax, ymax)
    poly4 = pad.add_polygon( [xpts4, ypts4], layer = 0)

    return 

def make_cavity_v3(fullfilepath,waveguide_width):
    """_summary_

    Args:
        gds_file_name (string): the file name of the cavity gds generated by waveguide solver
        waveguide_width (float with unit um): the end width of the tapping waveguide region
    Returns:
        _type_: _description_
    """
    
    # import the device from the desired gds file 
    cav = import_GDS(fullfilepath)
    # get the height and the width of the imported device 
    cavity_width = cav.ymax-cav.ymin 
    cavity_length = cav.xmax-cav.xmin
    xR = cav.xmax
    xL = cav.xmin
    cav.add_port('cav_wvg_port',midpoint=(xR,0),width = cavity_width, orientation= 0)
    
    trenchW = 6
    wvg_length = 5
    teth_length = 10
    teth_width = 1.3*cavity_width

    wvg_section = pg.rectangle((wvg_length,cavity_width))
    wvg_section.add_port('wvg_left',(0,cavity_width/2), orientation= 180)
    wvg_section.add_port('wvg_right',(wvg_length,cavity_width/2), orientation= 0)

    tether = make_tether_2(trenchW,teth_length,teth_width,cavity_width)

    # make the taper waveguide region
    tap_wvg = pg.taper(30,width1=cavity_width,width2=waveguide_width)
    tap_wvg.add_port('tap_wvg_left',(0,0), orientation= 180)

    E = Device('full_device')
    cav_ref = E.add_ref(cav)
    wvg_ref1 = E.add_ref(wvg_section)
    wvg_ref1.connect('wvg_left', destination=cav_ref.ports['cav_wvg_port'])

    teth_ref1 = E.add_ref(tether)
    teth_ref1.connect('tether_left', destination=wvg_ref1.ports['wvg_right'])

    tap_ref = E.add_ref(tap_wvg)
    tap_ref.connect('tap_wvg_left', destination=teth_ref1.ports['tether_right'])
    full_dev_width = 5 + abs(xL)+xR + teth_length+wvg_length + 30 + 10
    trench_rect = pg.rectangle((full_dev_width,trenchW)).move((xL,-trenchW/2))

    dev = pg.xor_diff(trench_rect,E)
    # dev = pg.boolean(dev,dev_text,'or')
    return dev

def make_array(dev,column_num,row_num,label_str):
    """this function makes an array of identical devices from the given phidl device onject dev 

    Args:
        dev (phidl device): the device we want to make an array of
        column_num (int): the number of columns for the array 
        row_num (int): the number of rows for the array 
        label_str (string): the string used to label the array 

    Returns:
        phidl device: the array made out of the device given 
    """
    mystyle = {(0.01, 0): {'fill': '#CC00FF','stroke': 'black'}}
    dev.write_svg('Nanobeam.svg', style=mystyle)

    obj_width = dev.ymax-dev.ymin #calculate total height of cell
    y_gap = 5 # added between each devices within an array 

    pad = Device()
    pad_ref = pad.add_array(dev, columns = column_num, rows = row_num, spacing = [0,obj_width+y_gap])
    pad.center = [0,0]

    # add label to pad
    label = Device()
    T = pg.text(label_str, layer = 0, size = 0.5)
    label << pg.invert(T, border = 0.15)
    label.align(elements = 'all', alignment = 'x')
    label.align(elements = 'all', alignment = 'y')
    label.center = [0,0]
    label.movex(pad.xmin-0.5)
    label.movey(pad.ymin-0.5)
    pad << label
    # add the alignment markers
    add_alignment_markers(pad)
    qp(pad) # debugging

    return pad 
    
def make_Grid(samples,grid_row,grid_col):
    """this function makes a grid using the arrays given

    Args:
        samples (list): list storing all the arrays 

    Returns:
        grid: the grid made out of the arrays given
    """
    
    grid_spacing_row = 10
    grid_spacing_col = 30
    
    Grid = pg.grid(samples,
            spacing = (grid_spacing_col,grid_spacing_row),
            separation = True,
            shape = (grid_row,grid_col),
            align_x = 'x',
            align_y = 'y',
            edge_x = 'x',
            edge_y = 'ymax')

    Grid.center = [0,0]
    qp(Grid)
    return Grid 

def gen_gds_list(gdsdir,gdsfiles):
    sample_list = []
    # list of minimum waveguide width we are sweeping through 
    w_list = [0.06,0.120,0.180] # in um
    for i in range(len(gdsfiles)):
        for j in range(len(w_list)):
            fullfilepath = gdsdir+gdsfiles[i]
            dev = make_cavity_v3(fullfilepath,w_list[j])
            nums = re.findall(r'[\d]+',gdsfiles[i])
            v_num = nums[0]
            a_num = nums[1]
            w_num = str(j+1)
            label_str = 'v' + v_num + 'a'+ a_num + 'w' + w_num
            pad = make_array(dev,column_num=1,row_num=3,label_str=label_str) # make 5 arrays at the same time
            sample_list.append(pad)
    return sample_list