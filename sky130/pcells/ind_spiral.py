'''Generate a simple asymmetric spiral inductor'''
import numpy as np
from math import ceil

import gdsfactory as gf
from gdsfactory.typings import Float2, LayerSpec
from sky130.pcells.via_generator import via_generator

@gf.cell
def ind_spiral(
    NUM_TURNS=3,
    NUM_SIDES=8,
    INNER_RADIUS=10,
    TRACK_SPACING = 1,
    TRACK_WIDTH = 2,
    VIA_FARM_LENGTH = 4,
    top_layer: LayerSpec = (71,20),
    bot_layer: LayerSpec = (70,20),
    via_layer: LayerSpec = (70,44)
)-> gf.Component:
    """Returns a spiral inductor

    .. plot::
      :include-source:

      import sky130

      c = sky130.pcells.ind_spiral(top_layer=(71,20),bot_layer = 70,20)
      c.plot()
    """
    assert NUM_SIDES>=4, f'{NUM_SIDES} is not greater than 3' # 3 works, but not pretty

    c = gf.Component()
    
    # TESTING
    # NUM_TURNS = 3 # number of turns
    # NUM_SIDES = 12 # number of sides
    # INNER_RADIUS = 5
    # TRACK_SPACING = 1 # spacing between adjacent spiral tracks
    # TRACK_WIDTH = 2 # track width
    


    # VIA_FARM_LENGTH = 3*TRACK_WIDTH # overlap between underpass and port1 track
    # VIA_FARM_LENGTH = TRACK_WIDTH
    # -------------------------------------------------------

    THETA = 2*np.pi/NUM_SIDES
    # Outer spiral
    v_outer = []
    v_inner = []

    rad_rate = (TRACK_SPACING+TRACK_WIDTH)/NUM_SIDES # in N angular steps, the radius grows by (S+W)
    LIM = int(NUM_SIDES*NUM_TURNS-NUM_SIDES/2+1) if NUM_SIDES%2==0 else int(NUM_SIDES*NUM_TURNS-(NUM_SIDES-1)/2) # keeps it one half turn away from the full N*t segments, so that the two ports are radially opposite one another (redues the shunt capacitance, which limits the SRF)
    for k in range(LIM):
        val_cos = np.cos(k*THETA)
        val_sin = np.sin(k*THETA)
        rad_inner = INNER_RADIUS+k*rad_rate
        rad_outer = rad_inner + TRACK_WIDTH
        v_inner.append((rad_inner*val_cos,rad_inner*val_sin))
        v_outer.append((rad_outer*val_cos,rad_outer*val_sin))

    if NUM_SIDES%2==1: # for odd sided spirals, we'll need an extra half step to reach the other port
        k += 0.5
        val_cos = np.cos(k*THETA)
        val_sin = np.sin(k*THETA)
        rad_inner = INNER_RADIUS+k*rad_rate
        rad_outer = rad_inner + TRACK_WIDTH
        v_inner.append((rad_inner*val_cos,rad_inner*val_sin))
        v_outer.append((rad_outer*val_cos,rad_outer*val_sin))

    vertices = []
    vertices.extend(v_inner)
    v_outer.reverse()
    vertices.extend(v_outer)
    
    # A small pad in the center
    #WIP
    dx = (TRACK_WIDTH/2)*np.tan(THETA/2)
    dy = TRACK_WIDTH/2
    # we apply this correction to the innermost vertex in the outer spiral
    temp_pt = vertices[-1]
    # p = (temp_pt[0]-dx, temp_pt[1]+dy)
    # vertices[-1] = p

    # Continue with the same tangent
    p = (temp_pt[0]-dx, -dy)
    vertices.append(p)

    p = (INNER_RADIUS-VIA_FARM_LENGTH, -dy)
    vertices.append(p)

    p = (INNER_RADIUS-VIA_FARM_LENGTH, dy)
    vertices.append(p)

    # move the first point along Y by dy of inner spiral to via farm pad edge
    p = (vertices[0][0]-dx,vertices[0][1]+TRACK_WIDTH/2)
    vertices[0] = p

    # --------------------------------------------------------
    x_vals = [v[0] for v in vertices]
    y_vals = [v[1] for v in vertices]

    # Outer Metal Tracks for ports
    # L = 3*W
    L = NUM_TURNS*(TRACK_SPACING+TRACK_WIDTH)
    
    # underpass
    pos_inner_track_edge = INNER_RADIUS-VIA_FARM_LENGTH
    len_underpass = VIA_FARM_LENGTH + INNER_RADIUS + NUM_TURNS*(TRACK_WIDTH+TRACK_SPACING)+TRACK_SPACING-pos_inner_track_edge
    c << gf.components.rectangle(size=(len_underpass,TRACK_WIDTH), centered=True, layer=bot_layer).movex(pos_inner_track_edge + 0.5*len_underpass) # port 1

    # port 1 track
    pos_P1_track_edge = pos_inner_track_edge + len_underpass - VIA_FARM_LENGTH
    len_port1_track = L+VIA_FARM_LENGTH
    c << gf.components.rectangle(size=(len_port1_track,TRACK_WIDTH), centered=True, layer=top_layer).movex(pos_P1_track_edge + 0.5*len_port1_track) # port 1

    # port 2 track
    pos_P2_track_egde = rad_inner
    len_port2_track = L+VIA_FARM_LENGTH
    c << gf.components.rectangle(size=(len_port2_track,TRACK_WIDTH), centered=True, layer=top_layer).movex(-(pos_P2_track_egde + 0.5*(L+VIA_FARM_LENGTH))) # port 2
  

    from shapely.geometry.polygon import Polygon
    c.add_polygon([tuple(x_vals),tuple(y_vals)], layer=top_layer)
    # plotVertices(vertices)

    # Adding the vias
    #FIXME: Via farm is not exactly centered on the track. Need to fix the offset.
    # via_drc_backoff = 0.2*W # backoff between via farm and r=track edges to avoid DRCs
    # via_wid = V-2*via_drc_backoff
    # via_len = W-2*via_drc_backoff
    # c << via_generator(width=via_wid,length=via_len,via_layer=via_layer).movex(pos_inner_track_edge+via_drc_backoff).movey(-via_len/2+via_drc_backoff)

    #TESTING:
    # via_drc_backoff = 0#0.25 # backoff between via farm and r=track edges to avoid DRCs
    via_wid = VIA_FARM_LENGTH # num of columns, X dimension of via farm
    via_len = TRACK_WIDTH # num of rows, Y dimension of via farm
    via_size=(0.17,0.17)
    via_spacing = (0.17,0.17)
    via_enclosure = (0.06,0.06)
    
    via_centering_X = via_wid%(via_size[0]+via_spacing[0])-via_spacing[0]
    via_centering_Y = via_len%(via_size[1]+via_spacing[1])-via_spacing[1]
    # print(f'Centering offsets are (dX,dY) = {(via_centering_X,via_centering_Y)}')

    vf = via_generator(width=via_wid,length=via_len,via_layer=via_layer, via_size=via_size,via_spacing=via_spacing, via_enclosure=via_enclosure)
    c << vf.movey(-via_len/2 + via_centering_Y/2).movex(-via_wid/2 + via_centering_X/2).movex(pos_inner_track_edge+via_wid/2)
    
    c << vf.movey(-via_len/2 + via_centering_Y/2).movex(-via_wid/2 + via_centering_X/2).movex(pos_P1_track_edge+via_wid/2)
    #.movex(pos_inner_track_edge+via_drc_backoff).movey(-via_len/2+via_drc_backoff)
    # c << gf.components.via(size=[5*0.17, 5*0.17], spacing=[0.17, 0.17], enclosure=1.0, layer=via_layer, bbox_offset=0)
    return c

def plotVertices(vertices):
    import matplotlib.pyplot as plt
    x_vals = [v[0] for v in vertices]
    y_vals = [v[1] for v in vertices]
    plt.plot(x_vals,y_vals,'-')
    plt.scatter(x_vals,y_vals)
    n = list(range(len(vertices)))
    for i, txt in enumerate(n):
        plt.annotate(txt, (x_vals[i], y_vals[i]))
    plt.grid()
    plt.show()

if __name__ == "__main__":

    # c = p_n_poly(p_poly_width= 5.73, p_poly_length=2)
    c = ind_spiral(
    NUM_TURNS=3,
    NUM_SIDES=8,
    INNER_RADIUS=10,
    TRACK_SPACING = 1,
    TRACK_WIDTH = 2,
    VIA_FARM_LENGTH = 4
    )
    # c.show(show_ports=True)

    c.plot()
    # scene = c.to_3d()
    # scene.show()
    
