import openvsp as vsp



def add_wing(**kwargs):
    '''This function accepts the arguments required to generate a wing.  Currently this function can accept:
    span, sweep, root_chord, tip_chord, x_loc, y_loc, z_loc, span_panel_count, chord_panel_count
    '''
    
    ## Create a new OpenVSP wing
    wid = vsp.AddGeom("WING","")

    ## Set parameters for the wing
    if 'span' in kwargs:
        span = kwargs['span']
        vsp.SetParmVal( wid, "Span", "XSec_1", span)
    if 'sweep' in kwargs:
        sweep = kwargs['sweep']
        vsp.SetParmVal(wid,"Sweep","XSec_1", sweep)
    if 'root_chord' in kwargs:
        root_chord = kwargs["root_chord"]
        vsp.SetParmVal( wid, "Root_Chord", "XSec_1", root_chord)
    if 'tip_chord' in kwargs:
        tip_chord = kwargs['tip_chord']
        vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", tip_chord )
    if 'x_loc' in kwargs:
        x_loc = kwargs['x_loc']
        vsp.SetParmVal( wid, "X_Rel_Location", "XForm", x_loc )
    if 'y_loc' in kwargs:
        y_loc = kwargs['y_loc']
        vsp.SetParmVal( wid, "Y_Rel_Location", "XForm", y_loc )  
    if 'z_loc' in kwargs:
        z_loc = kwargs['z_loc']
        vsp.SetParmVal( wid, "Z_Rel_Location", "XForm", z_loc )
    # if 'area' in kwargs:
    #     area = kwargs['area']
    #     vsp.SetParmVal(wid,"Area",'XSec_1',area)
    #     print("set area")
    # if 'taper' in kwargs:
    #     taper = kwargs['taper']
    #     vsp.SetParmVal(wid,"Taper",'XSec_1',taper)
    #     print("set taper")
    

    ## set panel density

    # set chordwise number of panels
    if "chord_panel_count" in kwargs:
        vsp.SetParmVal( wid, "Tess_W", "Shape", 25)
    # set spanwise number of panels
    if "span_panel_count" in kwargs:
        vsp.SetParmVal( wid, "SectTess_U", "XSec_1", 25)


    ## default settings
    # set inner cap to none
    vsp.SetParmVal( wid, "CapUMinOption", "EndCap", 0)
    # set outer cap to round
    vsp.SetParmVal( wid, "CapUMaxOption", "EndCap", 2)



    # Generate the geometry
    vsp.Update()
    return wid
    
    
def gen_multi_wing_geom(x_loc,y_loc,z_loc):
    '''This function generates a pair of wings with the main wing at the origin and the second wing at the points given in the input.
    This function returns the area of the wings for use in MachLine'''
    mainID = add_wing(span=10,sweep=0,root_chord=1,tip_chord=1)
    upperID = add_wing(span=10,sweep=0,root_chord=1,tip_chord=1,x_loc=x_loc,y_loc=y_loc,z_loc=z_loc)

    # vsp.WriteVSPFile("my_aircraft.vsp3")
    # Export the OpenVSP file
    vsp.WriteVSPFile("my_aircraft.vsp3")
    
    area_main = vsp.GetParmVal(mainID,"Area","XSec_1")
    area_upper = vsp.GetParmVal(upperID,"Area","XSec_1")
    area_t = area_main + area_upper
    # Save the geometry to a file (STL format)
    file_name = "multi_lifting_surface_{0}_{1}_{2}_".format(x_loc,y_loc,z_loc)
    file_name+=".stl"
    vsp.ExportFile(file_name,1,2,0,0)
    return area_t

if __name__ == "__main__":
    
    
    area = gen_multi_wing_geom(0,0,5)
    print("Area", area)
