import openvsp as vsp



def add_wing(**kwargs):
    '''This function accepts the arguments required to generate a wing.  Currently this function can accept:
    span, sweep, root_chord, tip_chord, x_loc, y_loc, z_loc, span_panel_count, chord_panel_count, airfoilType,ThickLoc
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
        chord_panel_count = kwargs["chord_panel_count"]
        vsp.SetParmVal( wid, "Tess_W", "Shape", chord_panel_count)
    # set spanwise number of panels
    if "span_panel_count" in kwargs:
        span_panel_count = kwargs["span_panel_count"]
        vsp.SetParmVal( wid, "SectTess_U", "XSec_1", span_panel_count)

    # set airfoil type
    if "airfoilType" in kwargs:
        airfoilType = kwargs["airfoilType"]
        xsec_surf = vsp.GetXSecSurf( wid, 0 )
        vsp.ChangeXSecShape( xsec_surf, 0, airfoilType )
        xsec_surf1 = vsp.GetXSecSurf( wid, 1 )
        vsp.ChangeXSecShape( xsec_surf1, 1, airfoilType )

        if airfoilType == 10:
            # wedge airfoil
            if "ThickLoc" in kwargs:
                thickLoc = kwargs["ThickLoc"]
                xsec_id1 = vsp.GetXSec(xsec_surf,0)
                thickLoc1 = vsp.GetXSecParm(xsec_id1,"ThickLoc")
                vsp.SetParmVal( thickLoc1, thickLoc)
                xsec_id2 = vsp.GetXSec(xsec_surf,1)
                thickLoc2 = vsp.GetXSecParm(xsec_id2,"ThickLoc")
                vsp.SetParmVal( thickLoc2, thickLoc)
        elif airfoilType == 7:
            # NACA 4 digit 
            if "NACA4" in kwargs:
                naca4 = kwargs["NACA4"]
                # calculate naca parameters from input 4 digit id
                toc = naca4%100
                camber_loc = (naca4-toc)/100%10
                camber = ((naca4-toc)/100-camber_loc)/10

                toc /= 100
                camber_loc /=100
                camber /=10

                # get section id for middle section
                xsec_id1 = vsp.GetXSec(xsec_surf,0)
                # get id for location of max camber and set param
                thickLoc1 = vsp.GetXSecParm(xsec_id1,"CamberLoc")
                vsp.SetParmVal( thickLoc1, camber_loc)
                # get id for thickness ratio and set param
                toc_id1 = vsp.GetXSecParm(xsec_id1,"ThickChord")
                vsp.SetParmVal( toc_id1, toc)
                # get id for camber amount and set param
                camber_id1 = vsp.GetXSecParm(xsec_id1,"Camber")
                vsp.SetParmVal( camber_id1, camber)
                

                # get section id for tip section
                xsec_id2 = vsp.GetXSec(xsec_surf,1)
                # get id for location of max camber and set param
                thickLoc2 = vsp.GetXSecParm(xsec_id2,"CamberLoc")
                vsp.SetParmVal( thickLoc2, camber_loc)
                # get id for thickness ratio and set param
                toc_id2 = vsp.GetXSecParm(xsec_id1,"ThickChord")
                vsp.SetParmVal( toc_id2, toc)
                # get id for camber amount and set param
                camber_id2 = vsp.GetXSecParm(xsec_id1,"Camber")
                vsp.SetParmVal( camber_id2, camber)
            

    ## default settings
    # set inner cap to none
    vsp.SetParmVal( wid, "CapUMinOption", "EndCap", 0)
    # set outer cap to round
    vsp.SetParmVal( wid, "CapUMaxOption", "EndCap", 2)



    # Generate the geometry
    vsp.Update()
    area = vsp.GetParmVal(wid,"Area","XSec_1")
    return wid, area
    
    
def gen_multi_wing_geom(x_loc,y_loc,z_loc):
    '''This function generates a pair of wings with the main wing at the origin and the second wing at the points given in the input.
    This function returns the area of the wings for use in MachLine'''
    mainID,area_main = add_wing(span=10,sweep=30,root_chord=4,tip_chord=1,airfoilType=10,ThickLoc=0.25,chord_panel_count=30,span_panel_count=40)
    upperID,area_upper = add_wing(span=10,sweep=30,root_chord=4,tip_chord=1,x_loc=x_loc,y_loc=y_loc,z_loc=z_loc,airfoilType=10,ThickLoc=0.25,chord_panel_count=30,span_panel_count=40)

    # vsp.WriteVSPFile("my_aircraft.vsp3")
    # Export the OpenVSP file
    vsp.WriteVSPFile("my_aircraft.vsp3")
    
    # area_main = vsp.GetParmVal(mainID,"Area","XSec_1")
    # area_upper = vsp.GetParmVal(upperID,"Area","XSec_1")
    area_t = area_main + area_upper
    # Save the geometry to a file (STL format)
    file_name = "multi_lifting_surface_{0}_{1}_{2}_".format(x_loc,y_loc,z_loc)
    file_name+=".stl"
    vsp.ExportFile(file_name,1,2,0,0)
    return area_t



if __name__ == "__main__":
    
    
    naca4 = 6510

    toc = naca4%100
    camber_loc = (naca4-toc)/100%10
    camber = ((naca4-toc)/100-camber_loc)/10

    toc /= 100
    camber_loc /=100
    camber /=10

    print(toc)
    print(camber_loc)
    print(camber)
