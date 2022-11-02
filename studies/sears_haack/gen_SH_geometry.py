import numpy as np

def _generate_SH(rmax, length, num_ax = 100, num_theta = 30):
    
    R_L_ratio = rmax/length
    # volume = ((3*np.pi**2)/16.)*rmax*rmax*length
    # Generate radius distribution
    x_range = np.linspace(0,length,num_ax)
    r_range = rmax*(4*(x_range/length)*(1 - (x_range/length)))**(3./4.)
    theta_range = np.linspace(0, 2*np.pi, num_theta + 1)
    
    # Initialize vertex position arrays
    x = np.ones((num_ax, num_theta))
    y = np.ones((num_ax, num_theta))
    z = np.ones((num_ax, num_theta))
    
    # generate vertex values
    for i in range(num_ax):
        x[i,:] = x[i,:]*x_range[i]
        if i == 0 or i == (num_ax - 1):
            y[i,:] = y[i,:]*0.0
            z[i,:] = z[i,:]*0.0
        else:
            for j in range(num_theta):
                y[i,j] = r_range[i]*np.sin(theta_range[j])
                z[i,j] = r_range[i]*np.cos(theta_range[j])
    
    # Initialize concatenated vertex array
    vertices = np.ones((2 + (num_theta*(num_ax - 2)),3))
    
    # generate concatenated vertex values
    for i in range(num_ax):
        if i == 0:
            vertices[0,:] = [x[0,0], y[0,0], z[0,0]]
        elif i == (num_ax - 1):
            vertices[-1,:] = [x[-1,0], y[-1,0], z[-1,0]]
        else:
            for j in range(num_theta):
                vert_index = (i - 1)*num_theta + j + 1
                vertices[vert_index, :] = [x[i,j], y[i,j], z[i,j]]\
    
    # Assign vertices to triangular elements
    nVerts = len(vertices[:,0])
    nTris = (num_ax - 3)*num_theta*2 + 2*num_theta
    
    # index of verts defining each tri
    tri_verts = np.ones((nTris,3))
    
    # index of the first vertex of each axial slice
    mid_strt_ind = np.arange(2, nVerts - 2 - num_theta, num_theta)
    # index of the first vertex of each following axial slice
    mid_end_ind = np.arange(2 + num_theta, nVerts - 2, num_theta)
    
    # assign nose and tail triangles
    for i in range(nTris):
        if i < num_theta:
            #nose triangles
            if (i+1) == num_theta:
                tri_verts[i] = [1,i+2,2]
            else:   
                tri_verts[i] = [1,i+2,i+3]
        elif i >= nTris - num_theta:
            # tail triangles, this is a mess but it works
            '''This is causing inverted CFD cells'''
            count = num_theta*(num_ax - 2) - 2
            if i == nTris - 1:
                tri_verts[i] = [i - count, nVerts, nVerts - num_theta]
            else:   
                tri_verts[i] = [i - count, nVerts,  i - count + 1]
    
    # Assign all ofther triangles
    for i in range(len(mid_strt_ind)):
        A = np.arange(mid_strt_ind[i], mid_strt_ind[i] + num_theta)
        B = np.arange(mid_end_ind[i], mid_end_ind[i] + num_theta)
        c = []
        for j in range(num_theta):
            c.append((A[j], B[j], B[(j + 1) % num_theta]))
            c.append((A[j], B[(j + 1) % num_theta], A[(j + 1) % num_theta]))
        C = np.asarray(c)
        tri_verts[i*(2*num_theta) + num_theta:(i + 1)*(2*num_theta) + num_theta, :] = C
    tri_verts = tri_verts. astype(int)
    
    #temporary component assignment
    comp_num = np.ones((nTris,1))
    
    return nVerts, nTris, vertices, tri_verts, comp_num, R_L_ratio

def _write_TRI(nVerts, nTris, vertices, tri_verts, comp_num, filename = 'test_file.tri'):
    
    '''write to .tri file'''
    
    with open(filename, 'w') as export_handle:
        
        # Write header
        print("{0:<18}{1:<18}".format(nVerts, nTris), file=export_handle)
            
        for vertex in vertices:
            print("{0:<18.10}{1:<18.10}{2:<18.10}".format(*vertex), file=export_handle)
        
        for tri_int in tri_verts:
            print("{0:<12}{1:<12}{2:<12}".format(*tri_int), file=export_handle)
            
        for comp_num in comp_num:
            # print("{0:<4d}".format(int(comp_num[0])), file=export_handle)
            print("{0:<4d}".format(int(1)), file=export_handle)
            
'''------------RUN SCRIPT------------'''
if __name__ == "__main__":

    num_ax = 80 # number of axial panels
    num_theta = 30 # number of radial panels
    length = 0.6096 # length of SH body
    rmax = length * 0.037879 # maximum radius of SH body
        
    nVerts, nTris, vertices, tri_verts, comp_num, r_l_ratio = _generate_SH(rmax, length, num_ax, num_theta)
    print('R over L ratio: ', r_l_ratio)

    _write_TRI(nVerts, nTris, vertices, tri_verts, comp_num, filename = 'studies/sears_haack/meshes/SH_{0}_{1}.tri'.format(num_ax, num_theta))
