import numpy as np
import pandas as pd

def pdb2dataframe(filename):
    #Load PA molecule from PDB file
    file = open(filename,'r')
    raw_lines = file.readlines()
    lines = [line.rstrip() for line in raw_lines]

    res = []
    bead = []
    x_coord = []
    y_coord = []
    z_coord = []
    for line in lines:
        split_line = line.split(' ')
        split_line = list(filter(None, split_line))
        if split_line[0] == 'ATOM':
            res.append(split_line[3])
            bead.append(split_line[2])
            x_coord.append(float(split_line[5]))
            y_coord.append(float(split_line[6]))
            z_coord.append(float(split_line[7]))

    PA_pdb = {'RES': res, 'BEAD_TYPE': bead, 'X':z_coord, 'Y':y_coord, 'Z': x_coord}
    df = pd.DataFrame(PA_pdb)
    return df

def pdb2fiber(input_pdb_filename, output_pdb_filename, molname = 'PA', topfilename = 'PA_box.top', Lx = 16, Ly = 16, Lz = 16, n_radial_cross_section = 16, stacking_distance = 3.5, offset = 4):

    #box_edge_dimension (nm)
    #stacking_distance (angstroms)
    #offset (angstroms)

    
    df = pdb2dataframe(input_pdb_filename)
    n_stacking = np.round(Lz/(stacking_distance / 10))
    
    #Load Initial Coord
    c = df[['X','Y','Z']].to_numpy()

    # Isolate PALM Tail
    df_pam = df.loc[(df['RES'] == 'PAM')]
    df_pam_sc = df_pam.loc[(df['BEAD_TYPE'] != 'BB')]
    df_pam_sc3 = df_pam.loc[(df['BEAD_TYPE'] == 'SC3')]

    #Center PA based 
    center_pt = df_pam_sc3[['X','Y','Z']]
    center = np.tile(center_pt,(len(c),1))
    c = c - center

    #Center PALM Tail to Origin
    w = df_pam_sc[['X','Y','Z']].to_numpy()
    center = np.tile(center_pt,(len(w),1))
    w_centered = w-center

    #Fit X,Y coordinates
    x = w_centered[:,0]
    y = w_centered[:,1]
    m, b = np.polyfit(x, y, 1)

    #Plot New Basis
    x = w_centered[:,0]
    y = w_centered[:,1]
    m, b = np.polyfit(x, y, 1)

    #Use 1D line to create change of basis matrix
    v1 = [-(m+b),1,0]
    v2 = [1, m+b, 0]
    v3 = [0,0,1]
    v = np.matrix([v1 / np.linalg.norm(v1),v2 / np.linalg.norm(v2),v3])

    #Transform full coordinates to new basis
    c = np.matrix(c)
    u = np.identity(3)
    d = np.linalg.inv(v)*u*c.transpose()
    d = d.transpose()

    #Offset PA from origin to avoid duplicate atom coordinates after rotation
    off_center_pt = [0,offset, 0]
    center = np.tile(off_center_pt,(len(d),1))
    d = d + center

    #Assign matrix for rotation (v) and stacking (z)
    v = d[:,(0,1)]
    z = d[:,2]

    #Initialize GRO file
    f= open(output_pdb_filename,"w+")
    f.write("GRowing Old MAkes el Chrono Sweat\n")
    n_beads = len(v)*n_radial_cross_section*n_stacking
    beads_per_pa = len(df)
    f.write(" %d\n"%n_beads)
    bead_number = 1
    res_number = 1

    for k in range(int(n_stacking)):
        for i in range(n_radial_cross_section):

            #Rotate PA molecule in XY-plane
            theta = 2*np.pi / n_radial_cross_section
            R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
            R = np.asmatrix(R)
            vo = R*(v.transpose())

            #Write GRO file
            for j in range(len(v)):
                x_coord = v[j,0]/10 + Lx/2
                y_coord = v[j,1]/10 + Ly/2
                z_coord = z[j]/10 + stacking_distance*k/10
                f.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n'%(res_number,df['RES'][j], df['BEAD_TYPE'][j], bead_number, x_coord,y_coord,z_coord))
                bead_number = bead_number + 1
                if j == len(v)-1:
                    res_number = res_number + 1
                elif df['BEAD_TYPE'][j+1] == 'BB':
                    res_number = res_number + 1
            v = vo.transpose()

    f.write('%f\t%f\t%f\n'%(Lx,Ly,Lz))
    f.close()
    
    
    #Update PDB File
    # Calculate actual nmol added
    actual_nmol = int(n_beads / beads_per_pa)

    
    # read <molname>.top
    with open('%s.top'%molname, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+='%s'%molname + ' '*(9-len(molname)) + '%s\n'%actual_nmol
            else:
                data+=line
        
    if topfilename==None:
        # Update number of PA molecules in <molname>.top file    
        with open('%s.top'%molname, 'w') as f:
            f.write(data)
    else:
        # Create <newtopfilename> file    
        with open('%s'%topfilename, 'w') as f:
            f.write(data)

def pdb2bilayer(input_pdb_filename, output_pdb_filename, molname = 'PA', topfilename = 'PA_box.top',Lx = 16, Ly = 16, Lz = 16, stacking_distance = 3.5, offset = 4):

    #box_edge_dimension (nm)
    #stacking_distance (angstroms)
    #offset (angstroms)

    
    df = pdb2dataframe(input_pdb_filename)
    PA_per_Xedge = np.round(Lx / (stacking_distance / 10))
    PA_per_Zedge = np.round(Lx / (stacking_distance / 10))

    
    #Load Initial Coord
    c = df[['X','Y','Z']].to_numpy()

    # Isolate PALM Tail
    df_pam = df.loc[(df['RES'] == 'PAM')]
    df_pam_sc = df_pam.loc[(df['BEAD_TYPE'] != 'BB')]
    df_pam_sc3 = df_pam.loc[(df['BEAD_TYPE'] == 'SC3')]

    #Center PA based 
    center_pt = df_pam_sc3[['X','Y','Z']]
    center = np.tile(center_pt,(len(c),1))
    c = c - center

    #Center PALM Tail to Origin
    w = df_pam_sc[['X','Y','Z']].to_numpy()
    center = np.tile(center_pt,(len(w),1))
    w_centered = w-center

    #Fit X,Y coordinates
    x = w_centered[:,0]
    y = w_centered[:,1]
    m, b = np.polyfit(x, y, 1)

    #Plot New Basis
    x = w_centered[:,0]
    y = w_centered[:,1]
    m, b = np.polyfit(x, y, 1)

    #Use 1D line to create change of basis matrix
    v1 = [-(m+b),1,0]
    v2 = [1, m+b, 0]
    v3 = [0,0,1]
    v = np.matrix([v1 / np.linalg.norm(v1),v2 / np.linalg.norm(v2),v3])

    #Transform full coordinates to new basis
    c = np.matrix(c)
    u = np.identity(3)
    d = np.linalg.inv(v)*u*c.transpose()
    d = d.transpose()

    #Offset PA from origin to avoid duplicate atom coordinates after rotation
    off_center_pt = [0,offset, 0]
    center = np.tile(off_center_pt,(len(d),1))
    d = d + center

    #Assign matrix for rotation (v) and stacking (z)
    v = d[:,(0,1)]
    z = d[:,2]

    #Initialize GRO file
    f= open(output_pdb_filename,"w+")
    f.write("GRowing Old MAkes el Chrono Sweat\n")
    n_beads = len(v)*PA_per_Xedge*PA_per_Zedge*2
    beads_per_pa = len(df)
    f.write(" %d\n"%n_beads)
    bead_number = 1
    res_number = 1

    for i in range(int(PA_per_Xedge)):
        for k in range(int(PA_per_Zedge)):
            for r in range(2):
                #Rotate PA molecule in XY-plane
                theta = np.pi
                R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
                R = np.asmatrix(R)
                vo = R*(v.transpose())
                #Write GRO file
                for j in range(len(v)):
                    x_coord = v[j,0]/10 + i*stacking_distance/10
                    y_coord = v[j,1]/10 + r*stacking_distance/10
                    z_coord = z[j]/10 + k*stacking_distance/10
                    f.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n'%(res_number,df['RES'][j], df['BEAD_TYPE'][j], bead_number, x_coord,y_coord,z_coord))
                    bead_number = bead_number + 1
                    if j == len(v)-1:
                        res_number = res_number + 1
                    elif df['BEAD_TYPE'][j+1] == 'BB':
                        res_number = res_number + 1
                v = vo.transpose()

    f.write('%f\t%f\t%f\n'%(Lx,Ly,Lz))
    f.close()
    
    #Update PDB File
    # Calculate actual nmol added
    actual_nmol = int(n_beads / beads_per_pa)

    
    # read <molname>.top
    with open('%s.top'%molname, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+='%s'%molname + ' '*(9-len(molname)) + '%s\n'%actual_nmol
            else:
                data+=line
        
    if topfilename==None:
        # Update number of PA molecules in <molname>.top file    
        with open('%s.top'%molname, 'w') as f:
            f.write(data)
    else:
        # Create <newtopfilename> file    
        with open('%s'%topfilename, 'w') as f:
            f.write(data)
def pdb2micelle(input_pdb_filename, output_pdb_filename, molname = 'PA', topfilename = 'PA_box.top',Lx = 16, Ly = 16, Lz = 16, n_radial_cross_section = 16, offset = 4):

    #box_edge_dimensions (nm)
    #offset (angstroms)
    
    #Error Check
    if np.mod(n_radial_cross_section,2) == 1:
        print('Error: Number of PAs per cross section must be an even number')
        return exit
    #Convert pdb file to pandas DataFrame
    df = pdb2dataframe(input_pdb_filename)
    
    #Load Initial Coord
    c = df[['X','Y','Z']].to_numpy()

    # Isolate PALM Tail
    df_pam = df.loc[(df['RES'] == 'PAM')]
    df_pam_sc = df_pam.loc[(df['BEAD_TYPE'] != 'BB')]
    df_pam_sc3 = df_pam.loc[(df['BEAD_TYPE'] == 'SC3')]

    #Center PA based 
    center_pt = df_pam_sc3[['X','Y','Z']]
    center = np.tile(center_pt,(len(c),1))
    c = c - center

    #Center PALM Tail to Origin
    w = df_pam_sc[['X','Y','Z']].to_numpy()
    center = np.tile(center_pt,(len(w),1))
    w_centered = w-center

    #Fit X,Y coordinates
    x = w_centered[:,0]
    y = w_centered[:,1]
    m, b = np.polyfit(x, y, 1)

    #Plot New Basis
    x = w_centered[:,0]
    y = w_centered[:,1]
    m, b = np.polyfit(x, y, 1)

    #Use 1D line to create change of basis matrix
    v1 = [-(m+b),1,0]
    v2 = [1, m+b, 0]
    v3 = [0,0,1]
    v = np.matrix([v1 / np.linalg.norm(v1),v2 / np.linalg.norm(v2),v3])

    #Transform full coordinates to new basis
    c = np.matrix(c)
    u = np.identity(3)
    d = np.linalg.inv(v)*u*c.transpose()
    d = d.transpose()

    #Offset PA from origin to avoid duplicate atom coordinates after rotation
    off_center_pt = [0,offset, 0]
    center = np.tile(off_center_pt,(len(d),1))
    coord = d + center

    #Initialize GRO file
    beads_per_pa = len(df)
    n_beads = 0
    bead_number = 1
    res_number = 1
    f= open(output_pdb_filename,"w+")
    f.write("GRowing Old MAkes el Chrono Sweat\n")
    f.write(" %d\n"%n_beads)

    #Calculate theta
    theta = 2*np.pi / n_radial_cross_section
    phi_range = np.arange(-n_radial_cross_section/4,n_radial_cross_section/4+1,1)
    for i in phi_range:
        for k in range(n_radial_cross_section):
            new_coord = rotate3D(coord, theta*i, 0, theta*k)
           
            write = 1
            if k > 0 and i == phi_range[0]:
                write = 0
            elif k > 0 and i == phi_range[-1]:
                write = 0
            
            if write == 1:
                #Write GRO file
                for j in range(len(coord)):
                    x_coord = new_coord[j,0]/10
                    y_coord = new_coord[j,1]/10
                    z_coord = new_coord[j,2]/10
                    f.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n'%(res_number,df['RES'][j], df['BEAD_TYPE'][j], bead_number, x_coord,y_coord,z_coord))
                    bead_number = bead_number + 1
                    if j == len(coord)-1:
                        res_number = res_number + 1
                    elif df['BEAD_TYPE'][j+1] == 'BB':
                        res_number = res_number + 1
                    n_beads += 1

    f.write('%f\t%f\t%f\n'%(Lx,Ly,Lz))
    f.close()
    
    #adjust number of atoms
    f = open(output_pdb_filename,'r')
    lines = f.readlines()
    lines[1] = " %d\n"%n_beads
    
    f = open(output_pdb_filename,'w')
    f.writelines(lines)
    f.close()
    
    #Update PDB File
    # Calculate actual nmol added
    actual_nmol = int(n_beads / beads_per_pa)

    
    # read <molname>.top
    with open('%s.top'%molname, 'r') as f:
        data = ''
        for line in f:
            if len(line.split())!=0 and line.split()[0]=='%s'%molname:
                data+='%s'%molname + ' '*(9-len(molname)) + '%s\n'%actual_nmol
            else:
                data+=line
        
    if topfilename==None:
        # Update number of PA molecules in <molname>.top file    
        with open('%s.top'%molname, 'w') as f:
            f.write(data)
    else:
        # Create <newtopfilename> file    
        with open('%s'%topfilename, 'w') as f:
            f.write(data)

def rotate3D(coord, Xtheta, Ytheta, Ztheta):
    #Rotate X (rotation in YZ plane)
    R = np.array([[1, 0,0],[0, np.cos(Xtheta), -np.sin(Xtheta)],[0, np.sin(Xtheta), np.cos(Xtheta)]])
    Rx = np.asmatrix(R)
    coord0 = Rx*(coord.transpose())
    coord = coord0.transpose()
    
    #Rotate Y (rotation in XZ plane)
    R = np.array([[np.cos(Ytheta), 0, np.sin(Ytheta)],[0,1,0], [-np.sin(Ytheta), 0, np.cos(Ytheta)]])
    Ry = np.asmatrix(R)
    coord0 = Ry*(coord.transpose())
    coord = coord0.transpose()
    
    #Rotate Z (rotation in XY plane)
    R = np.array([[np.cos(Ztheta), -np.sin(Ztheta), 0], [np.sin(Ztheta), np.cos(Ztheta), 0], [0, 0, 1]])
    Rz = np.asmatrix(R)
    coord0 = Rz*(coord.transpose())
    coord = coord0.transpose()
    
    return coord
    

    
    