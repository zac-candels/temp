import os, math, re, sys, time
import numpy as np
import struct
from matplotlib import pyplot as plt

start = time.perf_counter()

def read_params(path):
    params = {}
    # matches:   key   =   integer   (ignoring anything after a '#')
    pat = re.compile(r'^\s*([A-Za-z_]\w*)\s*=\s*([0-9]+)')
    with open(path) as f:
        for line in f:
            m = pat.match(line)
            if m:
                key, val = m.group(1), int(m.group(2))
                params[key] = val
    return params

# Usage
path = "."
path_input = path + "/input.txt"
params = read_params(path_input)
saveInterval = params["saveInterval"]
Num_steps   = params["timesteps"]



plt.close('all')

# Directory where data is stored
datadir = path + "/data/"

# Open the Header file
HeaderFile = open(datadir+"Header.mat", 'rb')

# Simulation size, struct.unpack reads 4 bytes from Header.mat and interprets them as an integer
LX = struct.unpack('=i', HeaderFile.read(4))[0]
# Read the next 4 bytes...
LY = struct.unpack('=i', HeaderFile.read(4))[0]
LZ = struct.unpack('=i', HeaderFile.read(4))[0]

# 2D or 3D
ndim = struct.unpack('=i', HeaderFile.read(4))[0]

# What time to start at, what time to end at, and the time increment
tstart = 0
tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

#pressure = np.zeros((tend, LX, LY, LZ))

# Slice in the z direction to plot
zslice=0

# Where to save the plots
outDirName = "figures"
os.system("mkdir -p %s"%outDirName)

#Num_steps = int(100)
x_n = 0
centroid_pos = []    
            
vel_arr = []
for t_idx in range(tstart, tend + 1, tinc):
    # File containing boundary ids
    file_name_bdy = os.path.join(datadir, "BoundaryLabels_t" + str(t_idx) + ".mat")
    FileSolid = open(file_name_bdy, 'rb')
    dat=FileSolid.read()



    file_name_phi = os.path.join(datadir, "OrderParameter_t" + str(t_idx) + ".mat")
    File_phi = open(file_name_phi, 'rb')
    dat = File_phi.read()

    # Fill a numpy array of dimensions (LX,LY,LZ) with the data from the file in the format '=i' (integer). (4*LY*LZ,4*LZ,4) are steps taken in bytes for each dimension. E.g in the z direction we move 4 bytes to the next z value, in the y direction we move 4 bytes * the number of z values to the next y value, etc.
    solid = np.ndarray((LX, LY, LZ), '=i', dat, 0, (4 * LY * LZ, 4 * LZ, 4))
    FileSolid.close()

    phi = np.ndarray((LX, LY, LZ), '=d', dat, 0, (8 * LY * LZ, 8 * LZ, 8))
    liquid = np.array(phi[:,:])
    # Set order parameter in the solid to 0.5 for visualisation
    # liquid[np.where(np.logical_or(solid == 1, solid == -1))[0], np.where(np.logical_or(solid == 1, solid == -1))[1], np.where(np.logical_or(solid == 1, solid == -1))[2]] = 0.
    # liquid[np.where(np.logical_or(solid == 3, solid == 2))[0], np.where(np.logical_or(solid == 3, solid == 2))[1], np.where(np.logical_or(solid == 3, solid == 2))[2]] = 0.
    File_phi.close()
    phi = liquid[:,:,0]

    file_name = os.path.join(datadir, "Velocity_t" + str(t_idx) + ".mat")
    FileV = open(file_name, 'rb')
    dat = FileV.read()
    v = np.ndarray((LX, LY, LZ, ndim), '=d', dat, 0, (ndim * 8 * LY * LZ, ndim * 8 * LZ, ndim * 8, 8))
    FileV.close()


    v_x = v[:, :, 0, 0]
    v_y = v[:, :, 0, 1]
    phi_mult_vel = np.array([0.0, 0.0])
    phi_droplet = 0
        
    for i in range(len(phi[:,0])):
        for j in range(len(phi[0,:])):
            v_mag2 = v_x[i,j]**2 + v_y[i,j]**2
            if phi[i,j] > 0.5: 
                # if( v_mag2 > 0.000 ):
                #     print("\n||v||^2 = ", v_mag2) 
                phi_droplet += phi[i,j] 
                phi_mult_vel += phi[i,j] * np.array([ v_x[i,j],  v_y[i,j] ]) 

    phi_droplet = np.asarray( phi_droplet )
    phi_mult_vel = np.asarray( phi_mult_vel )

    droplet_vel = phi_mult_vel / phi_droplet 
    vel_arr.append(droplet_vel[0])
    print(droplet_vel[0])

vel_arr = np.asarray(vel_arr)

print("vel avg = ", np.mean(vel_arr))
print("std = ", np.std(vel_arr))
                

            



        
        





