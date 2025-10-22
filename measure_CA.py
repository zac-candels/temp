import os, math, re, sys
import numpy as np
import struct
from matplotlib import pyplot as plt

plt.close('all')
def coord_k(k, ly, lz):
    """From a k value, determines its xk, yk, and zk."""
    xk = math.floor(k/(ly*lz))
    yk = math.floor((k - xk*lz*ly)/lz)
    zk = k - xk*lz*ly - yk*lz
    return xk, yk, zk

def interpolate_x_direction(pts_x, pts_phi):
    pts_x = np.asarray(pts_x)
    pts_phi = np.asarray(pts_phi)
    gas_phase_x, gas_phase_phi = [], []
    liquid_phase_x, liquid_phase_phi = [], []
    if len(pts_x) > 2:
        #sorted_indices = np.argsort(pts_phi)
        #pts_phi = pts_phi[sorted_indices]
        #pts_x = pts_x[sorted_indices]
        mask_gas = pts_phi < 0.5
        mask_liquid = pts_phi >= 0.5
        
        gas_phase_phi = pts_phi[mask_gas]
        gas_phase_x = pts_x[mask_gas]
        
        liquid_phase_phi = pts_phi[mask_liquid]
        liquid_phase_x = pts_x[mask_liquid]
        
        if len(liquid_phase_phi) == 0:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_x[closest_gas_index]
            # Nothing to interpolate
            return closest_gas_pt
        elif len(gas_phase_phi) == 0:
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_x[closest_liquid_index]
            return closest_liquid_pt
        else:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_x[closest_gas_index]
            
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_x[closest_liquid_index]
        
        x0, x1 = [closest_liquid_pt, closest_gas_pt]
        phi0, phi1\
            = [gas_phase_phi[closest_gas_index],liquid_phase_phi[closest_liquid_index]]
            
    if len(pts_x) == 2:
        x0, x1 = pts_x[0], pts_x[1]
        phi0, phi1 = pts_phi[0], pts_phi[1]
    
    m = ( phi1 - phi0 )\
        /( x1 - x0 )
    x_interpolated = x0 + (1/m)*(0.5 - phi0)
    
    return x_interpolated


def interpolate_y_direction(pts_y, pts_phi):
    pts_y = np.asarray(pts_y)
    pts_phi = np.asarray(pts_phi)
    gas_phase_y, gas_phase_phi = [], []
    liquid_phase_y, liquid_phase_phi = [], []
    if len(pts_y) > 2:
        #sorted_indices = np.argsort(pts_phi)
        #pts_phi = pts_phi[sorted_indices]
        #pts_x = pts_x[sorted_indices]
        mask_gas = pts_phi < 0.5
        mask_liquid = pts_phi >= 0.5
        
        gas_phase_phi = pts_phi[mask_gas]
        gas_phase_y = pts_y[mask_gas]
        
        liquid_phase_phi = pts_phi[mask_liquid]
        liquid_phase_y = pts_y[mask_liquid]
        
        if len(liquid_phase_phi) == 0:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_y[closest_gas_index]
            # Nothing to interpolate
            return closest_gas_pt
        elif len(gas_phase_phi) == 0:
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_y[closest_liquid_index]
            return closest_liquid_pt
        else:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_y[closest_gas_index]
            
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_y[closest_liquid_index]
        
        y0, y1 = [closest_liquid_pt, closest_gas_pt]
        phi0, phi1\
            = [gas_phase_phi[closest_gas_index],liquid_phase_phi[closest_liquid_index]]
            
    if len(pts_y) == 2:
        y0, y1 = pts_y[0], pts_y[1]
        phi0, phi1 = pts_phi[0], pts_phi[1]
    
    m = ( phi1 - phi0 )\
        /( y1 - y0 )
    y_interpolated = y0 + (1/m)*(0.5 - phi0)
    
    return y_interpolated



plt.close('all')

# Directory where data is stored
datadir = "./data/"

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

def coord_k(k, ly, lz):
    """From a k value, determines its xk, yk, and zk."""
    xk = math.floor(k/(ly*lz))
    yk = math.floor((k - xk*lz*ly)/lz)
    zk = k - xk*lz*ly - yk*lz
    return xk, yk, zk
#pressure = np.zeros((tend, LX, LY, LZ))

# Slice in the z direction to plot
zslice=0

# Where to save the plots
outDirName = "figures"
os.system("mkdir -p %s"%outDirName)

t = 0

pattern_bdy = re.compile(r"BoundaryLabels_t(\d+)\.mat")
pattern_phi = re.compile(r"OrderParameter_t(\d+)\.mat")

files = os.listdir(datadir)
max_n = -1
max_m = -1
target_file_bdy = None
target_file_phi = None

for filename_bdy in files:
    match_bdy = pattern_bdy.match(filename_bdy)
    if match_bdy:
        n = int(match_bdy.group(1))
        if n > max_n:
            max_n = n
            target_file_bdy = filename_bdy
            

for filename_phi in files:
    match_phi = pattern_phi.match(filename_phi)
    if match_phi:
        m = int(match_phi.group(1))
        if m > max_m:
            max_m = m
            target_file_phi = filename_phi          
            

# File containing boundary ids
file_name_bdy = os.path.join(datadir, target_file_bdy)
FileSolid = open(file_name_bdy, 'rb')
dat=FileSolid.read()



file_name_phi = os.path.join(datadir, target_file_phi)
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


full_interface_x = []
full_interface_y = []
phi_interface = []
    
for i in range(len(phi[:,0])):
    for j in range(len(phi[0,:])):
        if phi[i,j] > 0.3 and phi[i,j] < 0.7:    
            if j >= 25:
                full_interface_x.append(i)
                full_interface_y.append(j)    
                phi_interface.append(phi[i,j])
              
full_interface_x = np.asarray(full_interface_x)
full_interface_y = np.asarray(full_interface_y)
phi_interface = np.asarray(phi_interface)
mid_pt = np.median(full_interface_x)

mask = full_interface_x < mid_pt
left_interface_x = full_interface_x[mask]
left_interface_y = full_interface_y[mask]
left_interface_phi = phi_interface[mask]
    

plt.figure()
plt.plot(full_interface_x, full_interface_y, 'o', color='r')
plt.title('Full interface')
plt.savefig("./interface.png")


# plt.figure()
# plt.plot(left_interface_x, left_interface_y, 'o')
# plt.title('left interface')


# Note that this only gives the left half of interface.
left_interpolated_interface_x = []
left_interpolated_interface_y = []

tol = 1e-5

# Go through each point (x,y,phi) in the interface. If for any value of y
# there are multiple x values, interpolate the value of phi between them
# to get a single value. This new value of x will replace the two previous ones.
for i in range(len(left_interface_x)):
    x_i, y_i, phi_i = left_interface_x[i], left_interface_y[i], left_interface_phi[i]
    pts_x, pts_y, pts_phi = [x_i], y_i, [phi_i]
    for j in range(len(left_interface_x)):
        # for each fixed point P_i, loop through every other 
        # point P_j to see if there are other points that have the 
        # same y-value but different x-values. That is to say,
        # there may be multiple points at the same y-value. We collect all of these into an array,
        # select the two closest to 0.5 - but on opposite sides of 0.5 - and do the same routine as
        # mentioned above - that is to say, interpolate phi between them and return a single point.
        if j == i:
            continue
        x_j, y_j, phi_j\
            = left_interface_x[j], left_interface_y[j], left_interface_phi[j]
        if abs(y_j - y_i) < tol: # if point P_j is at the same y-level as point P_i
            if abs(x_i - x_j) > 2*tol: # if the points are at different x-levels, create arrays pts_x and pts_phi
                pts_x.append(x_j), pts_phi.append(phi_j)

    if len(pts_x) == 1: # If there is only one point at the given y-level, there's nothing further to do.
                        # Include this point in the new interface.
        left_interpolated_interface_x.append(x_i)
        left_interpolated_interface_y.append(y_i)
    else: # If there are multiple points at the same y-level, interpolate the values of phi between them and return a single point, x_new.
          # add x_new to the new interface.
        x_new = interpolate_x_direction(pts_x, pts_phi)
        left_interpolated_interface_x.append( x_new )
        left_interpolated_interface_y.append( y_i )


# Do the same procedure as above, but this time look for points that have the same y-value but different x-values.

# for i in range(len(left_interface_x)):
#     x_i, y_i, phi_i = left_interface_x[i], left_interface_y[i], left_interface_phi[i]
#     if y_i < 4:
#         print('pause')
#     pts_x, pts_y, pts_phi = x_i, [y_i], [phi_i]
#     for j in range(len(left_interface_x)):
#         # for each fixed point P_i, I have to loop through every other 
#         # point P_j to see if there are other points that have the 
#         # same x-value but different y-values.
#         if j == i:
#             continue
#         x_j, y_j, phi_j\
#             = left_interface_x[j], left_interface_y[j], left_interface_phi[j]
#         if abs(x_j - x_i) < tol: # if point P_j is at the same x-level as point P_i
#             if abs(y_i - y_j) > tol: # if the points are at different y-levels, create arrays pts_y and pts_phi
#                 pts_y.append(y_j), pts_phi.append(phi_j)

#     if len(pts_y) == 1: # If there is only one point at the given x-level, there's nothing further to do.
#                         # Include this point in the new interface.
#         left_interpolated_interface_x.append(x_i)
#         left_interpolated_interface_y.append(y_i)
        
#         # left_interface_y = np.delete(left_interface_y, i)
#         # left_interface_x = np.delete(left_interface_x, i)
#         # phi_interface = np.delete(left_interface_phi, i)
#     else: # If there are multiple points at the same x-level, interpolate the values of phi between them and return a single point, y_new.
#           # add x_new to the new interface.
#         y_new = interpolate_y_direction(pts_y, pts_phi)
#         left_interpolated_interface_x.append( x_i )
#         left_interpolated_interface_y.append( y_new )
# #    plt.plot(x_i,y_new,'o')
        

        
        
left_interpolated_interface_x = np.asarray(left_interpolated_interface_x)
left_interpolated_interface_y = np.asarray(left_interpolated_interface_y)
        
right_interpolated_interface_x = []
right_interpolated_interface_y = []
       
for i in range(len(left_interpolated_interface_x)):
    dist = abs(left_interpolated_interface_x[i] - mid_pt)
    right_interpolated_interface_x.append(mid_pt + dist)
    right_interpolated_interface_y.append(left_interpolated_interface_y[i])
    
# plt.figure()
# plt.plot(left_interpolated_interface_x,left_interpolated_interface_y,'o')
# plt.plot(right_interpolated_interface_x,right_interpolated_interface_y,'o')

interface_x = np.concatenate([left_interpolated_interface_x,
                                          right_interpolated_interface_x])
interface_y = np.concatenate([left_interpolated_interface_y,
                                          right_interpolated_interface_y])


interface = np.column_stack([interface_x, interface_y])
interface = np.unique(interface, axis=0)

interface_x = interface[:,0]
interface_y = interface[:,1]

# Fits a circle to the numerical data gathered above

x_m = np.mean(interface_x)
y_m = np.mean(left_interpolated_interface_y)

plt.figure()
plt.plot(interface_x, interface_y, 'o')
plt.title("Droplet interface")
plt.savefig("droplet_interface.png")

#%% Fit a circle to the extracted interface using least squares

from scipy      import optimize

method_2 = "leastsq"

def calc_R(xc, yc):
    #calculate the distance of each 2D points from the center (xc, yc) 
    return np.sqrt((interface_x-xc)**2 + (interface_y-yc)**2)

def f_2(c):
    #calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) 
    Ri = calc_R(*c)
    return Ri - Ri.mean()

center_estimate = x_m, y_m
center_2, ier = optimize.leastsq(f_2, center_estimate)

x_c, y_c = center_2
Ri       = calc_R(*center_2)
Radius       = Ri.mean()
residu  = sum((Ri - Radius)**2)

y_min = min(interface_y)

tol = 1
hydrophilicity = ""
if y_min < y_c:
    hydrophilicity = "Hydrophobic"
elif y_min > y_c:
    hydrophilicity = "Hydrophilic"
elif abs(y_min - y_c) < tol:
    hydrophilicity = "Neither"
    

x_min = min(interface_x)
x_max = max(interface_x)
# Equation of circle fitted to the extracted droplet interface.
# 
x_fit = np.linspace(x_c - Radius**2, x_c + Radius**2, 5000000)
y_fit_top = y_c + np.sqrt(Radius**2 - (x_fit-x_c)**2)
y_fit_bottom = y_c + -np.sqrt(Radius**2 - (x_fit-x_c)**2)


x_fit = x_fit[~np.isnan(y_fit_top)]
y_fit_top = y_fit_top[~np.isnan(y_fit_top)]
y_fit_bottom = y_fit_bottom[~np.isnan(y_fit_bottom)]

circle_top = np.column_stack( [x_fit, y_fit_top] )
circle_bottom = np.column_stack( [x_fit, y_fit_bottom] )

circle_top = circle_top[~np.isnan(circle_top).any(axis=1)]
circle_bottom = circle_bottom[~np.isnan(circle_bottom).any(axis=1)]
x_fit = circle_bottom[:, 0]

smallest_y_index = np.argmin(left_interpolated_interface_y)
target = left_interpolated_interface_y[smallest_y_index]

if hydrophilicity == "Hydrophilic": # ie \theta < 90
    print("\n\n hydrophilic")
    closest_index = np.argmin(np.abs(y_fit_top - target))
    closest_val_y = y_fit_top[closest_index]
    closest_val_x = - np.sqrt( Radius**2 - (closest_val_y - y_c)**2 ) + x_c
    deriv = -(closest_val_x - x_c)/np.sqrt( Radius**2 - (closest_val_x - x_c)**2 )
    CA_1 = 180*np.arctan(deriv)/np.pi
elif hydrophilicity == "Hydrophobic": # ie \theta > 90
    print("\n\n hydrophobic")
    closest_index = np.argmin(np.abs(y_fit_bottom - target))
    closest_val_y = y_fit_bottom[closest_index]
    closest_val_x = - np.sqrt( Radius**2 - (closest_val_y - y_c)**2 ) + x_c
    deriv = (closest_val_x - x_c)/np.sqrt( Radius**2 - (closest_val_x - x_c)**2 )
    CA_1 = 180 + 180*np.arctan(deriv)/np.pi
elif hydrophilicity == "Neither":
    if abs(Radius**2 - (closest_val_x - x_c)**2) < 2*tol:
        deriv = "undefined"
        CA_1 = 90

    
print("theta = ", CA_1)

# print('\n circle center at (', x_c, ', ', y_c, ')')
# print('R = ', Radius)
# print('differentiation occurs at (', closest_val_x, ', ', closest_val_y)

[fig, ax] = plt.subplots()
ax.plot(interface_x, interface_y, 'o')
ax.plot(circle_top[:,0], circle_top[:,1], 'r')
ax.plot(circle_bottom[:, 0], circle_bottom[:, 1], 'r')
ax.annotate('**', xy=(closest_val_x, np.min(interface_y)) )







