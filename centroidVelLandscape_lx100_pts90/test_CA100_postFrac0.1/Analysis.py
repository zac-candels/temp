import os, math, re, sys
import numpy as np
import struct
from matplotlib import pyplot as plt
import imageio.v2 as imageio  # for movie creation

# Directory where data is stored
datadir = "data/"

# Open the Header file
HeaderFile = open(datadir+"Header.mat", 'rb')

# Simulation size
LX = struct.unpack('=i', HeaderFile.read(4))[0]
LY = struct.unpack('=i', HeaderFile.read(4))[0]
LZ = struct.unpack('=i', HeaderFile.read(4))[0]

# Dimensionality
ndim = struct.unpack('=i', HeaderFile.read(4))[0]

# Time range
tstart = 0
tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

# Slice in the z direction to plot
zslice = int(LZ/2)

# Where to save the plots
outDirName = "figures"
os.makedirs(outDirName, exist_ok=True)

# Store the frame file paths for the movie
frame_files = []

for t in range(tstart, tend + 1, tinc):
    # Boundary file
    file_name = datadir + f"BoundaryLabels_t{t}.mat"
    with open(file_name, 'rb') as f:
        dat = f.read()
    solid = np.ndarray((LX, LY, LZ), '=i', dat, 0, (4 * LY * LZ, 4 * LZ, 4))

    # Order parameter
    file_name = datadir + f"OrderParameter_t{t}.mat"
    with open(file_name, 'rb') as f:
        dat = f.read()
    C = np.ndarray((LX, LY, LZ), '=d', dat, 0, (8 * LY * LZ, 8 * LZ, 8))
    liquid = np.array(C[:,:])
    liquid[np.logical_or(solid == 1, solid == -1)] = 0.5
    liquid[np.logical_or(solid == 3, solid == 2)] = 0.5

    # Viscous dissipation
    file_name = datadir + f"ViscousDissipation_t{t}.mat"
    with open(file_name, 'rb') as f:
        dat = f.read()
    D = np.ndarray((LX, LY, LZ), '=d', dat, 0, (8 * LY * LZ, 8 * LZ, 8))
    dissipation = np.array(D[:,:])
    dissipation[np.logical_or(solid == 1, solid == -1)] = 0
    dissipation[np.logical_or(solid == 3, solid == 2)] = 0

    # Velocity field
    file_name = datadir + f"Velocity_t{t}.mat"
    with open(file_name, 'rb') as f:
        dat = f.read()
    v = np.ndarray((LX, LY, LZ, ndim), '=d', dat, 0, (ndim * 8 * LY * LZ, ndim * 8 * LZ, ndim * 8, 8))
    print(f"max v at t={t}: ", np.amax(v))

    # Plot setup
    fig, ax = plt.subplots() #plt.subplots(2, 1, figsize=(6, 5))
    output = f"{outDirName}/component_plot_{t:012d}.png"

    rgbv = np.flip(liquid[:, :, zslice]).T
    rgbdiss = np.flip(dissipation[:, :, zslice]).T

    im = ax.imshow(np.fliplr(rgbv), interpolation='nearest', origin='upper')
    #im2 = ax[1].imshow(rgbdiss, interpolation='nearest', origin='upper')

    stepx, stepy = 1, 1
    X, Z = np.meshgrid(np.linspace(0, LX - 1, int(LX / stepx)), np.linspace(0, LY - 1, int(LY / stepy)))
    ax.quiver(X.T, Z.T,
                 np.flip(np.flip(-v[0:LX:stepx, 0:LY:stepy, zslice, 0], 0), 1),
                 np.flip(v[0:LX:stepx, 0:LY:stepy, zslice, 1]),
                 width=0.0008, headwidth=7.5, headlength=7.5)

    #fig.colorbar(im2, orientation="horizontal")
    plt.tight_layout()
    plt.savefig(output, dpi=400, format='png')
    plt.close(fig)

    frame_files.append(output)

# --- Create movie ---
import imageio_ffmpeg

movie_name = f"./simulation.mp4"
print(f"Creating movie {movie_name} ...")

# Ensure ffmpeg is available
fps = 10
writer = imageio.get_writer(movie_name, format='ffmpeg', mode='I', fps=fps)
for frame in frame_files:
    img = imageio.imread(frame)
    writer.append_data(img)
writer.close()

print("Movie created successfully!")