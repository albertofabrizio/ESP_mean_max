import numpy as np
import ase.io as aio
from ase.io import cube
import scipy.spatial as sps

# Conversion Factor
bohr2ang=0.529177

# Set Filename
file="Wheeler_mep.cube"

#Dictionary of Vdw radii in Angs
radii_vdw={"H":1.2,"C":1.7,"O":1.52,"N":1.55}


# Read Atomic coordinates
nat = 9
filename = "Wheeler.xyz"
xyzfile = aio.read(filename,":")
coords = np.zeros((len(xyzfile),nat,3),float)
natoms = np.zeros((len(xyzfile)),int)
atomic_symbols = []
atomic_valence = []
for i in range(len(xyzfile)):
    coords[i] = np.asarray(xyzfile[i].get_positions(),float)/bohr2ang
    atomic_symbols.append(xyzfile[i].get_chemical_symbols())
    atomic_valence.append(xyzfile[i].get_atomic_numbers())
    natoms[i] = int(len(atomic_symbols[i]))
iconf = 0
coordinates = coords[iconf]
atoms = atomic_symbols[iconf]
valences = atomic_valence[iconf]

center=np.mean(coordinates[:5],axis=0)
index=np.linalg.norm(coordinates-center,axis=1).argmax()

# radial cutoff is the sum of position of the furthest atom from the center and its vdw radii
rcut=np.linalg.norm(coordinates[index]+radii_vdw[atoms[index]]/bohr2ang-center)

# Read Cube
nside = {}
cubefile = open(file ,"r")
lines = cubefile.readlines()
nside[0] =int(lines[3].split()[0])
nside[1] =int(lines[4].split()[0])
nside[2] =int(lines[5].split()[0])
npoints = 1
for i in range(3):
    npoints *= nside[i]
dx = float(lines[3].split()[1])
dy = float(lines[4].split()[2])
dz = float(lines[5].split()[3])
origin = np.asarray(lines[2].split(),dtype=float)[1:4]
cubefile.close()
cubefile = aio.read(file)
#density_regular=cube.read_cube_data(file)[0]
density_regular = np.concatenate(np.concatenate(cube.read_cube_data(file)[0]))

#construct the grid

grid_regular=np.transpose(np.asarray(np.meshgrid(dx*np.asarray(range(nside[0])),dy*np.asarray(range(nside[1])), dz*np.asarray(range(nside[2])) ) ),(2,1,3,0))
grid_regular=grid_regular.reshape((npoints,3))
grid_regular+=origin

# find the indices of the field at 3.25 A over the molecule
layer=3.25/bohr2ang
mask_1=grid_regular[:,2]>=layer-dz
mask_2=grid_regular[:,2]<=layer+dz

d2list = np.linalg.norm(grid_regular[:,:2]-center[:2],axis=1)
mask_4 = d2list <= rcut
mask_3=mask_1*mask_2*mask_4
amplitudes=density_regular[mask_3]
print("ESP mean ",np.mean(amplitudes)*627.51)
print("ESP max ",np.max(amplitudes)*627.51)
