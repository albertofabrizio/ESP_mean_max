import numpy as np
import ase.io as aio
from ase.io import cube

# Conversion Factor
bohr2ang=0.529177

file="Wheeler_mep.cube"

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
density_regular=cube.read_cube_data(file)[0]

density_regular.shape
