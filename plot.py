import numpy as np
import matplotlib.pyplot as pp

# getting params
params = np.loadtxt('cmake-build-release/params.txt')
Ni = int(params[0])
Nj = int(params[1])
Nk = int(params[2])
L = params[3]
m = 4.245404
# load in data
halfi = int(Ni/2)
halfj = int(Nj/2)
data = np.loadtxt('cmake-build-release/data.txt')
data = data.reshape(Ni, Nj, Nk)

plane = data[halfi]
line = plane[halfj]

del data

data2 = np.loadtxt('cmake-build-release/data2.txt')
data2 = data2.reshape(Ni, Nj, Nk)

plane2 = data2[halfi]
line2 = plane2[halfj]

del data2


data3 = np.loadtxt('cmake-build-release/data3.txt')
data3 = data3.reshape(Ni, Nj, Nk)

plane3 = data3[halfi]
line3 = plane3[halfj]

del data3

# load in density
dens_plot = np.loadtxt('cmake-build-release/f.txt')
dens_plot = dens_plot.reshape(Ni, Nj, Nk)

# plot plane


x = np.linspace(-L / 2, L / 2, Nk)
pp.plot(x, line, 'b')
pp.plot(x, line2, 'r')
pp.plot(x, line3, 'g')
pp.show()

pp.pcolor(plane)
pp.show()
pp.pcolor(plane2)
pp.show()
pp.pcolor(plane3)
pp.show()

# plot density
dens_plane = dens_plot[int(Ni/2)]
dens_line = dens_plane[int(Nj/2)]
pp.plot(x, dens_line)
pp.show()
pp.pcolor(dens_plane)
pp.show()
