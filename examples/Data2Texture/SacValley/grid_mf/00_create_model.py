import flopy
import numpy as np
import pandas as pd

ncol        = 508
nrow        = 960
xllcorner   = 530500.000000000000
yllcorner   = 4195005.951999999583
cellsize    = 321.868800000000

griddepth = pd.read_csv("../griddepth.csv")
top = np.loadtxt("../DEM_SV/top_median_fifthMile.asc", skiprows=6)
top = np.where(top<=-999999, 2000, top)
nlay = len(griddepth)
botm = top[np.newaxis,:,:] - griddepth.zbot.values.reshape([nlay, 1, 1])

mfname = "mf"
mf  = flopy.modflow.Modflow(mfname, )
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delc=cellsize, delr=cellsize, top=top, botm=botm, )
dis.top.format.fortran = '(508F10.3)'
for i in range(nlay):
    dis.botm[i].format.fortran = '(508F10.3)'
bas = flopy.modflow.ModflowBas(mf, ibound=["../grid/ibound.dat",]*nlay, strt=0.0)
lpf = flopy.modflow.ModflowLpf(mf, )
oc  = flopy.modflow.ModflowOc(mf, )
mf.write_input()
