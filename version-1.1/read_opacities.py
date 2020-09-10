import numpy as np
import asciidata as asc
from scipy.interpolate import RectBivariateSpline


Opac_Dir = '/Users/bhargavvaidya/Simulations/Disk_setup/Makemake/Opacities/Mean_Opacity/Semenov_Opacity/'
Opac_Dir = '/home/vaidya/SG_RAD_PLUTO_Ver2010-12-07/Makemake/Opacities/Mean_Opacity/Semenov_Opacity/'
Opac_Dir = '/ptmp/mpia/bvaidya/Massive_Disk/SGRAD_PLUTO_Ver07122010/Makemake/Opacities/Mean_Opacity/Semenov_Opacity/'

RossFileName = 'Total_Rossland.out'
PlanckFileName = 'Total_Planck.out'
ConstantGasOpacity = 0.01

def SemenovMeanOpacity(Temp = None, Density = None,Type=None):
    # Type deteremines the mean opacity type, 0 : Rosseland , 1: Planck.

    if Type == 0 or Type == None:
        Whole_Table = np.array(asc.open(Opac_Dir+RossFileName))
    else:
        Whole_Table = np.array(asc.open(Opac_Dir+PlanckFileName))

    
    rho = Whole_Table[1:,0]
    T = Whole_Table[0,1:]
    OpMatrix = Whole_Table[1:,1:]


    LgInputTemp = np.log10(Temp)
    LgInputDens = np.log10(Density)

    pt1 = (np.abs(rho-LgInputDens)).argmin()
    pt2 = (np.abs(T-LgInputTemp)).argmin()
    
    if LgInputTemp > 3.87:
        interpval = ConstantGasOpacity
    else:
        outgrid = RectBivariateSpline(rho,T,OpMatrix,kx=1,ky=1)
        interpval = outgrid(rho[pt1],T[pt2])[0][0]

    return interpval

    
        



    
