import pyPLUTO as pp
import PhyConst as phc
import read_opacities as rop
import os
import numpy as np
import asciidata as asc
import itertools
from scipy import integrate


##############--DATA ANALYSIS CLASS--##########################
#
# This Class has all the functions for the analysing the Data.
# and getting fancy quantities from data.
# CALLED AFTER pyPLUTO.pload object is defined.
#
################################################################


def disk_reso(Data,ul=None):
    print " This pro requires to have unit length of simumation."
    print "By Default its taken to be 1 AU"
    if ul==None : ul=1.0
    
    Delr = Data.dx1.max()*ul
    Deltheta = Data.x1[Data.n1-1]*Data.dx2.min()*ul
    DelPhi = Data.x1[Data.n1-1]*np.sin(Data.x2[0.5*Data.n2])*Data.dx3.min()*ul
    Rat1 = Delr/Delr
    Rat2 = Deltheta/Delr
    Rat3 = DelPhi/Delr
    print '-----------ANALYSIS RESULTS----------------'
    print 'Delta R = %f in the units of scale length ul'%(Delr)
    print 'Delta theta = %f in the units of scale length ul'%(Deltheta)
    print 'Delta phi = %f in the units of scale length ul'%(DelPhi)
    print 'THUS the ratio is :%f:%f:%f'%(Rat1,Rat2,Rat3)
    print '--------------------------------------------'


def OpMatrix(Data, Type=None,**kwargs):
    if Type==None: Type=0
    RefVel2 = (phc.G*phc.Msun*kwargs.get('Mstar',10.0))/(phc.au*kwargs.get('ul',1.0))
    RefT = kwargs.get('mu',2.353)*(phc.mH/phc.kB)*RefVel2
    Tmatrix = RefT*(Data.pr/Data.rho)
    MatrixOp = np.zeros(Tmatrix.shape)
    for i in range(Data.n1):
        for j in range(Data.n2):
            for k in range(Data.n3):
                MatrixOp[i,j,k] = rop.SemenovMeanOpacity(Temp=Tmatrix[i,j,k],Density=Data.rho[i,j,k],Type=Type)
        print i
    

    return MatrixOp



def Analysis(filepath=None,info=None):
    Values = np.asarray(asc.open(filepath+'analysis.out'))
    [a,b]=Values.shape
    
    MyArr = []
    finfo = open(filepath+'analysis.info','r')
    for line in finfo.readlines():
        if line.find('column')>=0:
            MyArr.append(line.split())

    print len(MyArr)
    ana_dict= dict([('col'+ str(i+1),[Values[i],' '.join(MyArr[i][3:])]) for i in range(a-1)])
    
    
    
    return ana_dict

class Rad_Average(object):

    def Sigma(self,Data,**kwargs):
        Sig = np.zeros([Data.x1.shape][0])
        for i in range(Data.n1):
            if (Data.n3 == 1):
                Sig[i] = (Data.rho[i,:]*Data.x1[i]*np.sin(Data.x2)*Data.dx2).sum()
            else:
                Sig[i] = (Data.rho[i,:,10]*Data.x1[i]*np.sin(Data.x2)*Data.dx2).sum()
                
        PhySig = kwargs.get('urho',1.0e-8)*phc.au*kwargs.get('ul',1.0)
       
        return PhySig*Sig

    def Pressure(self,Data,**kwargs):
        RefVel = np.sqrt(phc.G*phc.Msun/phc.au)*np.sqrt(kwargs.get('Mstar',10.0)/kwargs.get('ul',1.0))
        Press = np.zeros([Data.x1.shape][0])
        for i in range(Data.n1):
            if (Data.n3 == 1):
                Press[i] = (Data.pr[i,:]*Data.x1[i]*np.sin(Data.x2)*Data.dx2).sum()
            else:
                Press[i] = (Data.pr[i,:,10]*Data.x1[i]*np.sin(Data.x2)*Data.dx2).sum()
            
        PhyPress = kwargs.get('urho',1.0e-8)*RefVel*RefVel*phc.au*kwargs.get('ul',1.0)
        return PhyPress*Press

    def Csound(self,Data,**kwargs):
        RefVel = np.sqrt(phc.G*phc.Msun/phc.au)*np.sqrt(kwargs.get('Mstar',10.0)/kwargs.get('ul',1.0))
        P = self.Pressure(Data,**kwargs)
        S = self.Sigma(Data,**kwargs)
        Cs = kwargs.get('Gammae',1.0001)*np.sqrt(P/S)
        return Cs*1.0e-5

    def Temp(self,Data,**kwargs):
        Cs = self.Csound(Data,**kwargs)*1.0e5
        mu = kwargs.get('mu',2.353)
        KELVIN = phc.kB/(mu*phc.mH)
        Temperature = Cs*Cs/KELVIN
        return Temperature

    def Omega(self,Data,**kwargs):
       RefVel = np.sqrt(phc.G*phc.Msun/phc.au)*np.sqrt(kwargs.get('Mstar',10.0)/kwargs.get('ul',1.0))
       OmInt = np.zeros([Data.x1.shape][0])
       for i in range(Data.n1):
           if (Data.n3 == 1):
               OmInt[i] = (Data.v3[i,Data.n2/2]/Data.x1[i])
           else:
               OmInt[i] = (Data.v3[i,Data.n2/2,10]/Data.x1[i])
           
       PhyOm = RefVel/(phc.au*kwargs.get('ul',1.0))
       return PhyOm*OmInt

    def ToomreQ(self,Data,**kwargs):
       Cs = self.Csound(Data,**kwargs)*1.0e5
       Om = self.Omega(Data,**kwargs)
       S = self.Sigma(Data,**kwargs)
       Q = (Cs*Om)/(2.0*np.pi*phc.G*S)
       #Q = np.zeros([Data.x1.shape][0])
       #for i in range(Data.n1):
       #    Q[i] = Cs[i]*Om[i]/(2.0*np.pi*phc.G*S[i]).sum()
       
       return Q

    def cooltime(self,Data,**kwargs):
        Cs = self.Csound(Data,**kwargs)*1.0e5
        S = self.Sigma(Data,**kwargs)
        Trad = self.Temp(Data,**kwargs)
        #tau = self.OpDepth(Data,**kwargs)
        Constants = (phc.kB/(kwargs.get('mu',2.353)*phc.mH*2.0*phc.sigma))*(1.0/(kwargs.get('Gammae',5./3.) - 1.0))
        ftau = 1.0#tau + 1.0/tau
        tcool = (ftau*Constants*S)/Trad**3

        return tcool

            
            
            
            
class Vol_Average(object):
    def Quantities(self, *args,**kwargs):
        RefVel = np.sqrt(phc.G*phc.Msun/phc.au)*np.sqrt(kwargs.get('Mstar',10.0)/kwargs.get('ul',1.0))
        nstep = args[0]+1
        wdir = args[1]
        Timeyr = np.zeros(nstep)
        Mdisk = np.zeros(nstep)
        dV = np.zeros(nstep)
        IntEdisk = np.zeros(nstep)
        Sigdisk = np.zeros(nstep)
        Csdisk = np.zeros(nstep)
        Omdisk = np.zeros(nstep)
        Qdisk = np.zeros(nstep)
        for ns in range(nstep):
            D = pp.pload(ns,w_dir=wdir)
            dV = np.zeros(D.rho.shape)
            Sigma = np.zeros(D.rho.shape)
            if (D.n3 == 1):
                dV = D.x1[:,np.newaxis]*D.x1[:,np.newaxis]*np.sin(D.x2[np.newaxis,:])*D.dx1[:,np.newaxis]*D.dx2[np.newaxis,:]*D.dx3
                Sigma = (D.rho*dV)*(D.rho*D.x1[:,np.newaxis]*np.sin(D.x2[np.newaxis,:])*D.dx2[np.newaxis,:])
            else:
                dV = D.x1[:,np.newaxis,np.newaxis]*D.x1[:,np.newaxis,np.newaxis]*np.sin(D.x2[np.newaxis,:,np.newaxis])*D.dx1[:,np.newaxis,np.newaxis]*D.dx2[np.newaxis,:,np.newaxis]*D.dx3[np.newaxis,np.newaxis,:]
                Sigma = (D.rho*dV)*(D.rho*D.x1[:,np.newaxis,np.newaxis]*np.sin(D.x2[np.newaxis,:,np.newaxis])*D.dx2[np.newaxis,:,np.newaxis])
            
            Csound = (D.rho*dV)*(np.sqrt(kwargs.get('Gammae',1.0001)*(D.pr/D.rho)))
            IntE = (Csound*Csound)/(kwargs.get('Gammae',1.0001)*((kwargs.get('Gammae',1.0001)-1.0)))
            Omega = (D.rho*dV)*(D.v3/D.x1[:,np.newaxis,np.newaxis])

            Timeyr[ns] = ns*2.0*np.pi*((kwargs.get('ul',1.0)*phc.au)/RefVel)/phc.yr
            Mdisk[ns] = ((D.rho*dV).sum())*((kwargs.get('urho',1.0e-8)*(kwargs.get('ul',1.0)*phc.au)**3)/phc.Msun)
            Sigdisk[ns] = (1.0/(D.rho*dV).sum())*(Sigma.sum())*(kwargs.get('urho',1.0e-8)*(kwargs.get('ul',1.0)*phc.au))
            Csdisk[ns] = (1.0/(D.rho*dV).sum())*(Csound.sum())*RefVel
            IntEdisk[ns] = (1.0/(D.rho*dV).sum())*(IntE.sum())*(RefVel**2)
            Omdisk[ns] = (1.0/(D.rho*dV).sum())*(Omega.sum())*(RefVel/(kwargs.get('ul',1.0)*phc.au))
            Qdisk[ns] = Csdisk[ns]*Omdisk[ns]/(2.0*np.pi*phc.G*Sigdisk[ns])

        return {'Time':Timeyr,'Mdisk':Mdisk,'Sigma':Sigdisk,'Csound':Csdisk*1.0e-5,'Omega':Omdisk,'ToomreQ':Qdisk, 'IntE':IntEdisk}

    



class Stress(object):
    def Reynold(self,Data,Gammae=1.0001):
        D = Data
        T = pp.Tools()
        dOmdr = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        VertAvg_Sg = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        VertAvg_P = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        VertAvg_Vr  = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        Delta_Vr  = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        VertAvg_Vphi  = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        Delta_Vphi  = np.zeros([D.rho.shape[0],D.rho.shape[2]])
        AziAvg_Vphi =  np.zeros(D.rho.shape[0])
        AziAvg_Sg =  np.zeros(D.rho.shape[0])
        AziAvg_Vr =  np.zeros(D.rho.shape[0])
        AziAvg_Fa = np.zeros(D.rho.shape[0])
        AziAvg_Fr = np.zeros(D.rho.shape[0])
        AziAvg_Deno = np.zeros(D.rho.shape[0])
        
        dV = D.x1[:,np.newaxis,np.newaxis]*D.x1[:,np.newaxis,np.newaxis]*np.sin(D.x2[np.newaxis,:,np.newaxis])*D.dx1[:,np.newaxis,np.newaxis]*D.dx2[np.newaxis,:,np.newaxis]*D.dx3[np.newaxis,np.newaxis,:]
        Mdisk = ((D.rho*dV).sum())
        Vr_avg = (D.v1*D.rho*dV).sum()/Mdisk
        Vp_avg = (D.v3*D.rho*dV).sum()/Mdisk

        for k in range(Data.n3):
            for i in range(Data.n1):
                VertAvg_Sg[i,k] = (D.rho[i,:,k]*D.x1[i]*np.sin(D.x2)*D.dx2).sum()
                VertAvg_P[i,k] = (D.pr[i,:,k]*D.x1[i]*np.sin(D.x2)*D.dx2).sum()
                VertAvg_Vr[i,k] = (1.0/VertAvg_Sg[i,k])*(D.rho[i,:,k]*(D.v1[i,:,k]-Vr_avg)*D.x1[i]*np.sin(D.x2)*D.dx2).sum()
                VertAvg_Vphi[i,k] = (1.0/VertAvg_Sg[i,k])*(D.rho[i,:,k]*(D.v3[i,:,k]-Vp_avg)*D.x1[i]*np.sin(D.x2)*D.dx2).sum()

        [r2d,phi2d]= np.meshgrid(D.x1,D.x3)
        r2d = r2d.T
        Fa = r2d*VertAvg_Sg*VertAvg_Vr*VertAvg_Vphi
        
        for i in range(Data.n1):
            AziAvg_Sg[i] = (VertAvg_Sg[i,:]*D.dx3).sum()
            AziAvg_Vr[i] = (VertAvg_Vr[i,:]*D.dx3).sum()
            AziAvg_Vphi[i] = (VertAvg_Vphi[i,:]*D.dx3).sum()
            Delta_Vr[i,:] =  VertAvg_Vr[i,:] - AziAvg_Vr[i]
            Delta_Vphi[i,:] =  VertAvg_Vr[i,:] - AziAvg_Vphi[i]

        Fr = r2d*VertAvg_Sg*Delta_Vr*Delta_Vphi
        Omega = VertAvg_Vphi/r2d
        Ciso2 = Gammae*(VertAvg_P/VertAvg_Sg)
        

        for k in range(D.n3):
            dOmdr[:,k] = (r2d[:,k]/Omega[:,k])*T.deriv(Omega[:,k],D.x1)

        Deno =  r2d*VertAvg_Sg*Ciso2*np.abs(dOmdr)

        for i in range(D.n1):
            AziAvg_Fa[i] = (Fa[i,:]*D.dx3).sum()
            AziAvg_Fr[i] = (Fr[i,:]*D.dx3).sum()
            AziAvg_Deno[i] = (Deno[i,:]*D.dx3).sum()

        Alpha_A = AziAvg_Fa/AziAvg_Deno
        Alpha_R = AziAvg_Fr/AziAvg_Deno
        
        return {'AlphaA': Alpha_A,'AlphaR': Alpha_R, 'Fa':AziAvg_Fa, 'Fr': AziAvg_Fr, 'Deno':AziAvg_Deno}
        
        

        
            
            

        

        return {'SgRphi':VertAvg_Sg, 'VrRphi':VertAvg_Vr, 'VpRphi':VertAvg_Vphi}
    
                

            
    

    
class Force_Ana(object):
    def Gravity(self,Data):
        [r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
        r2d=r2d.T
        z2d=z2d.T
        Tool = pp.Tools()
        Grav_force_dict = {}
        Grav_force_dict['G_r'] = 1.0/(r2d*r2d)
        Grav_force_dict['G_z'] = np.zeros(r2d.shape)

        return Grav_force_dict

    def Pressure(self,Data,phi=10):
        Tool = pp.Tools()
        Prgrad = Tool.Grad(Data.pr[:,:,phi],Data.x1,Data.x2,Data.dx1,Data.dx2,polar=True)
        Press_force_dict ={}
        Press_force_dict['Fp_r'] = -1.0*(Prgrad[:,:,0]/Data.rho[:,:,phi])
        Press_force_dict['Fp_th'] = -1.0*(Prgrad[:,:,1]/Data.rho[:,:,phi])

        return Press_force_dict

    def RadPressure(self,Data,phi=10,ul=1.0,urho=1.0e-8,Mstar=10.0):
        Tool = pp.Tools()
        R_GasC = phc.NA*phc.kB
        mu_HHe = 2.353
        RefLength = ul*phc.au
        RefDensity = urho
        RefVel = np.sqrt((phc.G*phc.Msun*Mstar)/RefLength)
        RefTemp = (RefVel**2)*(mu_HHe/R_GasC)
        RefEnergy = RefDensity*(RefLength**3)*(RefVel**2)

        radconstant = 4.0*(phc.sigma)/(phc.c)
        dimless_radconst = radconstant*(RefTemp**4.0)*(RefLength**3.0)/(RefEnergy)
        
        Radpr = (1.0/3.0)*dimless_radconst*((Data.Temp[:,:,phi]/RefTemp)**4.0)
        RadPrgrad = Tool.Grad(Radpr,Data.x1,Data.x2,Data.dx1,Data.dx2,polar=True)
        Rad_Press_force_dict ={}
        Rad_Press_force_dict['RFp_r'] = -1.0*(RadPrgrad[:,:,0]/Data.rho[:,:,phi])
        Rad_Press_force_dict['RFp_th'] = -1.0*(RadPrgrad[:,:,1]/Data.rho[:,:,phi])

        return Rad_Press_force_dict

    def Centrifugal(self,Data,phi=10):
        [r2d,z2d] = np.meshgrid(Data.x1,Data.x2)
        r2d=r2d.T
        z2d=z2d.T
        Centri_force_dict={}
        Centri_force_dict['Fcf_r'] = (Data.v3[:,:,10]*Data.v3[:,:,10])/r2d
        Centri_force_dict['Fcf_z'] = np.zeros(r2d.shape)
	
        return Centri_force_dict


        
    
    
    
