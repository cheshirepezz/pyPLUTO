import pyPLUTO as pp
import PhyConst as phc
import os
import numpy as np
import asciidata
from scipy import integrate





##############--DATA ANALYSIS CLASS--##########################
#
# This Class has all the functions for the analysing the Data.
# and getting fancy quantities from data.
# CALLED AFTER pyPLUTO.pload object is defined.
#
################################################################
class Force(object):
	def al_perp(self,Data):
		[r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T
		
		Tool = pp.Tools()
		
		GrdA3 = Tool.Grad(r2d*Data.A3,Data.x1,Data.x2,Data.dx1,Data.dx2)
		magGrdA3 = np.sqrt(GrdA3[:,:,0]**2 + GrdA3[:,:,1]**2)
		grda3_dict={}
		grda3_dict['GrdA3r']=GrdA3[:,:,0]
		grda3_dict['GrdA3z']=GrdA3[:,:,1]
		grda3_dict['magGrdA3']=magGrdA3

		return grda3_dict

	def al_para(self,Data):
		Br = Data.b1
		Bz = Data.b2
		Bpol = np.sqrt(Data.b1**2 + Data.b2**2)

		Bpara_dict = {}
		Bpara_dict['Br'] = Br
		Bpara_dict['Bz'] = Bz
		Bpara_dict['Bpol'] = Bpol

		return Bpara_dict
		
	def Gravity(self,Data):
		rg = 0.21
		zg = 0.21
		[r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T

		Tool = pp.Tools()
		
		gravdeno = ((r2d + rg)**2 + (z2d + zg)**2)**(1.5)
		gravphi = -1.0/(((r2d + rg)**2 + (z2d + zg)**2)**(0.5))
		gradphi = Tool.Grad(gravphi,Data.x1,Data.x2,Data.dx1,Data.dx2)

		grda3 = self.al_perp(Data)
		Bpara = self.al_para(Data)

		Grav_afl = (1.0/(Bpara['Bpol']))*(Bpara['Br']*gradphi[:,:,0] + Bpara['Bz']*gradphi[:,:,1])
		Grav_tfl = (1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*gradphi[:,:,0] + grda3['GrdA3z']*gradphi[:,:,1])
		
		
		Grav_force_dict = {}
		Grav_force_dict['G_r'] = (1.0*(r2d+rg))/(gravdeno)
		Grav_force_dict['G_z'] = (1.0*(z2d+zg))/(gravdeno)
		Grav_force_dict['Grav_tfl']=Data.rho*Grav_tfl
		Grav_force_dict['Grav_afl']=Data.rho*Grav_afl
		
		return Grav_force_dict
	def Pressure(self,Data):
		Tool = pp.Tools()
		Prgrad = Tool.Grad(Data.pr,Data.x1,Data.x2,Data.dx1,Data.dx2)

		grda3 = self.al_perp(Data)
		Press_tfl = (1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*Prgrad[:,:,0] + grda3['GrdA3z']*Prgrad[:,:,1])

		Bpara = self.al_para(Data)
		Press_afl = (1.0/(Bpara['Bpol']))*(Bpara['Br']*Prgrad[:,:,0] + Bpara['Bz']*Prgrad[:,:,1])
		
		Press_force_dict ={}
		Press_force_dict['Fp_r'] = -1.0*(Prgrad[:,:,0]/Data.rho)
		Press_force_dict['Fp_z'] = -1.0*(Prgrad[:,:,1]/Data.rho)
		Press_force_dict['Press_tfl']= Press_tfl
		Press_force_dict['Press_afl']= Press_afl
		
		return Press_force_dict
	def Centrifugal(self,Data):
		[r2d,z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T

		grda3 = self.al_perp(Data)
		Centri_tfl = (1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*((Data.v3*Data.v3)/r2d))

		Bpara = self.al_para(Data)
		Centri_afl = (1.0/(Bpara['Bpol']))*(Bpara['Br']*((Data.v3*Data.v3)/r2d))
		
		Centri_force_dict={}
		Centri_force_dict['Fcf_r'] = (Data.v3*Data.v3)/r2d
		Centri_force_dict['Fcf_z'] = np.zeros(r2d.shape)
		Centri_force_dict['Centri_tfl'] = Data.rho*Centri_tfl
		Centri_force_dict['Centri_afl'] = Data.rho*Centri_afl
		
		return Centri_force_dict

	def Mag_Pressure(self,Data):
		magpol =  np.sqrt(Data.b1**2 + Data.b2**2)
		magpr = 0.5*magpol**2
		bphipr = 0.5*Data.b3**2
		
		Tool = pp.Tools()
		Grdpolmagpr = Tool.Grad(magpr,Data.x1,Data.x2,Data.dx1,Data.dx2)
		Grdphimagpr = Tool.Grad(bphipr,Data.x1,Data.x2,Data.dx1,Data.dx2)

		grda3 = self.al_perp(Data)
		PolMagpr_tfl = (1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*Grdpolmagpr[:,:,0] + grda3['GrdA3z']*Grdpolmagpr[:,:,1])
		PhiMagpr_tfl = (1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*Grdphimagpr[:,:,0] + grda3['GrdA3z']*Grdphimagpr[:,:,1])
		Bpara = self.al_para(Data)
		PolMagpr_afl = (1.0/(Bpara['Bpol']))*(Bpara['Br']*Grdpolmagpr[:,:,0] + Bpara['Bz']*Grdpolmagpr[:,:,1])
		PhiMagpr_afl = (1.0/(Bpara['Bpol']))*(Bpara['Br']*Grdphimagpr[:,:,0] + Bpara['Bz']*Grdphimagpr[:,:,1])
		

		
		[r2d,z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T
	
		pinch = 2.0*(bphipr)/r2d
		pinch_tfl = (1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*pinch)
		pinch_afl = (1.0/(Bpara['Br']))*(Bpara['Bpol']*pinch)
		

		Magnetic_pressure_dict={}
		Magnetic_pressure_dict['bpolpr_tfl']=PolMagpr_tfl
		Magnetic_pressure_dict['bphipr_tfl']=PhiMagpr_tfl
		Magnetic_pressure_dict['pinch_tfl'] = pinch_tfl
		Magnetic_pressure_dict['bpolpr_afl']=PolMagpr_afl
		Magnetic_pressure_dict['bphipr_afl']=PhiMagpr_afl
		Magnetic_pressure_dict['pinch_afl'] = pinch_afl

		return Magnetic_pressure_dict
		

		
	def Lorentz(self,Data):	
		[r2d,z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T

		Tool = pp.Tools()
		grb1 = Tool.Grad(Data.b1,Data.x1,Data.x2,Data.dx1,Data.dx2)
		grb2 = Tool.Grad(Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2)
		grb3 = Tool.Grad(Data.b3,Data.x1,Data.x2,Data.dx1,Data.dx2)
		grI =  Tool.Grad(r2d*Data.b3,Data.x1,Data.x2,Data.dx1,Data.dx2) # This is gradient of current r*Bphi used to estimate Jz

		Jr  = -grb3[:,:,1]
		Jphi=  grb1[:,:,1] - grb2[:,:,0]
		Jz = (1.0/r2d)*grI[:,:,0]

		Loren_force_dict={}
		Loren_force_dict['Fl_r']= (Jphi*Data.b2 - Jz*Data.b3)/Data.rho
		Loren_force_dict['Fl_z']= (Jr*Data.b3 - Jphi*Data.b1)/Data.rho
		Loren_force_dict['Fl_phi']=(Jz*Data.b1 - Jr*Data.b2)/Data.rho
		return Loren_force_dict
	def Stellar_Rad(self,Data,**kwargs):
		[r2d,z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T
		Gravforce = self.Gravity(Data)

		Mstar = kwargs.get('Mstar',30.0)
		urho = kwargs.get('urho',5.0e-14)
		ul = kwargs.get('ul',1.0)
		Gammae = kwargs.get('Gammae',0.2369)
		Qo = kwargs.get('Qo',1400.0)
		Alpha = kwargs.get('Alpha',0.55)

		print '-----------------------------------------------'
		print 'xfl   : ',kwargs.get('xfl',5.0)
		print 'Alpha : ',Alpha
		print 'Gammae: ',Gammae
		print 'Qo    : ',Qo
		print 'ul    : ',ul
		print 'urho  : ',urho
		print 'Mstar : ',Mstar
		print '-----------------------------------------------'
		
		sigmae = 0.4
		clight = 3.0e10
		G = 6.67e-8
		Msun = 2.0e33
		AU=1.5e13
		uvel = np.sqrt((G*Mstar*Msun)/(ul*AU))
		Dless = uvel/(urho*ul*AU*sigmae*clight)
		Kpara = (Dless**(Alpha))*((Qo**(1.0-Alpha))/(1.0-Alpha))

		
		Tool = pp.Tools()
		grv1 = Tool.Grad(Data.v1,Data.x1,Data.x2,Data.dx1,Data.dx2)
		grv2 = Tool.Grad(Data.v2,Data.x1,Data.x2,Data.dx1,Data.dx2)
		DvrDr = np.abs(grv1[:,:,0])
		DvrDz = np.abs(grv1[:,:,1])
		DvzDr = np.abs(grv2[:,:,0])
		DvzDz = np.abs(grv2[:,:,1])
		xrat = z2d/r2d
		prf = 1.0/(1.0+xrat**2)
		
		
		dvdl = prf*(DvrDr + (xrat**2)*DvzDz + xrat*(DvrDz + DvzDr)) 
		Mt = Kpara*((1.0/Data.rho)*dvdl)**(Alpha)
		
		Rad_r = Mt*Gammae*Gravforce['G_r']
		Rad_z = Mt*Gammae*Gravforce['G_z']

		grda3 = self.al_perp(Data)
		StRad_tfl = Data.rho*(1.0/(grda3['magGrdA3']))*(grda3['GrdA3r']*Rad_r + grda3['GrdA3z']*Rad_z)

		Bpara = self.al_para(Data)
		StRad_afl = Data.rho*(1.0/(Bpara['Bpol']))*(Bpara['Br']*Rad_r + Bpara['Bz']*Rad_z)

		Rad_force_dict={'dvdl':dvdl,'Mt':Mt,'Fr_r':Rad_r,'Fr_z':Rad_z,'StRad_tfl':StRad_tfl, 'StRad_afl':StRad_afl}
		return Rad_force_dict

	def Disk_Rad(self,Data,**kwargs):
		Ld = asciidata.open(kwargs.get('file','/Users/bhargavvaidya/test_linediskrpcor_55.dat'))
		r2d = np.asarray(Ld[0]).reshape(516,1028)
		z2d = np.asarray(Ld[1]).reshape(516,1028)
		Srl = np.asarray(Ld[2]).reshape(516,1028)
		Svl = np.asarray(Ld[3]).reshape(516,1028)
		Mstar = kwargs.get('Mstar',30.0)
		urho = kwargs.get('urho',5.0e-14)
		ul = kwargs.get('ul',0.1)
		Gammae = kwargs.get('Gammae',0.2369)
		Zeta = kwargs.get('Zeta',0.4644)
		Lambda = kwargs.get('Lambda',0.4969)
		Qo = kwargs.get('Qo',1400.0)
		Alpha = kwargs.get('Alpha',0.55)

		print '-----------------------------------------------'
		print 'xfl   : ',kwargs.get('xfl',5.0)
		print 'Alpha : ',Alpha
		print 'Gammae: ',Gammae
		print 'Zeta  : ',Zeta
		print 'Lambda: ',Lambda
		print 'Qo    : ',Qo
		print 'ul    : ',ul
		print 'urho  : ',urho
		print 'Mstar : ',Mstar
		print '-----------------------------------------------'
		
		sigmae = 0.4
		clight = 3.0e10
		G = 6.67e-8
		Msun = 2.0e33
		AU=1.5e13
		uvel = np.sqrt((G*Mstar*Msun)/(ul*AU))
		Dless = uvel/(urho*ul*AU*sigmae*clight)
		prefactor = (3.0/2.0)*(1.0/np.pi)*Gammae*Zeta*Lambda
		Kpara = (Dless**(Alpha))*((Qo**(1.0-Alpha))/(1.0-Alpha))

		Tool = pp.Tools()
		grv2 = Tool.Grad(Data.v2,Data.x1,Data.x2,Data.dx1,Data.dx2)
		DvzDz = np.abs(grv2[:,:,1])

		dvdl = DvzDz 
		Disk_Mt = Kpara*((1.0/Data.rho)*dvdl)**(Alpha)
		
		Disk_Rad_r = Disk_Mt*prefactor*Srl[2:514,2:1026]
		Disk_Rad_z = Disk_Mt*prefactor*Svl[2:514,2:1026]

		DiskRad_force_dict={'d_dvdl':dvdl,'d_Mt':Disk_Mt,'d_Fr_r':Disk_Rad_r,'d_Fr_z':Disk_Rad_z}
		return DiskRad_force_dict

	def proj_force(self,Data,CompX,CompY,**kwargs):
	    Flr = CompX
	    Flz = CompY
	    magpol_p = np.sqrt(Flr*Flr + Flz*Flz)
	    phi = np.arctan2(Data.b2,Data.b1)
	    theta = np.zeros(phi.shape)
	    Flper_r_p=np.zeros(phi.shape)
	    Flper_z_p=np.zeros(phi.shape)
	    magFlper=np.zeros(phi.shape)
	    alpha = np.arctan2(Flz,Flr)
	    for i in range(phi.shape[0]):
		    for j in range(phi.shape[1]):
			if (alpha[i,j] > phi[i,j]):
			    theta[i,j] =  (alpha[i,j] - phi[i,j])
			    magFlper[i,j] = magpol_p[i,j]*np.sin(theta[i,j])
			    Flper_r_p[i,j] = -1.0*magFlper[i,j]*np.sin(phi[i,j])
			    Flper_z_p[i,j] = magFlper[i,j]*np.cos(phi[i,j])
			else:
			    theta[i,j] =  2.0*np.pi + (alpha[i,j] - phi[i,j])
			    magFlper[i,j] = magpol_p[i,j]*np.sin(theta[i,j])
			    Flper_r_p[i,j] = -1.0*magFlper[i,j]*np.sin(phi[i,j])
			    Flper_z_p[i,j] = magFlper[i,j]*np.cos(phi[i,j])


	    magFlpara =  magpol_p*np.cos(theta)
	    Flpara_r_p = magFlpara*np.cos(phi)
	    Flpara_z_p = magFlpara*np.sin(phi)
	    magFlpara_p = np.sqrt(Flpara_r_p*Flpara_r_p + Flpara_z_p*Flpara_z_p)
	    magFlper_p = np.sqrt(Flper_r_p*Flper_r_p + Flper_z_p*Flper_z_p)
	    Dummy = np.zeros(Data.rho.shape)
	    para_flvalues_dict = self.newFline_vals(Data,magFlpara_p,Dummy,**kwargs)
	    perp_flvalues_dict = self.newFline_vals(Data,magFlper_p,Dummy,**kwargs)
	    Qy = para_flvalues_dict['Qy']
	    para_flvalues=para_flvalues_dict['Fl_Val']
	    perp_flvalues=perp_flvalues_dict['Fl_Val']
	    
	    
	    

	    projection_dict={'Magnitude': magpol_p,'Al_Para_r':Flpara_r_p,'Al_Para_z':Flpara_z_p,'Al_Per_r':Flper_r_p,'Al_Per_z':Flper_z_p, 'para_flvalues': para_flvalues,'perp_flvalues': perp_flvalues,'Qy':Qy}

	    return projection_dict

	def newFline_vals(self,Data,CompX,CompY,**kwargs):
		Im = pp.Image()
		fl_dict=Im.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2,kwargs.get('xfl',5.0),0.001)
		Qx = fl_dict['qx']
		Qy = fl_dict['qy']

		Final_Values = np.zeros(shape=(2,len(Qx)))

		for i in range(len(Qx)):
			Final_Values[:,i] = Im.field_interp(CompX,CompY,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])

		Final_Values_dict={'Qy':Qy,'Fl_Val':Final_Values}
		return Final_Values_dict

		
		
		

	def Fline_values(self,Data,**kwargs):
		Im = pp.Image()
		fl_dict=Im.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2,kwargs.get('xfl',5.0),0.001)
		Qx = fl_dict['qx']
		Qy = fl_dict['qy']
		Gdict=self.Gravity(Data)
		Pdict=self.Pressure(Data)
		Cdict=self.Centrifugal(Data)
		Ldict=self.Lorentz(Data)
		MagPdict = self.Mag_Pressure(Data)
		StRdict=self.Stellar_Rad(Data,**kwargs)
		


	
				
		
		Dummy = np.zeros(Data.rho.shape)
		
		Fl_Gr = np.zeros(shape=(2,len(Qx)))
		tFl_Gr = np.zeros(shape=(2,len(Qx)))
		aFl_Gr = np.zeros(shape=(2,len(Qx)))
		
		Fl_Pr = np.zeros(shape=(2,len(Qx)))
		tFl_Pr = np.zeros(shape=(2,len(Qx)))
		aFl_Pr = np.zeros(shape=(2,len(Qx)))
		
		Fl_Cf = np.zeros(shape=(2,len(Qx)))
		tFl_Cf = np.zeros(shape=(2,len(Qx)))
		aFl_Cf = np.zeros(shape=(2,len(Qx)))
		
		Fl_Lf = np.zeros(shape=(2,len(Qx)))
		tFl_Lf1 = np.zeros(shape=(2,len(Qx)))
		tFl_Lf2 = np.zeros(shape=(2,len(Qx)))
		tFl_Pinch = np.zeros(shape=(2,len(Qx)))
		aFl_Lf1 = np.zeros(shape=(2,len(Qx)))
		aFl_Lf2 = np.zeros(shape=(2,len(Qx)))
		aFl_Pinch = np.zeros(shape=(2,len(Qx)))

		Fl_StRf = np.zeros(shape=(2,len(Qx)))
		tFl_StRf = np.zeros(shape=(2,len(Qx)))
		aFl_StRf = np.zeros(shape=(2,len(Qx)))

		Fl_DiskRf = np.zeros(shape=(2,len(Qx)))
		

		for i in range(len(Qx)):
			Fl_Gr[:,i] = Im.field_interp(Gdict['G_r'],Gdict['G_z'],Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			Fl_Pr[:,i] = Im.field_interp(Pdict['Fp_r'],Pdict['Fp_z'],Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			Fl_Cf[:,i] = Im.field_interp(Cdict['Fcf_r'],Cdict['Fcf_z'],Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			Fl_Lf[:,i] = Im.field_interp(Ldict['Fl_r'],Ldict['Fl_z'],Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			Fl_StRf[:,i] = Im.field_interp(StRdict['Fr_r'],StRdict['Fr_z'],Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])

			tFl_Gr[:,i] = Im.field_interp(Gdict['Grav_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			tFl_Pr[:,i] = Im.field_interp(Pdict['Press_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			tFl_Cf[:,i] = Im.field_interp(Cdict['Centri_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			tFl_Lf1[:,i] = Im.field_interp(MagPdict['bpolpr_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			tFl_Lf2[:,i] = Im.field_interp(MagPdict['bphipr_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			tFl_Pinch[:,i] = Im.field_interp(MagPdict['pinch_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			tFl_StRf[:,i] = Im.field_interp(StRdict['StRad_tfl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_Gr[:,i] = Im.field_interp(Gdict['Grav_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_Pr[:,i] = Im.field_interp(Pdict['Press_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_Cf[:,i] = Im.field_interp(Cdict['Centri_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_Lf1[:,i] = Im.field_interp(MagPdict['bpolpr_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_Lf2[:,i] = Im.field_interp(MagPdict['bphipr_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_Pinch[:,i] = Im.field_interp(MagPdict['pinch_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			aFl_StRf[:,i] = Im.field_interp(StRdict['StRad_afl'],Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			
		if kwargs.get('DiskRad',False) == True:
			DiskRdict = self.Disk_Rad(Data,**kwargs)
			for i in range(len(Qx)):
				Fl_DiskRf[:,i] = Im.field_interp(DiskRdict['d_Fr_r'],DiskRdict['d_Fr_z'],Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])

		var_fl = {}
		keys_list =['Qy','Fl_Gr','Fl_Pr','Fl_Cf','Fl_Lf','Fl_StRf','Fl_DiskRf','tFl_Gr','tFl_Pr','tFl_Cf','tFl_Lf1','tFl_Lf2','tFl_Pinch','tFl_StRf', 'aFl_Gr','aFl_Pr','aFl_Cf','aFl_Lf1','aFl_Lf2','aFl_Pinch','aFl_StRf']
		matrix_list =[Qy,Fl_Gr,Fl_Pr,Fl_Cf,Fl_Lf,Fl_StRf,Fl_DiskRf,tFl_Gr,tFl_Pr,tFl_Cf,tFl_Lf1,tFl_Lf2,tFl_Pinch,tFl_StRf,aFl_Gr,aFl_Pr,aFl_Cf,aFl_Lf1,aFl_Lf2,aFl_Pinch,aFl_StRf]

		for i in range(len(keys_list)):
			var_fl[keys_list[i]]=matrix_list[i]
		
		#print var_fl.keys()
		#var_fl={'Qy':Qy,'Fl_Gr':Fl_Gr,'Fl_Pr':Fl_Pr,'Fl_Cf':Fl_Cf,'Fl_Lf':Fl_Lf,'Fl_StRf':Fl_StRf,'Fl_DiskRf':Fl_DiskRf}
		return var_fl


	
	


class quantities(object):
	def Mflux(self,Data,**kwargs):
		G = 6.67e-8
		Msun = 2.0e33
		AU = 1.5e13
		year = 365.0*3600.0*24.0
		
		rho = Data.rho
		v2 = Data.v2
		v1 = Data.v1
		x1 = Data.x1
		x2 = Data.x2

		rin = kwargs.get('rin',50.0) 
		zin = kwargs.get('zin',150.0)
		ul = kwargs.get('ul',1.0)
		urho = kwargs.get('urho',5.0e-14)
		Mstar = kwargs.get('Mstar',30.0)
		uvel = np.sqrt((G*Mstar*Msun)/(ul*AU))

		rneed = np.abs(x1-rin).argmin()
		
		if kwargs.get('normalize',False)==True:
			zneed2 = np.abs(x2-(zin-0.5*rin)).argmin()
		else:
			zneed2 = np.abs(x2-0.001).argmin()
		
		zneed = np.abs(x2-zin).argmin()

		if kwargs.get('scale',False)==True:
			
			phx1 = x1*ul*AU
			phx2 = x2*ul*AU
			phv1 = v1*uvel
			phv2 = v2*uvel
			phrho= rho*urho
			
			Mfr = 2.0*np.pi*phx1[rneed]*integrate.trapz(phrho[rneed,zneed2:zneed]*phv1[rneed,zneed2:zneed],phx2[zneed2:zneed])
			Mfz = 2.0*np.pi*integrate.trapz(phrho[0:rneed,zneed]*phv2[0:rneed,zneed]*phx1[0:rneed],phx1[0:rneed])

			Mfr = (Mfr/Msun)*year
			Mfz = (Mfz/Msun)*year
			
			print "Radial Mass flux is", Mfr
			print "Vertical Mass Flux is",Mfz
			print "Vertical/Radial : ", Mfz/Mfr
		else:
			
			Mfr = 2.0*np.pi*x1[rneed]*integrate.trapz(rho[rneed,zneed2:zneed]*v1[rneed,zneed2:zneed],x2[zneed2:zneed])
			Mfz = 2.0*np.pi*integrate.trapz(rho[0:rneed,zneed]*v2[0:rneed,zneed]*x1[0:rneed],x1[0:rneed])

			print "Radial Mass flux is", Mfr
			print "Vertical Mass Flux is",Mfz
			print "Vertical/Radial :", Mfz/Mfr

		return {'Mfr':Mfr,'Mfz':Mfz}

	def Crit_points(self,Data,xfl=None):
		Im = pp.Image()
		if xfl is None: xfl = 5.0
		Vpol = np.sqrt(Data.v1**2 + Data.v2**2)
		Bpol = np.sqrt(Data.b1**2 + Data.b2**2)
		Btot = np.sqrt(Data.b1**2 + Data.b2**2 + Data.b3**2)
		Valfven = Bpol/np.sqrt(Data.rho)
		Vfast = Btot/np.sqrt(Data.rho)
		Varat = Vpol/Valfven
		Vfrat = Vpol/Vfast
		Dummy = np.zeros(shape=Valfven.shape)
		fldict = Im.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2,xfl,0.001)
		Qx = fldict['qx']
		Qy = fldict['qy']
		Fl_Varat = np.zeros(shape=(2,len(Qx)))
		Fl_Vfrat = np.zeros(shape=(2,len(Qx)))
		
		for i in range(len(Qx)):
			Fl_Varat[:,i] = Im.field_interp(Varat,Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])
			Fl_Vfrat[:,i] = Im.field_interp(Vfrat,Dummy,Data.x1,Data.x2,Data.dx1,Data.dx2,Qx[i],Qy[i])

		Val1 = [Qx[np.abs(Fl_Varat[0,:]-1.0).argmin()],Qy[np.abs(Fl_Varat[0,:]-1.0).argmin()]]
		Val2 = [Qx[np.abs(Fl_Vfrat[0,:]-1.0).argmin()],Qy[np.abs(Fl_Vfrat[0,:]-1.0).argmin()]]

	
		print '------------------------------------------------'
		print '[Qx,Qy] at Alfven :', Val1
		print '[Qx,Qy] at Fast   :', Val2
		print 'Opening Angle at Alfven :',90.0-(180.0/np.pi)*np.arctan(Val1[1]/(Val1[0]-xfl))
		print 'Opening Angle at Fast   :',90.0-(180.0/np.pi)*np.arctan(Val2[1]/(Val2[0]-xfl))
		print '------------------------------------------------'

		return [Val1,Val2,90.0-(180.0/np.pi)*np.arctan(Val1[1]/(Val1[0]-xfl)),90.0-(180.0/np.pi)*np.arctan(Val2[1]/(Val2[0]-xfl))]

		
	def Current(self,Data):
		D = Data
		[r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T
		Cur = r2d*D.b3
		return Cur

	def Magspeed(self,Data):
		magspeeddict={}
		Vpol = np.sqrt(Data.v1**2 + Data.v2**2)
		Bpol = np.sqrt(Data.b1**2 + Data.b2**2)
		Btot = np.sqrt(Data.b1**2 + Data.b2**2 + Data.b3**2)
		Alfv = Bpol/np.sqrt(Data.rho)
		Alftot = Btot/np.sqrt(Data.rho)
		Gam = 5.0/3.0
		cs = np.sqrt((Gam*Data.pr)/Data.rho)
		Magslow = np.sqrt(0.5*(Alftot**2 + cs**2) - 0.5*np.sqrt((Alftot**2 + cs**2)**2 - 4.0*(Alfv**2)*(cs**2)))
		magfast = np.sqrt(0.5*(Alftot**2 + cs**2) + 0.5*np.sqrt((Alftot**2 + cs**2)**2 - 4.0*(Alfv**2)*(cs**2)))
		magspeeddict={'fast':magfast,'slow':Magslow,'alfven':Alfv}
		return magspeeddict
		
		

	def Force_Multi(self,Data,**kwargs):
		D = Data
		[r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T
		T = pp.Tools()
		gradvr = T.Grad(D.v1,D.x1,D.x2,D.dx1,D.dx2)
		gradvz = T.Grad(D.v2,D.x1,D.x2,D.dx1,D.dx2)
		dvrdr = gradvr[:,:,0]
		dvrdz = gradvr[:,:,1]
		dvzdr = gradvz[:,:,0]
		dvzdz = gradvz[:,:,1]
		xrat = z2d/r2d
		dvdl = (1.0/(1.0 + xrat**2.0))*(np.abs(dvrdr) + xrat*(np.abs(dvrdz)+np.abs(dvzdr)) + (xrat**2.0)*(np.abs(dvzdz)))

		G = 6.67e-8
		Msun = 2.0e33
		AU = 1.5e13
		year = 365.0*3600.0*24.0
		sigmae = 0.4
		clight = 3.0e10

		Qo = 1400.0
		alp = kwargs.get('alpha',0.55)
		ul = kwargs.get('ul',1.0)
		urho = kwargs.get('urho',5.0e-14)
		Mstar = kwargs.get('Mstar',30.0)
		uvel = np.sqrt((G*Mstar*Msun)/(ul*AU))
		
		Kpara = (Qo**(1.0-alp))/(1.0-alp)
		Dless = uvel/(sigmae*clight*urho*ul*AU)
		codeval = (1.0/D.rho)*np.abs(dvdl)
		Mt = Kpara*(codeval*Dless)**(alp)

		return Mt

	def dvdl(self,Data,**kwargs):
		D = Data
		[r2d, z2d] = np.meshgrid(Data.x1,Data.x2)
		r2d=r2d.T
		z2d=z2d.T
		T = pp.Tools()
		gradvr = T.Grad(D.v1,D.x1,D.x2,D.dx1,D.dx2)
		gradvz = T.Grad(D.v2,D.x1,D.x2,D.dx1,D.dx2)
		dvrdr = gradvr[:,:,0]
		dvrdz = gradvr[:,:,1]
		dvzdr = gradvz[:,:,0]
		dvzdz = gradvz[:,:,1]
		xrat = z2d/r2d
		if kwargs.get('Disk',False) == True:
			dvdl = np.abs(dvzdz)
		else:
			dvdl = (1.0/(1.0 + xrat**2.0))*(np.abs(dvrdr) + xrat*(np.abs(dvrdz)+np.abs(dvzdr)) + (xrat**2.0)*(np.abs(dvzdz)))

		return dvdl

	def Lsob(self,Data,**kwargs):
		D = Data
		if kwargs.get('Disk',False) == True:
			DvDl = self.dvdl(D,Disk=True)
		else:
			DvDl = self.dvdl(D,Star=True)

		Gamma = 5.0/3.0
		Vpol = np.sqrt(D.v1**2 + D.v2**2)
		csound = np.sqrt(Gamma*(D.pr/D.rho))
		Lsob2D = csound/DvDl
		
		VpolAvg = np.zeros(shape=D.v1.shape)
		LsobAvg = np.zeros(shape=D.v1.shape)
		Area = np.zeros(shape=D.v1.shape)
		for j in range(D.n2):
			for i in range(D.n1):
				LsobAvg[i,j] = Lsob2D[i,j]*D.dx1[i]*D.dx2[j]
				Area[i,j] = D.dx1[i]*D.dx2[j]
				VpolAvg[i,j] = Vpol[i,j]*D.dx1[i]*D.dx2[j]
		
		LsAvg = np.sum(LsobAvg)
		ArAvg = np.sum(Area)
		VpAvg = np.sum(VpolAvg)

		Pol_Vel_Avg = VpAvg/ArAvg
		Sob_Length_Avg = LsAvg/ArAvg
		Growth_Rate = Sob_Length_Avg/Pol_Vel_Avg
		LineInst_Dict = {'VpAvg':Pol_Vel_Avg,'Omega': Growth_Rate,'Lsob':Sob_Length_Avg}
		
		return LineInst_Dict
						
	
		

	
		
		
	
