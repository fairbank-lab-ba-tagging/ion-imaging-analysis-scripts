# Script to read and process labview data
date=20180717
import utility as ut	
#------------------------------------------------------------------------------------------
def lv_ave(runtable,**kwargs): #obsolete
	
	savefile=kwargs.get('savefile',False)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import asciitable as asc
	import gen_reader as gr

	datafile=open(runtable,'r',)
	
	runnum,lvnum=[],[]
	for line in datafile:
		line=line.strip()
		column=line.split(',')
		runnum.append(float(column[0]))
		lvnum.append(column[1])
	datafile.close()
	runnum=np.array(runnum)
	
	# Do the average on the lv files
	avpow=[]
	for i in lvnum:
		lv_matrix=gr.reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+str(i)+'_laser_temp.txt',header=False,delimeter=';')
		avpow.append(float(np.mean(lv_matrix[:,1])))
	avpow=np.array(avpow)
	runtable=np.vstack((runnum,avpow))
	if savefile!=False:
		plt.figure('lv_ave')
		plt.clf()
		plt.plot(runnum,1000*avpow,'bo')
		plt.xlabel('Run')
		plt.ylabel('Laser Power (mW)')
		plt.savefig('/home/james/anaconda3/plots/'+str(date)+'/'+savefile+'.png')
	return runtable
	
#------------------------------------------------------------------------------------------

def get_avpow(lvnum,**kwargs): #Obsolete, use lv_avpow instead
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import numpy as np
	import pandas as pd
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=pd.read_csv('/home/james/anaconda3/data/'+str(date)+'/lv/'+str(lvnum)+'_0001_laser_temp.txt',sep='\s+')
	avpow=np.mean(np.array(lv_matrix)[:,-1])*T_R*flat_T #1.57 is power meter compensation
	
	return avpow
	
#------------------------------------------------------------------------------------------
#For single frame aquisitions 
def lv_energy(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	
	return lv_matrix[:,1]*T_R*flat_T
#------------------------------------------------------------------------------------------

def get_avtime(lvnum):
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	avpow=np.mean(lv_matrix[:,-1])
	
	return avpow
#------------------------------------------------------------------------------------------

def lv_time(lvnum):
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	
	if len(str(lvnum))==1:
		lvnum='000'+str(lvnum)
	if len(str(lvnum))==2:
		lvnum='00'+str(lvnum)
	if len(str(lvnum))==3:
		lvnum='0'+str(lvnum)
	lv_matrix=gr.reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+str(lvnum)+'_laser_temp.txt',header=False,delimeter=';')
	avpow=lv_matrix[:,-1]
	
	return avpow
	
#------------------------------------------------------------------------------------------
#Use this one for calculating energy (W*s)
#For scans
def lv_energy_cont(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	powsum=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		powsum[i]=np.sum(lv_matrix[:,2])*T_R*flat_T*ut.avediff(lv_matrix[:,0])
	
	return powsum

#------------------------------------------------------------------------------------------

def lv_time_cont(lvnum):
	
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	ontime=np.zeros((f,))
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		for j in lv_matrix[:,2]:
			if j >=1e-7:
				ontime[i]+=1
		ontime[i]=ontime[i]*ut.avediff(lv_matrix[:,0])
	
	return ontime
#------------------------------------------------------------------------------------------
def lv_temp_cont(lvnum,**kwargs):
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	temp=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		temp[i]=np.mean(lv_matrix[:,1])
	
	return temp	

#------------------------------------------------------------------------------------------
#for scans
def lv_avpow_cont(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	powave=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		power=[]
		for j in lv_matrix[:,1]:
			if j*T_R*op_flat >= 1E-6:
				power.append(j*T_R*op_flat)
		powave[i]=np.mean(power)
	
	return powave
	
#------------------------------------------------------------------------------------------

def lv_avpow(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	avpow=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		avpow[i]=np.mean(lv_matrix[:,-1])*T_R*flat_T
	
	return avpow
#------------------------------------------------------------------------------------------

def trig_check(trignum,**kwargs):
	
	ret_full=kwargs.get('return_full',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import gen_reader as gr
	
	trigdat=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/trig_'+str(trignum)+'.txt')
	
	# Calculate dt and find the accelerometer level during trigger window
	dt,acc_sig=[],[]
	for i in range(np.shape(trigdat)[0]):
		if i>0:
			dt.append(trigdat[i,0]-trigdat[i-1,0])
		if trigdat[i,2]==1:
			acc_sig_1.append(trigdat[i,1])
	
	acc_sig_ave=np.mean(acc_sig)
	dt_ave=np.mean(dt)
	
	plt.figure('trig')
	plt.clf()
	plt.plot(acc_sig,'b')
	plt.title(str(date)+' run'+str(trignum)+' Trigger Signal while Laser Shutter Open')
	
	if ret_full==True:
		return dt,acc_sig
	return dt_ave,acc_sig_ave
	
#------------------------------------------------------------------------------------------

def trig_check2(trignum,**kwargs):
	
	savepath=kwargs.get('savepath',False)
	
	import numpy as np
	import matplotlib.pyplot as plt
	import gen_reader as gr
	import utility as ut
	
	trigdat=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/trig_'+str(trignum)+'.txt')
	
	# Calculate dt and find the accelerometer level during trigger window
	acc_sig,acc_sig_1=[],[]
	j=0
	for i in range(np.shape(trigdat)[0]):
		if trigdat[i,2]==0:
			j=0
		if trigdat[i,2]==1:
			if j==0:
				acc_sig.append(np.mean(acc_sig_1))
			acc_sig_1.append(trigdat[i,1])
			j+=1
	
	plt.figure('trig2')
	plt.clf()
	plt.plot(acc_sig[1:],'bo',markersize=3)
	plt.title(str(date)+' run'+str(trignum)+' Average Trigger Signal while Laser Shutter Open')
	plt.xlabel('Trigger Number')
	plt.ylabel('Accelerometer Signal (V)')
	if savepath!=False:
		ut.create_dir(savepath)
		plt.savefig(savepath+'trig_'+str(trignum )+'_acc_sig.png')
	
	return acc_sig[1:]	

#------------------------------------------------------------------------------------------
def lv_power_read(lvnum,**kwargs):
	
	T_R=kwargs.get('T_R',1.)
	op_flat=kwargs.get('op_flat',True)
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob
	
	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	if op_flat==True:
		flat_T=.89
	if op_flat==False:
		flat_T=1.
	
	pows=[]
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		pts=len(lv_matrix[:,1])
		pows.append(lv_matrix[:,1]*T_R*flat_T)

		
	return pows

#-------------------------------------------------------------------------------------------
def lv_beamprofile(dist,**kwargs):
	
	xory=kwargs.get('xory','x')
	amp=kwargs.get('amp',.005)
	offset=kwargs.get('offset',0.15)
	width=kwargs.get('width',.1)
	guess=kwargs.get('guess',False)
	plotgauss=kwargs.get('gauss',False)
	stepsize=kwargs.get('stepsize',.005)
	
	import pandas as pd
	import matplotlib.pyplot as plt
	import numpy as np
	from scipy.optimize import curve_fit
	from scipy.special import erf
	from funk import gauss
	
	#z, offset, amplitude, slope
	def erfunc(x,amp,offset,width):
		 return amp/2*erf((x-offset)/(width*np.sqrt(2)))+amp/2
	
	#Import the data, conver pandas to data arrays for position of razor and power reading
	lv_matrix=pd.read_csv('/home/james/anaconda3/data/'+str(date)+'/lv/'+xory+'_profile_in_'+str(dist)+'cm.txt',sep=';')
	posdat=np.array(lv_matrix)[:,0]
	pwrdat=np.flip(np.array(lv_matrix)[:,1],0)
	
	#Do the fit, return the parameters (Don't need to use the covariances)
	params,covs=curve_fit(erfunc,posdat,pwrdat) 

	#Take derivative of data to make gaussian
	derivdat=np.diff(pwrdat,n=1)
	
	#Plot data and fitted function
	plt.figure()
	plt.plot(posdat,pwrdat,'r') #Plots original data
	plt.plot(posdat,erfunc(posdat,*params),'k') #Plots fit
	if guess==True:
		plt.plot(posdat,erfunc(posdat,amp,offset,width),'b') #Plots guesses
	if plotgauss==True:
		plt.plot(posdat,gauss(posdat,params[0]*stepsize/(params[2]*np.sqrt(2*np.pi)),params[1],params[2],'sigma'),'g') #Plot gaussian using parameters, needs to be normalized
		plt.plot(posdat[1::],derivdat,'b') #Plots derivative of data to compare to gaussian
	plt.show()
	
	fwhm=2*np.sqrt(2*np.log(2))*params[2]
	
	print('The '+xory+' FWHM is: '+str(25.4*fwhm)+' mm')
	print('The '+xory+' beam w is '+str(25.4*fwhm/(np.sqrt(2*np.log(2))))+ ' mm')
	return(params) #Amplitude, Offset, Width (sigma)

#-------------------------------------------------------------------------------------------

#Program to read out stage positions from scan lv files
#File should have xpos, ypos as the last two columns
def stagepos(lvnum, **kwargs):
	
	import matplotlib.pyplot as plt
	import numpy as np
	import gen_reader as gr
	import glob

	lvnum=str(lvnum).zfill(4)
	f=len(glob.glob('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_*_laser_temp.txt'))
	
	xpos=np.zeros(f)
	ypos=np.zeros(f)
	for i in range(f):
		frnum=str(i+1).zfill(4)
		lv_matrix=gr.wf_reader('/home/james/anaconda3/data/'+str(date)+'/lv/'+lvnum+'_'+frnum+'_laser_temp.txt',header=False)
		xpos[i]=np.mean(lv_matrix[:,-2])
		ypos[i]=np.mean(lv_matrix[:,-1])

	return xpos,ypos
