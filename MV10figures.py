from numpy import *
from pylab import *
from mpl_toolkits.basemap import NetCDFFile as ncopen
import matplotlib.numerix.ma as M
import cPickle as pickle
from contour_ranges import *
from matplotlib import cm
from matplotlib.font_manager import FontProperties
from scipy.interpolate import UnivariateSpline as spline
#from scipy.interpolate import SmoothBivariateSpline as bispline
#from scipy.interpolate.interpolate import spline
from scipy.interpolate.interpolate import interp1d
from matplotlib import colors

labels = ['(a)','(b)','(c)','(d)','(e)','(f)']
dp   = -5.e3         # defined by differencing increasing index
           # dlat also defined by differencing increasing index
om = 2.*pi/86400.
g = 9.8
def r(Ro):
    return sqrt(287.*60./Ro)/2./om

Ro = ['0.02','1.3','10.5']
#Ro = ['0.02']
#Ro = ['10.5']
#Ro = ['1.3']
#dir = '/Users/mitch/Data/earth2titan/'
dir = '/Volumes/Data/'
d = {}

for R in Ro:
	d[R] = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_T42/atmos_daily.nc','r')        
print d.keys()
lon = d[R].variables['lon'][:]
ni = len(lon)
lat = d[R].variables['lat'][:]
dlat = (lat[0]-lat[1])*pi/180.
nj = len(lat)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
f = 2.*om*sinphi
lev = d[R].variables['pfull'][:]/1000.
exner  = lev**(2./7.)
nk = len(lev)

## FIGURE 1
'''
# Annual average temperatures
delh = 0.2
crange=arange(0.7,1.,delh/4)
To = 285.
Tbar = To*(1.+delh/3.*(1.-3.*sinphi**2.))
kap = 2./7.
gam = 4.e-3
R = 287.
blah = exp(R*gam/g/(2-2*kap)*(1.-lev**(2-2*kap)))
Tstar = zeros((nk,nj),'d')
for j in range(nj):
    Tstar[:,j] = Tbar[j]*blah*exner
Tstar = where(Tstar<200.,200.,Tstar)

li=7
lf=10
panels = '13'
f = figure(figsize=(6,2))
m=0
delh = 0.21053
To = 285.
Tbar = To*(1.+delh/3.*(1.-3.*sinphi**2.))
kap = 2./7.
gam = 4.e-3
R = 287.
blah = exp(R*gam/g/(2-2*kap)*(1.-lev**(2-2*kap)))
Tstar = zeros((nk,nj),'d')
for j in range(nj):
    Tstar[:,j] = Tbar[j]*blah*exner
Tstar = where(Tstar<200.,200.,Tstar)
crange=arange(0.7,1.,delh/4)
for R in Ro:
	m+=1
	print R
	d = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_T42/atmos_monthly.nc','r')        
	teq = array(d.variables['teq'][0,:,:,0])
	tave = array(average(average(d.variables['temp'][li:lf],axis=0),axis=-1))
	subplot(panels+str(m))
	contour(lat,lev,teq/To,crange,colors='k',linestyles='dotted',linewidths=0.5)
	contour(lat,lev,tave/To,crange,colors='0.5')
	xlim(lat[0],lat[-1])
	xticks(arange(-90,100,30),fontsize=8)
	xlabel('latitude [deg]',fontsize=9)
	if m==1: ylabel(r'$p/p_s$',fontsize=10)
	ylim(lev[-1],lev[0])
	yticks(fontsize=8)
	title(r'$Ro_T$='+R,fontsize=11)
subplots_adjust(bottom=0.21,top=0.85)

show()
'''

## Figure 2
'''
# Streamfunction, winds
d = {}
sym = {}
for R in Ro:
    print R
    d[R] = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_T42/atmos_monthly.nc','r')
    sym[R] = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_symmetric_T42/atmos_monthly.nc','r')    
lon = d[R].variables['lon'][:]
ni = len(lon)
lat = d[R].variables['lat'][:]
dlat = (lat[0]-lat[1])*pi/180.
nj = len(lat)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
f = 2.*om*sinphi
lev = d[R].variables['pfull'][:]/1000.
exner  = lev**(2./7.)
nk = len(lev)
kb = argmin(abs(lev-0.4))
lat2 = zeros((nk,nj),'d')
f2   = lat2.copy()
lev2 = lat2.copy()
for j in range(nj):
    for k in range(nk):
        f2[k,j]   = f[j]
        lat2[k,j] = lat[j]
        lev2[k,j] = lev[k]

lis = {}
lis['0.02'] = [25]
lis['1.3'] = [25]
lis['10.5'] = [25]
lfs = {}
lfs['0.02'] = [36]
lfs['1.3'] = [36]
lfs['10.5'] = [36]

panels = '23'
f = figure(figsize=(6,3))
m=0
for i in range(len(lis['0.02'])):
    for R in Ro:
        m+=1
        li = lis[R][i]
        lf = lfs[R][i]
        a = r(float(R))
        y = a*sinphi
        v = average(average(d[R].variables['vcomp'][li:lf,:,:,:],axis=-1),axis=0)
        psi1 = zeros((nk+1,nj),'d')
        for k in range(len(lev)):
            psi1[k+1] = psi1[k]+v[k]*dp/g*2.*pi*a*cosphi
        psi = 0.5*(psi1[0:nk]+psi1[1:nk+1])
        u = average(average(d[R].variables['ucomp'][li:lf,:,:,:],axis=-1),axis=0)
        uave = u
        quad = panels+str(m)
        print quad
        ax = subplot(quad)
        cmap = cm.colors.ListedColormap(['r','b'])
        norm = cm.colors.BoundaryNorm([-1.e6,1.,1.e6],cmap.N)
        contour(lat,lev,u,U_range[R],colors='b',linewidths=0.75)
        contour(lat,lev,u,[0.],linewidths=1.5,colors='b')
        contourf(lat,lev,psi*1.e-9,Had_range[R],cmap=cm.binary,extend='both')
	max = M.max(psi)*1.e-9
	maxind = psi.argmax()
	kh,jh = unravel_index(maxind,shape(psi))
	if R=='10.5': blah = '%1.1f' %max
	else: blah = '%i' %max
	text(lat[jh],lev[kh],blah,color='w',fontsize=8)
        if mod(m,3)==1:
            ylabel(r'$p/p_s$')
            yticks(arange(0.2,0.9,.2),fontsize=8)
        else:
            yticks(arange(0.2,0.9,.2))
            yticklabels = ax.get_yticklabels()
            setp(yticklabels, visible=False)
        if m>=4:
            xlabel('latitude',fontsize=9)
            xticks(arange(-80,90,40),fontsize=8)
        else:
            xticks(arange(-80,90,40),fontsize=8)
            xticklabels = ax.get_xticklabels()
            setp(xticklabels, visible=False)
        if m<4: title(r'$Ro_T$ = '+R,fontsize=12)
        xlim(lat[0],lat[-1])
        ylim(lev[-1],lev[0])
        subplots_adjust(wspace=0.02,hspace=0.02,bottom=0.15)

for i in range(len(lis['0.02'])):
    for R in Ro:
        m+=1
        li = lis[R][i]
        lf = lfs[R][i]
        a = r(float(R))
        y = a*sinphi
        v = average(average(sym[R].variables['vcomp'][li:lf,:,:,:],axis=-1),axis=0)
        psi1 = zeros((nk+1,nj),'d')
        for k in range(len(lev)):
            psi1[k+1] = psi1[k]+v[k]*dp/g*2.*pi*a*cosphi
        psi1 = zeros((nk+1,nj),'d')
        for k in range(len(lev)):
            psi1[k+1] = psi1[k]+v[k]*dp/g*2.*pi*a*cosphi
        psi = 0.5*(psi1[0:nk]+psi1[1:nk+1])
        u = average(average(sym[R].variables['ucomp'][li:lf,:,:,:],axis=-1),axis=0)
        uave = u
        quad = panels+str(m)
        print quad
        ax = subplot(quad)
        cmap = cm.colors.ListedColormap(['r','b'])
        norm = cm.colors.BoundaryNorm([-1.e6,1.,1.e6],cmap.N)
        contour(lat,lev,u,U_range[R],colors='b',linewidths=0.75)
        contour(lat,lev,u,[0.],linewidths=1.5,colors='b')
        contourf(lat,lev,psi*1.e-9,Had_range[R],cmap=cm.binary,extend='both')
	max = M.max(psi)*1.e-9
	maxind = psi.argmax()
	kh,jh = unravel_index(maxind,shape(psi))
	if R=='10.5': blah = '%1.1f' %max
	else: blah = '%i' %max
	text(lat[jh],lev[kh],blah,color='w',fontsize=8)
        if mod(m,3)==1:
            ylabel(r'$p/p_s$')
            yticks(arange(0.2,0.9,.2),fontsize=8)
        else:
            yticks(arange(0.2,0.9,.2))
            yticklabels = ax.get_yticklabels()
            setp(yticklabels, visible=False)
        if m>=4:
            xlabel('latitude',fontsize=9)
            xticks(arange(-80,90,40),fontsize=8)
        else:
            xticks(arange(-80,90,40),fontsize=8)
            xticklabels = ax.get_xticklabels()
            setp(xticklabels, visible=False)
        if m<4: title(r'$Ro$ = '+R,fontsize=12)
        xlim(lat[0],lat[-1])
        ylim(lev[-1],lev[0])
        subplots_adjust(wspace=0.02,hspace=0.02,bottom=0.15)
        gcf().text(0.925,0.65,'3-dim', \
                   horizontalalignment='center', \
                   rotation=-90,\
                   fontproperties=FontProperties(size=12))
        gcf().text(0.925,0.29,'2-dim', \
                   horizontalalignment='center', \
                   rotation=-90,\
                   fontproperties=FontProperties(size=12))
show()
'''


## Figure 3
'''
# Local Rossby number
# -zeta/f is problematic
# instead use Ubar/(f*a*cosphi)

f2d = zeros((nk,nj),'d')
cos2d = zeros((nk,nj),'d')
for k in range(nk):
	f2d[k] = 2.*om*sinphi
	cos2d[k] = cosphi

figure(figsize=(6,2))
crange = arange(0.,3.01,0.5)
m = 0
quad = '13'
fmts = {}
fmts['0.02'] = '%1.2f'
fmts['1.3'] = '%i'
fmts['10.5'] = '%i'
for R in Ro:
	m+=1
	subplot(quad+str(m))
	a = r(float(R))
	#u = average(d[R].variables['ucomp'][-1],axis=-1)
	u = average(average(d[R].variables['ucomp'][25:],axis=-1),axis=0)

	Ros = u/abs(f2d)/a/cos2d
	Ros = M.array(Ros)
	Ros = M.masked_where(abs(f2d)<2.*om/10.,Ros)

	#contourf(lat,lev,u,U_range[R],cmap=cm.binary,extend='both')
	contourf(lat,lev,Ros,crange,cmap=cm.binary,extend='both')
	if m==3: colorbar(fraction=0.05,pad=0.0,shrink=0.85)
	cs = contour(lat,lev,Ros,Ro_range[R],colors='0.2')
	clabel(cs,Ro_labels[R],fmt=fmts[R],fontsize=10,rightside_up=True,inline_spacing=2)
	contour(lat,lev,u,[0.],colors='k')
	print M.maximum(Ro_range[R])
	title(r'$Ro_T = %s$' %(R),fontsize=11)
	yticks(fontsize=9)
	ylim(lev[-1],lev[0])
	xticks(arange(-90,100,45),fontsize=9)
	xlim(-90,90)
	if m==1: ylabel(r'$p/p_s$',fontsize=12)
	xlabel('latitude [deg]',fontsize=10)
subplots_adjust(bottom=0.2,top=0.82,right=0.85)
show()
'''

## FIGURE 4
'''
# Steady-state PV
d = {}
for R in Ro:
    print R
    d[R] = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_T42/atmos_monthly.nc','r')
lat = d[R].variables['lat'][:]
lev = d[R].variables['pfull'][:]/1000.
nj = len(lat)
nk = len(lev)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
exner  = lev**(2./7.)
f = 2.*om*sinphi
lat2 = zeros((nk,nj),'d')
f2   = lat2.copy()
lev2 = lat2.copy()
nth = 100.
thetanew = arange(250.,400.,150./nth)
PVnew = zeros((nth,nj),'d')
for j in range(nj):
    for k in range(nk):
        f2[k,j]   = f[j]
        lat2[k,j] = lat[j]
        lev2[k,j] = lev[k]
sinphi2 = sin(lat2*pi/180.)
cosphi2 = cos(lat2*pi/180.)

kbs = [argmin(abs(lev-0.2)),argmin(abs(lev-0.5))]

li = 25
lf = 36
figure(figsize=(6,3))
l=0
for kb in kbs:
	for R in Ro:
		print R
		l+=1
		Rd = 287.
		a = r(float(R))
		y = a*lat*pi/180.
		uave = average(average(d[R].variables['ucomp'][li:lf,:,:,:],axis=0),axis=-1)
		Tave = average(average(d[R].variables['temp'][li:lf,:,:,:],axis=0),axis=-1)
		thetaave = Tave.copy()
		print d[R].variables['time'][li]
		for j in range(nj):
			thetaave[:,j]=Tave[:,j]/exner
		dthetaavedy1 = zeros((nk,nj+1),'d')
		phi = lat*pi/180.
		for j in range(nj-1):
			dthetaavedy1[:,j+1] = (thetaave[:,j]-thetaave[:,j+1]) \
			/(y[j]-y[j+1])
		dthetaavedy1[:,0] = dthetaavedy1[:,1]
		dthetaavedy1[:,-1] = dthetaavedy1[:,nj-1]
		dthetaavedy = 0.5*(dthetaavedy1[:,0:nj]+dthetaavedy1[:,1:nj+1])
		dTavedp1 = zeros((nk+1,nj),'d')
		dthetaavedp1 = zeros((nk+1,nj),'d')
		duavedp1 = zeros((nk+1,nj),'d')
		for k in range(nk-1):
			dTavedp1[k+1,:] = (Tave[k,:]-Tave[k+1,:])/dp
			dthetaavedp1[k+1,:] = (thetaave[k,:]-thetaave[k+1,:])/dp
			duavedp1[k+1,:] = (uave[k,:]-uave[k+1,:])/dp
		dTavedp1[0] = dTavedp1[1]
		dTavedp1[-1] = dTavedp1[nk-1]
		dTavedp = 0.5*(dTavedp1[0:nk]+dTavedp1[1:nk+1])
		dthetaavedp1[0] = dthetaavedp1[1]
		dthetaavedp1[-1] = dthetaavedp1[nk-1]
		dthetaavedp = 0.5*(dthetaavedp1[0:nk]+dthetaavedp1[1:nk+1])
		duavedp = 0.5*(duavedp1[0:nk]+duavedp1[1:nk+1])
		vor = average(average(d[R].variables['vor'][li:lf,:,:,:],axis=-1),axis=0)
		PV = -g*( (f2+vor)*dthetaavedp + duavedp*dthetaavedy )
		quad = '23'+str(l)
		ax = subplot(quad)
		if l<4: 
			plot(lat,PV[kb]*1.e7,'k:')
			ysc = 7
		else:   
			plot(lat,PV[kb]*1.e8,'k:')
			ysc = 8
		plot(lat,uave[kb],'k')
		plot(lat,uave[kb]*0.,linewidth=0.5,color=(0.5,0.5,0.5))
		ymax = 60
		if mod(l,3)==1:
			ylabel(r'$\bar{q}\times10^%i, \ \bar{u}$' %(ysc),fontsize=12)
			yticks(arange(-40,50,20),fontsize=8)
		else:
			yticks(arange(-40,50,20),fontsize=8)
			yticklabels = ax.get_yticklabels()
			setp(yticklabels, visible=False)
		if l>=4:
			xlabel('latitude')
			xticks(arange(-80,90,40),fontsize=8)
		else:
			xticks(arange(-80,90,40),fontsize=8)
			xticklabels = ax.get_xticklabels()
			setp(xticklabels, visible=False)
		if l<4: title(r'$Ro_T$ = '+R,fontsize=12)
		xlim(lat[0],lat[-1])
		subplots_adjust(wspace=0.,hspace=0.,bottom=0.15)
		ylim(-ymax,ymax)
		gcf().text(0.925,0.64,'%i hPa' %(200), \
			horizontalalignment='center', \
			rotation=-90,\
			fontproperties=FontProperties(size=12))
		gcf().text(0.925,0.25,'%i hPa' %(500), \
			horizontalalignment='center', \
			rotation=-90,\
			fontproperties=FontProperties(size=12))
show()
'''

## Figure 5
'''
# Horizontal geopotential structure
Ros = ['0.02','1.3','10.5']
lev = 400.
datdir = '/Volumes/Data/'
g = 9.8
fprops = {}
fprops['fontsize']=8
crange = {'0.02':arange(0.9,0.99,.01),'1.3':arange(0.9,0.99,0.01),'10.5':arange(0.9,0.99,0.01)}
ticks = {'0.02':arange(0.9,0.99,.04),'1.3':arange(0.9,0.99,0.04),'10.5':arange(0.9,0.99,0.04)}

def file(Ro):
	return 'Ro%s_noseas_mlapse_nuv0.01_T42/atmos_daily.nc' %Ro
	
def geopot(T,p):
	gz1 = zeros((nk+1,nj,ni),'d')
	for k in range(nk-1,-1,-1):
		gz1[k+1] = gz1[k]+287.*T[k]*dp/p[k]
	gz1[0] = gz1[1]
	return 0.5*(gz1[0:nk]+gz1[1:nk+1])

f = figure(figsize=(6,2.2))
i=0
for Ro in Ros:
	i+=1
	d = ncopen(datdir+file(Ro),'r')
	T = d.variables['temp'][-1]
	p = d.variables['pfull'][:]
	khere = argmin(abs(p-lev))
	lat = d.variables['lat'][:]
	lon = d.variables['lon'][:]
	dp = p[1]-p[0]
	nk,nj,ni = shape(T)
	gz = geopot(T,p)
	#plot(gz[:,nj/2,ni/2],-p)
	#show()
	subplot('13%s' %i)
	con = contourf(lon,lat,gz[khere]/g,1.e3*crange[Ro],extend='both',cmap=cm.binary)
	if i==1: ylabel('latitude',fontsize=9)
	xlabel('longitude',fontsize=9)
	title(r'$Ro_{T}=%s$' %Ro,fontsize=11)
	yticks(arange(-80,90,20),fontsize=8)
	xticks(arange(0,365,60),fontsize=8)
	ylim(lat[0],lat[-1])
	xlim(lon[0],lon[-1])
	#clabel(con,ticks[Ro],**fprops)
	#cb = colorbar(orientation='horizontal',format='%.2f',ticks=ticks[Ro],pad=0.2)
	#axes(cb.ax)
	#xticks(fontsize=8)
colorbar(con,fraction=0.05)
subplots_adjust(bottom=0.2,top=0.85,right=0.85)
show()
'''

## Figure 6
'''
# Sponge runs, Streamfunction, winds
tags = ['T42','hilat_sym','lolat_sym']
titles = {}
titles['T42'] = 'Control'
titles['hilat_sym'] = 'Hi-lat sponge'
titles['lolat_sym'] = 'Lo-lat sponge'
for tag in tags:
    d[tag] = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_'+tag+'/atmos_daily.nc','r')
lon = d[tag].variables['lon'][:]
ni = len(lon)
lat = d[tag].variables['lat'][:]
dlat = (lat[0]-lat[1])*pi/180.
nj = len(lat)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
f = 2.*om*sinphi
lev = d[tag].variables['pfull'][:]/1000.
exner  = lev**(2./7.)
nk = len(lev)
kb = argmin(abs(lev-0.4))
lat2 = zeros((nk,nj),'d')
f2   = lat2.copy()
lev2 = lat2.copy()
for j in range(nj):
    for k in range(nk):
        f2[k,j]   = f[j]
        lat2[k,j] = lat[j]
        lev2[k,j] = lev[k]

lis = {}
lis['T42'] = [721]
lis['hilat_sym'] = [300]
lis['lolat_sym'] = [300]
lfs = {}
lfs['T42'] = [1080]
lfs['hilat_sym'] = [360]
lfs['lolat_sym'] = [360]

panels = '13'
f = figure(figsize=(6,2))
m=0
for i in range(len(lis['T42'])):
    for tag in tags:
        m+=1
        li = lis[tag][i]
        lf = lfs[tag][i]
        a = r(float(10.5))
        y = a*sinphi
        v = average(average(d[tag].variables['vcomp'][li:lf,:,:,:],axis=-1),axis=0)
        psi1 = zeros((nk+1,nj),'d')
        for k in range(len(lev)):
            psi1[k+1] = psi1[k]+v[k]*dp/g*2.*pi*a*cosphi
        psi = 0.5*(psi1[0:nk]+psi1[1:nk+1])
        u = average(average(d[tag].variables['ucomp'][li:lf,:,:,:],axis=-1),axis=0)
        uave = u
        quad = panels+str(m)
        print quad
        ax = subplot(quad)
        cmap = cm.colors.ListedColormap(['r','b'])
        norm = cm.colors.BoundaryNorm([-1.e6,1.,1.e6],cmap.N)
        contour(lat,lev,u,U_range[R],colors='b',linewidths=0.75)
        contour(lat,lev,u,[0.],linewidths=1.5,colors='b')
        contourf(lat,lev,psi*1.e-9,Had_range[R],cmap=cm.binary,extend='both')
	max = M.max(psi)*1.e-9
	maxind = psi.argmax()
	kh,jh = unravel_index(maxind,shape(psi))
	if R=='10.5': blah = '%1.1f' %max
	else: blah = '%i' %max
	text(lat[jh],lev[kh],blah,color='w',fontsize=8)
	title(titles[tag],fontsize=10)
        if mod(m,3)==1:
            ylabel(r'$p/p_s$')
            yticks(arange(0.2,0.9,.2))
        else:
            yticks(arange(0.2,0.9,.2))
            yticklabels = ax.get_yticklabels()
            setp(yticklabels, visible=False)
	xlabel('latitude',fontsize=9)
	xticks(arange(-80,90,40),fontsize=8)
        xlim(lat[0],lat[-1])
        ylim(lev[-1],lev[0])
colorbar(fraction=0.05)
subplots_adjust(wspace=0.02,hspace=0.02,bottom=0.22,top=0.83,right=0.85)
show()
'''

## Figure 7
'''
# EP fluxes
cmap = cm.Blues
kwargs = {}
kwargs['scale'] = 4.5e-5
kwargs['width'] = 0.0055
kwargs['headwidth'] = 5
kwargs['headlength'] = 4
kwargs['color'] = 'k'
figure(figsize=(6,2))
file = open('/Users/mitch/Scripts/Earth2Titan/0.02_Fluxes_days1000to1080.dic','r')
d = pickle.load(file)
file.close()
lat = d['lat']
lev = d['lev']
nj = len(lat)
nk = len(lev)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
phi = lat*pi/180.
exner  = lev**(2./7.)
f = 2.*om*sinphi
lat2 = zeros((nk,nj),'d')
f2   = lat2.copy()
lev2 = lat2.copy()
for j in range(nj):
    for k in range(nk):
        f2[k,j]   = f[j]
        lat2[k,j] = lat[j]
        lev2[k,j] = lev[k]
sinphi2 = sin(lat2*pi/180.)
cosphi2 = cos(lat2*pi/180.)

quad = '131'
R = '0.02'
key = R
a = r(float(R))
y = a*phi
thetaave = zeros((nk,nj),'d')
vpthetapbar = zeros((nk,nj),'d')
for j in range(nj):
    thetaave[:,j] = d['Tave'][:,j]/exner
    vpthetapbar[:,j] = d['vpTpbar'][:,j]/exner
    uave = d['uave']
rhoR2 = zeros((nk,nj),'d')
for j in range(nj):
    rhoR2[:,j] = 1.e5*lev/287./d['Tave'][:,j]
# Vertical gradients
duavedp1 = zeros((nk+1,nj),'d')
dthetaavedp1 = zeros((nk+1,nj),'d')
for k in range(nk-1):
    duavedp1[k+1]     = (uave[k]-uave[k+1])/dp
    dthetaavedp1[k+1] = (thetaave[k]-thetaave[k+1])/dp
# project thermal gradient up and down at boundaries
dthetaavedp1[0] = dthetaavedp1[1]
dthetaavedp1[-1] = dthetaavedp1[nk-1]
duavedp = 0.5*(duavedp1[0:nk]+duavedp1[1:nk+1])
dthetaavedp = 0.5*(dthetaavedp1[0:nk]+dthetaavedp1[1:nk+1])
# Horizontal shear
duavedy1 = zeros((nk,nj+1),'d')
for j in range(nj-1):
    duavedy1[:,j+1] = (uave[:,j]*cosphi[j]-uave[:,j+1]*cosphi[j+1])/dlat
duavedy = 0.5*(duavedy1[:,0:nj]+duavedy1[:,1:nj+1])/a/cosphi2
# EP flux (Vallis, p. 536-537)
Fj = zeros((nk,nj),'d')
Fk = Fj.copy()
for k in range(nk):
    Fj[k,:] = -rhoR2[k]*cosphi*d['upvpbar'][k,:] \
              +rhoR2[k]*cosphi*duavedp[k,:]*vpthetapbar[k,:] \
              /dthetaavedp[k,:]
    Fk[k,:] = -f/g*cosphi*vpthetapbar[k,:]/dthetaavedp[k,:] \
              +  1./g*cosphi*duavedy[k,:]*vpthetapbar[k,:]/dthetaavedp[k,:] #\
    #+  1./g*cosphi*d['upompbar'][k,:]
divFj1 = zeros((nk  ,nj+1),'d')
divFk1 = zeros((nk+1,nj  ),'d')
for j in range(nj-1):
    divFj1[:,j+1] = (Fj[:,j]*cosphi[j]-Fj[:,j+1]*cosphi[j+1]) \
                    /dlat
divFj = 0.5*(divFj1[:,0:nj]+divFj1[:,1:nj+1])/a/cosphi2
for k in range(nk-1):
    divFk1[k+1,:] = -rhoR2[k]*g*(Fk[k,:]-Fk[k+1,:])/dp
divFk = 0.5*(divFk1[0:nk,:]+divFk1[1:nk+1,:])
totdiv = (divFk+divFj)/rhoR2/cosphi2    
ax = subplot(quad)
print quad
# EP flux divergence
cf = contourf(lat,lev,totdiv*1.e5,divEP_range[R],extend='both',cmap=cmap)
quiver(lat2[:,::4],lev2[:,::4],Fj[:,::4]/a,Fk[:,::4]/1.e5,**kwargs)
ylim(lev[-1],lev[0])
xlim(lat[0],lat[-1])
xticks(arange(-80,90,40),fontsize=8)
yticks(fontsize=8)
xlabel('latitude',fontsize=8)
title(r'$Ro_T$=0.02, steady',fontsize=10)
ylabel(r'$p/p_s$',fontsize=10)

file = open('/Users/mitch/Scripts/Earth2Titan/10.5_Fluxes_days230to305.dic','r')
d = pickle.load(file)
file.close()
quad = '132'
R = '10.5'
key = R
a = r(float(R))
y = a*sinphi
thetaave = zeros((nk,nj),'d')
vpthetapbar = zeros((nk,nj),'d')
for j in range(nj):
    thetaave[:,j] = d['Tave'][:,j]/exner
    vpthetapbar[:,j] = d['vpTpbar'][:,j]/exner
    uave = d['uave']
rhoR2 = zeros((nk,nj),'d')
for j in range(nj):
    rhoR2[:,j] = 1.e5*lev/287./d['Tave'][:,j]
# Vertical gradients
duavedp1 = zeros((nk+1,nj),'d')
dthetaavedp1 = zeros((nk+1,nj),'d')
for k in range(nk-1):
    duavedp1[k+1]     = (uave[k]-uave[k+1])/dp
    dthetaavedp1[k+1] = (thetaave[k]-thetaave[k+1])/dp
# project thermal gradient up and down at boundaries
dthetaavedp1[0] = dthetaavedp1[1]
dthetaavedp1[-1] = dthetaavedp1[nk-1]
duavedp = 0.5*(duavedp1[0:nk]+duavedp1[1:nk+1])
dthetaavedp = 0.5*(dthetaavedp1[0:nk]+dthetaavedp1[1:nk+1])
# Horizontal shear
duavedy1 = zeros((nk,nj+1),'d')
for j in range(nj-1):
    duavedy1[:,j+1] = (uave[:,j]*cosphi[j]-uave[:,j+1]*cosphi[j+1])/dlat
duavedy = 0.5*(duavedy1[:,0:nj]+duavedy1[:,1:nj+1])/a/cosphi2
# EP flux (Vallis, p. 536-537)
Fj = zeros((nk,nj),'d')
Fk = Fj.copy()
for k in range(nk):
    Fj[k,:] = -rhoR2[k]*cosphi*d['upvpbar'][k,:] \
              +rhoR2[k]*cosphi*duavedp[k,:]*vpthetapbar[k,:] \
              /dthetaavedp[k,:]
    Fk[k,:] = -f/g*cosphi*vpthetapbar[k,:]/dthetaavedp[k,:] \
              +  1./g*cosphi*duavedy[k,:]*vpthetapbar[k,:]/dthetaavedp[k,:] #\
    #+  1./g*cosphi*d['upompbar'][k,:]
divFj1 = zeros((nk  ,nj+1),'d')
divFk1 = zeros((nk+1,nj  ),'d')
for j in range(nj-1):
    divFj1[:,j+1] = (Fj[:,j]*cosphi[j]-Fj[:,j+1]*cosphi[j+1]) \
                    /dlat
divFj = 0.5*(divFj1[:,0:nj]+divFj1[:,1:nj+1])/a/cosphi2
for k in range(nk-1):
    divFk1[k+1,:] = -rhoR2[k]*g*(Fk[k,:]-Fk[k+1,:])/dp
divFk = 0.5*(divFk1[0:nk,:]+divFk1[1:nk+1,:])
totdiv = (divFk+divFj)/rhoR2/cosphi2    
ax = subplot(quad)
print quad
# EP flux divergence
cf = contourf(lat,lev,totdiv*1.e5,divEP_range[R],cmap=cmap,extend='both')
quiver(lat2[:,::4],lev2[:,::4],Fj[:,::4]/a,Fk[:,::4]/1.e5,**kwargs)
ylim(lev[-1],lev[0])
xlim(lat[0],lat[-1])
xticks(arange(-80,90,40),fontsize=8)
yticks(fontsize=8)
xlabel('latitude',fontsize=8)
yticklabels = ax.get_yticklabels()
setp(yticklabels, visible=False)
title(r'$Ro_T$=10.5, spinup',fontsize=10)

file = open('/Users/mitch/Scripts/Earth2Titan/10.5_Fluxes_days1000to1080.dic','r')
d = pickle.load(file)
file.close()
quad = '133'
R = '10.5'
key = R
a = r(float(R))
y = a*sinphi
thetaave = zeros((nk,nj),'d')
vpthetapbar = zeros((nk,nj),'d')
for j in range(nj):
    thetaave[:,j] = d['Tave'][:,j]/exner
    vpthetapbar[:,j] = d['vpTpbar'][:,j]/exner
    uave = d['uave']
rhoR2 = zeros((nk,nj),'d')
for j in range(nj):
    rhoR2[:,j] = 1.e5*lev/287./d['Tave'][:,j]
# Vertical gradients
duavedp1 = zeros((nk+1,nj),'d')
dthetaavedp1 = zeros((nk+1,nj),'d')
for k in range(nk-1):
    duavedp1[k+1]     = (uave[k]-uave[k+1])/dp
    dthetaavedp1[k+1] = (thetaave[k]-thetaave[k+1])/dp
# project thermal gradient up and down at boundaries
dthetaavedp1[0] = dthetaavedp1[1]
dthetaavedp1[-1] = dthetaavedp1[nk-1]
duavedp = 0.5*(duavedp1[0:nk]+duavedp1[1:nk+1])
dthetaavedp = 0.5*(dthetaavedp1[0:nk]+dthetaavedp1[1:nk+1])
# Horizontal shear
duavedy1 = zeros((nk,nj+1),'d')
for j in range(nj-1):
    duavedy1[:,j+1] = (uave[:,j]*cosphi[j]-uave[:,j+1]*cosphi[j+1])/dlat
duavedy = 0.5*(duavedy1[:,0:nj]+duavedy1[:,1:nj+1])/a/cosphi2
# EP flux (Vallis, p. 536-537)
Fj = zeros((nk,nj),'d')
Fk = Fj.copy()
for k in range(nk):
    Fj[k,:] = -rhoR2[k]*cosphi*d['upvpbar'][k,:] \
              +rhoR2[k]*cosphi*duavedp[k,:]*vpthetapbar[k,:] \
              /dthetaavedp[k,:]
    Fk[k,:] = -f/g*cosphi*vpthetapbar[k,:]/dthetaavedp[k,:] \
              +  1./g*cosphi*duavedy[k,:]*vpthetapbar[k,:]/dthetaavedp[k,:] #\
    #+  1./g*cosphi*d['upompbar'][k,:]
divFj1 = zeros((nk  ,nj+1),'d')
divFk1 = zeros((nk+1,nj  ),'d')
for j in range(nj-1):
    divFj1[:,j+1] = (Fj[:,j]*cosphi[j]-Fj[:,j+1]*cosphi[j+1]) \
                    /dlat
divFj = 0.5*(divFj1[:,0:nj]+divFj1[:,1:nj+1])/a/cosphi2
for k in range(nk-1):
    divFk1[k+1,:] = -rhoR2[k]*g*(Fk[k,:]-Fk[k+1,:])/dp
divFk = 0.5*(divFk1[0:nk,:]+divFk1[1:nk+1,:])
totdiv = (divFk+divFj)/rhoR2/cosphi2    
ax = subplot(quad)
print quad
# EP flux divergence
cf = contourf(lat,lev,totdiv*1.e5,divEP_range[R],extend='both',cmap=cmap)
quiver(lat2[:,::4],lev2[:,::4],Fj[:,::4]/a,Fk[:,::4]/1.e5,**kwargs)
ylim(lev[-1],lev[0])
xlim(lat[0],lat[-1])
xticks(arange(-80,90,40),fontsize=8)
ax.yaxis.set_label_position('right')
ax.yaxis.set_ticks_position('right')
yticks(fontsize=8)
xlabel('latitude',fontsize=8)
title(r'$Ro_T$=10.5, steady',fontsize=10)
ylabel(r'$p/p_s$',fontsize=10,rotation=270)

subplots_adjust(bottom=0.2,top=0.81)
show()
'''

## Figure 8
# WaveRegimes3.pdf

## Figure 9
'''
# div(UpVp) 
file = dir+'ffts/JGR_fig_upvp_400mbar.dic'
days = ['1to360','250to360','1080to1200']
titles = {}
titles['1to360'] = 'Days 1-360'
titles['250to360'] = 'Days 250-360'
titles['1080to1200'] = 'Days 1080-1200'
d = open(file,'r')
data = pickle.load(d)
a = 280.e3

jh = argmin(abs(lat-0.))

crange=arange(0.05,2.1,0.2)
figure(figsize=(6,2))
panels = '13'
m = 0
for day in days:
	m+=1
	n = data['n'][day]
	freq = data['freq'][day]
	upr_vpr = data['up_vp'][day]
	freq*=86400.
	#var = sqrt(average((upr_vpr[:,:,:]/1.e9)**2.,axis=1,weights=cosphi))
	nfreq = len(freq)
	nlat = len(lat)
	nlon = len(n)
	dupvpdy1 = zeros((nfreq,nlat+1,nlon),'d')
	y = a*sinphi
	for j in range(nlat-1):
		dupvpdy1[:,j+1] = (upr_vpr[:,j+1]*cosphi[j+1]-upr_vpr[:,j]*cosphi[j])/(y[j+1]-y[j])
	dupvpdy = 0.5*(dupvpdy1[:,0:nlat]+dupvpdy1[:,1:nlat+1])
	var = dupvpdy[:,jh,:]
	dfreq = abs(freq[1]-freq[8])
	for l in range(len(freq)):
		if mod(l,10)==0: print l
		var[l,:] = average(var, \
					weights=exp(-((freq-freq[l])/dfreq)**2.), \
					axis=0)
	if m==1: max = M.max(var)
	crange = arange(-max/2.+max/50.,max/2.+max/10.,max/10.)
	blah = M.array(var)
	#upper = M.masked_where(blah<=max/100.,blah)
	#upper = M.masked_where(blah>=-max/100.,blah)
	lower = M.masked_where(blah>=-max/100.,blah)
	ax = subplot(panels+str(m))
	contourf(n,freq,var,crange,extend='both')
	#contourf(n,freq,upper,crange,extend='max')
	#contourf(n,freq,lower,crange,extend='min')
	#colorbar()
	xlabel('Zonal wavenumber',fontsize=9)
	if m==1: 
		ylabel('Frequency [1/days]',fontsize=9)
	else: 
		yticklabels = ax.get_yticklabels()
		setp(yticklabels, visible=False)
	title(titles[day],fontsize=10)
	xlim(0,4)
	xticks(arange(0,5,1),fontsize=8)
	yticks(fontsize=9)
	ylim(-2,2)
subplots_adjust(top=0.83,bottom=0.2)
show()
'''

## Figure 10
'''
# Phasespeed spectrum timelapse 400 hPa
radius = 0.0435
Ro ='10.5'
win = 20*8
a    = 6400.e3*radius
xmin = -10*radius/0.125
xmax = 60*radius/0.125
figdir = '/Users/mitch/Desktop/Geoff/movies/figures/phasespeed/'

crange = arange(-1000+5.,1001,250)

#li = 190*8
#lf = 200*8
li = 220*8
lf = 301*8
panels='22'

lev  = '400'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc'
#file = dir+'earth_1K_km_lapse_3hr/atmos_3hr.nc'
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])
km = argmin(abs(p-float(lev)))
u = array(d.variables['ucomp'][li:lf,km,:,:])
uave1 = average(u,axis=-1)
v = array(d.variables['vcomp'][li:lf,km,:,:])
temp = array(d.variables['temp'][li:lf,:,:,:])
lon = array(d.variables['lon'][:])
ni = len(lon)
lat = array(d.variables['lat'][:])
nj = len(lat)
time1 = array(d.variables['time'][li:lf])
nt   = len(time1)
lmax = nt/win
gz1 = zeros((nt,nk+1,len(lat),len(lon)),'d')
for k in range(nk-1,1,-1):
    dlnp = dp/p[k]
    gz1[:,k,:,:] = gz1[:,k,:,:] - 287.*temp[:,k,:,:]*dlnp
gz = 0.5*(gz1[:,0:nk,:,:]+gz1[:,1:nk+1,:,:])
gz = gz[:,km,:,:]
y = a*lat*pi/180.

dt = abs(time1[1]-time1[0])*86400.

# fft
nfreq= win
n    = fftfreq(ni,d=1./ni)
freq = fftfreq(nfreq,d=dt)
upr_upr = zeros((lmax,win,nj,ni),'d')
upr_vpr = upr_upr.copy()
gzpr = upr_upr.copy()
uave_sm = zeros((lmax,nj),'d')
uave_nosm = zeros((lmax),'d')
time_sm = zeros((lmax),'d')
for l in range(lmax):
    print l
    ufft = fft2(u[l*win:(l+1)*win],axes=(0,-1))
    vfft = fft2(v[l*win:(l+1)*win],axes=(0,-1))
    gzfft = fft2(gz[l*win:(l+1)*win],axes=(0,-1))
    time_sm[l] = average(time1[l*win:(l+1)*win])
    uave_sm[l,:] = average(uave1[l*win:(l+1)*win,:],axis=0)
    # power spectrum of velocity perturbations
    upr_upr[l] = 2.*real(ufft*conj(ufft))
    # power spectrum of geopotential perturbations
    gzpr[l] = 2.*real(gzfft*conj(gzfft))
    # horizontal eddy cospectrum (stress) term contributing to the zonal mean
    upr_vpr[l] = 2.*real(ufft*conj(vfft))

nsort = argsort(n)
fsort = argsort(freq)
n = n[nsort]
freq = freq[fsort]
for j in range(nj):
    print j
    for l in range(lmax):
	# Sort for plotting purposes
	for ll in range(nfreq):
	    upr_vpr[l,ll,j,:] = upr_vpr[l,ll,j,nsort]
	    upr_upr[l,ll,j,:] = upr_upr[l,ll,j,nsort]
	    gzpr[l,ll,j,:] = gzpr[l,ll,j,nsort]
	for i in range(len(n)):
	    upr_vpr[l,:,j,i] = upr_vpr[l,fsort,j,i]
	    upr_upr[l,:,j,i] = upr_upr[l,fsort,j,i]
	    gzpr[l,:,j,i] = gzpr[l,fsort,j,i]

dupr_vpr_dy1 = zeros((lmax,nfreq,nj+1,ni),'d')
for j in range(nj-1):
    dupr_vpr_dy1[:,:,j+1] = (upr_vpr[:,:,j+1]-upr_vpr[:,:,j])/(y[j+1]-y[j])
dupr_vpr_dy = 0.5*(dupr_vpr_dy1[:,:,0:nj]+dupr_vpr_dy1[:,:,1:nj+1])

cmin = -60*radius/0.125
cmax = 60*radius/0.125
dcp = (cmax-cmin)/nfreq
cnew = arange(cmin,cmax,dcp)  
upr_vpr_new = zeros((lmax,len(cnew),nj,ni),'d')
dupr_vpr_dy_new = zeros((lmax,len(cnew),nj,ni),'d')
upr_upr_new = zeros((lmax,len(cnew),nj,ni),'d')
gzpr_new = zeros((lmax,len(cnew),nj,ni),'d')
 
figure(figsize=(6,6))
m=0
fold = 2.*pi*freq
for l in range(lmax):
    m+=1
    for i in range(ni/2+1,ni,1):
	print i,n[i]
	for j in range(nj):
	    cold = -2.*pi*a*cos(lat[j]*pi/180.)*freq[:]/n[i]
	    fnew = -cnew*n[i]/a/cos(lat[j]*pi/180.)
	    func = spline(fold,upr_vpr[l,:,j,i],s=1.,k=5)
	    funcdy = spline(fold,dupr_vpr_dy[l,:,j,i],s=1.,k=5)
	    funcu = spline(fold,upr_upr[l,:,j,i],s=1.,k=5)
	    funcgz = spline(fold,gzpr[l,:,j,i],s=1.,k=5)
	    new = func(fnew)
	    newdy = funcdy(fnew)
	    newu = funcu(fnew)
	    newgz = funcgz(fnew)
	    upr_vpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),new,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    dupr_vpr_dy_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newdy,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    upr_upr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newu,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    gzpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newgz,0.)*n[i]/a/cos(lat[j]*pi/180.)
    #contourf(cnew,lat,transpose(sum(dupr_vpr_dy_new[l,:,:,:],axis=-1))*a,arange(-10000,10100,2000),extend='both')
    #colorbar()
    #show()

    quad = panels+str(m)
    ax = subplot(quad)
    contourf(cnew,lat,transpose(sum(upr_vpr_new[l,:,:,:],axis=-1)),crange,extend='both',cmap=cm.binary)
    #contourf(cnew,lat,transpose(sum(dupr_vpr_dy_new[l,:,:,:],axis=-1))*a,arange(-10,11,2),extend='both',cmap=cm.binary)
    #contourf(cnew,lat,transpose(sum(upr_upr_new[l,:,:,:],axis=-1)),cmap=cm.binary)#,crange,extend='both')
    #contourf(cnew,lat,transpose(sum(gzpr_new[l,:,:,:],axis=-1)))#,crange,extend='both',cmap=cm.binary)
    plot(uave_sm[l],lat,'b')
    xlim(xmin,xmax)
    ylim(lat[0],lat[-1])
    xticks(fontsize=8)
    if mod(m,2)==1:
	yticks(arange(-80,90,20),fontsize=8)
	ylabel('latitude',fontsize=9)
    else:
	yticks(arange(-40,50,20),fontsize=8)
	yticklabels = ax.get_yticklabels()
	setp(yticklabels, visible=False)
    if m>2:
	xlabel('speed [m/s]',fontsize=9)
    title('Day %s' %(time_sm[l]),fontsize=10)
gcf().text(0.05,0.925,r'$Ro_T$='+Ro, \
	   fontproperties=FontProperties(size=10))
subplots_adjust(bottom=0.17,top=0.87,wspace=0.05)
show()

'''


## Figure 11
'''
# Phasespeed spectrum timelapse 700 hPa
radius = 0.0435
Ro ='10.5'
win = 20*8
a    = 6400.e3*radius
xmin = -10*radius/0.125
xmax = 60*radius/0.125
figdir = '/Users/mitch/Desktop/Geoff/movies/figures/phasespeed/'

crange = arange(-1000+5.,1001,250)

#li = 190*8
#lf = 200*8
li = 220*8
lf = 301*8
panels='22'

lev  = '700'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc'
#file = dir+'earth_1K_km_lapse_3hr/atmos_3hr.nc'
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])
km = argmin(abs(p-float(lev)))
u = array(d.variables['ucomp'][li:lf,km,:,:])
uave1 = average(u,axis=-1)
v = array(d.variables['vcomp'][li:lf,km,:,:])
temp = array(d.variables['temp'][li:lf,:,:,:])
lon = array(d.variables['lon'][:])
ni = len(lon)
lat = array(d.variables['lat'][:])
nj = len(lat)
time1 = array(d.variables['time'][li:lf])
nt   = len(time1)
lmax = nt/win
gz1 = zeros((nt,nk+1,len(lat),len(lon)),'d')
for k in range(nk-1,1,-1):
    dlnp = dp/p[k]
    gz1[:,k,:,:] = gz1[:,k,:,:] - 287.*temp[:,k,:,:]*dlnp
gz = 0.5*(gz1[:,0:nk,:,:]+gz1[:,1:nk+1,:,:])
gz = gz[:,km,:,:]
y = a*lat*pi/180.

dt = abs(time1[1]-time1[0])*86400.

# fft
nfreq= win
n    = fftfreq(ni,d=1./ni)
freq = fftfreq(nfreq,d=dt)
upr_upr = zeros((lmax,win,nj,ni),'d')
upr_vpr = upr_upr.copy()
gzpr = upr_upr.copy()
uave_sm = zeros((lmax,nj),'d')
uave_nosm = zeros((lmax),'d')
time_sm = zeros((lmax),'d')
for l in range(lmax):
    print l
    ufft = fft2(u[l*win:(l+1)*win],axes=(0,-1))
    vfft = fft2(v[l*win:(l+1)*win],axes=(0,-1))
    gzfft = fft2(gz[l*win:(l+1)*win],axes=(0,-1))
    time_sm[l] = average(time1[l*win:(l+1)*win])
    uave_sm[l,:] = average(uave1[l*win:(l+1)*win,:],axis=0)
    # power spectrum of velocity perturbations
    upr_upr[l] = 2.*real(ufft*conj(ufft))
    # power spectrum of geopotential perturbations
    gzpr[l] = 2.*real(gzfft*conj(gzfft))
    # horizontal eddy cospectrum (stress) term contributing to the zonal mean
    upr_vpr[l] = 2.*real(ufft*conj(vfft))

nsort = argsort(n)
fsort = argsort(freq)
n = n[nsort]
freq = freq[fsort]
for j in range(nj):
    print j
    for l in range(lmax):
	# Sort for plotting purposes
	for ll in range(nfreq):
	    upr_vpr[l,ll,j,:] = upr_vpr[l,ll,j,nsort]
	    upr_upr[l,ll,j,:] = upr_upr[l,ll,j,nsort]
	    gzpr[l,ll,j,:] = gzpr[l,ll,j,nsort]
	for i in range(len(n)):
	    upr_vpr[l,:,j,i] = upr_vpr[l,fsort,j,i]
	    upr_upr[l,:,j,i] = upr_upr[l,fsort,j,i]
	    gzpr[l,:,j,i] = gzpr[l,fsort,j,i]

dupr_vpr_dy1 = zeros((lmax,nfreq,nj+1,ni),'d')
for j in range(nj-1):
    dupr_vpr_dy1[:,:,j+1] = (upr_vpr[:,:,j+1]-upr_vpr[:,:,j])/(y[j+1]-y[j])
dupr_vpr_dy = 0.5*(dupr_vpr_dy1[:,:,0:nj]+dupr_vpr_dy1[:,:,1:nj+1])

cmin = -60*radius/0.125
cmax = 60*radius/0.125
dcp = (cmax-cmin)/nfreq
cnew = arange(cmin,cmax,dcp)  
upr_vpr_new = zeros((lmax,len(cnew),nj,ni),'d')
dupr_vpr_dy_new = zeros((lmax,len(cnew),nj,ni),'d')
upr_upr_new = zeros((lmax,len(cnew),nj,ni),'d')
gzpr_new = zeros((lmax,len(cnew),nj,ni),'d')
 
figure(figsize=(6,6))
m=0
fold = 2.*pi*freq
for l in range(lmax):
    m+=1
    for i in range(ni/2+1,ni,1):
	print i,n[i]
	for j in range(nj):
	    cold = -2.*pi*a*cos(lat[j]*pi/180.)*freq[:]/n[i]
	    fnew = -cnew*n[i]/a/cos(lat[j]*pi/180.)
	    func = spline(fold,upr_vpr[l,:,j,i],s=1.,k=5)
	    funcdy = spline(fold,dupr_vpr_dy[l,:,j,i],s=1.,k=5)
	    funcu = spline(fold,upr_upr[l,:,j,i],s=1.,k=5)
	    funcgz = spline(fold,gzpr[l,:,j,i],s=1.,k=5)
	    new = func(fnew)
	    newdy = funcdy(fnew)
	    newu = funcu(fnew)
	    newgz = funcgz(fnew)
	    upr_vpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),new,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    dupr_vpr_dy_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newdy,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    upr_upr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newu,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    gzpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newgz,0.)*n[i]/a/cos(lat[j]*pi/180.)
    #contourf(cnew,lat,transpose(sum(dupr_vpr_dy_new[l,:,:,:],axis=-1))*a,arange(-10000,10100,2000),extend='both')
    #colorbar()
    #show()

    quad = panels+str(m)
    ax = subplot(quad)
    contourf(cnew,lat,transpose(sum(upr_vpr_new[l,:,:,:],axis=-1)),crange,extend='both',cmap=cm.binary)
    #contourf(cnew,lat,transpose(sum(dupr_vpr_dy_new[l,:,:,:],axis=-1))*a,arange(-10,11,2),extend='both',cmap=cm.binary)
    #contourf(cnew,lat,transpose(sum(upr_upr_new[l,:,:,:],axis=-1)),cmap=cm.binary)#,crange,extend='both')
    #contourf(cnew,lat,transpose(sum(gzpr_new[l,:,:,:],axis=-1)))#,crange,extend='both',cmap=cm.binary)
    plot(uave_sm[l],lat,'b')
    xlim(xmin,xmax)
    ylim(lat[0],lat[-1])
    xticks(fontsize=8)
    if mod(m,2)==1:
	yticks(arange(-80,90,20),fontsize=8)
	ylabel('latitude',fontsize=9)
    else:
	yticks(arange(-40,50,20),fontsize=8)
	yticklabels = ax.get_yticklabels()
	setp(yticklabels, visible=False)
    if m>2:
	xlabel('speed [m/s]',fontsize=9)
    title('Day %s' %(time_sm[l]),fontsize=10)
gcf().text(0.05,0.925,r'$Ro_T$='+Ro, \
	   fontproperties=FontProperties(size=10))
subplots_adjust(bottom=0.17,top=0.87,wspace=0.05)
show()
'''

## Figure 12

# Phase speed comparison at steady-state
radius = 0.0435
Ro ='10.5'
win = 20*8
a    = 6400.e3*radius
xmin = -10*radius/0.125
xmax = 60*radius/0.125

crange = arange(-1000+100.,1001,250)

li = 340*8
lf = 360*8
#li = 220*8
#lf = 301*8

lev  = '400'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr_yr4/atmos_3hr.nc'
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])
km = argmin(abs(p-float(lev)))
u = array(d.variables['ucomp'][li:lf,km,:,:])
uave1 = average(u,axis=-1)
v = array(d.variables['vcomp'][li:lf,km,:,:])
lon = array(d.variables['lon'][:])
ni = len(lon)
lat = array(d.variables['lat'][:])
nj = len(lat)
time1 = array(d.variables['time'][li:lf])
nt   = len(time1)
lmax = nt/win

dt = abs(time1[1]-time1[0])*86400.

# fft
nfreq= win
n    = fftfreq(ni,d=1./ni)
freq = fftfreq(nfreq,d=dt)
upr_vpr = zeros((lmax,win,nj,ni),'d')
uave_sm = zeros((nj),'d')
ufft = fft2(u[0:win],axes=(0,-1))
vfft = fft2(v[0:win],axes=(0,-1))
time_sm = average(time1[0:win])
uave_sm = average(uave1[0:win,:],axis=0)
# horizontal eddy cospectrum (stress) term contributing to the zonal mean
upr_vpr = 2.*real(ufft*conj(vfft))

nsort = argsort(n)
fsort = argsort(freq)
n = n[nsort]
freq = freq[fsort]
for j in range(nj):
	# Sort for plotting purposes
	for ll in range(nfreq):
	    upr_vpr[ll,j,:] = upr_vpr[ll,j,nsort]
	for i in range(len(n)):
	    upr_vpr[:,j,i] = upr_vpr[fsort,j,i]

cmin = -60*radius/0.125
cmax = 60*radius/0.125
dcp = (cmax-cmin)/nfreq
cnew = arange(cmin,cmax,dcp)  
upr_vpr_new1 = zeros((len(cnew),nj,ni),'d')

fold = 2.*pi*freq
for i in range(ni/2+1,ni,1):
	print i,n[i]
	for j in range(nj):
		cold = -2.*pi*a*cos(lat[j]*pi/180.)*freq[:]/n[i]
		fnew = -cnew*n[i]/a/cos(lat[j]*pi/180.)		
		func = interp1d(fold,upr_vpr[:,j,i])
		for ll in range(len(fnew)):
			try: upr_vpr_new1[ll,j,i] = func(fnew[ll])*n[i]/a/cos(lat[j]*pi/180.)
			except: pass

lev  = '700'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr_yr4/atmos_3hr.nc'
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])
km = argmin(abs(p-float(lev)))
u = array(d.variables['ucomp'][li:lf,km,:,:])
uave1 = average(u,axis=-1)
v = array(d.variables['vcomp'][li:lf,km,:,:])
lon = array(d.variables['lon'][:])
ni = len(lon)
lat = array(d.variables['lat'][:])
nj = len(lat)
time1 = array(d.variables['time'][li:lf])
nt   = len(time1)
lmax = nt/win

dt = abs(time1[1]-time1[0])*86400.

# fft
nfreq= win
n    = fftfreq(ni,d=1./ni)
freq = fftfreq(nfreq,d=dt)
upr_vpr = zeros((lmax,win,nj,ni),'d')
uave_sm = zeros((nj),'d')
ufft = fft2(u[0:win],axes=(0,-1))
vfft = fft2(v[0:win],axes=(0,-1))
time_sm = average(time1[0:win])
uave_sm = average(uave1[0:win,:],axis=0)
# horizontal eddy cospectrum (stress) term contributing to the zonal mean
upr_vpr = 2.*real(ufft*conj(vfft))

nsort = argsort(n)
fsort = argsort(freq)
n = n[nsort]
freq = freq[fsort]
for j in range(nj):
	# Sort for plotting purposes
	for ll in range(nfreq):
	    upr_vpr[ll,j,:] = upr_vpr[ll,j,nsort]
	for i in range(len(n)):
	    upr_vpr[:,j,i] = upr_vpr[fsort,j,i]

cmin = -60*radius/0.125
cmax = 60*radius/0.125
dcp = (cmax-cmin)/nfreq
cnew = arange(cmin,cmax,dcp)  
upr_vpr_new2 = zeros((len(cnew),nj,ni),'d')

fold = 2.*pi*freq
for i in range(ni/2+1,ni,1):
	print i,n[i]
	for j in range(nj):
		cold = -2.*pi*a*cos(lat[j]*pi/180.)*freq[:]/n[i]
		fnew = -cnew*n[i]/a/cos(lat[j]*pi/180.)		
		func = interp1d(fold,upr_vpr[:,j,i])
		for ll in range(len(fnew)):
			try: upr_vpr_new2[ll,j,i] = func(fnew[ll])*n[i]/a/cos(lat[j]*pi/180.)
			except: pass


figure(figsize=(6,3))
ax = subplot(121)
contourf(cnew,lat,transpose(sum(upr_vpr_new1[:,:,:],axis=-1)),crange,extend='both',cmap=cm.binary)
plot(uave_sm,lat,'b')
xlim(xmin,xmax)
ylim(lat[0],lat[-1])		
xticks(fontsize=8)
yticks(arange(-80,90,20),fontsize=8)
xlabel('speed [m/s]',fontsize=9)
title('400 hPa',fontsize=10)
ax = subplot(122)
contourf(cnew,lat,transpose(sum(upr_vpr_new2[:,:,:],axis=-1)),crange,extend='both',cmap=cm.binary)
plot(uave_sm,lat,'b')
xlim(xmin,xmax)
ylim(lat[0],lat[-1])		
xticks(fontsize=8)
yticks(arange(-80,90,20),fontsize=8)
xlabel('speed [m/s]',fontsize=9)
title('700 hPa',fontsize=10)
gcf().text(0.05,0.925,r'$Ro_T$='+Ro, \
	   fontproperties=FontProperties(size=10))
subplots_adjust(bottom=0.17,top=0.87)
show()



## Figure 13
'''
# Global wave timelapse
freq0 = -0.6
nh = 1
#Ro = '1.3'
Ro = '10.5'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc'
#file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_T42/atmos_daily.nc'
dt = 3.*3600.
a = 280.e3
g = 9.8
om = 2.*pi/86400.
d = ncopen(file,'r')
stepsperday = 8

levs = ['400','700','900']
if nh==1: crange = arange(-2.5,2.6,0.25)
if nh==2: crange = arange(-0.5,0.55,0.05)
crange=arange(-5,6,1)

tmin = 260*stepsperday
#nt = 6*8
nt = 20*stepsperday
tmax = tmin+1*nt
panels = '13'
npanels=3
nl = nt
figure(figsize=(5.5,2))

p = array(d.variables['pfull'][:])*1.e2
dp = (p[0]-p[1])
p0 = 1.e5
nk = len(p)
lat = array(d.variables['lat'][:])
dlat = (lat[0]-lat[1])*pi/180.
nj = len(lat)
cosphi = cos(lat*pi/180.)
y = a*sin(lat*pi/180.)
lon = array(d.variables['lon'][:])
ni = len(lon)
n    = fftfreq(len(lon),d=1./len(lon))
freq = fftfreq(nl,d=dt)
n3d = zeros((nl,len(lat),len(lon)),'d')
freq3d = n3d.copy()
f3d = zeros((nl,nj,ni),'d')
cosphi2d = zeros((len(lat),len(lon)),'d')
for j in range(nj):
    for l in range(nl):
	n3d[l,j,:] = n
	f3d[l,j,:] = 2.*om*sin(lat[j]*pi/180.)
for i in range(len(lon)):
    cosphi2d[:,i] = cosphi
    for j in range(len(lat)):
	freq3d[:,j,i] = freq


lmin = tmin
lmax = lmin+nt
for ll in range(npanels):
    km = argmin(abs(p/100.-float(levs[ll])))
    # Geopotential calculation
    T = array(d.variables['temp'][lmin:lmax,:,:,:])
    gz1 = zeros((nl,nk+1,len(lat),len(lon)),'d')
    for k in range(nk-1,1,-1):
	dlnp = dp/p[k]
	gz1[:,k,:,:] = gz1[:,k,:,:] - 287.*T[:,k,:,:]*dlnp
	gz = 0.5*(gz1[:,0:nk,:,:]+gz1[:,1:nk+1,:,:])
    gz = gz[:,km,:,:]
    u = array(d.variables['ucomp'][lmin:lmax,:,:,:])
    v = array(d.variables['vcomp'][lmin:lmax,km,:,:])
    vor = array(d.variables['vor'][lmin:lmax,km,:,:])
    

    # PV calculation
    thetaave = T.copy()
    for i in range(ni):
	for j in range(nj):
	    for l in range(nl):
		thetaave[l,:,j,i]*=1./exner
    dthetaavedy1 = zeros((nl,nk,nj+1,ni),'d')
    phi = lat*pi/180.
    for j in range(nj-1):
	dthetaavedy1[:,:,j+1] = (thetaave[:,:,j]-thetaave[:,:,j+1]) \
	    /a/(phi[j]-phi[j+1])
    dthetaavedy1[:,:,0] = dthetaavedy1[:,:,1]
    dthetaavedy1[:,:,-1] = dthetaavedy1[:,:,nj-1]
    dthetaavedy = 0.5*(dthetaavedy1[:,:,0:nj]+dthetaavedy1[:,:,1:nj+1])
    dthetaavedp1 = zeros((nl,nk+1,nj,ni),'d')
    duavedp1 = zeros((nl,nk+1,nj,ni),'d')
    for k in range(nk-1):
	dthetaavedp1[:,k+1,:] = (thetaave[:,k,:]-thetaave[:,k+1,:])/dp
	duavedp1[:,k+1,:] = (u[:,k,:]-u[:,k+1,:])/dp
    dthetaavedp1[:,0] = dthetaavedp1[:,1]
    dthetaavedp1[:,-1] = dthetaavedp1[:,nk-1]
    dthetaavedp = 0.5*(dthetaavedp1[:,0:nk]+dthetaavedp1[:,1:nk+1])
    duavedp1[:,0] = duavedp1[:,1]
    duavedp1[:,-1] = duavedp1[:,nk]
    duavedp = 0.5*(duavedp1[:,0:nk]+duavedp1[:,1:nk+1])
    PV = -g*( (f3d+vor)*dthetaavedp[:,km] + duavedp[:,km]*dthetaavedy[:,km] )

    # Filtering
    u = u[:,km]
    freq3d *=86400.
    freqmax = -0.2
    gzfft = fft2(gz,axes=(0,-1))
    #gzfft = where(freq3d==freq0,gzfft,0.)
    #gzfft = where(freq3d>freqmin,gzfft,0.)
    gzfft = where(freq3d<freqmax,gzfft,0.)
    gzfft = where(n3d==nh,gzfft,0.)
    gz = ifft2(gzfft,axes=(0,-1))

    PVfft = fft2(PV,axes=(0,-1))
    #vorfft = where(freq3d==freq0,vorfft,0.)
    #vorfft = where(freq3d>freqmin,vorfft,0.)
    PVfft = where(freq3d<freqmax,PVfft,0.)
    PVfft = where(n3d==nh,PVfft,0.)
    PV = ifft2(PVfft,axes=(0,-1))

    l = nt/2
    subplot(panels+str(ll+1))
    uave = average(u[l],axis=-1)
    contourf(lon,lat,gz[l],crange,extend='both',cmap=cm.binary)
    #contourf(lon,lat,gz[l],cmap=cm.binary)
    #colorbar()
    #contour(lon,lat,PV[l]*1.e7,arange(-1,1.1,.2),colors='k')
    #contour(lon,lat,PV[l],30,colors='k')
    #quiver(lon,lat,u2[l],v2[l])
    #contour(lon,lat,vor[l],15,colors='k')
    title(levs[ll]+' hPa',fontsize=10)
    plot(uave*5,lat,'b',linewidth=2)
    yticks(arange(-80,90,20),fontsize=8)
    xticks(arange(0,361,90),fontsize=8)
    xlim(0,360)
    ylim(lat[0],lat[-1])
    if ll+1==1: ylabel('latitude [deg]',fontsize=9)
    xlabel(r'longitude or speed$\times5$',fontsize=9)

subplots_adjust(top=0.85,bottom=0.23,wspace=0.24)
gcf().text(0.925,0.35,'Day %i' %((tmin+nl/2)/stepsperday), 
	horizontalalignment='center', \
	rotation=-90,\
	fontproperties=FontProperties(size=12))
show()
    #savefig('/Users/mitch/Desktop/Geoff/movies/figures/10.5/N1/%s.png' %(l))
'''

## Figures 14 & 15
'''
# Vertical/horizontal structure
freq0 = -0.25
lats = [0.,40.,80.]
levs = [400.,700.,900.]
stepsperday=8
lmin = 340*stepsperday
lmax = 360*stepsperday
nl = lmax-lmin
nh = 1
Ro = '10.5'
#file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr_yr4/atmos_3hr.nc'
#file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_lat90to45NS_sym/atmos_daily.nc'
dt = 3.*3600.
a = 280.e3
g = 9.8
om = 2.*pi/86400.
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])*1.e2
dp = (p[0]-p[1])
p0 = 1.e5
nk = len(p)
lat = array(d.variables['lat'][:])
dlat = (lat[0]-lat[1])*pi/180.
nj = len(lat)
cosphi = cos(lat*pi/180.)
y = a*sin(lat*pi/180.)
lon = array(d.variables['lon'][:])
ni = len(lon)
n    = fftfreq(len(lon),d=1./len(lon))
freq = fftfreq(nl,d=dt)
n4d = zeros((nl,len(p),len(lat),len(lon)),'d')
freq4d = n4d.copy()
for j in range(nj):
	for k in range(nk):	
		for l in range(nl):
			n4d[l,k,j,:] = n
for i in range(len(lon)):
	for j in range(nj):
		for k in range(nk):
			freq4d[:,k,j,i] = freq

# Geopotential calculation
T = array(d.variables['temp'][lmin:lmax,:,:,:])
time = array(d.variables['time'][lmin:lmax])
print time[-1]
gz1 = zeros((nl,nk+1,len(lat),len(lon)),'d')
for k in range(nk-1,1,-1):
    dlnp = dp/p[k]
    gz1[:,k,:,:] = gz1[:,k,:,:] - 287.*T[:,k,:,:]*dlnp
gz = 0.5*(gz1[:,0:nk,:,:]+gz1[:,1:nk+1,:,:])
freq4d *=86400.
print shape(gz)

for k in range(nk):
    print k
    gzfft = fft2(gz[:,k],axes=(0,-1))
    gzfft = where(freq4d[:,k]==freq0,gzfft,0.)
    gzfft = where(n4d[:,k]==nh,gzfft,0.)
    gz[:,k] = ifft2(gzfft,axes=(0,-1))
figure(figsize=(6,2))
crange = arange(-0.5,0.51,0.1)

# Vertical 
m = 0
for lat0 in lats:
	m+=1
	jm = argmin(abs(lat-lat0))
	print jm
	gzh = gz[:,:,jm,:]
	gzhave = average(gzh,axis=-1)
	gzhp = gzh.copy()
	for i in range(ni):
		gzhp[:,:,i]-=gzhave
	lev = p/1.e5
	subplot('13'+str(m))
	contourf(lon,lev,gzhp[nl/2]/9.8,crange,cmap=cm.binary,extend='both')
	ylim(lev[-1],lev[0])
	xlim(lon[0],lon[-1])
	yticks(fontsize=8)
	if m==1: ylabel(r'$p/p_s$',fontsize=10)
	xlabel('longitude [deg]',fontsize=9)
	xticks(arange(0,360,90),fontsize=8)
	title('%iN latitude' %(lat0),fontsize=12)
subplots_adjust(bottom=0.22,top=0.84,wspace=0.22)
gcf().text(0.925,0.35,'Day %i' %(time[nl/2]), 
	horizontalalignment='center', \
	rotation=-90,\
	fontproperties=FontProperties(size=12))
show()

# Horizontal
m = 0
for lev0 in levs:
	m+=1
	km = argmin(abs(p/100.-lev0))
	print km
	gzh = gz[:,km,:,:]
	lev = p/1.e5
	subplot('13'+str(m))
	gzhave = average(gzh,axis=-1)
	gzhp = gzh.copy()
	for i in range(ni):
		gzhp[:,:,i]-=gzhave
	contourf(lon,lat,gzhp[nl/2]/9.8,crange,cmap=cm.binary,extend='both')
	ylim(lat[0],lat[-1])
	xlim(lon[0],lon[-1])
	yticks(arange(-80,90,20),fontsize=8)
	if m==1: ylabel('latitude',fontsize=10)
	xlabel('longitude [deg]',fontsize=9)
	xticks(arange(0,360,90),fontsize=8)
	title('%i hPa' %(lev0),fontsize=12)
subplots_adjust(bottom=0.22,top=0.84,wspace=0.22)
gcf().text(0.925,0.35,'Day %i' %(time[nl/2]), 
	horizontalalignment='center', \
	rotation=-90,\
	fontproperties=FontProperties(size=12))
show()
'''

## Figure 16
'''
# PV spinup time sequence
R='10.5'
#R='1.3'
d = {}
d = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_T42/atmos_daily.nc','r')
lat = d.variables['lat'][:]
lev = d.variables['pfull'][:]/1000.
nj = len(lat)
nk = len(lev)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
exner  = lev**(2./7.)
f = 2.*om*sinphi
lat2 = zeros((nk,nj),'d')
f2   = lat2.copy()
lev2 = lat2.copy()
for j in range(nj):
    for k in range(nk):
        f2[k,j]   = f[j]
        lat2[k,j] = lat[j]
        lev2[k,j] = lev[k]
sinphi2 = sin(lat2*pi/180.)
cosphi2 = cos(lat2*pi/180.)

#days = [165,166,167,168,169,170]
days = [180,200,220,240,260,280]
kb = argmin(abs(lev-0.4))
panels = '23'

figure(figsize=(6,4))
l=0
for li in days:
    print li
    l+=1
    Rd = 287.
    a = r(float(R))
    y = a*sinphi
    uave = average(d.variables['ucomp'][li,:,:,:],axis=-1)
    vorave = average(d.variables['vor'][li,:,:,:],axis=-1)
    Tave = average(d.variables['temp'][li,:,:,:],axis=-1)
    thetaave = Tave.copy()
    print d.variables['time'][li]
    for j in range(nj):
	thetaave[:,j]=Tave[:,j]/exner
    dthetaavedy1 = zeros((nk,nj+1),'d')
    phi = lat*pi/180.
    for j in range(nj-1):
	dthetaavedy1[:,j+1] = (thetaave[:,j]-thetaave[:,j+1]) \
	    /a/(phi[j]-phi[j+1])
	    # /(y[j]-y[j+1])
    dthetaavedy1[:,0] = dthetaavedy1[:,1]
    dthetaavedy1[:,-1] = dthetaavedy1[:,nj]
    dthetaavedy = 0.5*(dthetaavedy1[:,0:nj]+dthetaavedy1[:,1:nj+1])
    dthetaavedp1 = zeros((nk+1,nj),'d')
    duavedp1 = zeros((nk+1,nj),'d')
    for k in range(nk-1):
	dthetaavedp1[k+1,:] = (thetaave[k,:]-thetaave[k+1,:])/dp
	duavedp1[k+1,:] = (uave[k,:]-uave[k+1,:])/dp
    dthetaavedp1[0] = dthetaavedp1[1]
    dthetaavedp1[-1] = dthetaavedp1[nk-1]
    dthetaavedp = 0.5*(dthetaavedp1[0:nk]+dthetaavedp1[1:nk+1])
    duavedp1[0] = duavedp1[1]
    duavedp1[-1] = duavedp1[nk]
    duavedp = 0.5*(duavedp1[0:nk]+duavedp1[1:nk+1])
    # PV = (f2+vorave)
    fV = -g*f2*dthetaavedp
    vV = -g*vorave*dthetaavedp
    uV = -g*duavedp*dthetaavedy
    PV = -g*( (f2+vorave)*dthetaavedp + duavedp*dthetaavedy)
    quad = panels+str(l)
    ax = subplot(quad)
    #plot(lat,fV[kb]*1.e8,'m',linewidth=0.5)
    #plot(lat,uV[kb]*1.e8,'b',linewidth=0.5)
    #plot(lat,vV[kb]*1.e8,'g',linewidth=0.5)
    plot(lat,PV[kb]*1.e8,'k:')
    plot(lat,uave[kb],'k')
    plot(lat,uave[kb]*0.,linewidth=0.5,color=(0.5,0.5,0.5))
    ymax = 40
    if mod(l,3)==1:
	# ylabel(r'$(\bar{f}+\bar{\zeta})\cdot2\times10^5$, $\bar{u}$',fontsize=12)
	ylabel(r'$\bar{q}\times10^8$, $\bar{u}$',fontsize=12)
	yticks(arange(-40,50,20),fontsize=8)
    else:
	yticklabels = ax.get_yticklabels()
	setp(yticklabels, visible=False)
    if l<4: 
	title('Day '+str(li),fontsize=10)
	xticks(arange(-80,90,40),fontsize=8)
	xticklabels = ax.get_xticklabels()
	setp(xticklabels, visible=False)	
	xlim(lat[0],lat[-1])
    else:
	title('Day '+str(li),fontsize=10)
	xlim(lat[0],lat[-1])
	xlabel('latitude',fontsize=9)
	xticks(arange(-80,90,40),fontsize=8)
    subplots_adjust(wspace=0.,bottom=0.15)
    ylim(-ymax,ymax)
show()
'''

## Figure 17
'''
# PV
d = {}
Ro = ['10.5']
for R in Ro:
    print R
    d[R] = ncopen(dir+'Ro'+R+'_noseas_mlapse_nuv0.01_T42/atmos_monthly.nc','r')
lat = d[R].variables['lat'][:]
lev = d[R].variables['pfull'][:]/1000.
nj = len(lat)
nk = len(lev)
sinphi = sin(lat*pi/180.)
cosphi = cos(lat*pi/180.)
exner  = lev**(2./7.)
f = 2.*om*sinphi
lat2 = zeros((nk,nj),'d')
f2   = lat2.copy()
lev2 = lat2.copy()
nth = 100
thetanew = arange(250.,400.,150./nth)
PVnew = zeros((nth,nj),'d')
for j in range(nj):
    for k in range(nk):
        f2[k,j]   = f[j]
        lat2[k,j] = lat[j]
        lev2[k,j] = lev[k]
sinphi2 = sin(lat2*pi/180.)
cosphi2 = cos(lat2*pi/180.)

kbs = [argmin(abs(lev-0.2)),argmin(abs(lev-0.5))]
#kbs = [argmin(abs(lev-0.4)),argmin(abs(lev-0.2))]
th1 = 320.
th2 = 290.
khs = [argmin(abs(thetanew-th1)),argmin(abs(thetanew-th2))]

li = 25
lf = 36
stepsperday=30
#figure(figsize=(6,3))
l=0
for kh in khs:
    print thetanew[kh]
    for R in Ro:
        l+=1
        Rd = 287.
        a = r(float(R))
        y = a*lat*pi/180.
        uave = average(average(d[R].variables['ucomp'][li:lf,:,:,:],axis=0),axis=-1)
        Tave = average(average(d[R].variables['temp'][li:lf,:,:,:],axis=0),axis=-1)
        thetaave = Tave.copy()
        for j in range(nj):
            thetaave[:,j]=Tave[:,j]/exner
        dthetaavedy1 = zeros((nk,nj+1),'d')
	phi = lat*pi/180.
        for j in range(nj-1):
            dthetaavedy1[:,j+1] = (thetaave[:,j]-thetaave[:,j+1]) \
                                  /(y[j]-y[j+1])
        dthetaavedy1[:,0] = dthetaavedy1[:,1]
        dthetaavedy1[:,-1] = dthetaavedy1[:,nj-1]
        dthetaavedy = 0.5*(dthetaavedy1[:,0:nj]+dthetaavedy1[:,1:nj+1])
        dTavedp1 = zeros((nk+1,nj),'d')
        dthetaavedp1 = zeros((nk+1,nj),'d')
        duavedp1 = zeros((nk+1,nj),'d')
        for k in range(nk-1):
            dTavedp1[k+1,:] = (Tave[k,:]-Tave[k+1,:])/dp
            dthetaavedp1[k+1,:] = (thetaave[k,:]-thetaave[k+1,:])/dp
            duavedp1[k+1,:] = (uave[k,:]-uave[k+1,:])/dp
        dTavedp1[0] = dTavedp1[1]
        dTavedp1[-1] = dTavedp1[nk-1]
        dTavedp = 0.5*(dTavedp1[0:nk]+dTavedp1[1:nk+1])
        dthetaavedp1[0] = dthetaavedp1[1]
        dthetaavedp1[-1] = dthetaavedp1[nk-1]
        dthetaavedp = 0.5*(dthetaavedp1[0:nk]+dthetaavedp1[1:nk+1])
        duavedp = 0.5*(duavedp1[0:nk]+duavedp1[1:nk+1])
        vor = average(average(d[R].variables['vor'][li:lf,:,:,:],axis=-1),axis=0)
        PV = -g*( (f2+vor)*dthetaavedp + duavedp*dthetaavedy )
	#PV*=1.e6

	#contourf(lat,lev,thetaave,arange(270,400,3),extend='both')
	#colorbar()
	#clabel(contour(lat,lev,thetaave,[th1,th2],colors='k'),fmt='%3i')
	#xticks(arange(-90,100,30))
	#ylim(lev[-1],lev[0])
	#xlim(lat[0],lat[-1])
	#show()

	dPVdy1 = zeros((nk,nj+1),'d')
	phi = lat*pi/180.
	for j in range(nj-1):
	    #dPVdy1[:,j+1] = (PVnew[:,j]-PVnew[:,j+1])/(y[j]-y[j+1])
	    dPVdy1[:,j+1] = (PV[:,j]-PV[:,j+1])/(phi[j]-phi[j+1])
	dPVdy1[:,0] = dPVdy1[:,1]
	dPVdy1[:,-1] = dPVdy1[:,nj]
	dPVdy = M.array(0.5*(dPVdy1[:,0:nj]+dPVdy1[:,1:nj+1]))
	figure(figsize=(3,3))
	
	# on isobars
	contourf(lat,lev,1.e6*dPVdy,arange(-0.5,0.51,.05),extend='both')
	#colorbar()
	contour(lat,lev,1.e6*dPVdy,[0.],colors='m')
	clabel(contour(lat,lev,thetaave,[295.,300.,305.],colors='k'),fmt='%3i')
	ylim(lev[-1],lev[0])
	xlabel('latitude [deg]',fontsize=10)
	ylabel(r'$p/p_s$',fontsize=12)
	title(r'Days %s to %s' %(li*30,(lf-1)*30),fontsize=12)
	xticks(arange(-90,100,30),fontsize=9)
	yticks(fontsize=9)
	xlim(lat[0],lat[-1])
	subplots_adjust(left=0.18,bottom=0.14,top=0.87)
	show()
	
	# on isentropes
	for j in range(nj):
		func = interp1d(thetaave[::-1,j],PV[::-1,j])
		for k in range(nth):
			try: PVnew[k,j] = func(thetanew[k])
			except: pass
 
	dPVdy1 = zeros((nth,nj+1),'d')
	phi = lat*pi/180.
	for j in range(nj-1):
	    #dPVdy1[:,j+1] = (PVnew[:,j]-PVnew[:,j+1])/(y[j]-y[j+1])
	    dPVdy1[:,j+1] = (PVnew[:,j]-PVnew[:,j+1])/(phi[j]-phi[j+1])
	dPVdy1[:,0] = dPVdy1[:,1]
	dPVdy1[:,-1] = dPVdy1[:,nj]
	dPVdy = M.array(0.5*(dPVdy1[:,0:nj]+dPVdy1[:,1:nj+1]))
	for j in range(nj):
	    dPVdy[:,j] = M.masked_where(thetanew<thetaave[-1,j],dPVdy[:,j])
	dPVdy*=1.e6
	contourf(lat,thetanew,dPVdy,arange(-0.5,0.51,.05),extend='both')
	contour(lat,thetanew,dPVdy,arange(2,63,10),colors='w')
	contour(lat,thetanew,dPVdy,[0.],colors='k')
	title('Days %i to %i' %(stepsperday*li,stepsperday*lf),fontsize=12)
	if R=='10.5':  ylim(280,320)
	if R=='0.02':  ylim(260,340)
	if R=='1.3':  ylim(265,340)
	xlim(lat[0],lat[-1])
	xlabel('latitude [deg]',fontsize=10)
	ylabel('Potential temperature [K]',fontsize=12)
	xticks(arange(-90,100,30),fontsize=9)
	yticks(fontsize=9)
	subplots_adjust(left=0.18,bottom=0.14,top=0.87)
	show()

	PVnew*=1.e6
	contourf(lat,thetanew,PVnew,arange(-0.5,0.51,.05),extend='both')
	#colorbar()
	contour(lat,thetanew,PVnew,arange(-20,21,1.),colors='w')
	contour(lat,thetanew,PVnew,[0.],colors='m')
	ylabel('Potential Temperature [K]')
	xlabel('latitude [deg]')
	title(r'Ro=%s PV, days %s to %s' %(R,li*30,(lf-1)*30))
	if R=='0.02': ylim(260,340)
	if R=='10.5': ylim(280,340)
	if R=='1.3': ylim(250,340)
	xticks(arange(-80,90,20))
	xlim(lat[0],lat[-1])
	show()
'''

## Figure 18 
# barotropic.pdf


## Figure 19
'''
# Ro=1.3 Phasespeed spectrum timelapse 400 hPa
radius = 0.125
Ro ='1.3'
win = 20*8
a    = 6400.e3*radius
xmin = -10*radius/0.125
xmax = 60*radius/0.125
figdir = '/Users/mitch/Desktop/Geoff/movies/figures/phasespeed/'

crange = arange(-1000+10.,1001,250)

li = 220*8
lf = 301*8
panels='22'

lev  = '400'
file = dir+'Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc'
#file = dir+'earth_1K_km_lapse_3hr/atmos_3hr.nc'
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])
km = argmin(abs(p-float(lev)))
u = array(d.variables['ucomp'][li:lf,km,:,:])
uave1 = average(u,axis=-1)
v = array(d.variables['vcomp'][li:lf,km,:,:])
temp = array(d.variables['temp'][li:lf,:,:,:])
lon = array(d.variables['lon'][:])
ni = len(lon)
lat = array(d.variables['lat'][:])
nj = len(lat)
time1 = array(d.variables['time'][li:lf])
nt   = len(time1)
lmax = nt/win
gz1 = zeros((nt,nk+1,len(lat),len(lon)),'d')
for k in range(nk-1,1,-1):
    dlnp = dp/p[k]
    gz1[:,k,:,:] = gz1[:,k,:,:] - 287.*temp[:,k,:,:]*dlnp
gz = 0.5*(gz1[:,0:nk,:,:]+gz1[:,1:nk+1,:,:])
gz = gz[:,km,:,:]
y = a*lat*pi/180.

dt = abs(time1[1]-time1[0])*86400.

# fft
nfreq= win
n    = fftfreq(ni,d=1./ni)
freq = fftfreq(nfreq,d=dt)
upr_upr = zeros((lmax,win,nj,ni),'d')
upr_vpr = upr_upr.copy()
gzpr = upr_upr.copy()
uave_sm = zeros((lmax,nj),'d')
uave_nosm = zeros((lmax),'d')
time_sm = zeros((lmax),'d')
for l in range(lmax):
    print l
    ufft = fft2(u[l*win:(l+1)*win],axes=(0,-1))
    vfft = fft2(v[l*win:(l+1)*win],axes=(0,-1))
    gzfft = fft2(gz[l*win:(l+1)*win],axes=(0,-1))
    time_sm[l] = average(time1[l*win:(l+1)*win])
    uave_sm[l,:] = average(uave1[l*win:(l+1)*win,:],axis=0)
    # power spectrum of velocity perturbations
    upr_upr[l] = 2.*real(ufft*conj(ufft))
    # power spectrum of geopotential perturbations
    gzpr[l] = 2.*real(gzfft*conj(gzfft))
    # horizontal eddy cospectrum (stress) term contributing to the zonal mean
    upr_vpr[l] = 2.*real(ufft*conj(vfft))

nsort = argsort(n)
fsort = argsort(freq)
n = n[nsort]
freq = freq[fsort]
for j in range(nj):
    print j
    for l in range(lmax):
	# Sort for plotting purposes
	for ll in range(nfreq):
	    upr_vpr[l,ll,j,:] = upr_vpr[l,ll,j,nsort]
	    upr_upr[l,ll,j,:] = upr_upr[l,ll,j,nsort]
	    gzpr[l,ll,j,:] = gzpr[l,ll,j,nsort]
	for i in range(len(n)):
	    upr_vpr[l,:,j,i] = upr_vpr[l,fsort,j,i]
	    upr_upr[l,:,j,i] = upr_upr[l,fsort,j,i]
	    gzpr[l,:,j,i] = gzpr[l,fsort,j,i]

dupr_vpr_dy1 = zeros((lmax,nfreq,nj+1,ni),'d')
for j in range(nj-1):
    dupr_vpr_dy1[:,:,j+1] = (upr_vpr[:,:,j+1]-upr_vpr[:,:,j])/(y[j+1]-y[j])
dupr_vpr_dy = 0.5*(dupr_vpr_dy1[:,:,0:nj]+dupr_vpr_dy1[:,:,1:nj+1])

cmin = -60*radius/0.125
cmax = 60*radius/0.125
dcp = (cmax-cmin)/nfreq
cnew = arange(cmin,cmax,dcp)  
upr_vpr_new = zeros((lmax,len(cnew),nj,ni),'d')
dupr_vpr_dy_new = zeros((lmax,len(cnew),nj,ni),'d')
upr_upr_new = zeros((lmax,len(cnew),nj,ni),'d')
gzpr_new = zeros((lmax,len(cnew),nj,ni),'d')
 
figure(figsize=(6,6))
m=0
fold = 2.*pi*freq
for l in range(lmax):
    m+=1
    for i in range(ni/2+1,ni,1):
	print i,n[i]
	for j in range(nj):
	    cold = -2.*pi*a*cos(lat[j]*pi/180.)*freq[:]/n[i]
	    fnew = -cnew*n[i]/a/cos(lat[j]*pi/180.)
	    func = spline(fold,upr_vpr[l,:,j,i],s=1.,k=5)
	    funcdy = spline(fold,dupr_vpr_dy[l,:,j,i],s=1.,k=5)
	    funcu = spline(fold,upr_upr[l,:,j,i],s=1.,k=5)
	    funcgz = spline(fold,gzpr[l,:,j,i],s=1.,k=5)
	    new = func(fnew)
	    newdy = funcdy(fnew)
	    newu = funcu(fnew)
	    newgz = funcgz(fnew)
	    upr_vpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),new,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    dupr_vpr_dy_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newdy,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    upr_upr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newu,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    gzpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newgz,0.)*n[i]/a/cos(lat[j]*pi/180.)
    #contourf(cnew,lat,transpose(sum(dupr_vpr_dy_new[l,:,:,:],axis=-1))*a,arange(-10000,10100,2000),extend='both')
    #colorbar()
    #show()

    quad = panels+str(m)
    ax = subplot(quad)
    contourf(cnew,lat,transpose(sum(upr_vpr_new[l,:,:,:],axis=-1)),crange,extend='both',cmap=cm.binary)
    #contourf(cnew,lat,transpose(sum(dupr_vpr_dy_new[l,:,:,:],axis=-1))*a,arange(-10,11,2),extend='both',cmap=cm.binary)
    #contourf(cnew,lat,transpose(sum(upr_upr_new[l,:,:,:],axis=-1)),cmap=cm.binary)#,crange,extend='both')
    #contourf(cnew,lat,transpose(sum(gzpr_new[l,:,:,:],axis=-1)))#,crange,extend='both',cmap=cm.binary)
    plot(uave_sm[l],lat,'b')
    xlim(xmin,xmax)
    ylim(lat[0],lat[-1])
    xticks(fontsize=8)
    if mod(m,2)==1:
	yticks(arange(-80,90,20),fontsize=8)
	ylabel('latitude',fontsize=9)
    else:
	yticks(arange(-40,50,20),fontsize=8)
	yticklabels = ax.get_yticklabels()
	setp(yticklabels, visible=False)
    if m>2:
	xlabel('speed [m/s]',fontsize=9)
    title('Day %s' %(time_sm[l]),fontsize=10)
gcf().text(0.05,0.925,r'$Ro_T$='+Ro, \
	   fontproperties=FontProperties(size=10))
subplots_adjust(bottom=0.17,top=0.87,wspace=0.05)
show()
'''
