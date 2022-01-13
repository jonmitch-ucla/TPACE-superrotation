from pylab import *
from numpy import *
from mpl_toolkits.basemap import NetCDFFile as ncopen
import matplotlib.numerix.ma as M
from matplotlib.font_manager import FontProperties as FP
from matplotlib import rc
from matplotlib import cm
from scipy.interpolate import UnivariateSpline as spline
from matplotlib.font_manager import FontProperties

rc('text', usetex=True, dvipnghack=True)
rc('font', **{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('font', size=8)

Ros = ['0.02','1.3','10.5']
als = ['0','0.25','0.75','10']
kas = ['','_ka200','_ka400_2','_ka800_2']
radii = {}
radii['10.5'] = 280.e3
radii['1.3'] = 280.e3*10.5/1.3
radii['0.02'] = 6400.e3
radii['0.07'] = 3400.e3

d = {}
om = 2.*pi/86400.
g = 9.8
R = 287.

def m(lat,u,a):
    uave = average(average(u,axis=0),axis=-1)
    uave *= cos(lat[None,:])/om/a
    print shape(uave)
    return uave + cos(lat[None,:])**2.

def psi(v,lev,lat):
    cosphi = cos(lat)
    nk = len(lev)
    nj = len(lat)
    psi1 = zeros((nk+1,nj),'d')
    p1 = zeros(nk+1,'d')
    dp = abs(lev[0]-lev[1])*100.
    for k in range(nk-1):
	psi1[k+1] = psi1[k]+v[k]*dp/g*2.*pi*a*cosphi
    psi1[-1] = 0.
    psi1[nk] = 0.
    return 0.5*(psi1[0:nk]+psi1[1:nk+1])

# Evolution of equatorial, 250hPa winds
'''
dir = '/Volumes/Data/SeasonalSuperrotation/'
lev  = '250'
figdir = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/'+lev+'hPa_equatorial_winds.pdf'
print figdir
m = 0
for Ro in Ros:
    a = radii[Ro]
    for al in als:
	m += 1
	subplot('33'+str(m))
	file = dir+'Ro'+Ro+'_al'+al+'_mlapse_T42/atmos_daily.nc'
	print file
	d = ncopen(file,'r')
	p = d.variables['pfull'][:]
	lat = d.variables['lat'][:]
	lon = d.variables['lon'][:]
	time = array(d.variables['time'][:])
	time -= time[0]
	kh = argmin(abs(p-float(lev)))
	jh = argmin(abs(lat))
	u = average(d.variables['ucomp'][:,kh,jh,:],axis=-1)
	plot(time,u)
	xlim(time[0],time[-1])
	ylim(-30,30)
	d.close()
savefig(figdir)
'''

# Annual mean winds
'''
figure(figsize=(5,4))
dir = '/Volumes/Data/SeasonalSuperrotation/'
lev  = '250'
figdir = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/'+lev+'hPa_mean_winds.pdf'
print figdir
ubar = zeros((3,4),'d')
Rosim = zeros((3,4),'d')
for i in range(3):
    Ro = Ros[i]
    for j in range(4):
	al = als[j]
	file = dir+'Ro'+Ro+'_al'+al+'_mlapse_T42/atmos_daily.nc'
	print file
	d = ncopen(file,'r')
	p = d.variables['pfull'][:]
	lat = d.variables['lat'][:]
	lon = d.variables['lon'][:]
	time = array(d.variables['time'][:])
	time -= time[0]
	kh = argmin(abs(p-float(lev)))
	jh = argmin(abs(lat))
	ubar[i,j] = average(d.variables['ucomp'][:,kh,jh,:])
x = array(als).astype('d')
y = array(Ros).astype('d')
ax = subplot(111)
cf = contourf(ubar)
contour(ubar,[0.],linewidths=1.5)
#ax.set_xscale('log')
#ax.set_yscale('log')
#xlim(x[0],x[-1])
#ylim(y[0],y[-1])
xticks([0,1,2,3],als,fontsize=10)
yticks([0,1,2],Ros,fontsize=10)
xlabel(r'$\alpha$',fontsize=14)
ylabel(r'$Ro_T$',fontsize=14)
colorbar(cf)
title('Annual mean 250hPa, equatorial zonal wind',fontsize=14)
savefig(figdir)
'''

# Max annual mean Rossby number

figure(figsize=(5,4))
dir = '/Users/mitch/data/SeasonalSuperrotation/'
#dir = '/Volumes/Data/SeasonalSuperrotation/'
lev  = '250'
figdir = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/'+lev+'hPa_max_Ro.pdf'
print figdir
Rosim = zeros((3,4),'d')
latmax = zeros((3,4),'d')
for i in range(3):
    Ro = Ros[i]
    a = radii[Ro]
    for j in range(4):
	al = als[j]
	file = dir+'Ro'+Ro+'_al'+al+'_mlapse_T42/atmos_daily.nc'
	print file
	d = ncopen(file,'r')
	p = d.variables['pfull'][:]
	lat = d.variables['lat'][:]
	lon = d.variables['lon'][:]
	time = array(d.variables['time'][:])
	time -= time[0]
	kh = argmin(abs(p-float(lev)))
	jh = argmin(abs(lat))
	vorzon = average(average(d.variables['vor'][:,kh,:,:],axis=0),axis=-1)
	vmax = vorzon[argmax(abs(vorzon))]
	latmax[i,j] = lat[argmax(abs(vorzon))]
	print latmax
	latmax*=pi/180.
	#Rosim[i,j] = abs(vmax)
	Rosim[i,j] = abs(vmax/(2.*om*sin(latmax[i,j])))
x = array(als).astype('d')
y = array(Ros).astype('d')
ax = subplot(111)
#cf = contourf(latmax)
cf = contourf(Rosim)
contour(Rosim,[1.],linewidths=1.5)
xticks([0,1,2,3],als,fontsize=10)
yticks([0,1,2],Ros,fontsize=10)
xlabel(r'$\alpha$',fontsize=14)
ylabel(r'$Ro_T$',fontsize=14)
colorbar(cf)
title('Annual mean 250hPa, max Ro',fontsize=14)
#savefig(figdir)


# Annual mean superrotation index
'''
figure(figsize=(5,4))
dir = '/Volumes/Data/SeasonalSuperrotation/'
lev  = '250'
figdir = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/'+lev+'hPa_mean_superrot_index.pdf'
print figdir
SI = zeros((3,4),'d')
for i in range(3):
    Ro = Ros[i]
    a = radii[Ro]
    for j in range(4):
	al = als[j]
	file = dir+'Ro'+Ro+'_al'+al+'_mlapse_T42/atmos_daily.nc'
	print file
	d = ncopen(file,'r')
	p = d.variables['pfull'][:]
	lat = d.variables['lat'][:]
	lon = d.variables['lon'][:]
	time = array(d.variables['time'][:])
	time -= time[0]
	kh = argmin(abs(p-float(lev)))
	jh = argmin(abs(lat))
	u = d.variables['ucomp'][:]
	l = m(lat*pi/180.,u,a)
	print shape(l)
	SI[i,j] = average(average(l,axis=0),weights=cos(lat*pi/180.))
x = array(als).astype('d')
y = array(Ros).astype('d')
ax = subplot(111)
cf = contourf(SI)
contour(SI,[1.],linewidths=1.5)
#ax.set_xscale('log')
#ax.set_yscale('log')
#xlim(x[0],x[-1])
#ylim(y[0],y[-1])
xticks([0,1,2,3],als,fontsize=10)
yticks([0,1,2],Ros,fontsize=10)
xlabel(r'$\alpha$',fontsize=14)
ylabel(r'$Ro_T$',fontsize=14)
colorbar(cf)
title('Annual mean superrotation index',fontsize=14)
savefig(figdir)
'''

# Max annual streamfunction
'''
figure(figsize=(5,4))
dir = '/Volumes/Data/SeasonalSuperrotation/'
figdir = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/max_streamfunction.pdf'
print figdir
psimax = zeros((3,4),'d')
for i in range(3):
    Ro = Ros[i]
    a = radii[Ro]
    for j in range(4):
	al = als[j]
	file = dir+'Ro'+Ro+'_al'+al+'_mlapse_T42/atmos_daily.nc'
	print file
	d = ncopen(file,'r')
	p = d.variables['pfull'][:]
	lat = d.variables['lat'][:]
	lon = d.variables['lon'][:]
	lev = d.variables['pfull'][:]
	psimax[i,j] = 0.
	for l in range(360):
	    v = average(d.variables['vcomp'][l,:,:,:],axis=-1)
	    tmp = psi(v,lev,lat)
	    tmp = tmp.flatten()
	    maxpsi = tmp[argmax(abs(tmp))]
	    if maxpsi>psimax[i,j]: psimax[i,j]=maxpsi
x = array(als).astype('d')
y = array(Ros).astype('d')
ax = subplot(111)
cf = contourf(psimax)
#ax.set_xscale('log')
#ax.set_yscale('log')
#xlim(x[0],x[-1])
#ylim(y[0],y[-1])
xticks([0,1,2,3],als,fontsize=10)
yticks([0,1,2],Ros,fontsize=10)
xlabel(r'$\alpha$',fontsize=14)
ylabel(r'$Ro_T$',fontsize=14)
colorbar(cf)
title('Max meridional mass flux',fontsize=14)
savefig(figdir)
'''

# Phasespeed spectrum timelapse 
'''
dir = '/Volumes/Data/SeasonalSuperrotation/'
radius = 0.0435
Ro ='10.5'
al = '0.75'
win = 20*8
a    = 6400.e3*radius
xmin = -10*radius/0.125
xmax = 40. 
#xmax = 60*radius/0.125
cmin = -10. 
cmax = 40.
crange = arange(-1000+5.,1001,250)

#li = 190*8
#lf = 200*8
li = 20*8 
lf = 180*8
panels='42'

lev  = '250'
file = dir+'Ro'+Ro+'_al'+al+'_mlapse_T42_3hr/atmos_3hr.nc'
figdir = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/'+Ro+'_'+al+'_'+lev+'hPa'+'_days'+str(li/8)+'to'+str(lf/8)+'.pdf'
print figdir
#file = dir+'earth_1K_km_lapse_3hr/atmos_3hr.nc'
d = ncopen(file,'r')

p = array(d.variables['pfull'][:])
nk = len(p)
dp = abs(p[0]-p[1])*100.
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
time1 -= time1[0]
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

dcp = (cmax-cmin)/nfreq
cnew = arange(cmin,cmax,dcp)  
upr_vpr_new = zeros((lmax,len(cnew),nj,ni),'d')
cold3d = zeros((len(cnew),nj),'d')
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
	    cold = abs(-2.*pi*a*cos(lat[j]*pi/180.)*freq[:]/n[i])
	    if l==0 and n[i]==1.0: cold3d[:,j] = cold
	    fnew = -cnew*n[i]/a/cos(lat[j]*pi/180.)
	    func = spline(fold,upr_vpr[l,:,j,i],s=1.,k=5)
	    #funcdy = spline(fold,dupr_vpr_dy[l,:,j,i],s=1.,k=5)
	    #funcu = spline(fold,upr_upr[l,:,j,i],s=1.,k=5)
	    #funcgz = spline(fold,gzpr[l,:,j,i],s=1.,k=5)
	    new = func(fnew)
	    #newdy = funcdy(fnew)
	    #newu = funcu(fnew)
	    #newgz = funcgz(fnew)
	    upr_vpr_new[l,:,j,i] = where(abs(cnew)<=cold.max(),new,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    #dupr_vpr_dy_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newdy,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    #upr_upr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newu,0.)*n[i]/a/cos(lat[j]*pi/180.)
	    #gzpr_new[l,:,j,i] = where(abs(cnew)<=abs(cold),newgz,0.)*n[i]/a/cos(lat[j]*pi/180.)
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
    text(30,55,'Day %i' %(time_sm[l]+1),fontsize=8)
gcf().text(0.05,0.925,r'$Ro_T=%s, \alpha=%s$' %(Ro,al), \
	   fontproperties=FontProperties(size=10))
subplots_adjust(bottom=0.17,top=0.87,wspace=0.05)
savefig(figdir)
#show()
'''

# Half-year-mean angular momentum and hadley cell
'''
i = 0
for Ro in Ros:
    a = radii[Ro]
    for al in als:
	i+=1
	panel = '22'+str(i)
	d = ncopen('Ro'+Ro+'_al'+al+'_monthly.nc','r') 
	lat = d.variables['lat'][:]
	lev = d.variables['pfull'][:]
	u = d.variables['ucomp'][31:]
	v = d.variables['vcomp'][31:]
	subplot(panel)
	contourf(lat,lev,m(lat*pi/180.,u,a))
	colorbar()
	contour(lat,lev,psi(v,lev,lat*pi/180.),10,linewidths=0.5,colors='k')
	title('Ro ='+Ro+', al = '+al)
	if i==1 or i==3: ylabel('pressure [mbar]')
	if i==3 or i==4: xlabel('latitude')
	xticks(arange(-90,100,30))
	xlim(lat[0],lat[-1])
	ylim(lev[-1],lev[0])
show()
'''

# Half-year-mean hadley cell
'''
i = 0
for Ro in Ros:
    a = radii[Ro]
    for ka in kas:
	i+=1
	panel = '22'+str(i)
	d = ncopen('Ro'+Ro+'_al10'+ka+'_monthly.nc','r') 
	lat = d.variables['lat'][:]
	lev = d.variables['pfull'][:]
	u = d.variables['ucomp'][31:]
	v = d.variables['vcomp'][31:]
	subplot(panel)
	contourf(lat,lev,m(lat*pi/180.,u,a))
	colorbar()
	hc = psi(v,lev,lat*pi/180.)*1.e-9
	contour(lat,lev,hc,10,linewidths=0.5,colors='k')
	max = M.max(hc)
	maxind = hc.argmax()
	kh,jh = unravel_index(maxind,shape(hc))
	if Ro=='10.5': blah = '%1.1f' %max
	else: blah = '%i' %max
	text(lat[jh],lev[kh],blah,color='w',fontsize=8)
	title('Ro ='+Ro+' , '+ka)
	if i==1 or i==3: ylabel('pressure [mbar]')
	if i==3 or i==4: xlabel('latitude')
	xticks(arange(-90,100,30))
	xlim(lat[0],lat[-1])
	ylim(lev[-1],lev[0])
show()
'''

# Flipbook of U' and Ubar
'''
panels = '42'
lev = '250'
d = ncopen('Ro10.5_al1_3hr.nc','r')
kh = argmin(abs(d.variables['pfull'][:]-float(lev)))
lat = d.variables['lat'][:]
lon = d.variables['lon'][:]
i = 0
for l in range(1310,1318):
    i+=1
    u = d.variables['ucomp'][l,kh]
    v = d.variables['vcomp'][l,kh]
    gz = d.variables['height'][l,kh]
    uave = u.mean(axis=-1)
    up = u-uave[:,None]
    vp = v-v.mean(axis=-1)[:,None]
    gzp = gz-gz.mean(axis=-1)[:,None]
    subplot(panels+str(i))
    contourf(lon,lat,gzp,arange(-5,5.1,.5))
    colorbar()
    plot(uave*5,lat)
    title(str(l))
show()
'''

# Seasonal momentum fluxes
# seasons defined in 60-day (pi/3) increments
'''
i = 0
lis = [660,780,840,960]
lfs = [780,840,960,1020]
label = ['Fall','Winter','Spring','Summer']
Ro = '1.3'
al = '0.01'
tag = ''
lev = '250'
a = radii[Ro]
for i in range(len(lis)):
    li = lis[i]
    lf = lfs[i]
    panel = '22'+str(i+1)
    print panel
    d = ncopen('Ro'+Ro+'_al'+al+tag+'_mlapse_T42/atmos_daily.nc','r') 
    #d = ncopen('/Volumes/Data/Ro1.3_noseas_mlapse_nuv0.01_T42/atmos_daily.nc','r') 
    lat = d.variables['lat'][:]
    lon = d.variables['lon'][:]
    p = d.variables['pfull'][:]
    kh = argmin(abs(p-float(lev)))
    u = d.variables['ucomp'][li:lf,kh]
    v = d.variables['vcomp'][li:lf,kh]
    ubar = u.mean(axis=-1).mean(axis=0)
    vbar = v.mean(axis=-1).mean(axis=0)
    figure(1)
    subplot(panel)
    plot(lat,ubar)
    xticks(arange(-90,100,30))
    xlim(lat[0],lat[-1])
    ylim(-40,70)
    title(label[i])
    if i==0 or i==2: ylabel('Zonal winds [m/s]')
    if i==2 or i==3: xlabel('Latitude')
    up = u - ubar[None,:,None]
    vp = v - vbar[None,:,None]
    upvp = up*vp
    upvpbar = upvp.mean(axis=-1)
    ubarvbar = ubar*vbar
    figure(2)
    subplot(panel)
    plot(lat,upvpbar.mean(axis=0)*cos(lat*pi/180.))
    plot(lat,ubarvbar*cos(lat*pi/180.))
    ylim(-10,10)
    xticks(arange(-90,100,30))
    xlim(lat[0],lat[-1])
    title(label[i])
    if i==0 or i==2: ylabel('Momentum flux [m^2/s^2]')
    if i==2 or i==3: xlabel('Latitude')
show()
'''

# A-seasonal momentum fluxes
'''
lev = '250'

save = True
outfile = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/aseasonal_'+lev+'_mom_flux_winds.pdf'

upy = {}
upy['10.5'] = 1
upy['1.3'] = 5
upy['0.02'] = 50
uy = {}
uy['10.5'] = (-20,59)
uy['1.3'] = (-20,59)
uy['0.02'] = (-20,59)
f6 = FP(size=6)
f8 = FP(size=8)
f10 = FP(size=10)
f12 = FP(size=12)

Ro = '0.02'
a = radii[Ro]
d = ncopen('/Volumes/Data/Ro'+Ro+'_noseas_mlapse_nuv0.01_T42/atmos_daily.nc','r') 
lat = d.variables['lat'][:]
nj = len(lat)
dphi = abs(lat[0]-lat[1])*pi/180.
lon = d.variables['lon'][:]
p = d.variables['pfull'][:]
kh = argmin(abs(p-float(lev)))
u = d.variables['ucomp'][720:,kh]
v = d.variables['vcomp'][720:,kh]
om = d.variables['omega'][720:,kh]
T = d.variables['temp'][720:,kh]
rho = 1./T.copy()
rho *= p[kh]*100./R
w = om.copy()/(-rho*g)
ubar = u.mean(axis=-1)
vbar = v.mean(axis=-1)
wbar = w.mean(axis=-1)
up = u-ubar[:,:,None]
vp = v-vbar[:,:,None]
wp = w-wbar[:,:,None]
upvp = up*vp
upwp = up*wp
ubarvbar = ubar*vbar
ubarwbar = ubar*wbar

f = figure(figsize=(6,3))

upvpbareq = upvp.mean(axis=0).mean(axis=-1) * cos(lat*pi/180.)
upwpbareq = upwp.mean(axis=0).mean(axis=-1) * cos(lat*pi/180.)
ubarvbareq = ubarvbar.mean(axis=0) * cos(lat*pi/180.)
ubarwbareq = ubarwbar.mean(axis=0) * cos(lat*pi/180.)
ubareq = ubar.mean(axis=0)
vbareq = vbar.mean(axis=0)
ax = subplot(231)
line = plot(lat,ubarvbareq,'k',label=r'$[\bar u \bar v ]$')
plot(lat,upvpbareq,'k:',label=r'$[ \overline{u^\prime v^\prime}]$')
#plot(lat,upwpbareq,'k--',label=r'$[ \overline{u^\prime w^\prime}]$')
plot(lat,0.*lat,'grey')
legend(prop=f6,loc=2)
ylim(-upy[Ro],upy[Ro])
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
title(Ro)
ylabel('Momentum flux [$m^2/s^2$]')
xticklabels = ax.get_xticklabels()
setp(xticklabels, visible=False)
ax = subplot(234)
plot(lat,ubareq,'k',label='Zonal')
plot(lat,vbareq*20,'k--',label='Mer. x 20')
plot(lat,0.*lat,'grey',linewidth=0.5)
legend(prop=f6,loc=2)
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
ylim(uy[Ro])
ylabel('Winds [$m/s$]')
xlabel('Latitude')

Ro = '1.3'
a = radii[Ro]
#d = ncopen('/Volumes/Data/Ro'+Ro+'_noseas_mlapse_nuv0.01_T42/atmos_daily.nc','r') 
d = ncopen('/Volumes/Data/Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc','r') 
u = d.variables['ucomp'][:,kh]
v = d.variables['vcomp'][:,kh]
ubar = u.mean(axis=-1)
vbar = v.mean(axis=-1)
up = u-ubar[:,:,None]
vp = v-vbar[:,:,None]
upvp = up*vp
ubarvbar = ubar*vbar
upvpbareq = upvp.mean(axis=0).mean(axis=-1) * cos(lat*pi/180.)
ubarvbareq = ubarvbar.mean(axis=0) * cos(lat*pi/180.)
ubareq = ubar.mean(axis=0)
vbareq = vbar.mean(axis=0)
ax = subplot(232)
plot(lat,ubarvbareq,'k',label=r'$[\bar u \bar v ]$')
plot(lat,upvpbareq,'k:',label=r'$[ \overline{u^\prime v^\prime}]$')
plot(lat,0.*lat,'grey')
ylim(-upy[Ro],upy[Ro])
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
title(Ro)
xticklabels = ax.get_xticklabels()
setp(xticklabels, visible=False)
ax = subplot(235)
plot(lat,ubareq,'k',label='Zonal')
plot(lat,vbareq*20,'k--',label='Mer. x 20')
plot(lat,0.*lat,'grey',linewidth=0.5)
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
ylim(uy[Ro])
xlabel('Latitude')

Ro = '10.5'
a = radii[Ro]
#d = ncopen('/Volumes/Data/Ro'+Ro+'_noseas_mlapse_nuv0.01_T42/atmos_daily.nc','r') 
d = ncopen('/Volumes/Data/Ro'+Ro+'_noseas_mlapse_nuv0.01_3hr/atmos_3hr.nc','r') 
u = d.variables['ucomp'][:,kh]
v = d.variables['vcomp'][:,kh]
ubar = u.mean(axis=-1)
vbar = v.mean(axis=-1)
up = u-ubar[:,:,None]
vp = v-vbar[:,:,None]
upvp = up*vp
ubarvbar = ubar*vbar
upvpbareq = upvp.mean(axis=0).mean(axis=-1) * cos(lat*pi/180.)
ubarvbareq = ubarvbar.mean(axis=0) * cos(lat*pi/180.)
ubareq = ubar.mean(axis=0)
vbareq = vbar.mean(axis=0)
ax = subplot(233)
plot(lat,ubarvbareq,'k',label=r'$[\bar u \bar v ]$')
plot(lat,upvpbareq,'k:',label=r'$[ \overline{u^\prime v^\prime}]$')
plot(lat,0.*lat,'grey')
ylim(-upy[Ro],upy[Ro])
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
title(Ro)
xticklabels = ax.get_xticklabels()
setp(xticklabels, visible=False)
ax = subplot(236)
plot(lat,ubareq,'k',label='Zonal')
plot(lat,vbareq*20,'k--',label='Mer. x 20')
plot(lat,0.*lat,'grey',linewidth=0.5)
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
ylim(uy[Ro])
xlabel('Latitude')

if save==True: savefig(outfile)
else: show()
'''

# Momentum fluxes, seasonal composites
'''
# 3-HOURLY DATA
# 60 -120 winter
# 120-240 spring
# 240-300 summer
# 0-60 & 300-360 fall

lev = '250'
Ro = '1.3'
al = '0.75'
a = radii[Ro]

save = True
outfile = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/'+Ro+'_'+al+'_'+lev+'_mom_flux_winds.pdf'

upy = {}
upy['10.5'] = 1
upy['1.3'] = 5
upy['0.02'] = 50
uy = {}
uy['10.5'] = (-20,59)
uy['1.3'] = (-20,59)
uy['0.02'] = (-20,59)
f6 = FP(size=6)
f8 = FP(size=8)
f10 = FP(size=10)
f12 = FP(size=12)

if Ro=='0.02': li = 720
else: li = 0
d = ncopen('/Volumes/Data/SeasonalSuperrotation/Ro'+Ro+'_al'+al+'_mlapse_T42_3hr/atmos_3hr.nc','r') 
lat = d.variables['lat'][:]
nj = len(lat)
dphi = abs(lat[0]-lat[1])*pi/180.
lon = d.variables['lon'][:]
p = d.variables['pfull'][:]
kh = argmin(abs(p-float(lev)))
u = d.variables['ucomp'][li:,kh]
v = d.variables['vcomp'][li:,kh]
om = d.variables['omega'][li:,kh]
T = d.variables['temp'][li:,kh]
rho = 1./T.copy()
rho *= p[kh]*100./R
w = om.copy()/(-rho*g)
ubar = u.mean(axis=-1)
vbar = v.mean(axis=-1)
wbar = w.mean(axis=-1)
up = u-ubar[:,:,None]
vp = v-vbar[:,:,None]
wp = w-wbar[:,:,None]
upvp = up*vp
upwp = up*wp
ubarvbar = ubar*vbar
ubarwbar = ubar*wbar

f = figure(figsize=(6,3))

# Annual
li,lf = 0,-1
upvpbareq = upvp[li:lf].mean(axis=0).mean(axis=-1) * cos(lat*pi/180.)
ubarvbareq = ubarvbar[li:lf].mean(axis=0) * cos(lat*pi/180.)
ubareq = ubar[li:lf].mean(axis=0)
vbareq = vbar[li:lf].mean(axis=0)
ax = subplot(231)
plot(lat,ubarvbareq,'k',label=r'$[\bar u \bar v ]$')
plot(lat,upvpbareq,'k:',label=r'$[ \overline{u^\prime v^\prime}]$')
plot(lat,0.*lat,'grey')
ylabel('Momentum flux [$m^2/s^2$]')
ylim(-upy[Ro],upy[Ro])
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
title('Annual')
xticklabels = ax.get_xticklabels()
setp(xticklabels, visible=False)
ax = subplot(234)
plot(lat,ubareq,'k',label='Zonal')
plot(lat,vbareq*20,'k--',label='Mer. x 20')
plot(lat,0.*lat,'grey',linewidth=0.5)
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
ylim(uy[Ro])
xlabel('Latitude')
ylabel('Winds [$m/s$]')

# Equinox
if Ro=='0.02': li,lf = 120,240
else: li,lf = 120*8,240*8
upvpbareq = upvp[li:lf].mean(axis=0).mean(axis=-1) * cos(lat*pi/180.)
ubarvbareq = ubarvbar[li:lf].mean(axis=0) * cos(lat*pi/180.)
ubareq = ubar[li:lf].mean(axis=0)
vbareq = vbar[li:lf].mean(axis=0)
ax = subplot(232)
plot(lat,ubarvbareq,'k',label=r'$[\bar u \bar v ]$')
plot(lat,upvpbareq,'k:',label=r'$[ \overline{u^\prime v^\prime}]$')
plot(lat,0.*lat,'grey')
ylim(-upy[Ro],upy[Ro])
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
title('Equinox')
yticklabels = ax.get_yticklabels()
setp(yticklabels, visible=False)
xticklabels = ax.get_xticklabels()
setp(xticklabels, visible=False)
ax = subplot(235)
plot(lat,ubareq,'k',label='Zonal')
plot(lat,vbareq*20,'k--',label='Mer. x 20')
plot(lat,0.*lat,'grey',linewidth=0.5)
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
ylim(uy[Ro])
xlabel('Latitude')
yticklabels = ax.get_yticklabels()
setp(yticklabels, visible=False)

# Solstice
if Ro=='0.02': li1,lf1,li2,lf2 = 120,180,240,300
else: li1,lf1,li2,lf2 = 120*8,180*8,240*8,300*8
upvpbarsols = 0.5*(upvp[li1:lf1].mean(axis=0).mean(axis=-1) \
	       +upvp[lf2:li2:-1].mean(axis=0).mean(axis=-1) ) * cos(lat*pi/180.)
ubarvbarsols = 0.5*(ubarvbar[li1:lf1].mean(axis=0) \
	       +ubarvbar[lf2:li2:-1].mean(axis=0) ) * cos(lat*pi/180.)
ubarsols = 0.5*(ubar[li1:lf1].mean(axis=0) \
	       +ubar[lf2:li2:-1].mean(axis=0) ) 
vbarsols = 0.5*(vbar[li1:lf1].mean(axis=0) \
		+vbar[lf2:li2:-1].mean(axis=0) ) 
ax = subplot(233)
plot(lat,ubarvbarsols,'k')
plot(lat,upvpbarsols,'k:')
plot(lat,0.*lat,'grey')
ylim(-upy[Ro],upy[Ro])
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
title('Solstice')
yticklabels = ax.get_yticklabels()
setp(yticklabels, visible=False)
xticklabels = ax.get_xticklabels()
setp(xticklabels, visible=False)
ax = subplot(236)
plot(lat,ubarsols,'k',label='Zonal')
plot(lat,vbarsols*20,'k--',label='Mer. x 20')
plot(lat,0.*lat,'grey',linewidth=0.5)
xticks(arange(-90,100,30))
xlim(lat[0],lat[-1])
ylim(uy[Ro])
xlabel('Latitude')
yticklabels = ax.get_yticklabels()
setp(yticklabels, visible=False)

subplots_adjust(wspace=0.02,hspace=0.02,bottom=0.15)

if save: savefig(outfile)
else: show()
'''

# Surface temperatures
'''
save = True
outfile = '/Users/mitch/Documents/Manuscripts/SeasonalE2T/alphas.pdf'

per   = 3.1104e7
delh = 60.
om  = 2.*pi/per
dt  = 0.005*per
nt = int(per/dt)
phi = arange(-pi/2.,pi/2.+1.e-3,pi/40.)
ny = len(phi)
time = arange(0.,per,dt)
phi1 = time.copy()
phi0 = time.copy()
phi0 = pi/2.*sin(om*time)
nt = len(time)
phi2d = array([l*0+phi for l in range(nt)])
phi02d = transpose(array([j*0+phi0 for j in range(ny)]))
time2d = transpose(array([j*0+time for j in range(ny)]))
T  = zeros((nt,ny),'d')
Te = T.copy()
Tbar = T.copy()
T1 = T.copy()

To = 285.
# Instantaneous forcing temperatures
T1 = To + delh/3. - delh*sin(phi2d)**2. + 2*delh*sin(om*time2d)*sin(phi2d)
# Annual average temperatures
Tbar = To + delh/3.*(1.-3.*sin(phi2d)**2.)

# Time-variable forcing for specified surface thermal inertia
def Te(al):
    tauf = 1./om/al
    print al,tauf
    return Tbar + 2*delh/tauf*sin(phi2d)*(sin(om*time2d)/tauf -    om*cos(om*time2d))/((1./tauf)**2. +   om**2.)

# Plot surface temperatures
f = figure(figsize=(6,1.5))
subplot(141)
#lines = contour(time/per,phi*180./pi,swapaxes(T1,0,1),arange(175,351,25),colors='k')
al = 10.
contourf(time/per,phi*180./pi,swapaxes(Te(al),0,1),arange(175,351,25),extend='both',cmap=cm.gray)
ylabel('latitude')
yticks(arange(-60,90,30))
ylim(-90,90)
title(r'$\alpha$ = %3.2f' %(al))
xlabel('time')

subplot(142)
#lines = contour(time/per,phi*180./pi,swapaxes(T1,0,1),arange(175,351,25),colors='k')
al = 0.75
contourf(time/per,phi*180./pi,swapaxes(Te(al),0,1),arange(175,351,25),extend='both',cmap=cm.gray)
yticks(arange(-60,90,30))
ylim(-90,90)
title(r'$\alpha$ = %3.2f' %(al))
xlabel('time')

subplot(143)
#lines = contour(time/per,phi*180./pi,swapaxes(T1,0,1),arange(175,351,25),colors='k')
al = 0.25
contourf(time/per,phi*180./pi,swapaxes(Te(al),0,1),arange(175,351,25),extend='both',cmap=cm.gray)
yticks(arange(-60,90,30))
ylim(-90,90)
title(r'$\alpha$ = %3.2f' %(al))
xlabel('time')

subplot(144)
#lines = contour(time/per,phi*180./pi,swapaxes(T1,0,1),arange(175,351,25),colors='k')
al = 1.e-19
contourf(time/per,phi*180./pi,swapaxes(Te(al),0,1),arange(175,351,25),extend='both',cmap=cm.gray)
xlabel('time')
yticks(arange(-60,90,30))
ylim(-90,90)
title(r'$\alpha$ = %3.2f' %(al))

subplots_adjust(bottom=0.23,top=0.8,wspace=.25)
if save: savefig(outfile)
else: show()
'''
