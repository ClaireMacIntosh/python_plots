# !/usr/bin/env python
#  -*- coding: utf-8 -*-
# Program to plot HC-TC trends.

import numpy as np
import netCDF4 as nCDF
import glob
import os.path
import datetime as datetime
import matplotlib.pyplot as plt
# #######################################################################
# 18 aug 2016
# making plots for chris to check somethign with another data provider
# READ SST monthly (1993-2010)
# READ SST monthly anomaly (1993-2010)
# Make Annual mean SST
# make 3 mn anomaly (DJf)
# Plot
#   Annual mean SST
#   Annual mean STDEV over 18 years
#   intermon anomaly stdev
#   quarterly DJF stdev
# ======================================================================
# def calc_ndays(yearin, monthin):

def make_anomaly(sst,nyr):
    # nyr=np.shape(sst[0,:,:])[0]/12.
    sst_ann1=np.reshape(sst,(12,nyr,720,1440),order='f')
    clim=np.nanmean(sst_ann1,axis=1)
    clim_rep=np.copy(clim)
    for year in range(1,nyr):
        clim_rep=np.append(clim_rep,clim, axis=0)
    anom=np.subtract(sst,clim_rep)
    return anom

def make_monthly_sst(filesin):
    sst_daily = [nCDF.Dataset(f) for f in filesin]
    sst_d = np.array([sst_daily[f].variables['sst'][:] for f in xrange(np.shape(filesin)[0])])
    sst_d[np.where(sst_d<0.)] = np.nan
    sst_mn = np.nanmean(sst_d,axis=0)
    return sst_mn

def makefilelist(ystart, ystop, path=''):
    filelist = []
    for yy in range(ystart, ystop+1):
        if yy % 4 != 0:
            filelist = []
            for doy in range(1, 365):
                f = os.path.join(path, '{Y}/mw.fusion.{Y}.{d}.v04.0.nc'.format(Y=yy, d="%03d" % doy))
                filelist.extend(glob.glob(f))
        elif yy % 4 == 0:
            for doy in range(1, 366):
                f = os.path.join(path, '{Y}/mw.fusion.{Y}.{d}.v04.0.nc'.format(Y=yy, d="%03d" % doy))
                filelist.extend(glob.glob(f))

    return filelist


# =======================================
# read SL CCI data ( to end 2014), reashape, remove annual cycle
# =============================================

year_start = 1998
year_stop = 2004
nyr=(year_stop-year_start+1)

dirin = "/glusterfs/surft/data/RSS-MW-OISST-v4.0_nc/"
# dirin = "/Users/clairemacintosh/Data/SST_daily_mw_OISST_v04.0/"


for year in range(year_start, year_stop + 1):
    print 'year = ', year
    file_list = makefilelist(year, year, dirin)
    doy_list = [int(filename[-12:-9]) for filename in file_list]
    month_list = np.array([(datetime.datetime(year, 1, 1) + datetime.timedelta(days - 1)).month for days in doy_list])

    for month in xrange(1,13,1):
        print 'month = ', month
        month_bounds = np.where(month_list==month)[0]
        sst_mn = np.expand_dims(make_monthly_sst(file_list[month_bounds[0]:(month_bounds[-1]+1)]),axis=0)
        if month == 1:
            sst_ann = sst_mn
        else:
            sst_ann = np.append(sst_ann,sst_mn,axis=0)

    if year == year_start:
        sst_all = sst_ann
    else:
        sst_all=np.append(sst_all,sst_ann,axis=0)

    print 'sst_files = ', np.shape(sst_all)

plt.imshow(sst_all[0])
plt.show()
plt.imshow(sst_all[10])
plt.show()

anom = make_anomaly(sst_all,nyr)

sst_ann1=np.reshape(sst_all,(12,nyr,720,1440),order='f')
sst_ann_avg=np.nanmean(sst_ann1,axis=0)

# Avg of the SST anom
sst_ann2=np.reshape(anom,(12,nyr,720,1440),order='f')
sst_ann_avg2=np.nanmean(sst_ann2,axis=0)
#stop
#anom=varin.variables['anomaly'][:]
anom_seas1=np.copy(anom[2:-10,:,:])
print np.shape(anom), np.shape(anom_seas1)
print anom_seas1[0:12,90,100]
anom_seas=np.reshape(anom_seas1,(3,(nyr-1)*4,720,1440),order='f')
print np.shape(anom_seas)
print anom_seas[:,0,90,100]
print np.nanmean(anom_seas[:,0,90,100])
anom_seas=np.nanmean(anom_seas, axis=0)
print np.shape(anom_seas)
print anom_seas[0,90,100]
djf1=np.reshape(anom_seas,(4,(nyr-1),720,1440),order='f')
djf=djf1[3,:,:,:]
jja=djf1[1,:,:,:]
print np.shape(djf)
plt.imshow(djf[0,:,:])
plt.show()
#=================================================
# make stdevs
#================================================
sst_ann_std=np.std(sst_ann_avg,axis=0)

anom_std=np.std(anom,axis=0)

djf_std=np.std(djf,axis=0)

jja_std=np.std(jja,axis=0)


print '---------------------'
print 'max sst ann stdev',np.nanmax(sst_ann_std)
print 'max anom stdev',np.nanmax(anom_std)
print 'max djf ann stdev',np.nanmax(djf_std)
print 'max jja ann stdev',np.nanmax(jja_std)
print '---------------------'
print 'min sst ann stdev',np.nanmin(sst_ann_std)
print 'min anom stdev',np.nanmin(anom_std)
print 'min djf ann stdev',np.nanmin(djf_std)
print 'min jja ann stdev',np.nanmin(jja_std)
print '----------------------'

print np.shape(sst_ann_std), np.shape(anom_std), np.shape(djf_std)
#stop
#======================================================
# plot
#=========================================================

colors=['#00ffff','#0033cc','#0066ff','#3399ff','#66ccff','#ccffff','#ffff00','#ffcc00','#ff9933','#ff6600','#ff3300','#ff0000']
xtv=[0,90*4,180*4,270*4,360*4]
xtl=['-180','-90','0','90','180']
ytv=[30*4,60*4,90*4,120*4,150*4]
ytl=['60','30','0','-30','-60']
'''
x=np.arange(0,360,1)
y=np.arange(0,180,1)
X,Y=np.meshgrid(x,y)
'''
fig=plt.figure(figsize=(11,8))
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
'''
left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 1.0   # the amount of height reserved for white space between subplots
'''
ax1=fig.add_subplot(321)
sst_avg=np.flipud(sst_avg) 
plt.imshow(sst_avg, interpolation='nearest',cmap=plt.get_cmap('Blues'))
plt.title('Time mean SST 1992-2015 \n from monthly input data',fontsize=8)
plt.colorbar(extend='both', shrink=0.8)
plt.xticks(xtv,xtl,fontsize=8)
plt.yticks(ytv,ytl,fontsize=8)
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.0)
'''
ax1=fig.add_subplot(322)
sst_ann_std=np.flipud(sst_ann_std)
clevs=[0.0,0.16,0.4,0.6,0.8,1.0,1.2,1.6,2]
colors=['#00ffff','#0033cc','#0066ff','#3399ff','#66ccff','#ccffff','#ffff00','#ffcc00','#ff9933','#ff6600','#ff3300','#ff0000']
vmin=0
vmax=2.4 
#plt.imshow(sst_ann_std,interpolation='nearest',vmin=vmin, vmax=vmax, cmap=cmap)
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='white')
# draw parallels and meridians.

m.drawparallels(np.arange(-90.,91.,30.),labels=[1,0,0,0],fontsize=10)
m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1],fontsize=10)
m.drawmapboundary(fill_color='white')
nlats = 180; nlons = 360; delta = 2.*np.pi/(nlons-1)
lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
lons = (-np.pi+delta*np.indices((nlats,nlons))[1,:,:])
x, y = m(lons*180./np.pi, lats*180./np.pi)

cs=plt.contourf(x,y,sst_ann_std,clevs,extend='both',colors=colors)
plt.title('STD of annual mean SST \n (n=18),1993-2010')
plt.colorbar(extend='both', ticks=clevs,shrink=0.8)

plt.xticks(xtv,xtl,fontsize=8)
plt.yticks(ytv,ytl,fontsize=8)

'''


ax1=fig.add_subplot(324)
clevs=np.arange(0,2.6,0.2)
vmin=0
vmax=2.4 
anom_std=np.flipud(anom_std) 
plt.imshow(anom_std, interpolation='nearest',vmin=vmin, vmax=vmax,cmap=plt.get_cmap('rainbow',12))
plt.title('STD of monthly SST anomaly \n (n=288), 1992-2015')
plt.colorbar(extend='both', ticks=clevs,shrink=0.8)
plt.xticks(xtv,xtl,fontsize=8)
plt.yticks(ytv,ytl,fontsize=8)
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.0)

ax1=fig.add_subplot(323)
sst_ann_std=np.flipud(sst_ann_std)
clevs=np.arange(0,2.6,0.2)
vmin=0
vmax=2.4 
plt.imshow(sst_ann_std,interpolation='nearest',vmin=vmin, vmax=vmax, cmap=plt.get_cmap('rainbow',10))
plt.title('STD of annual mean SST \n (n=24),1992-2015')
plt.colorbar(extend='both', ticks=clevs,shrink=0.8)
plt.xticks(xtv,xtl,fontsize=8)
plt.yticks(ytv,ytl,fontsize=8)
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.0)


ax1=fig.add_subplot(325)
clevs=np.arange(0,2.6,0.2)
djf_std=np.flipud(djf_std) 
plt.imshow(djf_std,interpolation='nearest', vmin=vmin, vmax=vmax,cmap=plt.get_cmap('rainbow',12))
plt.title('STD of annual DJF SST anomaly \n (n=23), Dec 1992- Feb 2015')
plt.colorbar(extend='both', ticks=clevs,shrink=0.8)
plt.xticks(xtv,xtl,fontsize=8)
plt.yticks(ytv,ytl,fontsize=8)
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.0)

ax1=fig.add_subplot(326)
jja_std=np.flipud(jja_std) 
plt.imshow(jja_std,interpolation='nearest', vmin=vmin, vmax=vmax,cmap=plt.get_cmap('rainbow',12))
plt.title('STD of annual JJA SST anomaly \n (n=23), Jun 1992- Aug 2015')
plt.colorbar(extend='both', ticks=clevs,shrink=0.8)
plt.xticks(xtv,xtl,fontsize=8)
plt.yticks(ytv,ytl,fontsize=8)
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.0)




plt.tight_layout()




plt.savefig('SSMi_SST_STD_for_Chris_Nov16.png')
plt.savefig('SSMi_SST_STD_for_Chris_Nov16.eps')
plt.show()


fig=plt.figure(figsize=(11,8))
ax1=fig.add_subplot(111)
#sst_ann_std=np.flipud(sst_ann_std)
clevs=[0.0,0.16,0.4,0.6,0.8,1.0,1.2,1.6,2]
colors=['#00ffff','#0033cc','#0066ff','#3399ff','#66ccff','#ccffff','#ffff00','#ffcc00','#ff9933','#ff6600','#ff3300','#ff0000']
colors=['#0033cc','#0066ff','#3399ff','#00ffff','#66ff66','#ffff00','#ff9933','#ff3300','#ff0000']

vmin=0
vmax=2.4 
#plt.imshow(sst_ann_std,interpolation='nearest',vmin=vmin, vmax=vmax, cmap=cmap)
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines() 
m.fillcontinents(color='white',lake_color='white')
# draw parallels and meridians.

m.drawparallels(np.arange(-90.,91.,30.),labels=[1,0,0,0],fontsize=10)
m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1],fontsize=10)
m.drawmapboundary(fill_color='white')
nlats = 180*4; nlons = 360*4; delta = 2.*np.pi/(nlons-1)
lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
lons = (-np.pi+delta*np.indices((nlats,nlons))[1,:,:])
x, y = m(lons*180./np.pi, lats*180./np.pi)

cs=plt.contourf(x,y,sst_ann_std,clevs,extend='both',colors=colors)
plt.title('STD of annual mean SST \n (n=24),1998-2015')
plt.colorbar(extend='both', ticks=clevs,shrink=0.8)

#plt.xticks(xtv,xtl,fontsize=8)
#plt.yticks(ytv,ytl,fontsize=8)


plt.savefig('SSMi_SST_STD_for_Chris_Nov16_newscale.png')
plt.savefig('SSMi_SST_STD_for_Chris_Nov16_newscale.eps')
plt.show()
