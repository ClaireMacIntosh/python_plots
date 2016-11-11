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
            for doy in range(1, 3):
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
year_stop = 1999


# dirin = "/glusterfs/surft/data/RSS-MW-OISST-v4.0_nc/"
dirin = "/Users/clairemacintosh/Data/SST_daily_mw_OISST_v04.0/"


for year in range(year_start, year_stop + 1):
    file_list = makefilelist(year, year, dirin)
    doy_list = [int(filename[-12:-9]) for filename in file_list]
    month_list = np.array([(datetime.datetime(year, 1, 1) + datetime.timedelta(days - 1)).month for days in doy_list])

    for month in xrange(1,2,1):
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



'''

