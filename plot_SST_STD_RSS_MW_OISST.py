#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Program to plot HC-TC trends.

import numpy as np
import netCDF4 as nCDF
import glob
import os.path
import matplotlib.pyplot as plt
########################################################################
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
#======================================================================
def calc_ndays(yearin, monthin):


def make_monthly(yearin, filesin):


def makefilelist(ystart, ystop, path=''):
    filelist = []
    for yy in range(ystart, ystop+1):
        if yy%4!=0:
            for doy in range(1,365):
                f = os.path.join(path,'{Y}/mw.fusion.{Y}.{d].v04.0.nc'.format(Y=yy, d=doy))
        elif yy%4==0:
            for doy in range(1,366):
                f = os.path.join(path,'{Y}/mw.fusion.{Y}.{d].v04.0.nc'.format(Y=yy, d=doy))
        filelist.extend(glob.glob(f))
    return filelist


#=======================================
# read SL CCI data ( to end 2014), reashape, remove annual cycle
#=============================================

year_start=2002
year_stop=2002
nmon = 12*(year_start - year_stop+1)
#year_stop=2013
sst=np.zeros((nmon, 720,1440))

dirin="/glusterfs/surft/data/RSS-MW-OISST-v4.0_nc/"
#NOTE!!!!!!!!!! DOESNT: ACCOUNT FOR LEAP YEARS AT THE OMENT!!!!!!!!!

for year in range(year_start, year_stop + 1):
    files = makefilelist(year, year, dirin)
    #read the year
    #make 12*monthly means
    make_monthly(year,files)
    #add to 

print files


