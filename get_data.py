import datetime as dt
import numpy as np
import pandas as pd
import glob
import seawater as sw
from netCDF4 import Dataset
from scipy.io import loadmat # only needed for get_WW_temporary

import sys
sys.path.append('/Users/zerickso/PythonCode')
from SB_support import readSB
from seawater_scattering import *

def min_max_filter(data,depth=None,fillen=7):
    # filter to separate out baseline and spikes using sequential running minimum and maximum filters
    # only use non-nan values
    # as in Briggs et al., DSRI, 2011 doi:10.1016/j.dsr.2011.07.007

    bad_vals = np.isnan(data)
    data = data[~bad_vals]
    if depth is not None:
        depth = depth[~bad_vals]

    window = np.floor(fillen/2).astype('int')
    tmp = np.array([np.min(data[i-window:i+window+1]) for i in range(window,len(data)-window)])
    baseline = np.array([np.max(tmp[i-window:i+window+1]) for i in range(window,len(tmp)-window)])
    spikes = data[window*2:-window*2]-baseline
    if depth is not None:
        depth_out = depth[window*2:-window*2]

    if depth is None:
        return baseline,spikes
    else:
        return baseline,spikes,depth_out

def median_filter(data,depth=None,fillen=7):
    bad_vals = np.isnan(data)
    data = data[~bad_vals]
    if depth is not None:
        depth = depth[~bad_vals]

    window = np.floor(fillen/2).astype('int')
    baseline = np.array([np.median(data[i:i+fillen]) for i in range(len(data)-fillen)])
    spikes = data[window:-window-1]-baseline
    if depth is not None:
        depth_out = depth[window:-window-1]

    if depth is None:
        return baseline,spikes
    else:
        return baseline,spikes,depth_out



def get_time(SB,refdate=dt.datetime(2018,1,1)):
    data = SB.data
    length = len(data['date'])
    years  = [int(str(data['date'][i])[0:4]) for i in range(length)]
    months = [int(str(data['date'][i])[4:6]) for i in range(length)]
    days   = [int(str(data['date'][i])[6:8]) for i in range(length)]
    hours  = [int(data['time'][i][0:2]) for i in range(length)]
    minutes= [int(data['time'][i][3:5]) for i in range(length)]
    seconds= [int(data['time'][i][6:8]) for i in range(length)]
    time = [dt.datetime(year=y,month=mo,day=d,hour=h,minute=mi,second=s) for y,mo,d,h,mi,s in zip(years,months,days,hours,minutes,seconds)]
    dtime = [t-refdate for t in time] # did this get put in?
    return np.array([t.days + t.seconds/86400 for t in dtime])

def get_starttime(SB,refdate=dt.datetime(2018,1,1)):
    year  = int(SB.headers['start_date'][0:4])
    month = int(SB.headers['start_date'][4:6])
    day   = int(SB.headers['start_date'][6:8])
    hour  = int(SB.headers['start_time'][0:2])
    minute= int(SB.headers['start_time'][3:5])
    second= int(SB.headers['start_time'][6:8])
    time = dt.datetime(year=year,month=month,day=day,hour=hour,minute=minute,second=second)
    dtime = time-refdate
    return dtime.days+dtime.seconds/86400

def get_startlat(SB):
    return float(SB.headers['north_latitude'][:-5])
def get_startlon(SB):
    return float(SB.headers['east_longitude'][:-5])

def get_FLBB_SR(folder,get_data=True,get_dives=None,verbose=False,remove_seawater=True):
    SB_files = np.sort(glob.glob('%s/EXPORTS_EXPORTSNP_CTDbinned_rosette_survey_*.sb'%folder))
    if get_dives is not None:
        # first get all stations
        all_stations = np.empty(len(SB_files))
        for i in range(len(SB_files)):
            SB = readSB(SB_files[i],no_warn=True)
            all_stations[i] = int(SB.headers['calibration_files'][7:10]) # this is more accurate than the 'station' variable
        # now re-make SB_files list to be one for each station (need to do it this way in case stations are duplicated)
        SB_files = np.array([SB_files[np.where(all_stations==station)[0][0]] for station in get_dives])

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        station = int(SB.headers['calibration_files'][7:10]) # this is more accurate than the 'station' variable
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (station,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            VSF_sw = betasw_HZP2019(700,142,T[-1],S[-1],sw.pres(z[-1],lat=50.5))[0].squeeze() if remove_seawater else 0 # arguments are wv, ang, T, S, p
            # VSF_dk = (np.array(SB.data['bbp'])-2.7788e-5)/(2*np.pi*1.17) # this is the median value from station 92, which was the dark cast
            val = np.array(SB.data['bbp'])/(2*np.pi*1.17)
            # convert using chi factor for 142deg angle used here (note this is different than that found in Zhang et al., 2021)
            VSF += (val-VSF_sw,)

    if get_data:
        if verbose:
            print('Output is time, lat, lon, station, T, S, z, VSF')
        return np.array(t),np.array(lat),np.array(lon),np.array(sta),np.array(T),np.array(S),np.array(z),np.array(VSF)
    else:
        if verbose:
            print('Output is time, lat, lon, station')
        return np.array(t),np.array(lat),np.array(lon),np.array(sta)

def get_FLBB_RR(folder,get_data=True,get_dives=None,verbose=False,remove_seawater=True):
    SB_files = np.sort(glob.glob('%s/EXPORTS_EXPORTSNP_CTDbinned_rosette_process_*.sb'%folder))
    if get_dives is not None:
        # first get all stations
        all_stations = np.empty(len(SB_files))
        for i in range(len(SB_files)):
            SB = readSB(SB_files[i],no_warn=True)
            all_stations[i] = int(SB.headers['calibration_files'][10:13]) # this is way more accurate than the 'station' variable
        # now re-make SB_files list to be one for each station (need to do it this way in case stations are duplicated)
        SB_files = np.array([SB_files[np.where(all_stations==station)[0][0]] for station in get_dives])

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        station = int(SB.headers['calibration_files'][10:13]) # this is more accurate than the 'station' variable
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (station,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            VSF_sw = betasw_HZP2019(700,142,T[-1],S[-1],sw.pres(z[-1],lat=50.5))[0].squeeze() if remove_seawater else 0 # arguments are wv, ang, T, S, p
            # VSF_dk = (np.array(SB.data['bbp'])-2.16497e-5)/(2*np.pi*1.17) # this is the median value from station 48, which was the dark cast
            val = np.array(SB.data['bbp'])/(2*np.pi*1.17)
            val = val * 5.98e-6 / 7.13e-6 # this is a conversion from the calibration scale factor that is in the SeaBASS R2 dataset - see Supplemental Info
            # convert using chi factor for 142deg angle used here (note this is different than that found in Zhang et al., 2021)
            VSF += (val-VSF_sw,)

    if get_data:
        if verbose:
            print('Output is time, lat, lon, station, T, S, z, VSF')
        return np.array(t),np.array(lat),np.array(lon),np.array(sta),np.array(T),np.array(S),np.array(z),np.array(VSF)
    else:
        if verbose:
            print('Output is time, lat, lon, station')
        return np.array(t),np.array(lat),np.array(lon),np.array(sta)

def get_FLBB_LF(filename,get_data=True,get_dives=None,verbose=False,remove_seawater=True):
    LF = readSB(filename,no_warn=True)
    LF_time = get_time(LF)
    LF_lon = np.array(LF.data['lon'])
    LF_lat = np.array(LF.data['lat'])

    LF_dep = np.array(sw.dpth(np.array(LF.data['pressure']),lat=50.5))

    if get_data:
        LF_temp = np.array(LF.data['wt'])
        LF_sal = np.array(LF.data['sal'])
        LF_VSF = np.array(LF.data['vsf700'])
        if remove_seawater:
            VSF_ss = betasw_wetlabs(700,20,140,41,LF_temp,LF_sal,np.array(LF.data['pressure']))
            LF_VSF -= VSF_ss

    begin_dive = np.where((np.diff(LF_time)>.01) & (np.diff(LF_time)<.4))[0]
    end_dive = np.array([np.where((np.diff(LF_dep[i:i+1000])<-1) & (LF_dep[i+1:i+1000]>150))[0][0]+i for i in begin_dive])

    if get_dives is None:
        inds = range(len(begin_dive))
    else:
        inds = get_dives

    LF_lon_dive = np.array(LF_lon[begin_dive][inds])
    LF_lat_dive = np.array(LF_lat[begin_dive][inds])
    LF_time_dive = np.array(LF_time[begin_dive][inds])

    if get_data:
        LF_temp_dive = np.array([LF_temp[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])])
        LF_sal_dive = np.array([LF_sal[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])])
        LF_VSF_dive = np.array([LF_VSF[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])])
        LF_dep_dive = np.array([LF_dep[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])])

        if verbose:
            print('Output is time, lat, lon, T, S, z, VSF')
        return LF_time_dive,LF_lat_dive,LF_lon_dive,LF_temp_dive,LF_sal_dive,LF_dep_dive,LF_VSF_dive

    else:
        if verbose:
            print('Output is time, lat, lon')
        return LF_time_dive,LF_lat_dive,LF_lon_dive

def get_WW(folder,get_data=True,get_dives=None,verbose=False,remove_seawater=False):
    SB_files = np.sort(glob.glob('%s/*.sb'%folder))[1::] # take out the first one, which is just floating at the surface

    if get_dives is not None:
        SB_files = SB_files[get_dives]

    lat = (); lon = (); t = ();
    if get_data:
        T = (); S = (); z = (); VSF = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            VSF_sw = betasw_HZP2019(700,124,T[-1],S[-1],sw.pres(z[-1],lat=50.5))[0].squeeze() if remove_seawater else 0 # arguments are wv, ang, T, S, p
            VSF += (np.array(SB.data['bbp700'])/(np.pi*2*1.02)-VSF_sw,) # convert using chi factor for 124deg angle used here

    if get_data:
        if verbose:
            print('Output is time, lat, lon, T, S, z, VSF')
        return np.array(t),np.array(lat),np.array(lon),np.array(T),np.array(S),np.array(z),np.array(VSF)
    else:
        if verbose:
            print('Output is time, lat, lon')
        return np.array(t),np.array(lat),np.array(lon)

def get_WW_temporary(folder,filename,get_data=True,get_dives=None,verbose=False,remove_seawater=False):
    # first get WW from SeaBASS
    if get_data:
        t,lat,lon,T,S,z,VSF = get_WW(folder,get_data=get_data,get_dives=get_dives,verbose=verbose,remove_seawater=remove_seawater)
    else:
        t,lat,lon = get_WW(folder,get_data=get_data,get_dives=get_dives,verbose=verbose,remove_seawater=remove_seawater)
    # then get data from the .mat file
    WW_matfile_data = loadmat(filename)
    WW_mat_lat = np.hstack((WW_matfile_data['epoch2'][0,0]['lat'][:,0],WW_matfile_data['epoch3'][0,0]['lat'][:,0]))
    iii = np.isfinite(WW_mat_lat)
    WW_mat_lat = WW_mat_lat[iii]
    if get_dives is None:
        get_dives = range(len(WW_mat_lat))
    WW_mat_lat = WW_mat_lat[get_dives]
    if verbose:
        print('Check if locations match?',np.array_equal(np.round(WW_mat_lat,4),lat))
    if get_data:
        WW_mat_temp = [T[0:len(zz)] for zz,T in zip(z,np.vstack((WW_matfile_data['epoch2'][0,0]['Tempave'],WW_matfile_data['epoch3'][0,0]['Tempave']))[iii][get_dives])]
        WW_mat_sal = [S[0:len(zz)] for zz,S in zip(z,np.vstack((WW_matfile_data['epoch2'][0,0]['Salave'],WW_matfile_data['epoch3'][0,0]['Salave']))[iii][get_dives])]

    if get_data:
        if verbose:
            print('Output is time, lat, lon, T, S, z, VSF')
        return t,lat,lon,np.array(WW_mat_temp),np.array(WW_mat_sal),z,VSF
    else:
        if verbose:
            print('Output is time, lat, lon')
        return t,lat,lon

def get_SG(file,wavelength=None,get_data=True,get_dives=None,verbose=False,remove_seawater=False):
    nc = Dataset(file,'r')

    inds = range(len(nc.variables['time'])) if get_dives is None else get_dives

    rd = dt.datetime(1970,1,1)
    sd = dt.datetime(2018,1,1)

    t = np.nanmedian(nc.variables['time'][inds],axis=1)
    t = [rd+dt.timedelta(seconds=tt)-sd for tt in t]
    t = np.array([tt.days+tt.seconds/86400 for tt in t])
    lat = np.nanmedian(nc.variables['lat'][inds],axis=1)
    lon = np.nanmedian(nc.variables['lon'][inds],axis=1)

    if get_data:
        T = nc.variables['T'][inds]
        S = nc.variables['S'][inds]
        z = np.repeat(np.atleast_2d(nc.variables['z'][:]),T.shape[0],axis=0)
        VSF_sw = betasw_HZP2019(700,124,T,S,sw.pres(z,lat=50.5))[0].squeeze() if remove_seawater else 0 # arguments are wv, ang, T, S, p
        VSF = nc.variables['baseline_%03dnm'%wavelength][inds]-VSF_sw

    if get_data:
        if verbose:
            print('Output is time, lat, lon, T, S, z, VSF')
        return t,lat,lon,T,S,z,VSF
    else:
        if verbose:
            print('Output is time, lat, lon')
        return t,lat,lon

def get_HS(folder,get_data=True,get_dives=None,verbose=False):
    SB_files = np.sort(glob.glob('%s/*.sb'%folder))

    if get_dives is not None:
        SB_files = SB_files[get_dives]

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = (); VSF_sd = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            p = np.array(SB.data['pressure'])
            z += (sw.dpth(p,lat=50.5),)
            VSF += (np.array(SB.data['bbp700'])/(2*np.pi*1.18),) # convert using chi factor for 142deg angle used here
            VSF_sd += (np.array(SB.data['bbp700_sd'])/(2*np.pi*1.18),) # convert using chi factor for 142deg angle used here

    if get_data:
        if verbose:
            print('Output is time, lat, lon, T, S, z, VSF, VSF_sd')
        return np.array(t),np.array(lat),np.array(lon),np.array(T),np.array(S),np.array(z),np.array(VSF),np.array(VSF_sd)
    else:
        if verbose:
            print('Output is time, lat, lon')
        return np.array(t),np.array(lat),np.array(lon)

def get_BB9_RR(folder,wavelength=None,get_data=True,get_dives=None,verbose=False):
    SB_files = np.sort(glob.glob('%s/*.sb'%folder))

    if get_dives is not None:
        SB_files = SB_files[get_dives]

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = (); VSF_sd = ();

    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            p = np.array(SB.data['pressure'])
            z += (sw.dpth(p,lat=50.5),)
            if wavelength==470:
                wvs = [440,488]
            elif wavelength==700:
                wvs = [650,715]
            else:
                print('Error: cannot have wavelength other than 470 or 700 right now')
            VSF1 = np.array(SB.data['bbp%03d'%wvs[0]])/(2*np.pi*1.076)
            VSF1_sd = np.array(SB.data['bbp%03d_sd'%wvs[0]])/(2*np.pi*1.076)
            VSF2 = np.array(SB.data['bbp%03d'%wvs[1]])/(2*np.pi*1.076)
            VSF2_sd = np.array(SB.data['bbp%03d_sd'%wvs[1]])/(2*np.pi*1.076)
            VSF += (((wvs[1]-wavelength)*VSF1 + (wavelength-wvs[0])*VSF2)/(wvs[1]-wvs[0]),)
            VSF_sd += (((wvs[1]-wavelength)*VSF1_sd + (wavelength-wvs[0])*VSF2_sd)/(wvs[1]-wvs[0]),)

    if get_data:
        if verbose:
            print('Output is time, lat, lon, T, S, z, VSF, VSF_sd')
        return np.array(t),np.array(lat),np.array(lon),np.array(T),np.array(S),np.array(z),np.array(VSF),np.array(VSF_sd)
    else:
        if verbose:
            print('Output is time, lat, lon')
        return np.array(t),np.array(lat),np.array(lon)
