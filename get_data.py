import datetime as dt
import numpy as np
import pandas as pd
import glob
import seawater as sw
import warnings
from netCDF4 import Dataset
from scipy.io import loadmat # only needed for get_WW_temporary

from SB_support import readSB
from seawater_scattering import *

END_DATE = 252 # last day is Sept. 9, which is yearday 251 (assuming Jan 1 is 0)

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
        return baseline#,spikes
    else:
        #return baseline,spikes,depth_out
        return baseline,depth_out

def get_time(SB,refdate=dt.datetime(2018,1,1)):
    data = SB.data
    length = len(data['date'])
    years  = [int(str(data['date'][i])[0:4]) for i in range(length)]
    months = [int(str(data['date'][i])[4:6]) for i in range(length)]
    days   = [int(str(data['date'][i])[6:8]) for i in range(length)]
    hours  = [int(data['time'][i][0:2]) for i in range(length)]
    minutes= [int(data['time'][i][3:5]) for i in range(length)]
    seconds= [int(data['time'][i][6:8]) for i in range(length)]
    time = [dt.datetime(year=y,month=mo,day=d,hour=h,minute=mi,second=s)
            for y,mo,d,h,mi,s in zip(years,months,days,hours,minutes,seconds)]
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

def get_FLBB_SR(files, get_data=True, get_dives=None, remove_seawater=True, VSF_dark=0):
    SB_files = np.sort(glob.glob(files))
    if get_dives is not None:
        # first get all stations
        all_stations = np.empty(len(SB_files))
        for i in range(len(SB_files)):
            SB = readSB(SB_files[i],no_warn=True)
            # 'calibration_file' is more accurate than 'station' var
            all_stations[i] = int(SB.headers['calibration_files'][7:10])
        # now re-make SB_files list to be one for each station
        # (need to do it this way in case stations are duplicated)
        SB_files = np.array([SB_files[np.where(all_stations==station)[0][0]]
                             for station in get_dives])

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        # 'calibration_file' is more accurate than 'station' var
        station = int(SB.headers['calibration_files'][7:10])
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (station,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            if remove_seawater:
                #VSF_sw = betasw_HZP2019(700,142,T[-1],S[-1],sw.pres(z[-1],lat=50.5))[0].squeeze()
                VSF_sw = betasw_wetlabs(700,20,142,42,T[-1],S[-1],sw.pres(z[-1],lat=50.5))
            else:
                VSF_sw = 0
            VSF += ((np.array(SB.data['bbp'])/(2*np.pi*1.17) - VSF_dark - VSF_sw),)
            # convert using chi factor for 142deg angle used here (1.17)
            # (note this is different than that found in Zhang et al., 2021)

    t = np.array(t,dtype='object')
    lat = np.array(lat,dtype='object')
    lon = np.array(lon,dtype='object')
    sta = np.array(sta,dtype='object')
    if get_data:
        T = np.array(T,dtype='object')
        S = np.array(S,dtype='object')
        z = np.array(z,dtype='object')
        VSF = np.array(VSF,dtype='object')
    output = (t,lat,lon,sta)
    if get_data:
        output += (T,S,z,VSF)

    return output

def get_FLBB_RR(files, get_data=True, get_dives=None, remove_seawater=True, VSF_dark=0):
    SB_files = np.sort(glob.glob(files))
    if get_dives is not None:
        # first get all stations
        all_stations = np.empty(len(SB_files))
        for i in range(len(SB_files)):
            SB = readSB(SB_files[i],no_warn=True)
            # 'calibration_files' is more accurate than 'station' var
            all_stations[i] = int(SB.headers['calibration_files'][10:13])
        # now re-make SB_files list to be one for each station
        # (need to do it this way in case stations are duplicated)
        SB_files = np.array([SB_files[np.where(all_stations==station)[0][0]]
                             for station in get_dives])

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        if get_starttime(SB) < 233:
            continue # Roger Revelle was optically suspect at beginning of cruise
        # 'calibration_files' is more accurate than 'station' var
        station = int(SB.headers['calibration_files'][10:13])
        if (station==45) | (station==48):
            continue # bad stations

        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (station,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            if remove_seawater:
                #VSF_sw = betasw_HZP2019(700,142,T[-1],S[-1],sw.pres(z[-1],lat=50.5))[0].squeeze()
                VSF_sw = betasw_wetlabs(700,20,142,42,T[-1],S[-1],sw.pres(z[-1],lat=50.5))
            else:
                VSF_sw = 0
            VSF += (((np.array(SB.data['bbp'])/(2*np.pi*1.17) * 5.98e-6 / 7.13e-6) - VSF_dark - VSF_sw),)
            # convert using chi factor for 142deg angle used here (1.17)
            # (note this is different than that found in Zhang et al., 2021)
            # conversion from the calibration scale factor in the SeaBASS R2 dataset
            # see Supplemental Info

    t = np.array(t,dtype='object')
    lat = np.array(lat,dtype='object')
    lon = np.array(lon,dtype='object')
    sta = np.array(sta,dtype='object')
    if get_data:
        T = np.array(T,dtype='object')
        S = np.array(S,dtype='object')
        z = np.array(z,dtype='object')
        VSF = np.array(VSF,dtype='object')
    output = (t,lat,lon,sta)
    if get_data:
        output += (T,S,z,VSF)

    return output

def get_FLBB_LF(filename, get_data=True, get_dives=None, remove_seawater=True, VSF_dark=0, remove_deep_vals=True):
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
            VSF_sw = betasw_wetlabs(700,20,140,41,LF_temp,LF_sal,np.array(LF.data['pressure']))
            LF_VSF -= VSF_sw

    begin_dive = np.where((np.diff(LF_time)>.01) & (np.diff(LF_time)<.4))[0]
    end_dive = np.array([np.where((np.diff(LF_dep[i:i+1000])<-1)
                                  & (LF_dep[i+1:i+1000]>150))[0][0]+i
                         for i in begin_dive])

    if get_dives is None:
        inds = np.arange(len(begin_dive))
    else:
        inds = get_dives

    LF_lon_dive = np.array(LF_lon[begin_dive][inds])
    LF_lat_dive = np.array(LF_lat[begin_dive][inds])
    LF_time_dive = np.array(LF_time[begin_dive][inds])

    if get_data:
        LF_temp_dive = np.array([LF_temp[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])],dtype='object')
        LF_sal_dive = np.array([LF_sal[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])],dtype='object')
        LF_VSF_dive = np.array([LF_VSF[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])],dtype='object')
        LF_dep_dive = np.array([LF_dep[b:e] for b,e in zip(begin_dive[inds],end_dive[inds])],dtype='object')

        if remove_deep_vals:
            z_inds = [(d<100) for d in LF_dep_dive]
            LF_temp_dive = np.array([T[i] for T,i in zip(LF_temp_dive,z_inds)],dtype='object')
            LF_sal_dive = np.array([S[i] for S,i in zip(LF_sal_dive,z_inds)],dtype='object')
            LF_VSF_dive = np.array([V[i] for V,i in zip(LF_VSF_dive,z_inds)],dtype='object')
            LF_dep_dive = np.array([z[i] for z,i in zip(LF_dep_dive,z_inds)],dtype='object')

    bad_inds = (LF_time_dive>END_DATE) # beyond time of the main EXPORTS experiment
    LF_time_dive = LF_time_dive[~bad_inds]
    LF_lat_dive = LF_lat_dive[~bad_inds]
    LF_lon_dive = LF_lon_dive[~bad_inds]
    inds = inds[~bad_inds]

    if get_data:
        LF_temp_dive = LF_temp_dive[~bad_inds]
        LF_sal_dive = LF_sal_dive[~bad_inds]
        LF_dep_dive = LF_dep_dive[~bad_inds]
        LF_VSF_dive = LF_VSF_dive[~bad_inds]

    output = (LF_time_dive,LF_lat_dive,LF_lon_dive,inds,)

    if get_data:
                output += (LF_temp_dive,LF_sal_dive,LF_dep_dive,LF_VSF_dive,)


    return output

def get_WW(files, get_data=True, get_dives=None, remove_seawater=False, VSF_dark=0):
    SB_files = np.sort(glob.glob(files))

    if get_dives is not None:
        SB_files = SB_files[get_dives.astype('int32')]

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (i,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            if remove_seawater:
                VSF_sw = betasw_wetlabs(700,20,124,41,T[-1],S[-1],sw.pres(z[-1],lat=50.5))
            else:
                VSF_sw = 0
            VSF += (np.array(SB.data['bbp700'])/(np.pi*2*1.1)-VSF_sw-VSF_dark,)
            # convert using chi factor for 124deg angle used here

    output = (np.array(t,dtype='object'),np.array(lat,dtype='object'),np.array(lon,dtype='object'),np.array(sta,dtype='object'))
    if get_data:
        output += (np.array(T,dtype='object'),np.array(S,dtype='object'),np.array(z,dtype='object'),np.array(VSF,dtype='object'),)

    return output


def get_SG(file,get_data=True,get_dives=None,remove_seawater=False,wavelength=700):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)

        nc = Dataset(file,'r')

        inds = np.arange(len(nc.variables['time'])) if get_dives is None else get_dives

        rd = dt.datetime(1970,1,1)
        sd = dt.datetime(2018,1,1)

        t = np.nanmedian(nc.variables['time'][inds],axis=1)
        t = [rd+dt.timedelta(seconds=tt)-sd for tt in t]
        t = np.array([tt.days+tt.seconds/86400 for tt in t])
        lat = np.nanmedian(nc.variables['lat'][inds],axis=1)
        lon = np.nanmedian(nc.variables['lon'][inds],axis=1)

        if get_data:
            T = np.ma.getdata(nc.variables['T'][inds])
            S = np.ma.getdata(nc.variables['S'][inds])
            z = np.repeat(np.atleast_2d(np.ma.getdata(nc.variables['z'][:])),T.shape[0],axis=0)
            if remove_seawater:
                VSF_sw = betasw_wetlabs(700,20,124,41,T,S,sw.pres(z,lat=50.5))
            else:
                VSF_sw = 0
            VSF = np.ma.getdata(nc.variables['baseline_%03dnm'%wavelength][inds])-VSF_sw

        bad_inds = (t>END_DATE) # beyond time of the main EXPORTS experiment
        t = t[~bad_inds]
        lat = lat[~bad_inds]
        lon = lon[~bad_inds]
        inds = inds[~bad_inds]

        if get_data:
            T = T[~bad_inds]
            S = S[~bad_inds]
            VSF = VSF[~bad_inds]
            z = z[~bad_inds]


    if get_data:
        return t,lat,lon,inds,T,S,z,VSF
    else:
        return t,lat,lon,inds

def get_HS(files,get_data=True,get_dives=None):
    SB_files = np.sort(glob.glob(files))

    if get_dives is not None:
        SB_files = SB_files[get_dives.astype('int')]

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = (); VSF_sd = ();
    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)
        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (i,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            p = np.array(SB.data['pressure'])
            z += (sw.dpth(p,lat=50.5),)
            # convert using chi factor for 142deg angle used here
            VSF += (np.array(SB.data['bbp700'])/(2*np.pi*1.18),)
            # convert using chi factor for 142deg angle used here
            VSF_sd += (np.array(SB.data['bbp700_sd'])/(2*np.pi*1.18),)

    output = (np.array(t,dtype='object'),np.array(lat,dtype='object'),np.array(lon,dtype='object'),np.array(sta,dtype='object'))

    if get_data:
        output += (np.array(T,dtype='object'),np.array(S,dtype='object'),np.array(z,dtype='object'),np.array(VSF,dtype='object'),np.array(VSF_sd,dtype='object'),)

    return output

def get_HS_day(files,get_data=True,get_dives=None):
    output = get_HS(files,get_data=get_data,get_dives=get_dives)

    t = output[0]
    good_inds = (t%1 < 0.1) | (t%1 > 0.8)
    output = [o[good_inds] for o in output]

    return output

def get_HS_night(files,get_data=True,get_dives=None):
    output = get_HS(files,get_data=get_data,get_dives=get_dives)

    t = output[0]
    good_inds = (t%1 > 0.2) & (t%1 < 0.6)
    output = [o[good_inds] for o in output]

    return output

def get_BB9_RR(files,get_data=True,get_dives=None,wavelength=700):
    SB_files = np.sort(glob.glob(files))
    sta = np.arange(len(SB_files))

    if get_dives is not None:
        SB_files = SB_files[get_dives.astype('int')]
        sta = sta[get_dives.astype('int')]

    lat = (); lon = (); t = ();
    if get_data:
        T = (); S = (); z = (); VSF = (); VSF_sd = ();

    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)

        if 'bbp650' not in SB.data.keys():
            # sometimes this isn't recorded for some reason...do not use these
            continue


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
            VSF_sd += (((wvs[1]-wavelength)*VSF1_sd +(wavelength-wvs[0])*VSF2_sd)/(wvs[1]-wvs[0]),)

    output = np.array(t,dtype='object'),np.array(lat,dtype='object'),np.array(lon,dtype='object'),sta
    if get_data:
        output += (np.array(T,dtype='object'),np.array(S,dtype='object'),
                   np.array(z,dtype='object'),np.array(VSF,dtype='object'),np.array(VSF_sd,dtype='object'),)

    return output

def get_BB9_SR(files,get_data=True,get_dives=None,wavelength=700):
    SB_files = np.sort(glob.glob(files))

    if get_dives is not None:
        SB_files = SB_files[get_dives.astype('int')]

    lat = (); lon = (); t = (); sta = ();
    if get_data:
        T = (); S = (); z = (); VSF = (); VSF_sd = ();

    for i in range(len(SB_files)):
        SB = readSB(SB_files[i],no_warn=True)

        if np.max(np.array(SB.data['depth']))<50: # do not use profiles shallower than 50 m
            continue

        lat += (get_startlat(SB),)
        lon += (get_startlon(SB),)
        t += (get_starttime(SB),)
        sta += (i,)

        if get_data:
            T += (np.array(SB.data['wt']),)
            S += (np.array(SB.data['sal']),)
            z += (np.array(SB.data['depth']),)
            if wavelength==470:
                wvs = [441,488]
            elif wavelength==700:
                wvs = [652,717]
            else:
                print('Error: cannot have wavelength other than 470 or 700 right now')
            VSF1 = np.array(SB.data['bbp%03d'%wvs[0]])/(2*np.pi*1.0772)
            VSF1_sd = np.array(SB.data['bbp%03d_sd'%wvs[0]])/(2*np.pi*1.0772)
            VSF2 = np.array(SB.data['bbp%03d'%wvs[1]])/(2*np.pi*1.0772)
            VSF2_sd = np.array(SB.data['bbp%03d_sd'%wvs[1]])/(2*np.pi*1.0772)
            VSF += (((wvs[1]-wavelength)*VSF1 + (wavelength-wvs[0])*VSF2)/(wvs[1]-wvs[0]),)
            VSF_sd += (((wvs[1]-wavelength)*VSF1_sd +
                        (wavelength-wvs[0])*VSF2_sd)/(wvs[1]-wvs[0]),)

    output = (np.array(t,dtype='object'),np.array(lat,dtype='object'),np.array(lon,dtype='object'),np.array(sta,dtype='object'))
    if get_data:
        output += (np.array(T,dtype='object'),np.array(S,dtype='object'),
                   np.array(z,dtype='object'),np.array(VSF,dtype='object'),np.array(VSF_sd,dtype='object'),)

    return output

def get_MCOMS_Argo(filename,get_data=True,get_dives=None):
    Argo_nc = Dataset(filename,'r')

    inds = np.arange(Argo_nc.dimensions['N_PROF'].size) if get_dives is None else get_dives

    inds = inds.astype('int')

    tt = [dt.datetime(1950,1,1)+dt.timedelta(days=t) - dt.datetime(2018,1,1)
          for t in Argo_nc.variables['JULD'][inds].filled(np.nan)]
    time = np.array([t.days + t.seconds/86400 for t in tt])

    # doesn't matter because no matches anyway
    #bad_inds = time>END_DATE
    #inds = inds[~bad_inds]
    #time = time[~bad_inds]

    lat = Argo_nc.variables['LATITUDE'][inds].filled(np.nan)
    lon = Argo_nc.variables['LONGITUDE'][inds].filled(np.nan)

    if get_data:
        depth = sw.dpth(Argo_nc.variables['PRES_ADJUSTED'][inds].filled(np.nan),lat=lat[0])
        temp = Argo_nc.variables['TEMP_ADJUSTED'][inds].filled(np.nan)
        sal = Argo_nc.variables['PSAL_ADJUSTED'][inds].filled(np.nan)
        VSF = Argo_nc.variables['BBP700_ADJUSTED'][inds].filled(np.nan)/(1.142*2*np.pi) # was changed from dividing by 1.1, need to re-run

    output = (time,lat,lon,inds,)
    if get_data:
        output += (temp,sal,depth,VSF,)

    return output
