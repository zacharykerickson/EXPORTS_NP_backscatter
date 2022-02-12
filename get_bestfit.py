import numpy as np
from scipy import odr
import seawater as sw
import warnings

def find_outliers(beta,x,y,dx,dy,inds=None):
    if len(beta)==1:
        beta = [1,beta[0]]
    xhat = (dy**2*x + dx**2*beta[0]*(y-beta[1]))/(dy**2+(dx*beta[0])**2)
    error = np.sqrt((xhat**2-2*xhat*x + x**2)/dx**2 + ((beta[0]*xhat)**2+2*beta[0]*xhat*(beta[1]-y)+(beta[1]-y)**2)/dy**2)
    error[x>xhat] *= -1
    if inds is None:
        outliers = np.abs(error/np.std(error))>3
    else:
        outliers = np.abs(error/np.std(error[inds]))>3
    return outliers

def offset_func(beta,x):
    return x+beta

def offset_model():
    return odr.Model(offset_func)

def find_nearby(times,lats,lons,TIME_THRESHOLD,DIST_THRESHOLD,
    find_closest=True,DDIST_DTIME_CONVERSION=2,return_distances=False,eliminate_duplicates=True):
    inds = ()
    for t,lat,lon in zip(times[0],lats[0],lons[0]):
        distances = np.array([sw.dist([lat,llat],[lon,llon])[0][0]
                              for llat,llon in zip(lats[1],lons[1])])
        dtimes = np.abs(times[1]-t)
        inds += (np.where((dtimes<=TIME_THRESHOLD) & (distances<=DIST_THRESHOLD))[0],)

    if find_closest:
        closest_ind = np.empty(shape=len(times[0]))
        for i in range(len(closest_ind)):
            if len(inds[i])==0:
                closest_ind[i] = np.nan
            elif len(inds[i])==1:
                closest_ind[i] = inds[i][0]
            else:
                ddist = np.array([sw.dist([lats[0][i],llat],[lons[0][i],llon])[0][0]
                                  for llat,llon in zip(lats[1][inds[i]],lons[1][inds[i]])])
                dtime = np.abs(times[0][i]-times[1][inds[i]])
                closest_ind[i] = inds[i][np.argmin(dtime+ddist*DDIST_DTIME_CONVERSION)]
        inds1 = np.arange(len(times[0]))
        inds1 = inds1[np.isfinite(closest_ind)]
        closest_ind = closest_ind[np.isfinite(closest_ind)].astype('int')
        ddist = np.array([sw.dist([lats[0][i],lats[1][j]],[lons[0][i],lons[1][j]])[0][0]
                          for i,j in zip(inds1,closest_ind)])
        dtime = np.abs(times[0][inds1]-times[1][closest_ind])
        # make sure no duplicate closest_inds
        if eliminate_duplicates:
            for i in np.unique(closest_ind):
                ii = np.where(closest_ind==i)[0]
                if len(ii)>1:
                    ind_to_keep = np.argmin(dtime[ii]+ddist[ii]*2)
                    ii = np.delete(ii,ind_to_keep)
                    closest_ind[ii] = -999
                    good_inds = closest_ind!=-999
                    inds1 = inds1[good_inds]
                    closest_ind = closest_ind[good_inds]
                    ddist = ddist[good_inds]
                    dtime = dtime[good_inds]
        dtime = dtime.astype('float64') # for some reason this is otherwise an 'object' array

        if return_distances:
            return inds1,closest_ind,ddist,dtime
        return inds1,closest_ind
    else:
        if return_distances:
            print('Sorry, cannot return distances without using setting find_closest to True')
        return inds

def interpolate_data(data,data_dep,deps,
                     return_uncert=False,ddep=None,minimum_uncert=0,average_type='mean',return_n=False):
    if ddep is None and return_uncert:
        print('Cannot return uncertainty without ddep set')
        return

    if average_type=='mean':
        avefunc = np.nanmean
    elif average_type=='median':
        avefunc = np.nanmedian

    if ddep is None: # just do a linear interpolation
        interp_data = np.array([np.interp(z,dat_z[np.isfinite(dat*dat_z)],
                                          dat[np.isfinite(dat*dat_z)],left=np.nan,right=np.nan)
                                for z,dat_z,dat in zip(deps,data_dep,data)])
        if return_n:
            print('Cannot return number of data points without ddep set')
        return interp_data
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            interp_data = np.array([[avefunc(dat[np.abs(dat_z-z)<ddep]) for z in zz]
                                for zz,dat_z,dat in zip(deps,data_dep,data)])
            if return_uncert:
                interp_uncert = np.array([[np.nanstd(dat[np.abs(dat_z-z)<ddep]) for z in zz]
                                      for zz,dat_z,dat in zip(deps,data_dep,data)])
                n = np.array([[np.sum((np.abs(dat_z-z)<ddep) & np.isfinite(dat)) for z in zz] for zz,dat_z,dat in zip(deps,data_dep,data)])

                interp_uncert = np.sqrt(interp_uncert**2 + minimum_uncert**2/n)

                if return_n:
                    return interp_data,interp_uncert,n
                else:
                    return interp_data,interp_uncert
            else:
                return interp_data
