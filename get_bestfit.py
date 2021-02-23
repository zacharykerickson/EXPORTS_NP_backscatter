import numpy as np
from scipy.optimize import curve_fit, minimize
import seawater as sw

def find_nearby(times,lats,lons,TIME_THRESHOLD,DIST_THRESHOLD,find_closest=True,DDIST_DTIME_CONVERSION=2,return_distances=False,eliminate_duplicates=True):
    inds = ()
    for t,lat,lon in zip(times[0],lats[0],lons[0]):
        distances = np.array([sw.dist([lat,llat],[lon,llon])[0][0] for llat,llon in zip(lats[1],lons[1])])
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
                ddist = np.array([sw.dist([lats[0][i],llat],[lons[0][i],llon])[0][0] for llat,llon in zip(lats[1][inds[i]],lons[1][inds[i]])])
                dtime = np.abs(times[0][i]-times[1][inds[i]])
                closest_ind[i] = inds[i][np.argmin(dtime+ddist*DDIST_DTIME_CONVERSION)]
        inds1 = np.arange(len(times[0]))
        inds1 = inds1[np.isfinite(closest_ind)]
        closest_ind = closest_ind[np.isfinite(closest_ind)].astype('int')
        ddist = np.array([sw.dist([lats[0][i],lats[1][j]],[lons[0][i],lons[1][j]])[0][0] for i,j in zip(inds1,closest_ind)])
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


        if return_distances:
            return inds1,closest_ind,ddist,dtime
        return inds1,closest_ind
    else:
        if return_distances:
            print('Sorry, cannot return distances right now without using setting find_closest to True')
        return inds

def interpolate_data(data,data_dep,deps,return_uncert=False,ddep=None,minimum_uncert=0,average_type='mean'):
    if ddep is None and return_uncert:
        print('Cannot return uncertainty without ddep set')
        return

    if average_type=='mean':
        avefunc = np.nanmean
    elif average_type=='median':
        avefunc = np.nanmedian

    if ddep is None: # just to a linear interpolation
        interp_data = np.array([np.interp(z,dat_z[np.isfinite(dat*dat_z)],dat[np.isfinite(dat*dat_z)],left=np.nan,right=np.nan) for z,dat_z,dat in zip(deps,data_dep,data)])
        return interp_data
    else:
        #interp_data = np.array([np.nanmean(dat[np.abs(dat_z-z)<ddep]) for z,dat_z,dat in zip(deps,data_dep,data)])
        interp_data = np.array([[avefunc(dat[np.abs(dat_z-z)<ddep]) for z in zz] for zz,dat_z,dat in zip(deps,data_dep,data)])
        if return_uncert:
            #interp_uncert = np.array([np.nanstd(dat[np.abs(dat_z-z)<ddep]) for z,dat_z,dat in zip(deps,data_dep,data)])
            interp_uncert = np.array([[np.nanstd(dat[np.abs(dat_z-z)<ddep]) for z in zz] for zz,dat_z,dat in zip(deps,data_dep,data)])
            interp_uncert[interp_uncert<minimum_uncert] = minimum_uncert

            return interp_data,interp_uncert
        else:
            return interp_data

def interpolate_data_depthdens(data,data_dep,data_dens,deps,dens,DEP_THRESHOLD,return_uncert=False,ddep=None,ddens=None,minimum_uncert=0):
    interp_data_depth = interpolate_data(data,data_dep,deps,return_uncert=return_uncert,ddep=ddep,minimum_uncert=minimum_uncert)
    interp_data_dens = interpolate_data(data,data_dens,dens,return_uncert=return_uncert,ddep=ddens,minimum_uncert=minimum_uncert)

    interp_data = ()
    if return_uncert:
        interp_uncert = ()
        for i in range(len(data)):
            inds = deps[i]<=DEP_THRESHOLD
            d = interp_data_dens[0][i]
            d[inds] = interp_data_depth[0][i][inds]
            u = interp_data_dens[1][i]
            u[inds] = interp_data_depth[1][i][inds]
            interp_data += (d,)
            interp_uncert += (u,)
        return interp_data,interp_uncert
    else:
        for i in range(len(data)):
            inds = deps[i]<=DEP_THRESHOLD
            d = interp_data_dens[i]
            d[inds] = interp_data_depth[i][inds]
            interp_data += (d,)
        return interp_data


def fit_data(params,x,y,dx,dy,verbose=False):
    if np.isscalar(dx):
        dx = dx*np.ones(len(x))
    if np.isscalar(dy):
        dy = dy*np.ones(len(y))

    if len(params)==1:
        func = ssr2
        init_cond = [0]
    elif len(params)==2:
        func = ssr
        init_cond = [1,0]

    fit_params = minimize(func,init_cond,args=(x,y,dx,dy),method='BFGS')
    error = func(fit_params.x,x,y,dx,dy)

    if verbose:
        print('Output is: best fit parameters, st. dev., squared error')
    return fit_params.x,np.sqrt(np.diag(fit_params.hess_inv)),error**2
