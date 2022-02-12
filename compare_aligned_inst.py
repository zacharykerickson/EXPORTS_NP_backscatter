## Import standard packages
import sys
from os import path
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw
import pandas as pd
from matplotlib import gridspec
from scipy import stats, odr

## Import local python files
import config as cf # For all of the configuration info
from seawater_scattering import betasw_HZP2019 # For seawater scattering
import get_data # For collecting instrument data
import get_bestfit # For data analysis and fitting procedures

if cf.SAVE_INFO and not path.exists(cf.SAVE_ALIGNED_INFO_FN):
    f = open(cf.SAVE_ALIGNED_INFO_FN,'a')
    f.write('inst1, inst2, num_profiles, n, n_rem_TS, n_rem_out, rsq, line_a, line_a_sd, line_a_95CI, line_b, line_b_sd, line_b_95CI, offset_b, offset_b_sd, offset_b_95CI\n')
    f.close()

## Get data for each of the instruments
inst1   = sys.argv[1]
label1  = cf.label[inst1]
sca1    = cf.sca[inst1]
chi1    = cf.chi[inst1]
unc_chi1= cf.unc_chi[inst1]
res_vsf1= cf.res_vsf[inst1]
unc_sf1 = cf.unc_sf[inst1]
unc_dk1 = cf.unc_dk[inst1]
files1  = cf.data_files[inst1]
func1   = cf.data_function[inst1]

inst2   = sys.argv[2]
label2  = cf.label[inst2]
sca2    = cf.sca[inst2]
chi2    = cf.chi[inst2]
unc_chi2= cf.unc_chi[inst2]
res_vsf2= cf.res_vsf[inst2]
unc_sf2 = cf.unc_sf[inst2]
unc_dk2 = cf.unc_dk[inst2]
files2  = cf.data_files[inst2]
func2   = cf.data_function[inst2]

if sca1==sca2:
    unc_chi1 = 0
    unc_chi2 = 0

figure_name = 'aligned_bbp_'+inst1+'_'+inst2+'.png'
yy = [np.min([cf.max_dep[inst1],cf.max_dep[inst2]]),0]

# get alignment statistics for each instrument
def get_params(df,inst):
    index = np.where((df['inst1']==inst) & (df['inst2']=='FLBBRR'))[0][0]
    n = df['n'][index]
    if (np.abs(1-df['line_a'][index])<df['line_a_95CI'][index]) & (inst != 'FLBBSR'):
        a = 1; sa = 0;
        if np.abs(df['offset_b'][index])<df['offset_b_95CI'][index]:
            b = 0; sb = 0;
        else:
            b = df['offset_b'][index]; sb = df['offset_b_sd'][index];
    else:
        a = df['line_a'][index]; sa = df['line_a_sd'][index];
        b = df['line_b'][index]; sb = df['line_b_sd'][index];
    return n,a,sa,b,sb

df = pd.read_csv(cf.SAVE_INFO_FN,skipinitialspace=True)
df = df.append(pd.DataFrame([['FLBBRR','FLBBRR',1,0,0,0,0,0]],
                       columns=('inst1','inst2','line_a','line_a_sd',
                                'line_b','line_b_sd','offset_b','offset_b_sd')),
                       ignore_index=True)
df = df.append(pd.DataFrame([['MCOMSBGC','FLBBRR',1,0,0,0,0,0]],
                       columns=('inst1','inst2','line_a','line_a_sd',
                                'line_b','line_b_sd','offset_b','offset_b_sd')),
                       ignore_index=True)
_,a_inst1,a_unc_inst1,b_inst1,b_unc_inst1 = get_params(df,inst1)
_,a_inst2,a_unc_inst2,b_inst2,b_unc_inst2 = get_params(df,inst2)

if cf.OUTPUT_INFO:
    print('For %s: Using alignment a=%.2f (%.2f) and b=%.1f (%.1f) x 1e-4 m-1'
          %(inst1,a_inst1,a_unc_inst1,b_inst1*1e4,b_unc_inst1*1e4))
    print('For %s: Using alignment a=%.2f (%.2f) and b=%.1f (%.1f) x 1e-4 m-1'
          %(inst2,a_inst2,a_unc_inst2,b_inst2*1e4,b_unc_inst2*1e4))
    print('')

## Standard variables
CENTRAL_LAT = 50.5
CENTRAL_LON = -145
DEG_LAT_TO_KM = sw.dist(CENTRAL_LAT+np.array([-.5,.5]),CENTRAL_LON,units='km')[0][0]
DEG_LON_TO_KM = sw.dist(CENTRAL_LAT,CENTRAL_LON+np.array([-.5,.5]),units='km')[0][0]
aspect = 1/np.cos(CENTRAL_LAT*np.pi/180)

## Find nearby indices
data1_loc = func1(files1,get_data=False)
data2_loc = func2(files2,get_data=False)

inds1,inds2,dist,dtime = get_bestfit.find_nearby((data1_loc[0],data2_loc[0]),
                                                 (data1_loc[1],data2_loc[1]),
                                                 (data1_loc[2],data2_loc[2]),
                                                 cf.TIME_THRES,cf.DIST_THRES,
                                                 return_distances=True)
if cf.OUTPUT_INFO:
    with np.printoptions(precision=2,suppress=True):
        print('%d matches found'%len(inds1))
        print('%s dives:'%label1,data1_loc[3][inds1])
        print('%s dives:'%label2,data2_loc[3][inds2])
        print('Distances (km)          ',dist)
        print('Time separations (hours)',dtime*24)
        print('')

if len(inds1)==0:
    if cf.OUTPUT_INFO:
        print('Quitting program because no matches were found between %s and %s.'%(inst1,inst2))
    exit()

## Get data from dives
data1 = func1(files1,get_dives=data1_loc[3][inds1])
data2 = func2(files2,get_dives=data2_loc[3][inds2])

## filter all data with median filters
good_inds_1 = [np.isfinite(np.array(T.astype('float32'))) & np.isfinite(np.array(S.astype('float32'))) & np.isfinite(np.array(b.astype('float32')))
               for T,S,b in zip(data1[4],data1[5],data1[7])]
if inst1 == 'BB2FLSG': # Seaglider VSF data is already filtered
    filtered_dep1 = [var[i] for var,i in zip(data1[6],good_inds_1)]
    filtered_VSF1 = [var[i] for var,i in zip(data1[7],good_inds_1)]
else:
    filtered_dep1 = [get_data.median_filter(var.astype('float32')[i]) for var,i in zip(data1[6],good_inds_1)]
    filtered_VSF1 = [get_data.median_filter(var.astype('float32')[i]) for var,i in zip(data1[7],good_inds_1)]
filtered_T1 = [np.interp(fd,d.astype('float32')[i],var.astype('float32')[i]) for fd,d,var,i in zip(filtered_dep1,data1[6],data1[4],good_inds_1)]
filtered_S1 = [np.interp(fd,d.astype('float32')[i],var.astype('float32')[i]) for fd,d,var,i in zip(filtered_dep1,data1[6],data1[5],good_inds_1)]
filtered_dens1 = [sw.pden(S,T,sw.pres(z,lat=CENTRAL_LAT))
                  for T,S,z in zip(filtered_T1,filtered_S1,filtered_dep1)]

good_inds_2 = [np.isfinite(T.astype('float32')) & np.isfinite(S.astype('float32')) & np.isfinite(b.astype('float32'))
               for T,S,b in zip(data2[4],data2[5],data2[7])]
if inst2 == 'BB2FLSG': # Seaglider VSF data is already filtered
    filtered_dep2 = [var[i] for var,i in zip(data2[6],good_inds_2)]
    filtered_VSF2 = [var[i] for var,i in zip(data2[7],good_inds_2)]
else:
    filtered_dep2 = [get_data.median_filter(var.astype('float32')[i]) for var,i in zip(data2[6],good_inds_2)]
    filtered_VSF2 = [get_data.median_filter(var.astype('float32')[i]) for var,i in zip(data2[7],good_inds_2)]
filtered_T2 = [np.interp(fd,d.astype('float32')[i],var.astype('float32')[i]) for fd,d,var,i in zip(filtered_dep2,data2[6],data2[4],good_inds_2)]
filtered_S2 = [np.interp(fd,d.astype('float32')[i],var.astype('float32')[i]) for fd,d,var,i in zip(filtered_dep2,data2[6],data2[5],good_inds_2)]
filtered_dens2 = [sw.pden(S,T,sw.pres(z,lat=CENTRAL_LAT))
                  for T,S,z in zip(filtered_T2,filtered_S2,filtered_dep2)]

## interpolate to density bins
dpyc = 0.1
interp_dens = np.arange(1024.05,1027.4,dpyc)
n = len(filtered_VSF1)
VSF_interp1,VSF_interp1_sd = get_bestfit.interpolate_data(filtered_VSF1,filtered_dens1,
                                                          np.meshgrid(interp_dens,np.arange(n))[0],
                                                          ddep=dpyc/2,return_uncert=True,
                                                          minimum_uncert=res_vsf1)
dep_interp1 = get_bestfit.interpolate_data(filtered_dep1,filtered_dens1,
                                           np.meshgrid(interp_dens,np.arange(n))[0])
T_interp1 = get_bestfit.interpolate_data(filtered_T1,filtered_dens1,
                                         np.meshgrid(interp_dens,np.arange(n))[0])
S_interp1 = get_bestfit.interpolate_data(filtered_S1,filtered_dens1,
                                         np.meshgrid(interp_dens,np.arange(n))[0])

VSF_interp2,VSF_interp2_sd = get_bestfit.interpolate_data(filtered_VSF2,filtered_dens2,
                                                          np.meshgrid(interp_dens,np.arange(n))[0],
                                                          ddep=dpyc/2,return_uncert=True,
                                                          minimum_uncert=res_vsf2)
dep_interp2 = get_bestfit.interpolate_data(filtered_dep2,filtered_dens2,
                                           np.meshgrid(interp_dens,np.arange(n))[0])
T_interp2 = get_bestfit.interpolate_data(filtered_T2,filtered_dens2,
                                         np.meshgrid(interp_dens,np.arange(n))[0])
S_interp2 = get_bestfit.interpolate_data(filtered_S2,filtered_dens2,
                                         np.meshgrid(interp_dens,np.arange(n))[0])

dsalt = np.abs(S_interp1-S_interp2).ravel()
dtemp = np.abs(T_interp1-T_interp2).ravel()

VSF_interp1_sw = np.array([betasw_HZP2019(cf.wv,sca1,T,S,sw.pres(z,lat=CENTRAL_LAT))[0].squeeze()
                           for T,S,z in zip(T_interp1,S_interp1,dep_interp1)])
VSF_interp2_sw = np.array([betasw_HZP2019(cf.wv,sca2,T,S,sw.pres(z,lat=CENTRAL_LAT))[0].squeeze()
                           for T,S,z in zip(T_interp2,S_interp2,dep_interp2)])

x = VSF_interp1.ravel()*2*np.pi*chi1 * a_inst1 + b_inst1
y = VSF_interp2.ravel()*2*np.pi*chi2 * a_inst2 + b_inst2

dxs = (unc_chi1/chi1*np.ones(shape=x.shape), # frac. uncert. on chi
       (VSF_interp1_sd/VSF_interp1).ravel(), # frac. uncert. based on vertical var. within the density bin
       cf.unc_sw*(VSF_interp1_sw/VSF_interp1).ravel()) # frac. uncert. of seawater contribution
dx_tot = np.sqrt(np.sum(np.array([x**2*d**2 for d in dxs]),axis=0))

dys = (unc_chi2/chi2*np.ones(shape=y.shape), # frac. uncert. on chi
       (VSF_interp2_sd/VSF_interp2).ravel(), # frac. uncert. based on vertical var. within the density bin
       cf.unc_sw*(VSF_interp2_sw/VSF_interp2).ravel()) # frac. uncert. of seawater contribution
dy_tot = np.sqrt(np.sum(np.array([y**2*d**2 for d in dys]),axis=0))

inds = np.isfinite(x*y*dx_tot*dy_tot) & (dsalt<cf.DSAL_THRES) & (dtemp<cf.DTEMP_THRES)
TSoutliers = np.isfinite(x*y*dx_tot*dy_tot) & ((dsalt>=cf.DSAL_THRES) | (dtemp>=cf.DTEMP_THRES))

## Do a first fit to find outliers
data = odr.RealData(x[inds], y[inds], sx=dx_tot[inds], sy=dy_tot[inds])
odr_line = odr.ODR(data, odr.unilinear, beta0=[1,0])
fit_line = odr_line.run()
odr_offset = odr.ODR(data, get_bestfit.offset_model(), beta0=[0])
fit_offset = odr_offset.run()
outliers = get_bestfit.find_outliers(fit_line.beta,x,y,dx_tot,dy_tot,inds) & inds
iinds = inds & ~outliers
r = stats.pearsonr(x[iinds],y[iinds])[0]
if cf.OUTPUT_INFO:
    print('n=%d (%d outliers removed), r^2=%.3f'%(np.sum(iinds),np.sum(outliers),r**2))

## Do second fit without the outliers
data = odr.RealData(x[iinds], y[iinds], sx=dx_tot[iinds], sy=dy_tot[iinds])
odr_line = odr.ODR(data, odr.unilinear, beta0=[1,0])
fit_line = odr_line.run()
odr_offset = odr.ODR(data, get_bestfit.offset_model(), beta0=[0])
fit_offset = odr_offset.run()
if cf.OUTPUT_INFO:
    print('Best-fit  line:  b=%+.3f (%.3f) x1e-4, a=%.3f (%.3f)'
          %(fit_line.beta[1]*1e4,fit_line.sd_beta[1]*1e4,fit_line.beta[0],fit_line.sd_beta[0]))
    print('Best-fit offset, b=%+.3f (%.3f) x1e-4'
          %(fit_offset.beta[0]*1e4,fit_offset.sd_beta[0]*1e4))
    print('')

# total error also contains measurement-independent errors in scaling factors and darks
a_line_total_unc = np.sqrt(fit_line.sd_beta[0]**2 +
                           fit_line.beta[0]**2 * (unc_sf1**2 + unc_sf2**2 + a_unc_inst1**2 + a_unc_inst2**2))
b_line_total_unc = np.sqrt(fit_line.sd_beta[1]**2 +
                           fit_line.beta[0]**2 * ((unc_dk1*2*np.pi*chi1)**2 + b_unc_inst1**2) +
                           fit_line.beta[1]**2 * (unc_sf2**2 + a_unc_inst2**2) +
                           ((unc_dk2*2*np.pi*chi2)**2 + b_unc_inst2**2))
b_offset_total_unc = np.sqrt(fit_offset.sd_beta[0]**2 + (unc_dk1*2*np.pi*chi1)**2 + b_unc_inst1**2 + (unc_dk2*2*np.pi*chi2)**2 + b_unc_inst2**2)
if cf.OUTPUT_INFO:
    print('Total error, adding in uncertainty in the calibration and alignment coefficients:')
    print('Best-fit  line:  b=%+.3f (%.3f) x1e-4, a=%.3f (%.3f)'
          %(fit_line.beta[1]*1e4,b_line_total_unc*1e4,fit_line.beta[0],a_line_total_unc))
    print('Best-fit offset, b=%+.3f (%.3f) x1e-4'
          %(fit_offset.beta[0]*1e4,b_offset_total_unc*1e4))
    print('')

df_line = np.sum(iinds)-2
df_offset = np.sum(iinds)-1
CI_95_line_a = stats.t.interval(0.95,df_line,scale=a_line_total_unc)[1]
CI_95_line_b = stats.t.interval(0.95,df_line,scale=b_line_total_unc)[1]
CI_95_offset_b = stats.t.interval(0.95,df_offset,scale=b_offset_total_unc)[1]

if cf.OUTPUT_INFO:
    print('95% Confidence intervals:')
    print('Best-fit line (df=%d): b=%+.1f (%.1f) x1e-4, a=%.2f (%.2f)'
          %(df_line,fit_line.beta[1]*1e4,CI_95_line_b*1e4,fit_line.beta[0],CI_95_line_a))
    print('Best-fit offset (df=%d), b=%+.1f (%.1f) x1e-4'
          %(df_offset,fit_offset.beta[0]*1e4,CI_95_offset_b*1e4))
    print('')


if cf.SAVE_INFO:
    f=open(cf.SAVE_ALIGNED_INFO_FN,'a')
    f.write('%s, %s, %d, %d, %d, %d, %.3f, %.4f, %.4f, %.4f, %.3e, %.3e, %.3e, %.3e, %.3e, %.3e\n'%
            (inst1,inst2,n,np.sum(iinds),np.sum(TSoutliers),np.sum(outliers),r**2,
             fit_line.beta[0],a_line_total_unc,CI_95_line_a,
             fit_line.beta[1],b_line_total_unc,CI_95_line_b,
             fit_offset.beta[0],b_offset_total_unc,CI_95_offset_b))
    f.close()


dists = np.tile(np.atleast_2d(dist).T,len(interp_dens)).ravel()
dtimes = np.tile(np.atleast_2d(dtime).T,len(interp_dens)).ravel()

## Plot figure
XX = np.array([0,0.003])
if cf.PLOT_FIG:
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(ncols=4, nrows=2, figure=fig)
    ax = [fig.add_subplot(gs[:, 0:2]),fig.add_subplot(gs[0,2]),fig.add_subplot(gs[0,3]),fig.add_subplot(gs[1,2]),fig.add_subplot(gs[1,3])]
    ax[0].errorbar(x[iinds],y[iinds],xerr=dx_tot[iinds],yerr=dy_tot[iinds],ls='none',
               label='Fitted data (n=%d)'%np.sum(iinds),clip_on=False,marker=',',color='k')
    ax[0].plot(XX,np.polyval(fit_line.beta,XX),color='k',
           label='best-fit line ($r^2$=%.2f):\na=%.2f (%.2f), b=%+.1f (%.1f) x10$^{-4}$ m$^{-1}$'
           %(r**2,fit_line.beta[0],a_line_total_unc,fit_line.beta[1]*1e4,b_line_total_unc*1e4))
    ax[0].plot(XX,XX+fit_offset.beta,color='k',ls='--',label='best-fit offset:\nb=%+.1f (%.1f) x10$^{-4}$ m$^{-1}$'
           %(fit_offset.beta[0]*1e4,b_offset_total_unc*1e4))
    ax[0].plot(XX,XX,label='1:1',color='k',ls=':')
    ax[0].scatter(x[TSoutliers],y[TSoutliers],c='0.5',zorder=-1,marker='.',label='Discarded due to T/S (n=%d)'%np.sum(TSoutliers))
    ax[0].scatter(x[outliers],y[outliers],c='r',zorder=-2,marker='x',s=50,label='Outliers (n=%d)'%np.sum(outliers))
    e = ax[1].errorbar(x[iinds],np.concatenate(dep_interp1)[iinds],xerr=dx_tot[iinds],
                   clip_on=False,ls='none',marker='.',alpha=0.5,label=label1)
    for b in e[2]:
        b.set_clip_on(False)
    e = ax[1].errorbar(y[iinds],np.concatenate(dep_interp1)[iinds],xerr=dy_tot[iinds],
                   clip_on=False,ls='none',marker='.',alpha=0.5,label=label2)
    for b in e[2]:
        b.set_clip_on(False)
    ax[2].plot(np.vstack((data1[2],data2[2])),np.vstack((data1[1],data2[1])),c='k',lw=1,zorder=-1)
    ax[2].scatter(data1[2],data1[1],label=label1,s=10)
    ax[2].scatter(data2[2],data2[1],label=label2,s=10)
    ax[3].scatter(dist,dtime*24,clip_on=False,marker='.',c='k')
    ax[4].scatter(dsalt[iinds],dtemp[iinds],marker='.',c='k',label='Fitted data')
    ax[4].scatter(dsalt[outliers],dtemp[outliers],marker='x',c='r',s=50,label='Outliers')

    ax[0].set_xlim(cf.xx)
    ax[0].set_ylim(cf.xx)
    ax[0].set_xlabel(label1+' $b_{bp}$ (m$^{-1}$)')
    ax[0].set_ylabel(label2+' $b_{bp}$ (m$^{-1}$)')
    ax[0].legend(loc='upper left',frameon=False)
    ax[0].set_aspect(1)
    ax[0].set_xticks(np.arange(0,0.0021,0.0005))
    ax[0].set_yticks(np.arange(0,0.0021,0.0005))
    ax[1].set_xlim(cf.xx)
    ax[1].set_xticks(np.arange(0,0.0025,0.001))
    ax[1].set_ylim(yy)
    ax[1].set_xlabel('$b_{bp}$ (m$^{-1}$)')
    ax[1].set_ylabel('Equivalent depth (m)')
    ax[2].set_xlabel('Longitude')
    ax[2].set_ylabel('Latitude')
    ax[2].set_xlim([-145.8,-144])
    ax[2].set_ylim([49.6,51.1])
    ax[2].set_aspect(aspect)
    ax[3].set_xlim([0,cf.DIST_THRES])
    ax[3].set_ylim([0,cf.TIME_THRES*24])
    ax[3].set_xlabel('Spatial separation (km)')
    ax[3].set_ylabel('Temporal separation (hr)')
    ax[4].set_xlim([0,cf.DSAL_THRES])
    ax[4].set_ylim([0,cf.DTEMP_THRES])
    ax[4].set_xlabel('Salinity difference (PSU)')
    ax[4].set_ylabel('Temperature difference (deg C)')
    #ax[4].legend(loc='lower right',frameon=False)

    for a,p in zip(ax,['a','b','c','d','e']):
        a.text(0,1.01,p,fontsize=16,fontweight='bold',transform=a.transAxes)

    plt.tight_layout()
    ax[1].legend(loc=0,title='%d profiles'%len(data1[0]))

    if cf.SAVE_FIG:
        plt.savefig(figure_name,dpi=300)
        if cf.OUTPUT_INFO:
            print('Figure saved to %s'%figure_name)
    else:
        plt.show()
