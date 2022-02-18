def weights_intrp_mtx_ops(odata):
    '''POINT: interpolates data from the atmospheric grid onto the ocean grid using a weights matrix | CALLS: | NOTES: this is much faster than interpolating with scipy! | TAKES : one array of size lats 266 x lons 256 of atmospheric data | GIVES: an array of NEMO-Salish size of the new data interpolated | USAGE: | ROLE IN CLUSTER PROJECT: | written by Michael Dunphy, adapted by tjšj, 2017'''
    from IPython import embed
    import netCDF4 as nc
    import scipy.interpolate as spi
    import scipy.sparse as sp
    import numpy as np
    '''takes a matrix of size lon(256), lat(266) and interpolates it onnto the NEMO gird
    size is 266 x 256'''
    weightsfile = '/home/mdunphy/MEOPAR/NEMO-forcing/grid/weights-gem2.5-ops.nc'
    with nc.Dataset(weightsfile) as f:
        s1 = f.variables['src01'][:]-1  # minus one for fortran-to-python indexing
        s2 = f.variables['src02'][:]-1
        s3 = f.variables['src03'][:]-1
        s4 = f.variables['src04'][:]-1
        w1 = f.variables['wgt01'][:]
        w2 = f.variables['wgt02'][:]
        w3 = f.variables['wgt03'][:]
        w4 = f.variables['wgt04'][:]
       
    NO = odata.size   # number of operational grid points
    NN = s1.size      # number of NEMO grid points
    
    # Build matrix
    n = np.array([x for x in range(0,NN)])
    M1 = sp.csr_matrix((w1.flatten(), (n, s1.flatten())), (NN,NO))
    M2 = sp.csr_matrix((w2.flatten(), (n, s2.flatten())), (NN,NO))
    M3 = sp.csr_matrix((w3.flatten(), (n, s3.flatten())), (NN,NO))
    M4 = sp.csr_matrix((w4.flatten(), (n, s4.flatten())), (NN,NO))
    M = M1+M2+M3+M4
    
        # Interpolate by matrix multiply - quite fast
    ndata = M*odata.flatten()

    # Reshape to NEMO shaped array
    ndata=ndata.reshape(s1.shape)
    
    return ndata

def weights_intrp_mtx_gemlam(odata):
    '''POINT: interpolates data from the atmospheric grid onto the ocean grid using a weights matrix | CALLS: | NOTES: this is much faster than interpolating with scipy! | TAKES : one array of size lats 266 x lons 256 of atmospheric data | GIVES: an array of NEMO-Salish size of the new data interpolated | USAGE: | ROLE IN CLUSTER PROJECT: | written by Michael Dunphy, adapted by tjšj, 2017'''
    from IPython import embed
    import netCDF4 as nc
    import scipy.interpolate as spi
    import scipy.sparse as sp
    import numpy as np
    '''takes a matrix of size lon(256), lat(266) and interpolates it onnto the NEMO gird
    size is 266 x 256'''
    weightsfile = '/data/tjarniko/MEOPAR/grid/weights-gem2.5-gemlam_201702.nc'
    with nc.Dataset(weightsfile) as f:
        s1 = f.variables['src01'][:]-1  # minus one for fortran-to-python indexing
        s2 = f.variables['src02'][:]-1
        s3 = f.variables['src03'][:]-1
        s4 = f.variables['src04'][:]-1
        w1 = f.variables['wgt01'][:]
        w2 = f.variables['wgt02'][:]
        w3 = f.variables['wgt03'][:]
        w4 = f.variables['wgt04'][:]
       
    NO = odata.size   # number of operational grid points
    NN = s1.size      # number of NEMO grid points
    
    # Build matrix
    n = np.array([x for x in range(0,NN)])
    M1 = sp.csr_matrix((w1.flatten(), (n, s1.flatten())), (NN,NO))
    M2 = sp.csr_matrix((w2.flatten(), (n, s2.flatten())), (NN,NO))
    M3 = sp.csr_matrix((w3.flatten(), (n, s3.flatten())), (NN,NO))
    M4 = sp.csr_matrix((w4.flatten(), (n, s4.flatten())), (NN,NO))
    M = M1+M2+M3+M4
    
        # Interpolate by matrix multiply - quite fast
    ndata = M*odata.flatten()

    # Reshape to NEMO shaped array
    ndata=ndata.reshape(s1.shape)
    
    return ndata

import arrow
import numpy as np
import glob



def get_windar(start,end):
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    thres = arrow.get('2014-09-12')

    arrow_array = []
    ncar = []

    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)
    #print(dayslen)
    ts = np.zeros(dayslen)
    #print(arrow_array)

    # ops_y2014m09d12.nc

    for i in range(0,dayslen):

        tdate = arrow_array[i][0]
        yy = tdate.format('YYYY')
        mm = tdate.format('MM')
        dd = tdate.format('DD')
        ymd = f'y{yy}m{mm}d{dd}'

        if tdate < thres:
            dirstring = '/results/forcing/atmospheric/GEM2.5/gemlam/'
            guess = 'gemlam_'+ ymd +'.nc' 
        else:
            dirstring = '/results/forcing/atmospheric/GEM2.5/operational/'
            guess = 'ops_'+ ymd +'.nc'

        #tdir = '/data/tjarniko/results/RIV_PIL/BASE/ncs/'

        w = glob.glob(dirstring+guess)
        w = w[0]
        ncar.append(w)
    
    return ncar

def yearwinds_de(spacing,stn_b,stn_e,year):
    '''POINT: for a given set of stations, produces 3 timeseries of wind forcing data (magnitude, stress, energy)  station as lists in a single netcdf file (one netcdf file per station) | CALLS: import_bathy, make_stns, filster_station_in_domain (from map_fxn), reads in netcdf files produced by produce_interpolated_wind_netcdf | NOTES: days start at 0, stations start at 0, for a spacing of 10 there should be 580 stations, the start and end station option is so that it can be run in parallel if necessary, directory into which nc files are deposited is hardcoded under stn_name | TAKES : spacing (has been decided as 10 for this project), year, start station, end station | GIVES: netcdf files that contain 3 windrelated timeseries for each station, in a given repo | USAGE: yearwinds_de(10,1,580,2016) | ROLE IN CLUSTER PROJECT: | written 2017, tjšj '''
    import sys
    sys.path.append('./extraction_scripts')
    
    import map_fxn as mf
    import time
    import xarray as xr
    import numpy as np

    print('Spacing between stations: ' + str(spacing))
    
    noday = 365
    if year == 2016:
        noday = 366
    
    print('no. days, from yearwinds_de fxn')
    print(noday)
    
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    print('number of stns: ' + str(len(d_stn_x)))
    
    for stn in range(stn_b,stn_e):
    
        print('station is: ' ,str(stn))
        print('x is :', d_stn_x[stn])
        print('y is :', d_stn_y[stn])
        
        ts_x = d_stn_x[stn]
        ts_y = d_stn_y[stn]
        
        daily_windmag = np.zeros(noday)
        daily_windstress = np.zeros(noday)
        daily_windenergy = np.zeros(noday)
        
        end_day = noday
        
        start = f'{str(year)}-01-01'
        end = f'{str(year)}-12-31'

        start_run = arrow.get(start)
        end_run = arrow.get(end)

        arrow_array = []
        ncar = []

        for r in arrow.Arrow.span_range('day', start_run, end_run):
            arrow_array.append(r)

        dayslen = len(arrow_array)
        #print(dayslen)
        ts = np.zeros(dayslen)
        #print(arrow_array)

        # ops_y2014m09d12.nc

        for i in range(0,dayslen):

            tdate = arrow_array[i][0]
            yy = tdate.format('YYYYMMDD')
            print(i)
#         for day in range(0,end_day):
#             print(day)
#             ##open the winds file
#             #dayfile
            
            filestart = f'./WINDFILES_interp/windint_{yy}.nc'
            
            #print(day)
            #print(filestart) 
            wf = xr.open_dataset(filestart)
            ## extract wm, ws, we for the station for the day
            ## attach it to the corresponding array above
            ## save the resulting 3 lists to a netcdf file as below
            
            #avg_daily_windmag_windstress_energy
            
            daily_we = wf.daily_avg_windenergy[ts_y,ts_x]
            daily_we = daily_we.values
            daily_ws = wf.daily_avg_windstress[ts_y,ts_x]
            daily_ws = daily_ws.values
            daily_wm = wf.daily_avg_windmag[ts_y,ts_x]
            daily_wm = daily_wm.values
            ##be gentle about OB1 errors!!!!
            

            
            daily_windmag[i] = daily_wm
            daily_windstress[i] = daily_ws
            daily_windenergy[i] = daily_we
            
            
            winds = xr.Dataset({'wind_mags':(['t'], daily_windmag), 'wind_stresses':(['t'], daily_windstress), 
                                'wind_energy':(['t'], daily_windenergy)})
            stn_name = '/data/tjarniko/MEOPAR/analysis_tereza/notebooks/CLUSTER_PAPER/CLEAN/NC_HINDCAST/' + str(year) + '/WIND_TS/stn_' + str(stn) + '_wind_data_sp' + str(spacing) + '.nc' 
            winds.to_netcdf(stn_name)
            
