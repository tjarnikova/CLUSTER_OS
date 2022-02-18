def yearhalo_de(spacing,stn_b,stn_e,year):
    
    import map_fxn as mf
    import time
    import xarray as xr
    import numpy as np
    print('Spacing between stations: ' + str(spacing))
    
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    
    no_stns = len(d_stn_x)
    monlist = ['jan','feb','mar','apr','may','jun','jul','aug','sep',
              'oct','nov','dec']

    if year == 2016: 
        daylist = [31,29,31,30,31,30,31,31,30,31,30,31]
        noday = 366  
    else:
        daylist = [31,28,31,30,31,30,31,31,30,31,30,31]
        noday = 365    

    
    print(year)
    print(noday)
        
    # a list of all the model outputs (tracers) for this year
    trace_list = mf.create_physlist(year)
    
    print("Number of stations:" + str(no_stns))
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])

    for stn in range(stn_b,stn_e):
    
        print('station is: ' ,str(stn))
        print('x is :', d_stn_x[stn])
        print('y is :', d_stn_y[stn])
        
        ts_x = d_stn_x[stn]
        ts_y = d_stn_y[stn]
        
        daily_halocline = np.zeros(noday)
        
        end_day = noday
        for day in range(1,end_day+1):
            if (day%20 ==0) :
                print(day)
            trc = trace_list[day-1]
            #print(trc)
            nemo = xr.open_mfdataset(trc)
            halo = halo_de(nemo,grid,ts_x,ts_y)
            q = isinstance(halo, (np.ndarray))
            if (q == True) :
                print('found an array!')
                daily_halocline[day-1] = halo[0]
            else:
                daily_halocline[day-1] = halo

        halo = xr.Dataset({'halocline_depth':(['t'], daily_halocline)})
        stn_name = '/data/tjarniko/MEOPAR/analysis-tereza/notebooks/CLUSTER_201905/NC_HINDCAST/' +  str(year)+ '/HALO_TS/stn_' + str(stn)  + 'halo_depth_sp' + str(spacing)+ '.nc'
        halo.to_netcdf(stn_name)
        

def halo_de(nemo,grid,ts_x,ts_y):
    import numpy as np
    
    #retreive distances between grid cells ("distance between this grid cell and the one before it")

    delt_z = grid.variables['e3t'][0,:,ts_y:ts_y+1,ts_x:ts_x+1]
    #remove the stupid singletons
    delt_z = delt_z.values
    delt_z = np.squeeze(delt_z)
    
    sal = nemo.vosaline[0,:,ts_y:ts_y+1,ts_x:ts_x+1]
    depth = nemo.deptht
    dv = depth.values[:]
    sv = sal.values[:]
    sv = np.squeeze(sv)

    #for this grid cell, where do we no longer have salinity values?
    bottom = np.where(sv ==0)
    bottom = bottom[0]
    bottom = bottom[0]

    fxnl_sal = sv[0:bottom]
    fxnl_e3t = delt_z[0:bottom]
    fxnl_depth = dv[0:bottom]
    
    sal_gradient = np.zeros_like(fxnl_depth)

    #calculating gradient
    for i in range(1,bottom):
        sal_gradient[i] = (fxnl_sal[i]-fxnl_sal[i-1])/fxnl_e3t[i]
        
    mg = max(sal_gradient)
    #print(mg)

    halocline = np.where(sal_gradient == mg)
    halocline = np.squeeze(halocline)

    halocline = fxnl_depth[halocline]

    
    return halocline
