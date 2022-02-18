def yearved_de(spacing,stn_b,stn_e,year):
    
    import netCDF4 as nc
    import map_fxn as mf
    import time
    import xarray as xr
    import numpy as np
    print('Spacing between stations: ' + str(spacing))
    
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    grid2 = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    
    no_stns = len(d_stn_x)
    monlist = ['jan','feb','mar','apr','may','jun','jul','aug','sep',
              'oct','nov','dec']
    if year == 2016: 

        daylist = [31,29,31,30,31,30,31,31,30,31,30,31]
        noday = 366
        print(noday)
    if year != 2016:
        daylist = [31,28,31,30,31,30,31,31,30,31,30,31]
        noday = 365
        print(noday)
        

    # a list of all the model outputs (tracers) for this year    
    trace_list = mf.create_physlist_W(year)
    
    print("Number of stations:" + str(no_stns))
    
    for stn in range(stn_b,stn_e):
    
        print('station is: ' ,str(stn))
        print('x is :', d_stn_x[stn])
        print('y is :', d_stn_y[stn])
        
        ts_x = d_stn_x[stn]
        ts_y = d_stn_y[stn]
        
        daily_ved = np.zeros(noday)
        
        end_day = noday
        for day in range(1,end_day+1):

            if (day%20 ==0) :
                print(day)
            trc = trace_list[day-1]
            nemo = nc.Dataset(trc)
            #print(nemo)
            #w = np.array(nemo.variables['vert_eddy_diff'])
            ved  = nemo['vert_eddy_diff'][0,:,ts_y,ts_x]
            where0 = (np.where(ved == 0))
            #first out of bounds cell
            stop = where0[0][1]
            #get only applicable veds
            veds = ved[0:stop]
            #get applicable depths
            depths = grid2['e3w_0'][0,0:stop,ts_y,ts_x]

            #get total depth to fially divide by 
            totdepth = np.sum(depths)
            #weight veds by depths - ie size of cell
            veds_depths = veds*depths
            #get avg
            avg_ved = np.sum(veds_depths/totdepth)

            
            daily_ved[day-1] = avg_ved
    
        ved = xr.Dataset({'daily_ved':(['t'], daily_ved)})
        stn_name = '/data/tjarniko/MEOPAR/analysis_tereza/notebooks/CLUSTER_PAPER/CLEAN/NC_HINDCAST/' + str(year) + '/VED_TS/stn_' + str(stn)  + 'avg_ved_sp' + str(spacing)+ '.nc'
        ved.to_netcdf(stn_name)
