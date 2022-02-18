def bio_de(spacing, stn_start, stn_end, year):
    import arrow
    import time
    import xarray as xr
    import numpy as np
    import map_fxn as mf
    
    if year == 2016:
        noday = 366
    else:
        noday = 365

    day_start = 0
    day_end = noday
    print('Spacing between stations: ' + str(spacing))
    print('First day: ' + str(day_start))
    print('Last day: ' + str(day_end))
    print('testing reload')
    print(year)

    dirname = '/data/tjarniko/MEOPAR/analysis-tereza/notebooks/CLUSTER_201905/NC_HINDCAST/' + str(year) +'/BIO_TS/'
    
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    
    no_stns = len(d_stn_x)
    

    # a list of all the model outputs (tracers) for this year
    trace_list = mf.create_tracelist(year)
    
    print("Number of stations:" + str(no_stns))
    
    stn_names = list()
    for i in range(0,no_stns):
        stn_name = 'stn_' + str(i)
        stn_names.append(stn_name)
        #print(stn_name)
    

    for n in range(stn_start,stn_end):
        print('stn is: ' + str(n))
        start_time = time.time()
        
        MYRI_ar = np.zeros(noday)
        MICZ_ar = np.zeros(noday)
        PHY2_ar = np.zeros(noday)
        PHY_ar = np.zeros(noday)
        for i in range(day_start,day_end):
            if(i%20 == 0):
                print('day is:' + str(i))
            b = i
            e = i+1
            start_time = time.time()
            trc = (trace_list[i])
            nemo = xr.open_mfdataset(trc)
            t1 = b
            t2 = e
            x1 = d_stn_x[n]
            #print(x1)
            x2 = d_stn_x[n]+1
            y1 = d_stn_y[n]
            #print(y1)
            y2 = d_stn_y[n]+1
            z1 = 0
            z2 = 39
            vol,spacetime, MYRI_d, MICZ_d, PHY2_d, PHY_d = int_phyt_fxn(grid,nemo,t1,t2,x1,x2,y1,y2,z1,z2)
            
            
            MYRI_ar[i] = sum(np.squeeze(MYRI_d))
            MICZ_ar[i] = sum(np.squeeze(MICZ_d))
            PHY2_ar[i] = sum(np.squeeze(PHY2_d))
            PHY_ar[i] = sum(np.squeeze(PHY_d))
            
        jan = xr.Dataset({'MYRI': (['t'],  MYRI_ar),'MICZ': (['t'],  MICZ_ar)
         ,'PHY2': (['t'],  PHY2_ar),'PHY': (['t'],  PHY_ar)})
            

        stn_name = dirname + stn_names[n] + '_sp' + str(spacing)+ '.nc'
            #print(stn_name)
        jan.to_netcdf(stn_name)
            

def int_phyt_fxn(grid,d_set,t1,t2,x1,x2,y1,y2,z1,z2):
    import numpy as np



    delt_y = grid.variables['e1t'][0,y1:y2,x1:x2]
    delt_x = grid.variables['e2t'][0,y1:y2,x1:x2]
    delt_z = grid.variables['e3t'][0,z1:z2,y1:y2,x1:x2]
    #print(delt_y)
    #print(delt_x)
    #print('defining deltas, areas, and volumes')
    d_y2 = delt_y[0]
    d_y = d_y2[0]
    #print(d_y)
    #d_x = delt_x.values[:,:]
    d_x2 = delt_x[0]
    d_x = d_x2[0]
    #print(d_x)
    
    d_z = delt_z
    #print(d_z)
    area = d_y*d_x
    #print(area)
    ar_z = np.empty_like(d_z)
    ar_z[:,:,:] = area
    vol = ar_z * d_z
    
    #print('defining biology model slice')
    l_MYRI = d_set.ciliates[:,z1:z2,y1:y2,x1:x2]    
    l_MICZ = d_set.microzooplankton[:,z1:z2,y1:y2,x1:x2]
    l_PHY2 = d_set.flagellates[:,z1:z2,y1:y2,x1:x2]
    l_PHY = d_set.diatoms[:,z1:z2,y1:y2,x1:x2]
    
    spacetime = np.empty_like(l_PHY.values)

    spacetime[:,:,:,:] = vol
    #print('multiplying out volumes by phytoplankton')
    MYRI_d = spacetime*l_MYRI.values
    MICZ_d = spacetime*l_MICZ.values
    PHY2_d = spacetime*l_PHY2.values
    PHY_d  = spacetime*l_PHY.values
    
    return vol,spacetime, MYRI_d, MICZ_d, PHY2_d, PHY_d
