
def freshwater_retrieval_4m(stn_b,stn_e,threshold,year):
    import netCDF4 as nc 
    import bathy_fxn as bf
    import xarray as xr
    import numpy as np
    import map_fxn as mf
    
    dirst = '/data/tjarniko/MEOPAR/analysis-tereza/notebooks/CLUSTER_201905/NC_HINDCAST/' + str(year)
    monlist = ['jan','feb','mar','apr','may','jun','jul','aug','sep',
      'oct','nov','dec']
    if year == 2016: 
        daylist = [31,29,31,30,31,30,31,31,30,31,30,31]
        noday = 366
    else:
        daylist = [31,28,31,30,31,30,31,31,30,31,30,31]
        noday = 365    

    tracelist = mf.create_physlist(year)
    
    print(year)
    print(noday)
    
    #load station data:
    stn_nc = '/data/tjarniko/MEOPAR/at3/notebooks/CLUSTER/verze2pt0/ncfiles/stn_data_sp10.nc'
    
    stn_dat = nc.Dataset(stn_nc)
    inds = np.array(stn_dat['index'])
    stn_xs = np.array(stn_dat['stn_x'])
    stn_ys = np.array(stn_dat['stn_y'])
    depths = np.array(stn_dat['depth'])
    
    #load bathymetry data:
    bathy = nc.Dataset('//data/tjarniko/MEOPAR/grid/bathymetry_201702.nc', 'r')
    print('bathy loaded')   
    #load grid:
    grid = nc.Dataset('//data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    
    #establish which cells in the domain are shallower/deeper than given threshold
    depth_bin = bf.threshold_value(bathy, threshold)
    
    trc = tracelist[0]
    #print(trc)
    
    #retrieve depths
    trd = nc.Dataset(trc)
    depth_list = np.array(trd.variables['deptht'])
    #print(depth_list)
    #what cell in z-direction corresponds to threshold depth?
    depth_ind = bf.find_closest_cell_to_given_depth(trc,threshold)
    
    
    print('are you still with me computer')
    for s in range(stn_b,stn_e):
        
        #retrieve station data for this station
        stn = inds[s]
        stn_x = stn_xs[s]
        stn_y = stn_ys[s]
        depth = depths[s]
        
        print('stn')
        print(s)
        #deep or shallow stn? shallow = 1, deep = 2
        DOS = 0
        if (depth >= threshold):
            DOS = 2
        if (depth < threshold):
            DOS = 1
        #print('status')
        #print(DOS)
        
        #depths of cells:
        #print(grid)
        dz = np.array(grid.variables['e3t_0'])
        dz = np.squeeze(dz)
        
        dz = dz[:,stn_y,stn_x]
        #print(dz)
        zeds_to_mult = dz[0:4]
        
        fwi_ar = np.zeros(noday)
        fwi_pm_ar = np.zeros(noday)
        
        for d in range(0,noday):
            #print(d)
            if (d%50 == 0):
                print(d)
            trc = tracelist[d]
            #print(trc)
            trd = nc.Dataset(trc)
            
            sal = np.array(trd.variables['vosaline'])
            sal = np.squeeze(sal)

            #print(sal.shape)
            
            #establish salinity threshold
            #if the station is deep enough, it's simple - its just the salinity at that depth
            
            if (DOS ==2):
                sal_at_depth = sal[depth_ind,stn_y,stn_x]
                #salinity profile 
                sal_prof = sal[0:4,stn_y,stn_x]                
                #multiply
                differences = np.zeros_like(sal_prof)
                differences = -(sal_prof-sal_at_depth)
                # if things are saltier we are just not going to worry about it
                differences[differences <0] =0
                fwi = zeds_to_mult * differences
                fwi = sum(fwi)
                
                
                
            if (DOS ==1):
                x_close, y_close, closest = bf.closest_deep_station_with_array(stn_y,stn_x,depth_bin, sal)
                sal_at_depth = sal[depth_ind,y_close,x_close]
                
                #where does this stn end?
                sal_prof_pr = sal[:,stn_y,stn_x]
                #where_ends = np.where(sal_prof_pr == 0)
                #where_ends = np.squeeze(where_ends)
                #if where_ends.shape == ():
                #    where_ends = where_ends
                #else:
                #    where_ends = where_ends[0]
                sal_prof = sal_prof_pr[0:4]                
                differences = np.zeros_like(sal_prof)
                differences = -(sal_prof-sal_at_depth)
                # if things are saltier we are just not going to worry about it
                differences[differences <0] =0
                zeds_to_mult = dz[0:4]
                fwi = zeds_to_mult * differences
                fwi = sum(fwi)
                
            fwi_ar[d] = fwi

        stn_name = dirst + '/FWI_TS/stn_' + str(s) + '_fwi4m_data_sp10_threshold' + str(threshold) + '.nc' 
        print(stn_name)
          #     jan = xr.Dataset({'MYRI': (['t'],  MYRI_ar),'MICZ': (['t'],  MICZ_ar)
         #,'PHY2': (['t'],  PHY2_ar),'PHY': (['t'],  PHY_ar)})
            
            
        fresh = xr.Dataset({'freshwater_index':(['t'], fwi_ar) })
        print(fresh)
        fresh.to_netcdf(stn_name)
  
    return
