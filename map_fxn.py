import xarray as xr
import glob
import numpy as np
import arrow

def import_bathy(filepath):
    filepath = '/data/tjarniko/MEOPAR/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = xr.open_dataset(filepath)
    return grid

def make_stns(spacing):
    xs =[]
    ys =[]
    count = 0
    stn_x = []
    stn_y = []
    for i in range (0,200):
        x = 2+i*spacing
        if(x>397):
            break
        else:
            xs.append(x)
        
    for j in range(0,400):
        y = 2+j*spacing
        if(y>897):
            break
        else:
            ys.append(y)
            
    for a in range(0,len(xs)):
        for b in range(0,len(ys)):
            ts_x = xs[a]
            ts_y = ys[b]
            stn_x.append(ts_x)
            stn_y.append(ts_y)
            
    return stn_x, stn_y

def filter_stn_in_domain(stn_x,stn_y,fmask):
    d_stn_x = []
    d_stn_y = []
    for s in range(0,len(stn_x)):
        x = stn_x[s]
        y = stn_y[s]
        stn = fmask.values[y-1:y+2,x-1:x+2]
        #if there are no zeroes in s
        if((0 in stn) == False):
            d_stn_x.append(x)
            d_stn_y.append(y)
            
    return d_stn_x, d_stn_y

def create_physlist(year):
    import arrow
    start = str(year) + '-01-01'
    end = str(year) + '-12-31'

    start_run = arrow.get(start)
    end_run = arrow.get(end)

    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    date_array = []

    for i in range(0,len(arrow_array)):    
        q = arrow_array[i][0]  
        numdat = q.format('YYYYMMDD').lower() 
        date_array.append(numdat)

    respath = '/data/tjarniko/avg/'
    trace_list = []
    for i in range(0,len(date_array)):
        
        Day_OI = date_array[i]
        print(Day_OI)
        tracers = 'SalishSea_1*' + Day_OI + '_grid_T.nc'
        
        ptt = respath + tracers
        trace = glob.glob(ptt)
        trace = (trace[0])
        print(trace)
        trace_list.append(trace)
    
    return trace_list

def create_tracelist(year):
    import arrow
    start = str(year) + '-01-01'
    end = str(year) + '-12-31'

    start_run = arrow.get(start)
    end_run = arrow.get(end)

    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    date_array = []

    for i in range(0,len(arrow_array)):    
        q = arrow_array[i][0]  
        numdat = q.format('YYYYMMDD').lower() 
        date_array.append(numdat)

    respath = '/data/tjarniko/avg/'
    trace_list = []
    for i in range(0,len(date_array)):
        
        Day_OI = date_array[i]
        print(Day_OI)
        tracers = 'SalishSea_1*' + Day_OI + '_ptrc_T.nc'
        
        ptt = respath + tracers
        trace = glob.glob(ptt)
        trace = (trace[0])
        print(trace)
        trace_list.append(trace)
    
    return trace_list

def create_physlist_W(year):
    import arrow
    start = str(year) + '-01-01'
    end = str(year) + '-12-31'

    start_run = arrow.get(start)
    end_run = arrow.get(end)

    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    date_array = []

    for i in range(0,len(arrow_array)):    
        q = arrow_array[i][0]  
        numdat = q.format('YYYYMMDD').lower() 
        date_array.append(numdat)

    respath = '/data/tjarniko/avg/'
    trace_list = []
    for i in range(0,len(date_array)):
        
        Day_OI = date_array[i]
        print(Day_OI)
        tracers = 'SalishSea_1*' + Day_OI + '_grid_W.nc'
        
        ptt = respath + tracers
        trace = glob.glob(ptt)
        trace = (trace[0])
        print(trace)
        trace_list.append(trace)
    
    return trace_list


