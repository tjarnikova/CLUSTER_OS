code repository for figures and tables in Jarnikova et al 2022
https://os.copernicus.org/preprints/os-2021-66/

Each figure is in its own notebook, format _FF01_map.ipynb or similar. 
Each figure has a load_data section that should rely on the folder CP_DATA


Workflow for clustering:
- 1) extract signals for clustering using the 5 functions in extract_sigs. to recreate this on a given filesystem would need to change paths. 
        - for the this purposes of the paper repository the data is as accessed as .ncs  in ./CP_DATA/NC_HINDCAST
        - the DATA folder is provided at the repository xx
        - on salish it's originally stored: f'/data/tjarniko/MEOPAR/analysistereza/notebooks/CLUSTER_201905/NC_HINDCAST/{year}/BIO_TS/'
- 2) extracted signals are then used in the figures. 

Data requirements 
Fig 1 (map): 
    - 1 day of model output, on salish cluster xx-not rerun or available 
    - bathy and mesh files are in ./CP_DATA
    
Fig 2 (example signals):
    - extracted yearly station signals are in ./CP_DATA/NC_HINDCAST 
    
Fig 3 (