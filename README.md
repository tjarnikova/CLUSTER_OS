code repository for figures and tables in Jarnikova et al 2022
https://os.copernicus.org/preprints/os-2021-66/

> Each figure is in its own notebook, format _FF01_map.ipynb or similar. 
> Each figure has a load_data section that should rely only on the folder CP_DATA, which holds the extracted signals
> The extraction codes shown in ./extract_sigs were run on salish at the University of British Columbia, to recreate this on a different filesystem would need to change paths, the code is provided for transparency
> Extracted signals (at UBC) are at:
f'/data/tjarniko/MEOPAR/analysistereza/notebooks/CLUSTER_201905/NC_HINDCAST/{year}/BIO_TS/'
> They are also copied to CP_DATA/NC_HINDCAST, which is provided in the code repository along with model code: xx FRDR 

Workflow for clustering:
- 1) EXTRACT > ./extract_sigs
(see above)

- 2) AGGLOMERATE, CLUSTER, VISUALIZE > eg _FF03
extracted signals are agglomerated and clustered in each of the figure notebooks (eg _FF03 for freshwater)


Data requirements 
Fig 1 (map): 
    - 1 day of model output, on salish cluster xx-not rerun or available 
    - bathy and mesh files are in ./CP_DATA
    
Fig 2 (example signals):
    - extracted yearly station signals are in ./CP_DATA/NC_HINDCAST/
    
Fig 3,4 (FWI, HALO):
    - extracted yearly station signals are in ./CP_DATA/NC_HINDCAST/
    - the datamats, linkmats, signalmats are in ./CP_DATA/
