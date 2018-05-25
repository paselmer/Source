PRO camal_hdf5_nrb_reader, filename

@/cpl3/CAMAL/Source/L1A/camal_hdf5_nrb_common

; Get the file ID
file_id = H5F_OPEN(filename)

; Read in the meta data (#records, structure)
dataset_id = H5D_OPEN(file_id,'meta')
meta = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,meta

; Read in the PGain ()
dataset_id = H5D_OPEN(file_id,'PGain')
PGain = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,PGain

; Read in the nav data (#records, structure)
dataset_id = H5D_OPEN(file_id,'nav')
nav = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,nav

; Read in the laserspot data ( #records x 2[lat,lon] )
dataset_id = H5D_OPEN(file_id,'laserspot')
laserspot = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,laserspot

; Read in the bin altitude data ( #records x 2[lat,lon] )
dataset_id = H5D_OPEN(file_id,'bin_alt_array')
bin_alt_array = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,bin_alt_array

; Read in the number of fixed frame bins ( )
dataset_id = H5D_OPEN(file_id,'num_ff_bins')
num_ff_bins = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,num_ff_bins

; Read in the number of records (  )
dataset_id = H5D_OPEN(file_id,'num_recs')
num_recs = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,num_recs

; Read in the off-nadir angle ( #records, radians )
dataset_id = H5D_OPEN(file_id,'ONA')
ONA = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,ONA

; Read in the DEM elevation at nadir ( #records, meters )
dataset_id = H5D_OPEN(file_id,'DEM_nadir')
DEM_nadir = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,DEM_nadir

; Read in the DEM surface type at nadir ( #records, meters )
dataset_id = H5D_OPEN(file_id,'DEM_nadir_surftype')
DEM_nadir_surftype = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,DEM_nadir_surftype

; Read in the DEM elevation at laserspot ( #records, meters )
dataset_id = H5D_OPEN(file_id,'DEM_laserspot')
DEM_laserspot = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,DEM_laserspot

; Read in the DEM surface type at laserspot ( #records, meters )
dataset_id = H5D_OPEN(file_id,'DEM_laserspot_surftype')
DEM_laserspot_surftype = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,DEM_laserspot_surftype

; Read in the energy monitor data ( #wavelengths[355,532,1064] x #records )
dataset_id = H5D_OPEN(file_id, 'EM')
EM = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,EM

; Read in the NRB data (#chans x #records x #bins, float64 array)
dataset_id = H5D_OPEN(file_id,'nrb')
nrb = H5D_READ(dataset_id)
H5D_CLOSE, dataset_id
help,nrb

end
