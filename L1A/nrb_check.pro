; Main program to play around with L1A files and make sure data look okay.
; Basically, for testing purposes. Don't take anything here too seriously.

filename = '/cpl3/CAMAL/Data/CAMAL_17/20171206/L1/NRB_CAMAL_17_20171206_nav.hdf5'
dset = '20171206'
outdir = '/cpl3/CAMAL/Analysis/CAMAL_17/'+dset+'/'

@/cpl3/CAMAL/Source/L1A/camal_hdf5_nrb_common
camal_hdf5_nrb_reader, filename
STOP


;**************** Convert laserspot lat/lon from radians to degrees ****************

laserspot = laserspot * (180.0/!dpi)


;**************** Plot the lat/lon on a map *******************

r0 = 0
r1 = num_recs-1
colors=plot_colors()
window,3,xsize=750,ysize=750
marg=0.1
map_set,/mercator,/continents,/grid,label=1, $ ;plot lat/lon on a map
    limit=[min(nav[r0:r1].lat)-marg,min(nav[r0:r1].lon)-marg,max(nav[r0:r1].lat)+marg,max(nav[r0:r1].lon)+marg],charsize=1.5
plots,nav[r0:r1].lon,nav[r0:r1].lat,color=colors[15],psym=5
plots,laserspot[1,r0:r1],laserspot[0,r0:r1],color=colors[1],psym=1
;a=tvrd(true=1)
;write_png,outdir+'map_'+dset+'.png',a

end
