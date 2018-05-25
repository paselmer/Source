; This main program compares the lidar and DEM surface altitudes using the 
; CAMAL NRB file as input.

dset = '20070809'
tit_tag = 'CPL 19 Oct 17'
filename = '/cpl3/ACEPOL-17/L1/NRB_ACEPOL-17_19oct17_cls.hdf5'
outdir = '/cpl3/ACEPOL-17/analysis/'

@/cpl3/CAMAL/Source/L1A/camal_hdf5_nrb_common
cpl_hdf5_nrb_reader, filename


; Only use 1 wavelength to save memory. NRB files are huge at raw rez.

nrb = transpose( nrb[*,*,1] )
nrb = nrb[*,0:821] ;if you go too far below 0.0 alt, you'll reach no-data zone
bin_alt_array = bin_alt_array[0:821]

; Make a curtain plot of the NRB
STOP
labeled_image,nrb,1e9,findgen(num_recs)/100.0,bin_alt_array,'recs','alt','t',20,20,1300,900,70,50,70,50,'layerZ', $
                   20.0,'/home/selmer/li.png',0,0,/roundxnumbers,/roundynumbers,/png
STOP

; Initialize some analysis variables

num_recs = num_recs[0]
surface_alt = fltarr(num_recs)-9999.9
dem_lid_dif = fltarr(num_recs)-9999.9
margin = 600.0 ; margin for surface altitude detection
gthreshfact = 6.5
gwidthinbins = 3
backgroundalt = -300.0
mingthresh = 1e3

; Loop through every NRB profile, scanning for lidar surface altitude

window,0,xsize=750,ysize=1000
colors=plot_colors()

for i = 0UL, num_recs-1 do begin
 
    ; This means the ocean
    if DEM_laserspot[i] eq -9999.0 then DEM_laserspot[i] = 0.0
 
    find_surface_alt,nrb,i,bin_alt_array,surface_alt,DEM_laserspot[i], $
        margin,gthreshfact,gwidthinbins,backgroundalt,Lsearch_limit,mingthresh
	
    dem_lid_dif[i] = DEM_laserspot[i] - surface_alt[i]	
    
    ;if (dem_lid_dif[i] gt 500.0) and (surface_alt[i] gt -9990.0) then begin
    ;    print,i,surface_alt[i],DEM_laserspot[i],DEM_laserspot[i]+margin
    ;    plot,nrb[i,*],bin_alt_array,yrange=[-200,18e3],xrange=[-0.5e9,2e9],charsize=1.5
    ;        oplot,[0,2e9],[surface_alt[i],surface_alt[i]],color=colors[1]
    ;    pause
    ;endif

endfor

; Save for possible later comparison.r analysis
laserspot_nav = laserspot
nav_nav = nav
ONA_nav = ONA

valid_locs = where(surface_alt gt -9999.0,nv)
if (nv gt 0) then begin

    print,'There are ',nv,' valid points.'
    surface_alt = surface_alt[valid_locs]
    dem_lid_dif = dem_lid_dif[valid_locs]
    DEM_laserspot = DEM_laserspot[valid_locs]
    roll = nav[valid_locs].roll
    ONA = ONA[valid_locs]

endif

; Print out stats of the comparison

avg = mean(dem_lid_dif)
med = median(dem_lid_dif)
mn  = min(dem_lid_dif)
mx  = max(dem_lid_dif)
sdv = stddev(dem_lid_dif)

print,'************************************************************'
print,'Here are the DEM-lidar comparison stats. All units meters...'
print,'Min:     ',mn
print,'Mean:    ',avg
print,'Median:  ',med
print,'Max:     ',mx
print,'Std Dev: ',sdv
print,'************************************************************'

window,0,xsize=1400,ysize=1105
!p.multi = [0, 1, 3]
chsz=3.0
clrind=7
plot,surface_alt,charsize=chsz,xstyle=1,title='Lidar (w) & DEM (g) surf. evel., '+tit_tag,ytitle='elevation (m)'
oplot,DEM_laserspot,color=colors[1]
plot,dem_lid_dif,charsize=chsz,psym=1,xstyle=1,title='DEM minus Lidar surf. elev.',ytitle='diff. (m)',ystyle=1,yrange=[-1000,1500]
xyouts,n_elements(dem_lid_dif)/2,1200,'Min:     '+strtrim(string(mn),2),charsize=1.5,color=colors[clrind]
xyouts,n_elements(dem_lid_dif)/2,1090,'Mean:    '+strtrim(string(avg),2),charsize=1.5,color=colors[clrind]
xyouts,n_elements(dem_lid_dif)/2,980, 'Median:  '+strtrim(string(med),2),charsize=1.5,color=colors[clrind]
xyouts,n_elements(dem_lid_dif)/2,870, 'Max:     '+strtrim(string(mx),2),charsize=1.5,color=colors[clrind]
xyouts,n_elements(dem_lid_dif)/2,760, 'Std Dev: '+strtrim(string(sdv),2),charsize=1.5,color=colors[clrind]
plot,ONA*(180.0/!dpi),charsize=chsz,xstyle=1,title='Off-nadir angle',ytitle='ONA (deg)',xtitle='valid records'
a=tvrd(true=1)
write_png,outdir+'dem_lid_surf_elev_summary_'+dset+'.png',a
STOP

!P.MULTI=0
make_histogram_plot,dem_lid_dif,30.0,'DEM, lidar surf. elev. diff. '+tit_tag,'meters'
a=tvrd(true=1)
write_png,outdir+'dem_lid_dif_dist_'+dset+'.png',a
STOP

;plot,ONA,dem_lid_dif,psym=1
;plot,DEM_laserspot,dem_lid_dif,psym=1
;plot,surface_alt,dem_lid_dif,psym=1

;----------------------------------------------------------------------------------------------
; BELOW HERE IS CODE TO COMPARE NRB FILES PRODUCED WITH DIFFERENT NAV SOURCES -----------------
; ---------------------------------------------------------------------------------------------
;----------------------------------------------------------------------------------------------
filename = '/cpl3/CAMAL/Data/CAMAL_17/'+dset+'/L1/NRB_CAMAL_17_'+dset+'_gps.hdf5'
camal_hdf5_nrb_reader, filename

; Plot headings
window,0
plot,nav_nav.heading,charsize=1.5;,xrange=[1.0e4,2.0e4]
oplot,nav.heading,color=colors[1]

window,1
; Plot laserspots
plot,laserspot_nav[0,*],charsize=1.5,yrange=[0.55,0.65];,xrange=[5.5e4,6.0e4]
oplot,laserspot[0,*],color=colors[1]

window,2
plot,ONA_nav,charsize=1.5
oplot,ONA,color=colors[1]

STOP
window,3
plot,laserspot_nav[0,*]*(180.0/!dpi),charsize=1.5,yrange=[30,40],xrange=[1.0e4,2.0e4]
oplot,laserspot[0,*]*(180.0/!dpi),color=colors[1]

plot,nav_nav.drift-nav.drift,psym=1
STOP
plot,nav_nav.roll-nav.roll,psym=1
STOP
plot,nav_nav.pitch-nav.pitch,psym=1
STOP
plot,nav_nav.heading-nav.heading,psym=1

r0 = 4.0e4
r1 = 4.1e4;num_recs-1
colors=plot_colors()
window,3,xsize=750,ysize=750
marg=0.1
map_set,/mercator,/continents,/grid,label=1, $ ;plot lat/lon on a map
    limit=[min(nav[r0:r1].lat)-marg,min(nav[r0:r1].lon)-marg,max(nav[r0:r1].lat)+marg,max(nav[r0:r1].lon)+marg],charsize=1.5
plots,nav[r0:r1].lon,nav[r0:r1].lat,color=colors[15],psym=5
plots,laserspot[1,r0:r1]*(180.0/!dpi),laserspot[0,r0:r1]*(180.0/!dpi),color=colors[1],psym=1
plots,laserspot_nav[1,r0:r1]*(180.0/!dpi),laserspot_nav[0,r0:r1]*(180.0/!dpi),color=colors[10],psym=1
;a=tvrd(true=1)
;write_png,outdir+'map_'+dset+'.png',a


end


