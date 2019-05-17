
; >>>>>>>>>>>>>>>>>>>>>>>>>>>>> LOAD NEW >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
proj_dir = '/cpl3/Seac4rs/'
wlstr = '532'
datestr = '16aug13'
sortie='13954'
tit_tag = 'CPL '+datestr+' '+wlstr
filename = proj_dir+'L1/NRB_Seac4rs_'+datestr+'_iwg1.hdf5'
outdir = '/cpl3/Seac4rs/analysis/'
alt1 = -0.5e3
alt2 = 20e3
r1 = 0
r2 = -99 ;-99 for all
maxscale =1e9
!PATH = '/cpl/dhlavka/Cpl/Source/:' + !PATH

@/cpl3/CAMAL/Source/L1A/camal_hdf5_nrb_common
cpl_hdf5_nrb_reader, filename
if r2 eq 11000 then r2 = num_recs[0]-1

; Only use 1 wavelength to save memory. NRB files are huge at raw rez.
if wlstr eq '355' then nrb = transpose( nrb[*,r1:r2,0] )
if wlstr eq '532' then nrb = transpose( nrb[*,r1:r2,1] )
if wlstr eq '1064' then nrb = transpose( nrb[*,r1:r2,2]+nrb[*,r1:r2,3] )
vlocs = where((bin_alt_array ge alt1) and (bin_alt_array le alt2))
nrb = nrb[*,vlocs] ;if you go too far below 0.0 alt, you'll reach no-data zone
bin_alt_array = bin_alt_array[vlocs]
; Save for avg'd prof plot comparison
newnrb = nrb
newz = bin_alt_array

; Make a curtain plot of the NRB

labeled_image,nrb,maxscale,findgen(n_elements(nrb[*,0]))/100.0,bin_alt_array,'recs','alt (m)',tit_tag+' - new NRB (Python, HDF5)',20,20,1300,900,70,50,70,50,'layerZ', $
                   20.0,outdir+'new_'+wlstr+'_NRB_'+datestr+'.png',0,0,/roundxnumbers,/roundynumbers,/png

STOP

; >>>>>>>>>>>>>>>>>>>>>>>>>>>>> LOAD OLD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

NRB_file = proj_dir+'Analy_complex/NRB_'+sortie+'_'+datestr+'.xdr'

fnum_nrb= 10
fnum_cal= 11
openr, fnum_nrb, NRB_file, /XDR

nwl= 4
eof_nrb= 0           ; end of file flag for nrb file
vsmo= 1   ; Number of vertical bins to average to smooth lidar signals (1,3,5)
hsmo= 1   ; Number of horizontal bins to average to smooth lidar signals (1,3,5)

; Point to the end of NRB_FIL and read the total number of profiles contained
; in the file, then re-position to the beginning.
nrb_fstat= fstat(fnum_nrb)
proftot= -1L
point_lun, fnum_nrb, nrb_fstat.size-4L
readu, fnum_nrb, proftot
point_lun, fnum_nrb, 0
print, 'Total number of profiles in NRB file= ', proftot

; Read the first short record of NRB_FIL (start and end annulus angles per chan)
; to determine whether nrb data has 4 or 10 channels (76 bytes).
sortie= 0L & jday= 0L & year= 0L & num_vbins= 0L & vres= 0.0
nwl= 0L & met_assign= 0L & polrat= 0.0
ann_angle_beg= fltarr(5)
ann_angle_end= fltarr(5)
ann_angle_bye= 0.0
readu, fnum_nrb, sortie,jday,year,num_vbins,vres,nwl,met_assign, $
                 ann_angle_beg,ann_angle_end,ann_angle_bye,polrat
print, 'Short rec= ', sortie,jday,year,vres,met_assign,ann_angle_beg(2),polrat

; Develop begining and ending time stamps for the period of interest.
fhr= -1L & fmin= -1L & fsec= -1L
print, "Day of year at file start= ", jday
sjday= jday & shr= -1L & smin= -1L & ssec= -1L
ejday= jday & ehr= 99L & emin= 99L & esec= 99L
if (shr lt fhr) then sjday= jday+1L
if (ehr lt fhr) then ejday= jday+1L
shr2= shr
if (sjday ne jday) then shr2= shr+24L
stime= shr2*3600L + smin*60L + ssec
ehr2= ehr
if (ejday ne jday) then ehr2= ehr+24L
etime= ehr2*3600L + emin*60L + esec
ftime= fhr*3600L + fmin*60L + fsec
;  Modify following paragraph appropiately.
sec_adv= 0L    ; for flight segments with large or several laser off periods (default)
if (sec_adv lt 0L) then sec_adv= 0L
xx_sec= 25000L ; for flight segments with large or several laser off periods (default)

month= -1
day= -1
jday_beg= DJUL_DAY(year,month,day,shr,smin,ssec,fix(sjday))
jday_end= DJUL_DAY(year,month,day,ehr,emin,esec,fix(ejday))
print, 'User-input start and end jday for period= ', jday_beg,jday_end

if (ann_angle_beg(2) eq 0.0 and ann_angle_end(2) eq 0.0) then begin
  byt_adv= sec_adv*33248L + 76L       ; 4 CHAN record
endif else begin
  byt_adv= sec_adv*76496L + 76L       ; 10 CHAN record
endelse
point_lun, fnum_nrb, byt_adv          ; Data file is now pointed to at or just
                                      ; before period of interest.

profcnt= sec_adv

; Develop height array for for CPL frame (meters).
ht_m= fltarr(num_vbins)
for bin= 0,num_vbins-1 do begin
  ht_m(bin)= bin_ht(bin,vres)
endfor


; Read 4 channel nrb data file framed by jday_beg and jday_end and optionally
; smooth nrb data.  This is the original configuration.
print, 'Reading 4-channel data.'
NRB_RD4_BLOCK_PAS, num_vbins,nwl,eof_nrb,proftot,profcnt,frame_top, $  
           jday_beg,jday_end,xx_sec,vsmo,hsmo, $
           fnum_nrb,djday,saturate,bground,energy,nrb_smo,depol_smo, $
           sdev_smo,gnd_ht,gnd_bin,nlay,lay_type,lay_top,lay_bot,  $
           lay_top_bin,lay_bot_bin,hr,minu,sec,heading,pitch,roll, $
           plnht,lat,lon,sun_el,sun_az,maxlay,nrb,depol_rat,nrb_sdev,nchan, $
           0,jday_beg,jday_end

if wlstr eq '355' then nrb = transpose( reform( nrb_smo[*,0,r1:r2] ) )
if wlstr eq '532' then nrb = transpose( reform( nrb_smo[*,1,r1:r2] ) )
if wlstr eq '1064' then nrb = transpose( reform( nrb_smo[*,2,r1:r2] ) )
vlocs = where((ht_m ge alt1) and (ht_m le alt2))
nrb = nrb[*,vlocs]
ht_m = ht_m[vlocs]
; Save for avg'd prof plot comparison
oldnrb = nrb
oldz = ht_m

labeled_image,nrb,maxscale,findgen(n_elements(nrb[*,0]))/100.0,ht_m,'recs','alt (m)',tit_tag+' - old NRB (IDL, XDR)',20,20,1300,900,70,50,70,50,'layerZ', $
                   20.0,outdir+'old_'+wlstr+'_NRB_'+datestr+'.png',0,0,/roundxnumbers,/roundynumbers,/png



; Compare average profile of new vs. old
new = mean(newnrb,dimension=1)
old = mean(oldnrb,dimension=1)
window,2,xsize=800,ysize=1000
colors=plot_colors()
plot,new,newz,color=0,background=colors[160],charsize=2.0,thick=2.0,$
    title=datestr,xtitle='NRB',ytitle='alt (m)',yrange=[alt1,alt2],ystyle=1, $
    xrange=[0,1e8],xstyle=1
oplot,old,oldz,color=colors[3],thick=2.0
a=tvrd()
write_png,'avgd_prof_plot.png',a

end
