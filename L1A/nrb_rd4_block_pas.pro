pro NRB_RD4_BLOCK_PAS, num_vbins,nwl,eof_nrb,proftot,profcnt,frame_top, $
               jday_beg,jday_end,xx_sec,vsmo,hsmo, $
               fnum_nrb,djday,saturate,bground,energy,nrb_smo,depol_smo, $
               sdev_smo,gnd_ht,gnd_bin,nlay,lay_type,lay_top,lay_bot,  $
               lay_top_bin,lay_bot_bin,hr,minu,sec,heading,pitch,roll, $
               plnht,lat,lon,sun_el,sun_az,maxlay,nrb,depol_rat,nrb_sdev,nchan, $
	       exclude,jday1,jday2

; This subroutine reads the 4 channel nrb data file framed by jday_beg and
; jday_end, and will
; optionally smooth nrb data to specifications.  This version of NRB_RD4_SPECS
; is designed for reading one user-specified block of data.
;
; Patrick Selmer modified code so that all background values are now saved in [11/10/16]
; the variable "bground" as an array.
;
; Patrick Selmer modified code so that it now will invalidate all data between  [1/4/18]
; new input parameters, "jday1" and "jday2" if "exclude" is 1.
; This was first created for the ACEPOL campaign.

nchan= 4

; Define the size, type, and meaning of the nrb read parameters.
sectot= 0L              ;number of seconds in aggregate
secseq= 0L              ;sequence order of current second in aggregate
rawcnt= 0L              ;shotcounter of fifth (1-10) raw record in second
jday_dec= 0.0D          ;decimal julian day
satu_ht= fltarr(nchan)  ;height (m) where saturation first occurs per channel
                        ;(standard sequence: 355,532,1064par,1064per)
bkgrd= fltarr(nchan)    ;background count value per channel (standard seq)
emon= fltarr(nwl)       ;energy monitor per wavelength (355,532,1064)
reserv1= 0L             ;unused reserve
navrec= bytarr(256)     ;navigation data
frame_top= 0.0          ;top height (meters) of CPL frame
frame_bot= 0.0          ;bottom height (meters) of CPL frame
frame_vres= 0.0         ;vertical resolution of CPL frame
frame_nbins= 0L         ;number of vertical bins in CPL frame
nbs= fltarr(num_vbins,nchan) ;normalized relative backscatter profile per chan,
                             ;seq= 355,532,1064,unused
depol1064= fltarr(num_vbins) ;depolarization ratio profile for 1064 nm
nbs_sdev= fltarr(num_vbins,nchan) ;nrb standard deviation profile per wvlength,
                                  ;depol standard dev profile in last chan pos
grd_ht= 0.0             ;ground height location (meters amsl) {-999.0}
grdht_ind_er2= 0L       ;ground height location index (ER-2 reference) {-999}
grdht_ind_frame= 0L     ;ground height location index (CPL frame) {-999}
maxlay= 10L             ;maximum num. of layers allowed in layer location prog.
numlay= 0L              ;number of layers actually detected for current second
ltype= bytarr(maxlay)   ;type code of layer:  0= dummy
                                             ;1= PBL
                                             ;2= elevated aerosol
                                             ;3= cloud
                                             ;4= indeterminate
reserv2= bytarr(2)      ;unused reserve
lay_topht= fltarr(maxlay) ;all layer top heights found in current profile (m)
lay_botht= fltarr(maxlay) ;all layer bottom heights found in current profile (m)
top_ind_er2= fltarr(maxlay) ;all layer top bin ids from ER-2 reference
bot_ind_er2= fltarr(maxlay) ;all layer bottom bin ids from ER-2 reference
top_ind_frame= lonarr(maxlay) ;all layer top bin ids from CPL frame
bot_ind_frame= lonarr(maxlay) ;all layer bottom bin ids from CPL frame

print, 'Number of NRB profiles read= ', profcnt

; Read all records for current time block and store all necessary variables in
; memory.
  djday= dblarr(xx_sec)
  saturate= fltarr(nchan,xx_sec)
  bground= fltarr(4,xx_sec)  ;added 11/10/2016 by P.S.
  energy= 0
  nrb= fltarr(num_vbins,nchan,xx_sec)
  depol_rat= fltarr(num_vbins,xx_sec)
  ;nrb_sdev= fltarr(num_vbins,nchan,xx_sec)
  gnd_ht= fltarr(xx_sec)
  gnd_bin= 0
  nlay= lonarr(xx_sec)
  lay_type= intarr(maxlay,xx_sec)
  lay_top= fltarr(maxlay,xx_sec)
  lay_bot= fltarr(maxlay,xx_sec)
  lay_top_bin= 0
  lay_bot_bin= 0
  hr= intarr(xx_sec)
  minu= intarr(xx_sec)
  sec= intarr(xx_sec)
  heading= fltarr(xx_sec)
  pitch= fltarr(xx_sec)
  roll= fltarr(xx_sec)
  plnht= fltarr(xx_sec)
  lat= fltarr(xx_sec)
  lon= fltarr(xx_sec)
  sun_el= fltarr(xx_sec)
  sun_az= fltarr(xx_sec)
  nrb_sdev=0.0

  i= -1L

for itot= 0,xx_sec-1 do begin

; Call subroutine NRB_RD4 to read one long record of nrb data.
  NRB_RD4, eof_nrb,fnum_nrb,proftot,profcnt, $
         sectot,secseq,rawcnt,jday_dec,satu_ht,bkgrd,emon,reserv1,navrec, $
         frame_top,frame_bot,frame_vres,frame_nbins,nbs,depol1064,nbs_sdev, $
         grd_ht,grdht_ind_er2,grdht_ind_frame,maxlay,numlay,ltype, $
         reserv2,lay_topht,lay_botht,top_ind_er2,bot_ind_er2,top_ind_frame, $
         bot_ind_frame,nav_jday,nav_hr,nav_min,nav_sec,nav_heading,nav_pitch, $
         nav_roll,nav_gspd,nav_plnht_gps,nav_lat_gps,nav_lon_gps, $
         nav_sun_elv,nav_sun_az

  if (itot eq 0) then jday_dec0= jday_dec

  if (eof_nrb eq 1) then begin
    xx_sec= i+1L
    print, 'NRB File at EOF, xx_sec= ', xx_sec
    goto, fin_proc
  endif

; Test to see if profile is in correct time block.
  if (jday_dec lt jday_beg) then goto, readfile
  if (jday_dec gt jday_end+2) then goto, readfile
  if (jday_dec gt jday_end) then begin
    print, 'Current time beyond segment end...stopping read! ', jday_dec
	 goto, fin_proc
  endif

  ; This if block set data values to "invalids." This was first created
  ; to invalidate chunks of ACEPOL flights.
  ; If criteria for invalidation not met, values populated without
  ; inteference.
  if exclude and (jday_dec ge jday1) and (jday_dec le jday2) then begin
      i= i+1L
      djday(i)= jday_dec
      saturate(*,i)= satu_ht
      bground(*,i)= -999.999
      energy= emon
      nrb(*,*,i)= (nbs*0.0) -999.99
      depol_rat(*,i)= (depol1064*0.0)
      ;nrb_sdev(*,*,i)= nbs_sdev
      gnd_ht(i)= -999.9
      gnd_bin= -999
      nlay(i)= -999
      lay_type(*,i)= -999
      lay_top(*,i)= (lay_topht*0.0) - 999.9
      lay_bot(*,i)= (lay_botht*0.0) - 999.9
      lay_top_bin= (top_ind_frame*0.0) - 999
      lay_bot_bin= bot_ind_frame
      hr(i)= nav_hr
      minu(i)= nav_min
      sec(i)= nav_sec
      heading(i)= nav_heading
      pitch(i)= nav_pitch
      roll(i)= nav_roll
      plnht(i)= nav_plnht_gps
      lat(i)= nav_lat_gps
      lon(i)= nav_lon_gps
      sun_el(i)= nav_sun_elv
      sun_az(i)= nav_sun_az
  endif else begin
      i= i+1L
      djday(i)= jday_dec
      saturate(*,i)= satu_ht
      bground(*,i)= bkgrd
      energy= emon
      nrb(*,*,i)= nbs
      depol_rat(*,i)= depol1064
      ;nrb_sdev(*,*,i)= nbs_sdev
      gnd_ht(i)= grd_ht
      gnd_bin= grdht_ind_frame
      nlay(i)= numlay
      lay_type(*,i)= fix(ltype)
      lay_top(*,i)= lay_topht
      lay_bot(*,i)= lay_botht
      lay_top_bin= top_ind_frame
      lay_bot_bin= bot_ind_frame
      hr(i)= nav_hr
      minu(i)= nav_min
      sec(i)= nav_sec
      heading(i)= nav_heading
      pitch(i)= nav_pitch
      roll(i)= nav_roll
      plnht(i)= nav_plnht_gps
      lat(i)= nav_lat_gps
      lon(i)= nav_lon_gps
      sun_el(i)= nav_sun_elv
      sun_az(i)= nav_sun_az
  endelse  

  readfile:
endfor

fin_proc:
xx_sec= i+1L
print, 'Number of profiles in segment= ', xx_sec
;print, 'First jday_dec read= ', jday_dec0
;print, 'Last jday_dec read= ', jday_dec
if (xx_sec le 0L) then goto, fin_rout

; Re-set the dimensions of the segment arrays in case xx_sec has changed.
djday= djday(0:xx_sec-1)
saturate= saturate(*,0:xx_sec-1)
;bground= bground(*,0:xx_sec-1)
;energy= energy(*,0:xx_sec-1)
nrb= nrb(*,*,0:xx_sec-1)
depol_rat= depol_rat(*,0:xx_sec-1)
;nrb_sdev= nrb_sdev(*,*,0:xx_sec-1)
gnd_ht= gnd_ht(0:xx_sec-1)
;gnd_bin= gnd_bin(0:xx_sec-1)
nlay= nlay(0:xx_sec-1)
lay_type= lay_type(*,0:xx_sec-1)
lay_top= lay_top(*,0:xx_sec-1)
lay_bot= lay_bot(*,0:xx_sec-1)
;lay_top_bin= lay_top_bin(*,0:xx_sec-1)
;lay_bot_bin= lay_bot_bin(*,0:xx_sec-1)
hr= hr(0:xx_sec-1)
minu= minu(0:xx_sec-1)
sec= sec(0:xx_sec-1)
heading= heading(0:xx_sec-1)
pitch= pitch(0:xx_sec-1)
roll= roll(0:xx_sec-1)
plnht= plnht(0:xx_sec-1)
lat= lat(0:xx_sec-1)
lon= lon(0:xx_sec-1)
sun_el= sun_el(0:xx_sec-1)
sun_az= sun_az(0:xx_sec-1)

  nrb_smo=nrb
  nrb=0
  depol_smo=depol_rat
  depol_rat=0
  sdev_smo=0.0

; Call Subroutine HORIZ_AVG to do a running average of the data with time as
; specified by hsmo.
;sdev_flg= 1   ; Do perform sdev and depol smooth analysis
;HORIZ_AVG, xx_sec,hsmo,nrb,depol_rat,nrb_sdev,nrb_h,depol_h,sdev_h, $
;           num_vbins,nchan,sdev_flg

; Call Subroutine VERTI_AVG to do a running average of the data with height as
; specified by vsmo.
;VERTI_AVG, num_vbins,xx_sec,vsmo,nrb_h,depol_h,sdev_h,nrb_smo,depol_smo, $
;           sdev_smo,nchan,sdev_flg

; Free up unused memory.
;nrb_h=0 & depol_h=0 & sdev_h=0

fin_rout:
return
end


