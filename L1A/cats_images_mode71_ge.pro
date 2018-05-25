;;;;;;;;;;;;;;;;; cats_images_mode71_ge.pro ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;       I INTEND FOR THIS TO BE A USEFUL DESCRIPTION (10/20/16)
;
;	Reads in CATS L1B Data for Mode 7.1 and creates all neccessary 
;       files to plot curtain plots in Google Earth. These files are
;       saved in a new, appropriately named directory. The following
;       four files are produced:
;            - A PNG curtain plot image
;            - A KML file
;            - A COLLADA Model (.DAE) file
;            - A PNG image showing the ground track
;            - A KML file that plots the points of 
;              the ground track
;
;       The zmult variable controls the height of the image displayed
;       in Google Earth.
;
;       How's it done?
;
;       A KML file references a COLLADA model file, which have DAE extensions.
;       The COLLADA model file describes how a PNG image of the data should
;       be displayed in Google Earth. The most important aspect of the COLLADA
;       file to understand is the 3 sets of numbers.
;
;       The first set falls under a <geometery> block. It describes the coordinates
;       of the curtain in order of x,y,z. The curtain's coordinates are "model"
;       coordinates, with the origin being 0 lon,0 lat,0 alt. All the numbers
;       are simply distances from the origin. The PNG image is a series of rows
;       and columns of data. As of 10/20/16, the COLLADA file is hard-coded to
;       break the image into 274 columns (horizontal direction) and 30 rows
;       (vertical direction). Every 274 lines of model coordinate triplets 
;       represent 274 columns of the image, proceeding from the first through
;       last rows.
;
;       The second set of numbers is under a <source> block and describes
;       the normalized coordinates (x,y) of the PNG image with the origin
;       being at the bottom left. Each line corresponds to a line of the
;       model coordinates (the first set of numbers).
;
;       The third and final set of numbers, which is under a <triangles> block
;       is a bit trickier to explain. Every 3 numbers in a row describes a triangle.
;       There are 12 numbers in a row therefore each row describes 4 triangle.
;       Imagine four points on the image: A - bottom left, B - bottom right,
;       C - top left, D - top right. The four triangles from left to right in
;       each row  are ABD, DBA, ADC, CDA. The pattern repeats. The numbers themselves
;       should always be integers and are indexes of the "arrays" created by
;       the first two sets of numbers.
;
;       A KML file must reference this model in order to display it in Google
;       Earth.
;       
;
;
;	Inputs:
;		Via config file
;
;	Notes:
;
;		Pressure units in Pa
;		Temperature units in K
;		All Height units in output structures are in km asl
;		Extinction units in km-1
;		Backscatter units in (km sr)-1
;
;	Version:
;
;		Version ?. Last edited by Patrick Selmer Nov. 28, 2016
;					
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO cats_images_mode71_ge,config_path

!PATH= '/Data/Code/Dev/dhlavka:' + !PATH		        
@/Data/Code/Dev/dhlavka/cats_l1b_v2_06_parms_common    ;common for L1B data
opath = '/Data/Code/Dev/selmer/testing_images/' ;files will be stored under custom dir under opath

depolbscutoff = 3.0e-4
;zmult = 22.0 ;commented out. This is a decent value to use.

ae=6378.137d0;Wikipedia world geodetic system #84;WGS84;See 05Oct11 [Equatorial radius of the Earth]
b=6356.752d0;See GLAS atbd Laser footprint location(geolocation) and surface profiles.(Schutz,2002)
f=(ae-b)/ae
esquared=2.0d0*f-f^2
eccentricity=sqrt(esquared)
dre=1.0d0*6371.22d0;r(km) earth volume equivalent sphere


;******************************* Read in data from config file *****************************************
dum=''
openr,clun,config_path+'cats_images_mode71_ge.config',/get_lun
readf,clun,dum ;skip first header line
readf,clun,dum ;skip second header line
readf,clun,dum
strings = strsplit(dum,' ',/extract)
file_path=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
file_name=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
dialog_path=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
opath=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
absmax=strings[2]
absmax=float(absmax)
readf,clun,dum
strings = strsplit(dum,' ',/extract)
xtotsize=strings[2]
xtotsize=long(xtotsize)
readf,clun,dum
strings = strsplit(dum,' ',/extract)
ytotsize=strings[2]
ytotsize=long(ytotsize)
readf,clun,dum
strings = strsplit(dum,' ',/extract)
tit1=strings[2:n_elements(strings)-1]
tit1=strjoin(tit1,' ')
readf,clun,dum
strings = strsplit(dum,' ',/extract)
tit2=strings[2:n_elements(strings)-1]
tit2=strjoin(tit2,' ')
readf,clun,dum
strings = strsplit(dum,' ',/extract)
num_x_axes=strings[2]
num_x_axes=fix(num_x_axes)
readf,clun,dum
strings = strsplit(dum,' ',/extract)
wavelength=strings[2]
wavelength=long(wavelength)
readf,clun,dum
strings = strsplit(dum,' ',/extract)
rec_sum=strings[2]
rec_sum=long(rec_sum)
readf,clun,dum
strings = strsplit(dum,' ',/extract)
timelimit1=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
timelimit2=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
alt1=float(strings[2])
readf,clun,dum
strings = strsplit(dum,' ',/extract)
alt2=float(strings[2])
readf,clun,dum
strings = strsplit(dum,' ',/extract)
ifov=strings[2]
readf,clun,dum
strings = strsplit(dum,' ',/extract)
zmult=float(strings[2])


;******************************* Read in and Manipulate CATS Data *****************************************

l1b_path = '/Data/Dev/Std/L1/L1b/Products/'  ;;'/Data/Dev/Std/L1/L1b/Products/processed/'
fnlen = 61;length of filename

if (file_path ne 'x') then begin
    infile = file_path + file_name
endif else begin
    infile = dialog_pickfile(path=dialog_path)
endelse


gname = strmid(infile,fnlen,57,/reverse_offset) ;granule name w/no path
type='FFOV'
plotdate='11 Feb. 2015'
rawplotdate=strmid(gname,26,10)
subs=strsplit(rawplotdate,'-',/extract)
yr=subs[0]
dy=subs[2]
months = ['January','February','March','April','May','June','July','August','September','October','November','December']
CASE subs[1] OF
  '01':mn=months[0]
  '02':mn=months[1]
  '03':mn=months[2]
  '04':mn=months[3]
  '05':mn=months[4]
  '06':mn=months[5]
  '07':mn=months[6]
  '08':mn=months[7]
  '09':mn=months[8]
  '10':mn=months[9]
  '11':mn=months[10]
  '12':mn=months[11]
ENDCASE
plotdate = string(dy)+' '+mn+' '+string(yr)

;Read in L1B Data
CATS_L1B_HDF5_v2_06_READER, infile			;Sub-Routine to read in L1B data

numrecs=Number_Profiles
;Compute hr, min, sec for each profile
hr=intarr(numrecs)
minu=intarr(numrecs)
secs=intarr(numrecs)
timestrarr = strarr(numrecs)
for r = 0,numrecs-1 do begin
  DJUL_DAY_INV_V2, Profile_UTC_Time(r),hour,min,sec
  hr(r)=hour
  minu(r)=min
  secs(r)=sec
  shr=strmid(hr(r),6,2)
  if hr(r) lt 10 then shr='0'+strmid(shr,1,1)
  smin=strmid(minu(r),6,2)
  if minu(r) lt 10 then smin='0'+strmid(smin,1,1)
  ssec=strmid(secs(r),6,2)
  if secs(r) lt 10 then ssec='0'+strmid(ssec,1,1)
  timestrarr[r] = shr+':'+smin+':'+ssec
endfor

;Subset arrays here
srec=where(timestrarr eq timelimit1,ntl1)
if (ntl1 eq 0) then begin
    srec=0
endif else begin
    srec=srec[0]
endelse
erec=where(timestrarr eq timelimit2,ntl2)
if (ntl2 eq 0) then begin
    erec=numrecs-1
endif else begin
    erec=erec[0]
endelse
altitude = bin_altitude_array
zlocs = where( (altitude lt alt2) and (altitude gt alt1) , nz1)
altitude = altitude[zlocs]*zmult

shr=strmid(hr(srec),6,2)
if hr(srec) lt 10 then shr='0'+strmid(shr,1,1)
smin=strmid(minu(srec),6,2)
if minu(srec) lt 10 then smin='0'+strmid(smin,1,1)
ssec=strmid(secs(srec),6,2)
if secs(srec) lt 10 then ssec='0'+strmid(ssec,1,1)
ehr=strmid(hr(erec),6,2)
if hr(erec) lt 10 then ehr='0'+strmid(ehr,1,1)
emin=strmid(minu(erec),6,2)
if minu(erec) lt 10 then emin='0'+strmid(emin,1,1)
esec=strmid(secs(erec),6,2)
if secs(erec) lt 10 then esec='0'+strmid(esec,1,1)
timerange=shr+':'+smin+':'+ssec+' - '+ehr+':'+emin+':'+esec+' UTC'
starttimelabel = '_' + strjoin( strsplit(strmid(timerange,0,8),':',/extract),'-')
createtimelabel = '_' + strjoin( strsplit(strmid(systime(),11,8),':',/extract),'-') 
; Define custom output directory for all the output files.
; Will be under opath
custom_dir = opath+'Google_Earth'+'_'+rawplotdate+starttimelabel+createtimelabel+'/'
command = 'mkdir '+custom_dir
spawn,command

nprof=lindgen(erec-srec+1)+srec
nprof1=[srec,erec]
;latitude=[ISS_Latitude(srec),ISS_Latitude(erec)]
;lat = interpol(Latitude,nprof1,nprof)
if (ifov eq 'r') then begin
  lat = CATS_Right_FOV_Latitude[srec:erec] ;ISS_Latitude[srec:erec]
  lon = CATS_Right_FOV_Longitude[srec:erec];ISS_Longitude[srec:erec]
endif else begin
  lat = CATS_Left_FOV_Latitude[srec:erec] ;ISS_Latitude[srec:erec]
  lon = CATS_Left_FOV_Longitude[srec:erec];ISS_Longitude[srec:erec]
endelse
timestrarr = timestrarr[srec:erec]

nrecs=erec-srec+1
nbin=nz1 ;Number_Bins
cats_atb_532_rfov = transpose(Total_Attenuated_Backscatter532_Right_FOV[zlocs,srec:erec])
cats_atb_1064_rfov = transpose(Total_Attenuated_Backscatter1064_Right_FOV[zlocs,srec:erec])
cats_atb_1064perp_rfov = transpose(Perpendicular_Attenuated_Backscatter1064_Right_FOV[zlocs,srec:erec])
cats_atb_1064par_rfov = transpose(Total_Attenuated_Backscatter1064_Right_FOV[zlocs,srec:erec]-Perpendicular_Attenuated_Backscatter1064_Right_FOV[zlocs,srec:erec])
cats_atb_532_lfov = transpose(Total_Attenuated_Backscatter532_Left_FOV[zlocs,srec:erec])
cats_atb_1064_lfov = transpose(Total_Attenuated_Backscatter1064_Left_FOV[zlocs,srec:erec])
cats_atb_1064perp_lfov = transpose(Perpendicular_Attenuated_Backscatter1064_Left_FOV[zlocs,srec:erec])
cats_atb_1064par_lfov = transpose(Total_Attenuated_Backscatter1064_Left_FOV[zlocs,srec:erec]-Perpendicular_Attenuated_Backscatter1064_Left_FOV[zlocs,srec:erec])

avrecs=nrecs/rec_sum[0]
cats_bks_532_rfov=dblarr(avrecs,nbin)
cats_bks_1064_rfov=dblarr(avrecs,nbin)
cats_depol_532_rfov=dblarr(avrecs,nbin)
cats_depol_1064_rfov=dblarr(avrecs,nbin)
cats_bks_532_lfov=dblarr(avrecs,nbin)
cats_bks_1064_lfov=dblarr(avrecs,nbin)
cats_depol_532_lfov=dblarr(avrecs,nbin)
cats_depol_1064_lfov=dblarr(avrecs,nbin)
data = dblarr(avrecs,nbin) ;array to pass to plot_image program
temp = dblarr(nz1)
cratiorfov = dblarr(avrecs,nbin) 
cratiolfov = dblarr(avrecs,nbin)
latitude=fltarr(avrecs)
longitude=fltarr(avrecs)
avgtimestrarr=strarr(avrecs)
CTRS = fltarr(avrecs,nbin,3) ;[xc,yc,zc]
modelx = fltarr(avrecs,nbin) ;;
modely = fltarr(avrecs,nbin)
modelz = fltarr(avrecs,nbin)
for r=0,avrecs-1 do begin
   begrec=float(r)*float(rec_sum)
   endrec=begrec+float(rec_sum)-1.
   if r eq avrecs-1 then begin
      latitude(r)=mean(lat(begrec:nrecs-1))
      longitude(r)=mean(lon(begrec:nrecs-1))
      avgtimestrarr(r)=timestrarr(nrecs-1)
      for b=0,nbin-1 do begin
	cats_bks_532_rfov(r,b)=mean(cats_atb_532_rfov(begrec:nrecs-1,b))
	cats_bks_1064_rfov(r,b)=mean(cats_atb_1064_rfov(begrec:nrecs-1,b))
	cats_bks_1064perp_rfov=mean(cats_atb_1064perp_rfov(begrec:nrecs-1,b))
	cats_bks_1064par_rfov=mean(cats_atb_1064par_rfov(begrec:nrecs-1,b))
	cats_bks_532_lfov(r,b)=mean(cats_atb_532_lfov(begrec:nrecs-1,b))
	cats_bks_1064_lfov(r,b)=mean(cats_atb_1064_lfov(begrec:nrecs-1,b))
	cats_bks_1064perp_lfov=mean(cats_atb_1064perp_lfov(begrec:nrecs-1,b))
	cats_bks_1064par_lfov=mean(cats_atb_1064par_lfov(begrec:nrecs-1,b))
	if cats_bks_532_rfov(r,b) ge 0.002 then begin	 	 
	  cats_depol_1064_rfov(r,b)=cats_bks_1064perp_rfov / cats_bks_1064par_rfov
	  cratiorfov(r,b) = cats_bks_1064_rfov(r,b) / cats_bks_532_rfov(r,b)
	endif 
	if cats_bks_532_lfov(r,b) ge 0.002 then begin	 	 
	  cats_depol_1064_lfov(r,b)=cats_bks_1064perp_lfov / cats_bks_1064par_lfov
	  cratiolfov(r,b) = cats_bks_1064_lfov(r,b) / cats_bks_532_lfov(r,b)
	endif 
      endfor
    endif else begin
    latitude(r) = mean(lat(begrec:endrec))
    longitude(r) = mean(lon(begrec:endrec))
    avgtimestrarr(r) = timestrarr( (begrec+endrec)/2 )
    for b=0,nbin-1 do begin
          ; Convert [lon,lat,alt] to CTRS coordinates
      CTRS[r,b,*] = xyzg(ae,esquared,altitude[b],latitude[r],longitude[r]) ;returns in km
      transform_coordinates,CTRS[r,b,0]-ae,CTRS[r,b,1],CTRS[r,b,2],90.0,90.0,0.0,phi,xm,ym,zm
      modelx[r,b] = xm ;should hopefully be in model coordinates
      modely[r,b] = ym ;should hopefully be in model coordinates
      modelz[r,b] = zm ;should hopefully be in model coordinates  
      cats_bks_532_rfov(r,b)=mean(cats_atb_532_rfov(begrec:endrec,b))
      cats_bks_1064_rfov(r,b)=mean(cats_atb_1064_rfov(begrec:endrec,b))
      cats_bks_1064perp_rfov=mean(cats_atb_1064perp_rfov(begrec:endrec,b))
      cats_bks_1064par_rfov=mean(cats_atb_1064par_rfov(begrec:endrec,b))
      cats_bks_532_lfov(r,b)=mean(cats_atb_532_lfov(begrec:endrec,b))
      cats_bks_1064_lfov(r,b)=mean(cats_atb_1064_lfov(begrec:endrec,b))
      cats_bks_1064perp_lfov=mean(cats_atb_1064perp_lfov(begrec:endrec,b))
      cats_bks_1064par_lfov=mean(cats_atb_1064par_lfov(begrec:endrec,b))
      if cats_bks_532_rfov(r,b) ge 0.002 then begin
         cats_depol_1064_rfov(r,b)=(cats_bks_1064perp_rfov)/ cats_bks_1064par_rfov
         cratiorfov(r,b) = cats_bks_1064_rfov(r,b) / cats_bks_532_rfov(r,b)	 
      endif
      if cats_bks_532_lfov(r,b) ge 0.002 then begin
         cats_depol_1064_lfov(r,b)=(cats_bks_1064perp_lfov)/ cats_bks_1064par_lfov
         cratiolfov(r,b) = cats_bks_1064_lfov(r,b) / cats_bks_532_lfov(r,b)	 
      endif
    endfor
  endelse
  IF (ifov eq 'r') THEN BEGIN
  CASE wavelength[0] OF
       1064: temp[0l:nz1-1l] = cats_bks_1064_rfov[r,*] ;!!!!!! SET CURTAINPLOT VARIABLE HERE <<<<<<<<<<<<<<<<<<
        532: temp[0l:nz1-1l] = cats_bks_532_rfov[r,*]
          0: temp[0l:nz1-1l] = cats_depol_1064_rfov[r,*]
	  1: temp[0l:nz1-1l] = cats_depol_532_rfov[r,*]
  ENDCASE
  ENDIF ELSE BEGIN
  CASE wavelength[0] OF
       1064: temp[0l:nz1-1l] = cats_bks_1064_lfov[r,*] ;!!!!!! SET CURTAINPLOT VARIABLE HERE <<<<<<<<<<<<<<<<<<
        532: temp[0l:nz1-1l] = cats_bks_532_lfov[r,*]
          0: temp[0l:nz1-1l] = cats_depol_1064_lfov[r,*]
	  1: temp[0l:nz1-1l] = cats_depol_532_lfov[r,*]
  ENDCASE
  ENDELSE
  ;temp = reverse(temp)
  data[r,*] = temp
endfor  
       
rec_num = findgen(avrecs)
t1=avgtimestrarr[0]
t1=strsplit(t1,':',/extract)
t1=strjoin(t1,'-')
t2=avgtimestrarr[r-1]
t2=strsplit(t2,':',/extract)
t2=strjoin(t2,'-')
timeframe=gname
strput,timeframe,t1,37
strput,timeframe,t2,46

;prototype image plot
x=latitude
y=reverse(altitude)
if (num_x_axes[0] gt 1) then begin
   ;POSITION=[leftAxisLoc, bottomAxisLoc, rightAxesLoc, topAxesLoc] ;use this for multi-axes
   posvec=[0.09,0.20-(ytotsize*5e-5),0.81,0.86+(ytotsize*5e-5)] 
endif else begin
   posvec=[0.09,0.14,0.81,0.86] 
endelse
tit1  = tit1+': '+plotdate
CASE num_x_axes OF
     3: xtit = 'Latitude (degrees), Longitude (degrees), UTC (HH:MM:SS)' 
     2: xtit = 'Latitude (degrees), Longitude (degrees)' 
     1: xtit = 'Latitude (degrees)' 
ENDCASE
IF (ifov eq 'r') THEN BEGIN
CASE wavelength[0] OF
     1064: data_label = '_1064_rfov' 
      532: data_label = '_532_rfov'
        0: data_label = '_1064depol_rfov'
	1: data_label = '_532depol_rfov'
ENDCASE
ENDIF ELSE BEGIN
CASE wavelength[0] OF
     1064: data_label = '_1064_lfov' 
      532: data_label = '_532_lfov'
        0: data_label = '_1064depol_lfov'
	1: data_label = '_532depol_lfov'
ENDCASE
ENDELSE
ytit = 'Altitude (km)'
cbtit = 'Backscatter (km!E-1!Nsr!E-1!N)'
if (wavelength eq 0) then cbtit = 'depol ratio'
filename = custom_dir+timeframe+data_label+'.png'
x2 = strtrim(string(longitude,format='(f7.2)'),2)
x3 = avgtimestrarr
;plot_image_v0,data,x,y,absmax[0],posvec,xtotsize[0],ytotsize[0],tit1,tit2,xtit,ytit,cbtit,filename,num_x_axes[0],x2,x3
data = rotate(data,7) ;put into correct orientation
modelx = rotate(modelx,7) ;put into correct orientation
modely = rotate(modely,7) ;put into correct orientation
modelz = rotate(modelz,7) ;put into correct orientation
bare_image,data,absmax,xtotsize[0],ytotsize[0],39,filename

; Create image mesh array for COLLADA model for GE
nxi = 274L ;xtotsize[0]
nyi = 30L ;ytotsize[0] 
x_img_mesh = findgen(nxi)/float(nxi)
y_img_mesh = findgen(nyi)/float(nyi)

; Interpolate model coordinates array to mesh grid size
modelx = congrid(modelx,nxi,nyi,/interp) 
modely = congrid(modely,nxi,nyi,/interp)
modelz = congrid(modelz,nxi,nyi,/interp)

; ****************************** The following code writes the COLLADA model file *******************************************

COLLADA_base_filename = 'COLLADA_'+rawplotdate+starttimelabel+createtimelabel+'.dae'
COLLADA_file = custom_dir+COLLADA_base_filename
openw,mlun,COLLADA_file,/get_lun

accessor_count = nxi*nyi
count1 = nxi*nyi*3L
count2 = nxi*nyi*2L
count3 = (nxi-1L)*(nyi-1L)*4L

openr,chunk1,'COLLADA_chunk1.txt',/get_lun
finfo = file_info('COLLADA_chunk1.txt')
text = bytarr(finfo.size)
readu,chunk1,text
writeu,mlun,text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk1
free_lun,chunk1
printf,mlun,'<init_from>'+timeframe+data_label+'.png'+'</init_from>'+string([13b,10b]) ;variable COLLADA line [PNG file name]

openr,chunk2,'COLLADA_chunk2.txt',/get_lun
finfo = file_info('COLLADA_chunk2.txt')
text = bytarr(finfo.size)
readu,chunk2,text
writeu,mlun,text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk2
free_lun,chunk2
printf,mlun,'<float_array id="mesh1-geometry-position-array" count="'+strtrim(string(count1),2)+'">'+string([13b,10b]) ;variable COLLADA line [model coordinates]

for j=0L, nyi-1 do begin
for k=0L, nxi-1 do begin

    printf,mlun,strtrim(string(modelx[k,j]),2)+' '+strtrim(string(modely[k,j]),2)+' '+strtrim(string(modelz[k,j]),2)
    
endfor
endfor

openr,chunk3,'COLLADA_chunk3.txt',/get_lun
finfo = file_info('COLLADA_chunk3.txt')
text = bytarr(finfo.size)
readu,chunk3,text
writeu,mlun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk3
free_lun,chunk3

printf,mlun,'<accessor source="#mesh1-geometry-position-array" count="'+strtrim(string(accessor_count),2)+'" stride="3">'+string([13b,10b])

openr,chunk4,'COLLADA_chunk4.txt',/get_lun
finfo = file_info('COLLADA_chunk4.txt')
text = bytarr(finfo.size)
readu,chunk4,text
writeu,mlun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk4
free_lun,chunk4

printf,mlun,'<float_array id="mesh1-geometry-uv-array" count="'+strtrim(string(count2),2)+'">'+string([13b,10b])

for j=0L, nyi-1 do begin
for k=0L, nxi-1 do begin;

    printf,mlun,strtrim(string(x_img_mesh[k]),2)+' '+strtrim(string(y_img_mesh[j]),2)
    
endfor
endfor

openr,chunk5,'COLLADA_chunk5.txt',/get_lun
finfo = file_info('COLLADA_chunk5.txt')
text = bytarr(finfo.size)
readu,chunk5,text
writeu,mlun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk5
free_lun,chunk5

printf,mlun,'<accessor source="#mesh1-geometry-uv-array" count="'+strtrim(string(accessor_count),2)+'" stride="2">'+string([13b,10b])

openr,chunk6,'COLLADA_chunk6.txt',/get_lun
finfo = file_info('COLLADA_chunk6.txt')
text = bytarr(finfo.size)
readu,chunk6,text
writeu,mlun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk6
free_lun,chunk6

printf,mlun,'<triangles material="calipso_data" count="'+strtrim(string(count3),2)+'">'+string([13b,10b])

openr,chunk6b,'COLLADA_chunk6b.txt',/get_lun
finfo = file_info('COLLADA_chunk6b.txt')
text = bytarr(finfo.size)
readu,chunk6b,text
writeu,mlun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk6b
free_lun,chunk6b

; Print <p> under <triangles> kml block...
i1 = 0L
i2 = 1L
i3 = nxi+1L
i4 = nxi+1L
i5 = 1L
i6 = 0L
i7 = 0L
i8 = nxi+1L
i9 = nxi
i10 = nxi
i11 = nxi+1L
i12 = 0L
for k=1L, (nxi-1L)*(nyi-1L) do begin;IMPORTANT - k must start at 1!
    
    printf,mlun,strtrim(string(i1),2)+' '+strtrim(string(i2),2)+' '+ $
                strtrim(string(i3),2)+' '+strtrim(string(i4),2)+' '+ $
                strtrim(string(i5),2)+' '+strtrim(string(i6),2)+' '+ $
                strtrim(string(i7),2)+' '+strtrim(string(i8),2)+' '+ $
                strtrim(string(i9),2)+' '+strtrim(string(i10),2)+' '+ $
                strtrim(string(i11),2)+' '+strtrim(string(i12),2)
    if ((k mod (nxi-1) eq 0) and (k ne 0)) then begin
        i1=i1+2L 
	i2=i2+2L
	i3=i3+2L
	i4=i4+2L
	i5=i5+2L
	i6=i6+2L
	i7=i7+2L
	i8=i8+2L
	i9=i9+2L
	i10=i10+2L
	i11=i11+2L
	i12=i12+2L
    endif else begin
        i1=i1+1L 
	i2=i2+1L
	i3=i3+1L
	i4=i4+1L
	i5=i5+1L
	i6=i6+1L
	i7=i7+1L
	i8=i8+1L
	i9=i9+1L
	i10=i10+1L
	i11=i11+1L
	i12=i12+1L
    endelse		
    
endfor

print,'accessor count = ',nxi*nyi
print,'count1 = ',nxi*nyi*3L
print,'count2 = ',nxi*nyi*2L
print,'count3 = ',(nxi-1L)*(nyi-1L)*4L

openr,chunk8,'COLLADA_chunk8.txt',/get_lun
finfo = file_info('COLLADA_chunk8.txt')
text = bytarr(finfo.size)
readu,chunk8,text
writeu,mlun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk8
free_lun,chunk8

;END WRITING OF COLLADA FILE


; ****************************** The following code writes the KML file that calls the COLLADA file *******************************************

openw,klun,custom_dir+'KML_'+rawplotdate+starttimelabel+createtimelabel+'.kml',/get_lun

openr,chunk1,'KML_chunk1.txt',/get_lun
finfo = file_info('KML_chunk1.txt')
text = bytarr(finfo.size)
readu,chunk1,text
writeu,klun,text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk1
free_lun,chunk1

printf,klun,'<name>CATS'+createtimelabel+'</name>'+string([13b,10b])
printf,klun,'<description>CATS_'+rawplotdate+starttimelabel+'</description>'+string([13b,10b])
printf,klun,'<LookAt>'+string([13b,10b])
printf,klun,'<longitude>'+strtrim(string(longitude[avrecs/2]),2)+'</longitude>'
printf,klun,'<latitude>'+strtrim(string(latitude[avrecs/2]),2)+'</latitude>'

openr,chunk2,'KML_chunk2.txt',/get_lun
finfo = file_info('KML_chunk2.txt')
text = bytarr(finfo.size)
readu,chunk2,text
writeu,klun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk2
free_lun,chunk2

printf,klun,'<href>'+COLLADA_base_filename+'</href>'
printf,klun,'</Link>'
printf,klun,'<ResourceMap>'
printf,klun,'<Alias>'
printf,klun,'<targetHref>'+timeframe+data_label+'.png'+'</targetHref>'
printf,klun,'<sourceHref>'+timeframe+data_label+'.png'+'</sourceHref>'

openr,chunk3,'KML_chunk3.txt',/get_lun
finfo = file_info('KML_chunk3.txt')
text = bytarr(finfo.size)
readu,chunk3,text
writeu,klun,[13b,10b],text,[13b,10b] ;Windows carriage return @ end! Is only 10b for UNIX
close,chunk3
free_lun,chunk3


;END WRITING OF FIRST CURTAIN-PLOT KML FILE

;**************** Produce a World Map of the data *******************
colors=plot_colors()
window,3,xsize=600,ysize=500
map_set,/mercator,/continents,/advance,/grid,label=1 ;plot lat/lon on a map
plots,longitude,latitude,color=colors[15]
a=tvrd(true=1)
write_png,custom_dir+'map_'+timeframe+'.png',a

close,/all
free_lun,clun
free_lun,mlun
free_lun,klun


;**************** Create KML file that plots ground track in Google Earth using placemarks *******************

;The code after this point is to write the lat/lon points of the selected data
;into a KML format that Google Earth can read. (11/25/15)

openw,txt,custom_dir+'kml_latlon_'+timeframe+'.kml',/get_lun

h01='<?xml version="1.0" encoding="UTF-8"?>'
h02='<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">'
h03='<Document>'
h04='    <name>'+timeframe+'.kmz</name>'
h05='    <Style id="sn_grn-pushpin">'
h06='        <IconStyle>'
h07='            <scale>1.2</scale>'
h08='            <Icon>'
h09='                <href>http://maps.google.com/mapfiles/kml/pushpin/grn-pushpin.png</href>'
h10='            </Icon>'
h11='        </IconStyle>'
h12='    </Style>'
h13='    <Style id="sh_grn-pushpin">'
h14='        <IconStyle>'
h15='            <scale>1.2</scale>'
h16='            <Icon>'
h17='                <href>http://maps.google.com/mapfiles/kml/pushpin/grn-pushpin.png</href>'
h18='            </Icon>'
h19='        </IconStyle>'
h20='    </Style>'
h21='    <StyleMap id="this-is-the-pushpin">'
h22='        <Pair>'
h23='            <key>normal</key>'
h24='            <styleUrl>#sn_grn-pushpin</styleUrl>' ;CAN CHANGE COLOR HERE <------
h25='        </Pair>'
h26='        <Pair>' 
h27='           <key>highlight</key>'
h28='           <styleUrl>#sh_grn-pushpin</styleUrl>' ;CAN CHANGE COLOR HERE <------
h29='        </Pair>'
h30='    </StyleMap>' 
printf,txt,h01
printf,txt,h02
printf,txt,h03
printf,txt,h04
printf,txt,h05
printf,txt,h06
printf,txt,h07
printf,txt,h08
printf,txt,h09
printf,txt,h10
printf,txt,h11
printf,txt,h12
printf,txt,h13
printf,txt,h14
printf,txt,h15
printf,txt,h16
printf,txt,h17
printf,txt,h18
printf,txt,h19
printf,txt,h20
printf,txt,h21
printf,txt,h22
printf,txt,h23
printf,txt,h24
printf,txt,h25
printf,txt,h26
printf,txt,h27
printf,txt,h28
printf,txt,h29
printf,txt,h30

l01='<Placemark>'
l02='   <LookAt>'
l03='       <longitude>' & l03b='</longitude>'
l04='       <latitude>' & l04b='</latitude>'
l05='       <altitude>0</altitude>'
l06='       <heading>-1.840268938694791e-010</heading>'
l07='       <tilt>0</tilt>'
l08='       <range>1000.000177907706</range>'
l09='       <gx:altitudeMode>relativeToSeaFloor</gx:altitudeMode>'
l10='   </LookAt>'
l11='   <styleUrl>#this-is-the-pushpin</styleUrl>' 
l12='   <Point>'
l13='       <gx:drawOrder>1</gx:drawOrder>'
l14='       <coordinates>' & l14b=',0</coordinates>'
l15='   </Point>'
l16='   <description>' & l16b='</description>'
l17='</Placemark>'

for i=0,r-1 do begin
  lats=strtrim(string(latitude[i]),2)
  lons=strtrim(string(longitude[i]),2)
  printf,txt,l01
  printf,txt,l02
  printf,txt,l03+lons+l03b
  printf,txt,l04+lats+l04b
  printf,txt,l05
  printf,txt,l06
  printf,txt,l07
  printf,txt,l08
  printf,txt,l09
  printf,txt,l10
  printf,txt,l11
  printf,txt,l12
  printf,txt,l13
  printf,txt,l14+lons+','+lats+l14b
  printf,txt,l15
  printf,txt,l16+lons+','+lats+l16b
  printf,txt,l17
endfor

t01='</Document>'
t02='</kml>'
printf,txt,t01
printf,txt,t02
close,txt
free_lun,txt

STOP

end
