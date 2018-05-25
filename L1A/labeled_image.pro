PRO labeled_image,bs,absmax,x,z,xtit,ytit,title,nxticks,nyticks,xresize,yresize,xmargin,ymargin,xpad,ypad,layerZ,dx,outfilename,xlab,ylab,$
                  roundxnumbers=roundxnumbers,roundynumbers=roundynumbers,png=png,layers=layers,dist=dist
;Written by Patrick Selmer - Mod: 8/4/2014
;
;DESCRIPTION:
;This procedure injests a two dimensional array of values
;and related information then produces an image using the
;RAINBOW+WHITE IDL color table. This program was designed
;to produce an image to the screen of a Windows 7 machine.
;
;ARGUMENTS:
;bs - The 2D array used to produce the image
;     *MUST be 2 dimensional floating point array
;     *MUST have the same number of elements in its first
;           dimension as x
;     *MUST have the same number of elemets in its second
;           dimension as z
;absmax - The number to which to scale the data contained
;         within bs
;         *MUST be a scalar floating point number
;x - 1D array of horizontal dimension
;    *MUST be a 1 dimensional floating point array
;z - 1D array of vertical dimension
;xtit - The title of the x-axis
;       *MUST be a string (not an array of strings)
;ytit - The title of the y-axis
;       *MUST be a string (not an array of strings)
;title - The title places across the top of the image
;        *MUST be a string (not an array of strings)
;nxticks - The number of desirable xtick marks
;          *MUST be an integer
;nyticks - The number of desirable ytick marks
;          *MUST be an integer
;xresize - The width of the area in device coordinates
;          which you would like the actual bs image to
;          occupy
;          *MUST be an integer
;yresize - The height of the area in device coordinates
;          which you would like the actual bs image to
;          occupy
;          *MUST be an integer
;xmargin - The left indent of the image
;          *MUST be an integer
;ymargin - The bottom margin of the image
;          *MUST be an integer
;xpad - Defines the room to the right of the image
;       *MUST be an integer
;ypad - Defines the room above the image
;       *MUST be an integer
;
;layerZ - 3D array which defines the heights of the tops
;         and bottoms of atmospheric layers. This is
;         plotted overtop the image of the bs array only
;         if the layers keyword is set. 1st dimension
;         represents either the top [0] or bottom [1] of
;         an atmospheric layer detected. The 2nd dimesion
;         represents the layer number. The 3rd dimension
;         corresponds to the x ordinate (nx elements).
;         *MUST be 3 dimensional floating point array
;
;dx - The unit spacing between x array elements. This isn't
;     used unless the dist keyword is set. If the dist
;     keyword isn't set then a dummy placeholder value could
;     be put into this value
;     *MUST be a floating point number if the dist keyword is
;     set
;
;xlab - If you would like custom x-axis string labels, make this
;       a string array with the same number of elements as x. If
;       not, then set this equal to a 0 of 16 bit integer type
;
;ylab - If you would like custom y-axis string labels, make this
;       a string array with the same number of elements as x. If
;       not, then set this equal to a 0 of 16 bit integer type
;
;outfilename - String containing full path name with extension
;              of png image file that is produced.
;
;KEYWORDS:
;roundxnumbers - Set this keyword if you want nice round
;                numbers to label the x tick marks
;roundynumbers - Set this keyword if you want nice round
;                numbers to label the y tick marks
;png - Set this keyword to save a png of the image in the
;      the current directory
;
;layers - Set this keyword to plot the layerZ array over
;         top of the bs image
;
;dist - Tells the program to make an x ordinate array
;       filled with "distance from starting point"
;       values. Invokes use of dx variable.
;
;NOTES:
;Program not coded to handle units less than 1. All vertical
;arrays (bs[x,*] and z) must be

;load array of color indicies [red = 14]
colors = plot_colors()

;Get data types of axis label input variables
xdatatype = size(xlab,/type)
ydatatype = size(ylab,/type)

;Get the dimensions...
nx=n_elements(x)
nz=n_elements(z)
nlayers=n_elements(layerZ[0,*,0])

;Instead of using non-linear x ordinates (such as lat/lon), make simple distance from origin array
if keyword_set(dist) then begin
	x=findgen(nx)
	x=x*dx
endif

;/reverse array so that image is in correct orientation\
temp=fltarr(nz)
for i=0l,nx-1l do begin
	temp[0l:nz-1l]=bs[i,*]
	temp=reverse(temp)
	bs[i,*]=temp
endfor


;/****** FORM IMAGE *******\
;loadct,file='/Data/Code/Dev/selmer/colors.palm.tbl' ;load color table (IDL rainbow=39)
loadct,39
COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr
bsnorm = bytscl(bs,min=0.0,max=absmax) ;scale/normalize
rc=bindgen(nx,nz)
gc=bindgen(nx,nz)
bc=bindgen(nx,nz)
tvlct,R_orig,G_orig,B_orig ;load into the color translation table

for i=0l,nx-1l do begin ;find rgb values that correspond to variable values
 for j=0,nz-1 do begin
  rc[i,j]=R_orig[bsnorm[i,j]] ;red component of the image
  gc[i,j]=G_orig[bsnorm[i,j]] ;green component of the image
  bc[i,j]=B_orig[bsnorm[i,j]] ;blue component of the image
 endfor
endfor

xticklen=10 ;x-axis tick length in device coordinates
yticklen=10 ;x-axis tick length in device coordinates
window,1,xsize=xresize+xmargin+xpad,ysize=yresize+ymargin+ypad
rc=congrid(rc,xresize,yresize)
gc=congrid(gc,xresize,yresize)
bc=congrid(bc,xresize,yresize)
bsimage = [ [[rc]], [[gc]], [[bc]] ]
tv,bsimage,xmargin,ymargin,true=3 ;*produce actual graphic image on screen*

;/*********** LABEL IMAGE ************\
;draw box around image
plots,[[xmargin,ymargin],[xresize+xmargin,ymargin]],/device
plots,[[xmargin,ymargin],[xmargin,ymargin+yresize]],/device
plots,[[xmargin,ymargin+yresize],[xmargin+xresize,ymargin+yresize]],/device
plots,[[xmargin+xresize,ymargin+yresize],[xmargin+xresize,ymargin]],/device

;/convert x and z into device coordinates\
cz2d = (max(z)-min(z))/double(yresize) ;conversion factor (z)
cx2d = (max(x)-min(x))/double(xresize) ;conversion factor (x)
AreAnyXNeg = where(x lt 0.0,countnegx)
AreAnyZNeg = where(z lt 0.0,countnegz)
zdev = z/cz2d
xdev = x/cx2d
;convert layers heights into device coordinates
layerZdev = layerZ/cz2d
;figure out the maximum number of x decimal spaces for label
for l=0,5 do begin
     diff = (max(x)-min(x))-10.0d0^double(l)
     if (diff lt 0.0d0) then break
endfor
CASE l OF
    0: maxdecimalplacesX = 4;less than one
    1: maxdecimalplacesX = 3;one's place
    2: maxdecimalplacesX = 2;ten's place
    3: maxdecimalplacesX = 1;hundred's place
    4: maxdecimalplacesX = 0;thousand's place
    5: maxdecimalplacesX = 0;ten thousand's place
ENDCASE
;figuring out the maximum number of z decimal spaces for label
for m=0,5 do begin
     ;it't the range that determines this
     diff = (max(z)-min(z))-10.0d0^double(m)
     if (diff lt 0.0d0) then break
endfor
CASE m OF
  	0: maxdecimalplacesZ = 4;less than one
    1: maxdecimalplacesZ = 3;one's place
    2: maxdecimalplacesZ = 2;ten's place
    3: maxdecimalplacesZ = 1;hundred's place
    4: maxdecimalplacesZ = 0;thousand's place
    5: maxdecimalplacesZ = 0;ten thousand's place
ENDCASE

;Things are a bit more complex if there are neg values
;since device coordinates are never negative
convert2negx=0.0d0
if (countnegx gt 0) then begin
   convert2negx = abs(min(xdev))
   xdev = xdev + convert2negx
endif else begin
   xdev = xdev - abs(min(xdev))
endelse
convert2negz=0.0d0
if (countnegz gt 0) then begin
   convert2negz = abs(min(zdev))
   zdev = zdev + convert2negz
   layerZdev = layerZdev + convert2negz
endif else begin
   zdev = zdev - abs(min(zdev))
endelse

;Make arrays with device coordinate locations of tick marks
deltaZticks = (max(zdev)-min(zdev))/nyticks
deltaXticks = (max(xdev)-min(xdev))/nxticks
Zticks = dblarr(nz)*0.0d0
Xticks = dblarr(nx)*0.0d0
Xticklabels = strarr(nx)
Zticklabels = strarr(nz)
xlab2 = strarr(nx)
ylab2 = strarr(nz)

;find center values of both x/z ordinate arrays
xcenter = (max(x) + min(x))/2.0
zcenter = (max(z) + min(z))/2.0
xrevflag=0
if (x[0] gt x[nx-1]) then xrevflag = 1
zrevflag=0
if (z[0] gt z[nz-1]) then zrevflag = 1


; ------- X Axis Labeling Prep ---------
if keyword_set(roundxnumbers) then begin

  ;l is available to work with in this code chunk
  if (l gt 0) then begin
    l=l-1
    topbound = double(fix(max(x)/10.0d0^double(l))+1)*10.0d0^double(l)
    if (countnegz gt 0) then begin
      botbound = double(fix(min(x)/10.0d0^double(l))-1)*10.0d0^double(l)
    endif else begin
      botbound = 0.000d0
    endelse
  endif else begin
    print,'Values are less than 1, program not coded to handle. Stopping...'
    print,'Try running again without the roundxnumbers keyword'
    stop ;code more if you really want...
  endelse
  deltaXticks = ((topbound-botbound)/cx2d)/double(nxticks)
  for i=0,nxticks-1 do begin
    o=0 ;tells Xticklables to include an extra char for neg sign
    ;need to bump up bottom tick by diff btwn abs min and min u wnt btm tick 2b
    Xticks[i] = double(i)*deltaXticks+double(xmargin)+( double(botbound)-double(min(x)) )/double(cx2d)
    ;convert from dev to X coordinates for human-interpretable label
    ;STOP
    if (countnegx gt 0) then begin
      convertback2originalX = (Xticks[i]-double(xmargin))*cx2d-double(abs(min(x)))
    endif else begin
      convertback2originalX = (Xticks[i]-double(xmargin))*cx2d+double(abs(min(x)))
    endelse
    ;correct if x ordinate is high-to-low descending order from left to right
    if (xrevflag eq 1) then begin
       d = xcenter - convertback2originalX
       convertback2originalX = convertback2originalX + 2.0*d
    endif
    ;check for number that is extremly small but not exactly zero due to computational artifacts
    if (abs(convertback2originalX*100.0d0) lt 1.0d0) then convertback2originalX = 0.0d0
    ;ensure consistent decimal places for x-labels
    for k=0,5 do begin
      diff = convertback2originalX-10.0d0^double(k)
      if (diff lt 0.0d0) then break
    endfor
    if (xdatatype eq 7) then begin
      difx = x-convertback2originalX
      mindifx = min(abs(difx),mnloc)
      xlab2[i] = xlab[mnloc]
    endif
    Xticklabels[i] = string( convertback2originalX )
    Xticklabels[i] = strcompress(Xticklabels[i],/remove_all)
    if (convertback2originalX lt 0.0) then o=1
    CASE k OF
        0: Xticklabels[i] = strmid(Xticklabels[i],0,2+o+maxdecimalplacesX);less than one
    	1: Xticklabels[i] = strmid(Xticklabels[i],0,2+o+maxdecimalplacesX);one's place
    	2: Xticklabels[i] = strmid(Xticklabels[i],0,3+o+maxdecimalplacesX);ten's place
    	3: Xticklabels[i] = strmid(Xticklabels[i],0,4+o+maxdecimalplacesX);hundred's place
    	4: Xticklabels[i] = strmid(Xticklabels[i],0,5+o+maxdecimalplacesX);thousand's place
       	5: Xticklabels[i] = strmid(Xticklabels[i],0,6+o+maxdecimalplacesX);ten thousand's place
    ENDCASE
  endfor
  Xticks = Xticks[ where(Xticks ne 0.0d0) ]
  goodlocs = where(Xticklabels ne '')
  Xticklabels = Xticklabels[ goodlocs ]
  xlab2 = xlab2[ goodlocs ]
endif else begin

  ;Make arrays with device coordinate locations of tick marks
  deltaZticks = (max(zdev)-min(zdev))/nyticks
  deltaXticks = (max(xdev)-min(xdev))/nxticks
  Zticks = dblarr(nz)*0.0d0
  Xticks = dblarr(nx)*0.0d0
  Xticklabels = strarr(nx)
  Zticklabels = strarr(nz)
  i=0
  repeat begin
      o=0 ;tells Xticklables to include an extra char for neg sign
      Xticks[i] = double(i)*deltaXticks+xmargin
      if (countnegx gt 0) then begin
        convertback2originalX = (Xticks[i]-double(xmargin))*cx2d - double(abs(min(x)))
      endif else begin
        convertback2originalX = (Xticks[i]-double(xmargin))*cx2d + double(abs(min(x)))
      endelse
      ;correct if x ordinate is high-to-low descending order from left to right
      if (xrevflag eq 1) then begin
         d = xcenter - convertback2originalX
         convertback2originalX = convertback2originalX + 2.0*d
      endif
      ;ensure consistent decimal places for x-labels
      for k=0,5 do begin
        diff = convertback2originalX-10.0d0^double(k)
        if (diff lt 0.0d0) then break
      endfor
      if (xdatatype eq 7) then begin
        difx = x-convertback2originalX
        mindifx = min(abs(difx),mnloc)
        xlab2[i] = xlab[mnloc]
      endif
      Xticklabels[i] = string( convertback2originalX )
      Xticklabels[i] = strcompress(Xticklabels[i],/remove_all)
      if (convertback2originalX lt 0.0) then o=1
      CASE k OF
    	0: Xticklabels[i] = strmid(Xticklabels[i],0,2+o+maxdecimalplacesX);less than one
    	1: Xticklabels[i] = strmid(Xticklabels[i],0,2+o+maxdecimalplacesX);one's place
    	2: Xticklabels[i] = strmid(Xticklabels[i],0,3+o+maxdecimalplacesX);ten's place
    	3: Xticklabels[i] = strmid(Xticklabels[i],0,4+o+maxdecimalplacesX);hundred's place
    	4: Xticklabels[i] = strmid(Xticklabels[i],0,5+o+maxdecimalplacesX);thousand's place
       	5: Xticklabels[i] = strmid(Xticklabels[i],0,6+o+maxdecimalplacesX);ten thousand's place
      ENDCASE
      i=i+1
  endrep until ((Xticks[i-1] ge max(xdev)+xmargin) or (i ge nx))
  Xticks = Xticks[ where(Xticks ne 0.0d0) ]
  goodlocs = where(Xticklabels ne '')
  Xticklabels = Xticklabels[ goodlocs ]
  xlab2 = xlab2[ goodlocs ]

endelse

; ------- Y (Z) Axis Labeling Prep ---------
;if keyword roundynumbers is set...
if keyword_set(roundynumbers) then begin

  ;m is available to work with in this code chunk
  if (m gt 0) then begin
    m=m-1
    topbound = double(fix(max(z)/10.0d0^double(m))+1)*10.0d0^double(m)
    if (countnegz gt 0) then begin
      botbound = double(fix(min(z)/10.0d0^double(m))-1)*10.0d0^double(m)
    endif else begin
      botbound = 0.000d0
    endelse
  endif else begin
    print,'Values are less than 1, program not coded to handle. Stopping...'
    print,'Try running again without the roundynumbers keyword'
    stop ;code more if you really want...
  endelse
  ;delta Z tick marks in device coordinates (conv magn not vect)
  deltaZticks = ((topbound-botbound)/cz2d)/double(nyticks)
  for i=0,nyticks-1 do begin
    o=0 ;tells Zticklables to include an extra char for neg sign
    ;need to bump up bottom tick by diff btwn abs min and min u wnt btm tick 2b
    Zticks[i] = double(i)*deltaZticks+double(ymargin)+( double(botbound)-double(min(z)) )/double(cz2d)
    ;convert from dev to Z coordinates for human-interpretable label
    if (countnegz gt 0) then begin
      convertback2originalZ = (Zticks[i]-double(ymargin))*cz2d-double(abs(min(z)))
    endif else begin
      convertback2originalZ = (Zticks[i]-double(ymargin))*cz2d+double(abs(min(z)))
    endelse
    ;;correct if z ordinate is high-to-low descending order from left to right
    ;if (zrevflag eq 1) then begin
       ;d = zcenter - convertback2originalZ
       ;convertback2originalZ = convertback2originalZ + 2.0*d
    ;endif
    ;check for number that is extremly small but not exactly zero due to computational artifacts
    if (abs(convertback2originalZ*100.0d0) lt 1.0d0) then convertback2originalZ = 0.0d0
    ;ensure consistent decimal places for x-labels
    for k=0,5 do begin
      diff = convertback2originalZ-10.0d0^double(k)
      if (diff lt 0.0d0) then break
    endfor
    if (ydatatype eq 7) then begin
      dify = z-convertback2originalZ
      mindify = min(abs(dify),mnloc)
      ylab2[i] = ylab[mnloc]
    endif
    Zticklabels[i] = string( convertback2originalZ )
    Zticklabels[i] = strcompress(Zticklabels[i],/remove_all)
    if (convertback2originalZ lt 0.0) then o=1
    CASE k OF
        0: Zticklabels[i] = strmid(Zticklabels[i],0,2+o+maxdecimalplacesZ);less than one
    	1: Zticklabels[i] = strmid(Zticklabels[i],0,2+o+maxdecimalplacesZ);one's place
    	2: Zticklabels[i] = strmid(Zticklabels[i],0,3+o+maxdecimalplacesZ);ten's place
    	3: Zticklabels[i] = strmid(Zticklabels[i],0,4+o+maxdecimalplacesZ);hundred's place
    	4: Zticklabels[i] = strmid(Zticklabels[i],0,5+o+maxdecimalplacesZ);thousand's place
       	5: Zticklabels[i] = strmid(Zticklabels[i],0,6+o+maxdecimalplacesZ);ten thousand's place
    ENDCASE
  endfor
  Zticks = Zticks[ where(Zticks ne 0.0d0) ]
  goodlocs = where(Zticklabels ne '')
  Zticklabels = Zticklabels[ goodlocs ]
  ylab2 = ylab2[ goodlocs ]

;if keyword roundynumbers is NOT set...
endif else begin

  j=0
  repeat begin
      o=0 ;tells Zticklables to include an extra char for neg sign
      Zticks[j] = double(j)*deltaZticks+ymargin
      convertback2originalZ = (Zticks[j]-ymargin)*cz2d- double(abs(min(z)))
      ;convert from dev to Z coordinates for human-interpretable label
      if (countnegz gt 0) then begin
        convertback2originalZ = (Zticks[i]-double(ymargin))*cz2d-double(abs(min(z)))
      endif else begin
        convertback2originalZ = (Zticks[i]-double(ymargin))*cz2d+double(abs(min(z)))
      endelse
      ;;correct if z ordinate is high-to-low descending order from left to right
      ;if (zrevflag eq 1) then begin
         ;d = zcenter - convertback2originalZ
         ;convertback2originalZ = convertback2originalZ + 2.0*d
      ;endif
      ;ensure consistent decimal places for z-labels
      for k=0,5 do begin
          diff = convertback2originalZ-10.0d0^double(k)
          if (diff lt 0.0d0) then break
      endfor
      if (ydatatype eq 7) then begin
        dify = z-convertback2originalZ
        mindify = min(abs(dify),mnloc)
        ylab2[i] = ylab[mnloc]
      endif
      Zticklabels[j] = string( convertback2originalZ )
      Zticklabels[j] = strcompress(Zticklabels[j],/remove_all)
      if (convertback2originalZ lt 0.0) then o=1
      CASE k OF
    	  0: Zticklabels[j] = strmid(Zticklabels[j],0,2+o+maxdecimalplacesZ);less than one
    	  1: Zticklabels[j] = strmid(Zticklabels[j],0,2+o+maxdecimalplacesZ);one's place
    	  2: Zticklabels[j] = strmid(Zticklabels[j],0,3+o+maxdecimalplacesZ);ten's place
    	  3: Zticklabels[j] = strmid(Zticklabels[j],0,4+o+maxdecimalplacesZ);hundred's place
    	  4: Zticklabels[j] = strmid(Zticklabels[j],0,5+o+maxdecimalplacesZ);thousand's place
    	  5: Zticklabels[j] = strmid(Zticklabels[j],0,6+o+maxdecimalplacesZ);ten thousand's place
      ENDCASE
      j=j+1
  endrep until ((Zticks[j-1] ge max(zdev)+ymargin) or (j ge nz))
  Zticks = Zticks[ where(Zticks ne 0.0d0) ]
  goodlocs = where(Zticklabels ne '')
  Zticklabels = Zticklabels[ goodlocs ]
  ylab2 = ylab2[ goodlocs ]

endelse






;plot x tick marks (in device coordinates)
if (xdatatype eq 7) then Xticklabels=xlab2 ;USE custom user axis labels if type eq "string" (which = 7)
if (ydatatype eq 7) then Zticklabels=ylab2
for i=0,n_elements(Xticks)-1 do begin
   if ((Xticks[i]-xmargin le xresize) and (Xticks[i] ge xmargin)) then begin ;avoid stray tick mark away from box
      plots,[ [Xticks[i],ymargin],[Xticks[i],ymargin-xticklen] ],/device
      xyouts,Xticks[i],ymargin-2*xticklen,Xticklabels[i],alignment=0.5,/device
   endif
endfor

;plot y tick marks (in device coordinates)
for j=0,n_elements(Zticks)-1 do begin
   if ((Zticks[j]-ymargin le yresize) and (Zticks[j] ge ymargin )) then begin ;avoid stray tick mark away from box
      plots,[ [xmargin,Zticks[j]],[xmargin-yticklen,Zticks[j]] ],/device
      xyouts,xmargin-1.2*yticklen,Zticks[j]-3,Zticklabels[j],alignment=1.0,/device
   endif
endfor

;plot layer heights over top if layers keyword is set
if keyword_set(layers) then begin

    if (xrevflag eq 1) then begin
       temp=fltarr(nx)
       for i=0,1 do begin
       	for j=0,nlayers-1 do begin
       	    temp[0:nx-1]=layerZdev[i,j,*]
       	    temp=reverse(temp)
       	    layerZdev[i,j,*]=temp
       	endfor
       endfor
    endif

	for i=0,nlayers-1 do begin
		for j=0,nx-1 do begin
		  ;don't plot points outside plot area
		  if ( (layerZdev[0,i,j] lt yresize) and (layerZdev[0,i,j] gt 0) ) then begin
		    plots,xdev[j]+xmargin,layerZdev[0,i,j]+ymargin,psym=3,color=colors[14],/device
		  endif
		  ;don't plot points outside plot area
		  if ( (layerZdev[1,i,j] lt yresize) and (layerZdev[1,i,j] gt 0) ) then begin
		    plots,xdev[j]+xmargin,layerZdev[1,i,j]+ymargin,psym=3,color=colors[14],/device
		  endif
		endfor
	endfor

endif

;Write title of image
xyouts,xmargin,ypad/3+ymargin+yresize,title,alignment=0.0,/device,charsize=1.5,charthick=2

;Label the x axis
xyouts,xresize/2+xmargin,ymargin-3*ymargin/4,xtit,alignment=0.5,/device,charsize=1.5,charthick=1

;Label the y axis
xyouts,xmargin-3*xmargin/4,yresize/2+ymargin,ytit,alignment=0.5,/device,orientation=90,charsize=1.5,charthick=1

;Add a color bar
;Currently hardcoded to put colorbar into right margin (uses space allotted by xpad)
totalx = float(xresize+xmargin+xpad)
totaly = float(yresize+ymargin+ypad)
x0 = xresize/2+xmargin
xextent = xresize/2
x1 = xextent + x0
y0 = totaly - ypad/2
y1 = totaly - ypad/4
x0 = float(x0)/float(totalx)
y0 = float(y0)/float(totaly)
x1 = float(x1)/float(totalx)
y1 = float(y1)/float(totaly)
pos=[x0, y0, x1, y1] ;defines both size and position of color bar
minpos = min(pos)
if (minpos le 0) then print, 'Hey operator, you probably need to make "xpad" larger...'
colorbar2,maxrange=absmax,position=pos,divisions=4,charsize=1.0
;Save a PNG file if png keyword is set
if keyword_set(png) then begin
    a=tvrd(true=1)
    write_png,outfilename,a
endif


end
