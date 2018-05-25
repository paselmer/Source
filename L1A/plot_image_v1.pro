PRO plot_image_v1,data,x,y,absmax,posvec,xtotsize,ytotsize,tit1,tit2,xtit,ytit,cbtit,filename,num_x_axes,x2,x3,ct

;DESCRIPTION:
;This code is essentially the successor to my labeled_image code. This code uses the built-in IDL plot
;procedure to create a plotting area, then fits an image inside the plotting area. The axes labels created
;by "PLOT" should exactly correspond to the the image within. In this version, a color bar is automatically
;generated and there is no way of turning that off (without screwing with the code of course ;-) ). The code
;automatically saves the plot to a .png file at the location specified by the "filename" variable.
;
;Color table number made to be an input for v1 [7/29/2016]

;INPUTS:
;data       -> 2-d array (float preferred) with first dimension corresponding to the horizontal direction
;              and the second dimension corresponding to the vertical direction
;x          -> 1-d array corresponding to horizontal units of data
;y          -> 1-d array corresponding to the vertical units of data
;absmax     -> The highest value of the color bar (scale)
;posvec     -> a vector in the style of IDL's "POSITION" which control the plotting area within the graphics window
;xtotsize   -> The total x-size of the graphics window
;ytotsize   -> The total y-size of the graphics window
;tit1       -> The primary title of the plot
;tit2       -> The secondary title of the plot
;xtit       -> The title of the x-axis
;ytit       -> The title of the y-axis
;cbtit      -> The title of the color bar
;filename   -> The name of the output png file
;num_x_axes -> Corresponds to the number of labels to the x-axis. This was conceived with labeling lat/lon/time all
;              on the x-axis at once (for lidar data).
;x2         -> if num_x_axes is set > 1 this is the second layer horizontal units of data
;x3         -> if num_x_axes is set > 2 this is the third layer horizontal units of data
;ct         -> The index of the IDL color table to be used
;
;OUTPUTS:
;(none)

;set some tunable parameters here
tl = -0.03 ;ticklength
bigchar = 2.0
littlechar = 1.5
littlecharthick = 1.5
cbchar = littlechar;/1.5

;Get the dimensions...
nx = n_elements(x)
ny = n_elements(y)

;create graphics window
window,1,xsize=xtotsize,ysize=ytotsize

;create an empty plotting area with all labels except color bar labels
colors=plot_colors() ;load an array of color indicies
fakex=findgen(nx)
nxticks=6
plot,fakex,y,position=posvec,xtick_get=v,/NODATA,xstyle=1,xticks=nxticks,xtickinterval=float(nx-1)/float(nxticks-1)
v = ulong(v)
!x.tickname = strtrim(string(x[v],format='(f6.2)'),2) ;*** AT SOME POINT CHANGE THIS LINE SO X VAR IS AUTO FORMATTED
plot,fakex,y,position=posvec,/NODATA,ticklen=tl,charsize=littlechar,color=0,background=colors[160], $
     ytitle=ytit,xthick=2,ythick=2,yrange=[min(y),max(y)],ystyle=1,charthick=littlecharthick,xstyle=1,xticks=nxticks,$
     xtickinterval=float(nx-1)/float((nxticks-1))
;Write out the title
xn = (posvec[0]+posvec[2])/2.0
yn = 1.0 - (20.0/ytotsize)
xyouts,xn,yn,tit1,charsize=bigchar,color=0,/NORMAL,alignment=0.5,charthick=2.0
xn = (posvec[0]+posvec[2])/2.0
yn = (1.0 + posvec[3])/2.0 - (1.0 - posvec[3])/6.0 
xyouts,xn,yn,tit2,charsize=littlechar,color=0,/NORMAL,alignment=0.5,charthick=littlecharthick
if (num_x_axes le 1) then begin
  ;Write out the x title
  xn = (posvec[0]+posvec[2])/2.0
  yn = posvec[1]/2.0 - posvec[1]/4.0
  xyouts,xn,yn,xtit,charsize=littlechar,color=0,/NORMAL,alignment=0.5,charthick=littlecharthick
endif else begin
  nv = n_elements(v)
  delx = (posvec[2] - posvec[0])/(nv-1)
  xa = findgen(nv)*delx + posvec[0]
  yn0 = posvec[1] ;- posvec[1]*(40.0/100.0)
  for i=2,num_x_axes do begin
      for j=0, nv-1 do begin
          yn = yn0 - ( abs(tl) + 1.05*((float(!d.y_ch_size)/float(ytotsize))*float(i+1)) + float(i-2)*(4.0/float(ytotsize)) ) 
          if (i eq 2) then xyouts,xa[j],yn,x2[v[j]],color=0,/NORMAL,alignment=0.5,charsize=littlechar,charthick=littlecharthick
          if (i eq 3) then xyouts,xa[j],yn,x3[v[j]],color=0,/NORMAL,alignment=0.5,charsize=littlechar,charthick=littlecharthick
      endfor
  endfor
  xn = (posvec[0]+posvec[2])/2.0
  yn = yn0 - ( abs(tl) + 1.05*((float(!d.y_ch_size)/float(ytotsize))*float(i+1)) + float(i-2)*(10.0/float(ytotsize)) ) 
  if (ytotsize lt 600) then yn = float(!d.y_ch_size)/float(ytotsize)  ;added 12/18/15
  xyouts,xn,yn,xtit,charsize=littlechar,color=0,/NORMAL,alignment=0.5,charthick=littlecharthick
endelse


;/*** FORM IMAGE ***\
loadct,ct
COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr
datanorm = bytscl(data,min=0.0,max=absmax) ;scale/normalize
rc=bindgen(nx,ny)
gc=bindgen(nx,ny)
bc=bindgen(nx,ny)
tvlct, R_orig, G_orig, B_orig ;load into the color translation table

;find rgb values that correspond to variable values
rc = R_orig[datanorm]
gc = G_orig[datanorm]
bc = B_orig[datanorm]

;compute posvec in terms of device coordinates
left_devc = posvec[0]*xtotsize
bot_devc = posvec[1]*ytotsize
right_devc = posvec[2]*xtotsize
top_devc = posvec[3]*ytotsize

;compute margins in terms of device coordinates
left_marg =  left_devc + 1 ;plus 1 is important
bot_marg = bot_devc + 1
right_marg = xtotsize - right_devc
top_marg = ytotsize - top_devc

;compute size of image on plot
xresize = right_devc - left_devc - 1 ;minus 2 is important
yresize = top_devc - bot_devc ;- 1

;Write image into plotting area
rc=congrid(rc,xresize,yresize)
gc=congrid(gc,xresize,yresize)
bc=congrid(bc,xresize,yresize)
data_image = [ [[rc]], [[gc]], [[bc]] ]
tv,data_image,left_marg,bot_marg,true=3 ;*produce actual graphic image on screen*

;Add a labeled color bar
;Currently hardcode to put a color bar in the right margin
xc0 = posvec[2] + (1.0-posvec[2])/4.0
yc0 = posvec[1]
xc1 = posvec[2] + 3.0*(1.0-posvec[2])/8.0
yc1 = posvec[3]
cbpos = [xc0,yc0,xc1,yc1]
colorbar2,maxrange=absmax,position=cbpos,division=4,charsize=cbchar,/vertical,/right,color=0
xn = 1.0-(1.0-posvec[2])/8.0
yn = (posvec[3]+posvec[1])/2.0
xyouts,xn,yn,cbtit,orientation=90,alignment=0.5,charsize=littlechar,color=0,/NORMAL

;save the image to a png file
a=tvrd(true=1)
write_png,filename,a


end
