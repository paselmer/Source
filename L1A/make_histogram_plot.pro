PRO make_histogram_plot,array,bsize,tit,xtit,a,relative=relative,zbuff=zbuff

; PURPOSE: 
; To create a histogram bar plot of an input array of data
;
; DESCRIPTION: 
; This code uses the IDL histogram function to produce a 
; frequency distribution of data which is then plotted.
; The "relative" keyword causes the procedure to plot
; relative frequency rather than absolute frequency
; 

freq = histogram(array,binsize=bsize,locations=lower_bounds)
if keyword_set(relative) then begin
   ytit = 'Relative Frequency'
   freq = freq / total(freq)
endif else begin
   ytit = 'Frequency'
endelse   
nb = n_elements(freq)

if (nb eq 0) then STOP, 'NO ELEMENTS IN HISTOGRAM! CHECK INPUT ARRAY.'

upper_bounds = lower_bounds + bsize


;******* Create actual histogram plot *******

xtotsize = 800
ytotsize = 700

if keyword_set(zbuff) then begin
   ; Plot graphic in z-buffer so it doesn't appear on screen
   myDevice = !D.NAME ; get device in use
   set_plot, 'Z'
   Erase
   Device, Set_Resolution=[xtotsize,ytotsize], Set_Pixel_Depth=24, Decomposed=1
endif else begin
   window,xsize=xtotsize,ysize=ytotsize
endelse

ncolors=6
cindex_array=lindgen(ncolors)
colors=lindgen(3,ncolors)
colors[0:2,0]=[0,0,255]     ;blue
colors[0:2,1]=[255,0,0]    ;red
colors[0:2,2]=[0,139,0]    ;green 4
colors[0:2,3]=[198,99,231] ;lavender
colors[0:2,4]=[210,105,30]  ;chocolate
colors[0:2,5]=[255,255,255] ;white
for i=0,ncolors-1 do begin
	cindex_array[i] = colors[0,i] + 256L * (colors[1,i] + 256L * colors[2,i])
endfor
colors=cindex_array

plot,upper_bounds,freq,xrange=[lower_bounds[0],upper_bounds[nb-1]],xtitle=xtit,$
     ytitle=ytit,/NODATA,color=0,background=colors[5],charsize=1.5,ticklen=-0.03, $
     xtick_get=v,ymargin=[5,5],xmargin=[12,3]

for j=0, nb-1 do begin
    polyfill,[lower_bounds[j],lower_bounds[j],upper_bounds[j],upper_bounds[j]],[0,freq[j],freq[j],0],color=colors[0]
    plots,[lower_bounds[j],lower_bounds[j]],[0,freq[j]],color=0
    plots,[upper_bounds[j],upper_bounds[j]],[0,freq[j]],color=0
    plots,[lower_bounds[j],upper_bounds[j]],[freq[j],freq[j]],color=0
endfor

;Manually place title using XYOUTS procedure
x = (!X.Window[1] - !X.Window[0]) / 2. + !X.Window[0]
y=0.95
xyouts, x, y, tit, /Normal, Alignment=0.5, Charsize=2.0, color=0
a = 0 ;nothing if 'zbuff' keyword not set

if keyword_set(zbuff) then begin
    a=tvrd(true=1)
    Set_Plot, myDevice
endif    

end
