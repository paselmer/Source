PRO find_surface_alt,tabs,i,z,surface_alt,guessalt,margin,gthreshfact,gwidthinbins,backgroundalt,Lsearch_limit,mingthresh,remove=remove

;Written by Patrick Selmer
;Last modified - 07/18/16

;PURPOSE:
;This code's purpose is to identify the altitude of the Earth's surface
;in CATS lidar data.

;DESCRIPTION:
;The search for the ground is limited to the domain bounded by the lowest altitude
;in the attenuated backscatter profile, tabs, and the altitude defined by guessalt+margin.
;The median of the region of the profile assumed to be background (defined to be
;the region below backgroundalt) is calculated. The background nosie, standard deviation
;of tabs, is computed. A threshold is computed using the median plus
;gthreshfact times the background noise. The tabs profile is searched low alt to high alt
;,bin by bin. The first bin where tabs is greater than the threshold is assumed to
;be the start of the ground signal. The search will continue until either the
;top of the search region is reached or the number of consecutive bins greater
;than the threshold is equal to gwidthinbins. The altitude of the maximum tabs
;of the bins selected as ground signal is considered the surface altitude.
;The surface altitude is stored in its appropriate profile index in the
;surface_alt array.

;VARIABLES:
;*All altitudes specified in meters*
;tabs - Input. Attenuated backscatter array (nprofs x nbins)
;i - Input. profile number index
;z - Input. altitude array (nbins)
;surface_alt - Output. Array where computed surface altitudes are stored
;guessalt - An input parameter which is the guess for the surface altitude
;margin - Input. The altitude slop above your guess. Sort of the positive margin of error
;gthreshfact - Input. the ground threshold factor. applied to background
;gwidthinbins - Input. the maximum number of bins you think the ground return can span
;backgroundalt - Input. the altitude below which you are confident the ground return does
;                exist
;Lsearch_limit - Output. Do not search for layers at or below this altitude.
;mingthresh - Input. Min value threshold is allowed to be.

;KEYWORDS:
;remove - this sets the R array equal to zero within bins deemed "surface return"

;NOTES:
;
; 7/5/16:
;Lsearch_limit is now set to be the highest of the consecutive bins above the threshold. 



searchtop = guessalt + margin ;top altitude of the ground search region

blwground = where(z lt backgroundalt,count)

select_bins = lonarr(gwidthinbins)


IF (count ne 0) THEN BEGIN ;simple error check. in theory this check should never be failed.

 ; set a default layer search limit
 Lsearch_limit = 0.0

 med = median(tabs[i,blwground])
 background = stddev(tabs[i,blwground])
 thresh = med + gthreshfact*background
 ori_thresh = thresh ;DELETE
 if (thresh lt mingthresh) then thresh = mingthresh
 ;IF ( (I GE 28704) AND (I LE 28717) ) THEN BEGIN
    ;WINDOW,0,XSIZE=700,YSIZE=850
    ;plot,tabs[i,*],z,charsize=1.5,xrange=[-0.01,0.1]
    ;PRINT,MED,BACKGROUND,THRESH
    ;stop
 ;ENDIF

 gtcount=0 ; counts bins where tabs is greater than the threshold
 j = n_elements(z)-1 ;assumes highest element is lowest alt
 while ((z[j] lt searchtop) and (gtcount lt gwidthinbins)) do begin
       if (tabs[i,j] gt thresh) then begin
           select_bins[gtcount] = j
           gtcount=gtcount+1
       endif else if (gtcount gt 0) then begin
           j=j-1
           break
       endif
       j=j-1 ;search upward
 endwhile
 
 if (gtcount gt 0) then begin 
    select_bins = select_bins[0:gtcount-1] ;make sure only as long as has to be
    Lsearch_limit = z[select_bins[gtcount-1]]
 endif
 ;MAXtabs=0.0
 if (gtcount eq 1) then begin
    maxtabs = max(tabs[i,select_bins[0:gtcount-1]],relbin)
    bin_of_maxtabs = select_bins[relbin];where(tabs[i,*] eq maxtabs)
    surface_alt[i] = z[bin_of_maxtabs[0]]
    ;if (tabs[i,bin_of_maxtabs[0]] lt 0.009) then begin 
    ;   print,tabs[i,bin_of_maxtabs[0]],surface_alt[i],i
    ;   COLORS=PLOT_COLORS()
    ;   PLOT,tabs[I,WHERE(Z LE SEARCHTOP)],Z[WHERE(Z LE SEARCHTOP)]
    ;   PRINT,I,THRESH,MAXtabs,SURFACE_ALT[I]
    ;   PRINT,MED,BACKGROUND,'gt count = 1'
    ;   stop
    ;endif
    if (keyword_set(remove) eq 1) then begin
       tabs[i,select_bins[gtcount-1]:select_bins[0]] = 0.0
    endif
 endif
 if (gtcount gt 1) then begin
    maxtabs = max(tabs[i,select_bins[0:gtcount-1]],relbin)
    bin_of_maxtabs = select_bins[relbin];where(tabs[i,*] eq maxtabs)
    ;calculate the second derivates - slope of slopes - code could be optimized
    diffofdiffs = fltarr(gtcount)
    for u=0, gtcount-1 do begin
       if (select_bins[u]+1 lt n_elements(z)) then diffofdiffs[u] = $
       (tabs[i,select_bins[u]]-tabs[i,select_bins[u]+1]) - (tabs[i,select_bins[u]-1]-tabs[i,select_bins[u]])
    endfor
    maxdod = max(diffofdiffs,relbin)
    bin_of_maxdod = select_bins[relbin]
    if ((bin_of_maxtabs eq bin_of_maxdod) or (gtcount lt gwidthinbins)) then begin
       surface_alt[i] = z[bin_of_maxtabs[0]]
       ;if (tabs[i,bin_of_maxtabs[0]] lt 0.009) then begin 
       ;   print,tabs[i,bin_of_maxtabs[0]],surface_alt[i],i
       ;   COLORS=PLOT_COLORS()
       ;   PLOT,tabs[I,WHERE(Z LE SEARCHTOP)],Z[WHERE(Z LE SEARCHTOP)]
       ;   PRINT,I,THRESH,MAXtabs,SURFACE_ALT[I]
       ;   PRINT,MED,BACKGROUND,'gt count > 1'
       ;   stop
       ;endif
    endif
    ;COLORS=PLOT_COLORS()
    ;PLOT,tabs[I,WHERE(Z LE SEARCHTOP)],Z[WHERE(Z LE SEARCHTOP)]
    ;OPLOT,  [0,30] ,  [surface_alt[i],surface_alt[i]]  , color=colors[1]
    ;PRINT,I,THRESH,MAXtabs,SURFACE_ALT[I]
    ;PRINT,MED,BACKGROUND
    ;STOP
 endif


 ;COLORS=PLOT_COLORS()
 ;IF (I GE 14267) THEN BEGIN
    ;PLOT,tabs[I,WHERE(Z LE SEARCHTOP)],Z[WHERE(Z LE SEARCHTOP)],psym=1
    ;OPLOT,  [0,30] ,  [surface_alt[i],surface_alt[i]]  , color=colors[1]
    ;PRINT,gtcount,I,THRESH,MAXtabs,SURFACE_ALT[I]
    ;PRINT,MED,BACKGROUND
    ;STOP
 ;ENDIF

ENDIF

end
