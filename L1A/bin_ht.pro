function BIN_HT, bin_no,vres

; This function calculates the height (amsl) of the input data bin assignment
; for CPL data.  Height will be given in the same units as the vertical
; resolution (vres), ie meters or kilometers.  The first bin, bin 0, is at the
; top of the CPL frame (near 22000 meters above msl).  Height calculations
; refer to the center of the bin.  There are a total of 900 bins when
; the resolution is near 30 meters.

; Allow for vres approximations.
if (vres ge 29.0 and vres le 30.2) then vres= 29.980
if (vres gt 30.2 and vres le 31.0) then vres= 30.787
if (vres ge 0.0290 and vres le 0.0302) then vres= 0.02998
if (vres gt 0.0302 and vres le 0.031) then vres= 0.030787

if (vres eq 15.0) then sbin= 21997.5 else $
if (vres eq 0.015) then sbin= 21.9975 else $
if (vres eq 29.980) then sbin= 22305.12 else $
if (vres eq 0.02998) then sbin= 22.30512 else $
if (vres eq 30.787) then sbin= 22566.871 else $
if (vres eq 0.030787) then sbin= 22.566871 else begin $
  bin_elv= -999.9
  print, "Problems ingesting vertical resolution in BIN_HT...", vres
  goto, end_proc
endelse

; Calculate height of bin assignment (msl).
bin_elv= sbin - (bin_no*vres)

end_proc:

return, bin_elv
end
;-----------------------------------------------------------------------------
