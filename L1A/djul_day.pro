function DJUL_DAY, year,month,day,hour,minute,sec,jday
;     This function calculates the julian day in double precision decimal
;     from either calander day (jday must be input as -1)
;     or integer day of year.

daynum=0.0D
ldpmon = intarr(12,2)
ldpmon(*,0) = [0,31,59,90,120,151,181,212,243,273,304,334]
ldpmon(*,1) = [0,31,60,91,121,152,182,213,244,274,305,335]        

if (jday(0) eq -1) then begin
  k=0
  if ((year mod 4) eq 0) then k=1
  jday=day+ldpmon((month-1),k)
endif  

daynum=(double(hour)/24.0D)+(double(minute)/1440.0D)+(double(sec)/86400.0D)+jday

return, daynum
end
;-----------------------------------------------------------------------------
