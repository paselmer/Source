in_dir = '/cpl3/CAMAL/Config/dttables/'

dttable = in_dir + 'dttable_camal_chan1_27238-022318.xdr'

dtc = fltarr(2701) ; array to hold corrected counts

openr,lun,dttable,/get_lun,/xdr

readu,lun,dtc

close,/all

end

