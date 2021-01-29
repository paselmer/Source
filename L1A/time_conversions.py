import numpy as np
import datetime as DT
import pdb


def is_leap(yr):
    """
    Simple function to determine if input year is a leap year
    according to the Gregorian calendar. Enter a year as an interger.
    Example: 2017
    Returns a 1 if it is a leap year, 0 if not.
    """
    
    lp = 0 # leap year? 1=yes, 0=no
    
    if yr % 4 == 0:
        if ((yr % 100 == 0) and (yr % 400 == 0)):
            lp = 1
        elif ((yr % 100 == 0) and (yr % 400 != 0)):
            lp = 0
        else:
            lp = 1
    
    return(lp)


def GPSWeek2UTC(GPSW,GPSms):
    """Function to convert GPS week into UTC"""
    
    # !!!!!!!!!!! READ !!!!!!!!!!!!!!!!!!!
    # Since I found a code online that did the exact same thing,
    # I've ceased work on this code as of 11/16/17. I may want to complete
    # this code at some point in the future, because its conversion does
    # not rely on fancy Python libraries like the one I found online. 
    
    # Initialize some variables
    
    DinY = 365
    ldpmon = np.zeros((2,12))
    ldpmon[0,:] = np.array([0,31,59,90,120,151,181,212,243,273,304,334])  #normal
    ldpmon[1,:] = np.array([0,31,60,91,121,152,182,213,244,274,305,335])  #leap
    months = np.array([1,2,3,4,5,6,7,8,9,10,11,12])
    leap = 0 # will change to 1 if leap year determined
    
    # The anchor to perform conversion
    # Converted on https://www.labsat.co.uk/index/php/en/gps-time-calculator
    
    anchorY = 2017
    anchorW = 1951 # = 28 May 2017, 00:00:00 UTC
    anchorDoY = 148 # May 28's day of year (JDay, as it were)
    anchorD2E = DinY - anchorDoY
    leaps_since_anchorY = 0
    
    # Compute total # of days since beginning of anchor year
    
    mindays = (GPSW - anchorW) * 7.0
    ndays = mindays + int( GPSms/(1e3*86400.0) )
    # The # of input secs minus the number of secs in the integer # of days
    GPSsecleft = (GPSms/1e3) - float(int( GPSms/(1e3*86400.0) ))*86400.0
    ndays_since_anchorY = ndays+(anchorDoY-1)
    
    # Compute the total number of years since beginning of anchor year
    
    years = int(ndays_since_anchorY/DinY)
    year = int(years + anchorY) #just make sure it's int for next step
    
    # Compute the number of leap years since the anchor year
    
    if years >= 1:
        print("More than 1 year since anchor year")
        for testyear in range(anchorY,year):
            testleap = is_leap(testyear)
            leaps_since_anchorY += testleap
    
    # Determine if current year is a leap year
    
    leap = is_leap(year)
    
    # Use the remaining days to figure out the month
    
    # This DoY, or "day of year" starts at zero. Zero being the first day.
    DoY = ndays_since_anchorY - (years*DinY + leaps_since_anchorY) 
    if DoY > 0:
        month_sel_mask = ldpmon[leap,:] <= DoY
        month = max(months[month_sel_mask])
    else:
        print("Retroactively correcting for intervening leap years.")
        year = year - 1
        leap = is_leap(year)
        DoY = (DinY+leap) + DoY
        pdb.set_trace()
    
    print(month,leaps_since_anchorY)
    return(None)
    

def weeksecondstoutc(gpsweek,gpsseconds,leapseconds):
    """ I, Patrick Selmer, did not write this. I found it on Github at the
        following URL:
        https://gist.github.com/jeremiahajohnson/eca97484db88bcf6b124
    """
    
    import datetime, calendar
    datetimeformat = "%Y-%m-%d %H:%M:%S"
    epoch = datetime.datetime.strptime("1980-01-06 00:00:00",datetimeformat)
    elapsed = datetime.timedelta(days=(gpsweek*7),seconds=(gpsseconds+leapseconds))
    #print(type(epoch+elapsed))
    #return datetime.datetime.strftime(epoch + elapsed,datetimeformat)
    return(epoch+elapsed)


def delta_datetime_vector(dt_vect):
    """ This function will take a numpy 1D array of datetime.datetime objects
        and return a 1D float64 array of the seconds between records.
        Returned array will have same size as input array; therefore first
        element of returned array will always be equal to zero.
    """
    
    if len(dt_vect.shape) > 1:
        print("Array must be 1D!")
        print("Input array has dimensions (shape) of ",len(dt_vect))
        pdb.set_trace()

    nr = dt_vect.shape[0]
    del_t = np.zeros(nr,dtype=DT.datetime)
    del_t[1:] = dt_vect[1:] - dt_vect[:-1]
    del_t[0] = del_t[1]*0.0
    
    del_secs = np.zeros(nr,dtype=np.float64)
    for k in range(0,nr): del_secs[k] = del_t[k].total_seconds()
    
    return del_secs
    

def gpstime(gpstime, msec, direction):
    """ Direct translation of Steve Palm's gpstime IDL procedure used
        in the CATS-ISS L1A process.
    
        Converts GPS time to Year,Month,Day,HH:MM:SS, and Fractional Day
        
        Output list:
        [hour, minute, second, month, day, year, fractional_day]
    """
    
    #************
    #
    # Routine to convert Day, Month, Year, HOur, Minute, Second to GPS time (direction = 1)
    # Routine to convert GPS time to Day, Month, Year, Hour, Minute, Second (direction = 2)
    # The Year must be >= 2000 and <= 2020
    #
    #************
    
    # The GPS time defined by "anchor" is 1/1/2009 at 00:00:00 UTC
    anchor = 914803215
    anchor_year = 2009
    leap_seconds_since_2009 = 0
    if gpstime > 1025136013: leap_seconds_since_2009 = 1 # June 30, 2012, 23:59:58 UTC
    if gpstime > 1119744014: leap_seconds_since_2009 = 2 # June 30, 2015, 23:59:58 UTC
    if gpstime > 1167264015: leap_seconds_since_2009 = 3 # December 31, 2016, 23:59:58 UTC
    hour = -1
    minute = -1
    second = -1
    month = -1
    day = -1
    year = -1
    fractional_day = -999.99
    null_list = [hour, minute, second, month, day, year, fractional_day]
    
    # define days in a year for years 2000 to 2020
    days_year = [366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366]

    # define days in months
    days_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    if direction == 1: # UTC time to GPS time
        
        if (year < 2000) or (year > 2020): 
            print('Returning nothing. Year out-of-bounds for CATS.')
            return null_list
            
        fraction_day = 0.0
        
        seconds_in_day = 24 * 3600
        year_delta = year - anchor_year
        
        if year >= anchor_year:
            i1 = anchor_year - 2000
            i2 = year - 2000 - 1
            
            sum1 = 0
            for i in range(i1,i2+1):
                sum1 = sum1 + (seconds_in_day * days_year[i])
                
            sum2 = sum2 * seconds_in_day
            sum3 = (day - 1) * seconds_in_day
            sum4 = hour * 3600 + minute*60 + second
            
            gpstime = anchor + sum1 + sum2 + sum3 + sum4
            
            leap_year = 0
            if (year % 4) == 0: leap_year = 1
            print('leap_year = ', leap_year)
            if leap_year and (month > 2): gpstime = gpstime + seconds_in_day
            
            if (year == 2012) and (month > 6): gpstime += 1 # leap second on 6/30 23:59:60  2012
            if (year > 2012): gpstime += 1   # leap second 6/30 23:59:60  2015
            if (year == 2015) and (month > 6): gpstime += 1 # leap second on 6/30 23:59:60  2015
            if (year > 2015): gpstime += 1  # leap second on 6/30 23:59:60  2015
            
        else:
            
            i1 = year - 2000 + 1
            i2 = anchor_year - 2000 - 1
            
            sum1 = 0
            for i in range(i1,i2+1):
                sum1 = sum1 + (seconds_in_day * days_year[i])
                
            i1 = month
            i2 = 11
            sum2 = 0
            for i in range(i1,i2+1):
                sum2 = sum2 + days_month[i]
                
            sum2 = sum2 * seconds_in_day
            days = days_month[month-1] - day
            sum3 = days * seconds_in_day
            sum4 = hour * 3600 + minute * 60 + second
            sum4 = seconds_in_day - sum4
            
            gpstime = anchor - sum1 - sum2 - sum3 - sum4
            if (year == 2015) and (month > 6): gpstime += 1 # leap second on 6/30 23:59:60  2015
            if (year > 2015): gpstime += 1  # leap second on 6/30 23:59:60  2015
            
            
    if direction == 2: # GPS time to year, month, day, fraction day
        
        tdelta = gpstime - anchor - leap_seconds_since_2009  # seconds from 1/1/2009, 00:00:00
        hours = int(tdelta / 3600)
        seconds = int(tdelta - (hours * 3600))
        days = int(hours / 24)    # total days since 1/1/2009, 00:00:00
        
        if days > 3660:
            print('Problem in gpstime. gpstime = ', gpstime)
            print('days = ', days)
            year = -1
            return null_list
            
        i = 9
        sum = 0
        while (sum <= days) and  (i < 20):
            sum = sum + days_year[i]
            i += 1
            
        # print('days, sum, i = ', days, sum, i)
        
        hour = int( hours - (days * 24) )
        years = int( days / 365 )
        y2 = float(days) / 365.140
        years = int(y2)
        
        # print('years = ', years , y2)
        day = days_year[i-1] - sum + days + 1 # julian day
        
        #print(sum,days,days_year[i-1])
        #print('day = ',day)
        
        minute = int( seconds / 60 )
        second = seconds - ( minute * 60 )
        
        fractional_day = (float(hour)*3600.0 + float(seconds) + float(msec)/1000.0) / 86400.0
        #print(hour,seconds,msec,fraction_day)
        
        year = 2009 + years
        if (year == 2012) or (year == 2016) or (year == 2020): days_month[1]=29
        
        i = 0
        day_sum = 0
        while (day_sum < day):
          day_sum = day_sum + days_month[i]
          i += 1
          
        #print(i,day_sum,days_month[i-1])
        
        if (i > 0): day = day - (day_sum - days_month[i-1])
        if (i == 0): i = 1
        month = i
        
        if (month == 3) and (day == 0):
            day = 29
            month -= 1
            
    return [hour, minute, second, msec, month, day, year, fractional_day]
                        
            

    
    
    
# Commented out code below is for testing of above conversion algorithms 

### Test "gpstime" function ###
# 914803215 = 1/1/2009  00:00:00 UTC
# 1072569616 = 1/1/2014
# 1261872016 = 1/1/2020
# 1119744000 = 2015-06-30T23-59-44
# 1041836952 = 1/10/2013, 07:08:56
# 1188251359 = 8/31/2017, 21:49:01
# 1113096818 = 4/15/2015, 01:33:22
#GPS_time = 1113096818 #Raw_CATS.Data_Coarse_Time
#msec = 0            #long(Raw_CATS.Fine_Time) * 20L				; This is micro seconds
#msec = int(msec / 1000)	  #This is milliseconds
#r = gpstime(GPS_time, msec, 2)
#print(r)
    
#ans = weeksecondstoutc(1811,164196.732,16) ## --> '2014-09-22 21:36:52'
#ans = weeksecondstoutc(2086,259200.0,16) ## --> '2020-01-01 00:00:00'
#print(ans)
#print(type(ans))

#something = GPSWeek2UTC(2086,259200000) # 01 Jan 2020, 00:00:00 UTC
#something = GPSWeek2UTC(1968,248000000)
#something = GPSWeek2UTC(1982,248000000)
    
    
    
