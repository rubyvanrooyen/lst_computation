#!/usr/bin/python
"""
Computation of Sidereal times from date/time provided by user

Author: R van Rooyen
Date: December 2013
Revisions:
2014-01-02 : R van Rooyen -- Wraps datetime.datetime class to extend functionality to include LST
"""

import datetime
import time
import numpy


# Class SiderealTime wraps the datetime.datetime object and extend functionality to include LST dates
# normal datetime object attributes + lst variates
class SiderealTime():
  def __init__(self, *args, **kwargs):
    if len(args) > 0: self.date_object = datetime.datetime(*args, **kwargs)
    else:             self.date_object = datetime.datetime.utcfromtimestamp(self.utc2ut1(time.time()))
    self.s_day = 23.9344695917 # hours
    self.m_day = 24.065709825  # hours
    self.m_angle = 360./self.s_day  # degrees

  def __getattr__(self,name):
    return getattr(self.date_object, name)


  ##Astronomical Calendar (Explanatory Supplement) P 604
  # Converting between Gregorian Calendar Date and Julian Day
  def ymd2jd (self, year, month, day):
    """
    Converts a calender date (year, month, day) to a Julian day number

    @param  self   This object
    @param  year
    @param  month
    @param  day

    @return julian_day [integer day value reference to noon]
    """
    julian_day = int((1461*(year + 4800 + int((month - 14)/12.)))/4.) \
               + int((367*(month - 2 - 12*int((month - 14)/12.)))/12.) \
               - int((3*int((year + 4900 + int((month - 14)/12.))/100.))/4.) \
               + day - 32075 # 12h based
    return julian_day

  ##ftp://maia.usno.navy.mil/ser7/ser7.dat
  # Universal Time (UT1)
  def utc2ut1(self, UTC_timestamp):
    """
    Approximate dT = UTC-UT1 [seconds]

    @param  self          This object
    @param  UTC_timestamp Linux timestamp number of seconds since unix epoch as seconds in UTC

    @return dT correction [seconds]
    """
    date_obj = datetime.datetime.utcfromtimestamp(UTC_timestamp)
    JD = self.ymd2jd(date_obj.year, date_obj.month, date_obj.day) # days -- starting at noon
    MJD = JD - 2400000.5 # days -- starting at midnight
    t = 2000.0 + (MJD - 51544.03) / 365.2422
    # t is the date in Besselian day fraction
    UT2_minus_UT1 = 0.022*numpy.sin(2*numpy.pi*t) - 0.012*numpy.cos(2*numpy.pi*t) - 0.006*numpy.sin(4*numpy.pi*t) + 0.007*numpy.cos(4*numpy.pi*t)
    # DUT1= (UT1-UTC) = -0.0614 - 0.00107 (MJD - 56625) - (UT2-UT1)
    DUT1 = -0.0614 - 0.00107*(MJD - 56625) - UT2_minus_UT1 # seconds
    return (UTC_timestamp+DUT1)

  ##Astronomical almanac B3-B10
  def eqeq(self, epoch=2000):
    """
    Equation of the equinoxes: Eq = GMST - GAST [seconds]

    @param  self   This object
    @param  epoch  Julian epoch, default assumed to be 2000 Jan 1d, 12h UT1

    @return eqeq   Earth rotation correction to mean sidereal [seconds]
    """
    # Assume epoch 2000 January 1d 12h UT1
    epoch = self.ymd2jd(epoch,1,1)
    # GMST at 0h UT
    JD = self.ymd2jd(self.date_object.year, 1, 0)    # beginning of the year on day 0 at 12h (noon)
    Du = (JD - 0.5) - epoch # at 0h (midnight)
    D=self.date_object.timetuple().tm_yday
    H = (self.date_object.hour*3600+self.date_object.minute*60+(self.date_object.second+self.date_object.microsecond/1e6))/3600.
    T=((Du*24*3600+67)/3600./24. + D + 24./self.s_day*H)/36525. #UT
    e = numpy.deg2rad((23.4393 - 0.0000004*T)%360)         # obliquity
    L = numpy.deg2rad((280.4665 + 36000.7698*T)%360)       # mean longitude of the Sun
    L1 = numpy.deg2rad((218.3165 + 481267.8813*T)%360)     # mean longitude of the Moon
    Om = numpy.deg2rad((125.04452 - 1934.136261*T)%360)    # longitude of the ascending node of the Moon
    delta_phi = -17.2*numpy.sin(Om) -1.32*numpy.sin(2*L) -0.23*numpy.sin(2*L1) +0.21*numpy.sin(2*Om) # nutation in longitude
    eqeq = delta_phi*numpy.cos(e)/15. # seconds
    return eqeq

  def gmst(self, epoch=2000): # default epoch is 2000, but can be changed using epoch optional argument
    """
    Greenwich Mean Sidereal Time [hours]

    @param  self   This object
    @param  epoch  Julian epoch, default assumed to be 2000 Jan 1d, 12h UT1

    @return gmst   Greenwich mean sidereal time [hours]
    """
    # Assume epoch 2000 January 1d 12h UT1
    epoch = self.ymd2jd(epoch,1,1)
    # GMST at 0h UT
    JD = self.ymd2jd(self.date_object.year, 1, 0)    # beginning of the year on day 0 at 12h (noon)
    Du = (JD - 0.5) - epoch # at 0h (midnight)
    D=self.date_object.timetuple().tm_yday
    H = (self.date_object.hour*3600+self.date_object.minute*60+(self.date_object.second+self.date_object.microsecond/1e6))/3600.
    ERA = 24.*((0.7790572732640 + \
                0.00273781191135448*(Du*self.m_day%self.m_day)*self.m_angle) + \
                0.00273781191135448*D + (1./self.s_day)*H) %24.
    T=((Du*24*3600+67)/3600./24. + D + 24./self.s_day*H)/36525. #UT
    part = (0.014506 + 4612.156534*T + 1.3915817*T*T - 0.00000044*T**3 - 0.000029956*T**4 - 3.68e-8*T**5)/3600./15.
    GMST = (ERA+part)%24
    return GMST

  def gast(self, epoch=2000):
    """
    Greenwich Apparent Sidereal Time [hours]

    @param  self   This object
    @param  epoch  Julian epoch, default assumed to be 2000 Jan 1d, 12h UT1

    @return gast   Greenwich apparent sidereal time [hours]
    """
    eqeq = self.eqeq(epoch=epoch)
    print "Equation of the equinoxes correction = ", eqeq
    return (self.gmst(epoch=epoch) + eqeq/3600.)

  # Assume East longitude positive and epoch 2000 for gmt
  def lmst(self, longitude, epoch=2000):
    """
    Local Mean Sidereal Time [hours]

    @param  self      This object
    @param  longitude Longitude of location (floating point value with East as possitive)
    @param  epoch     Julian epoch, default assumed to be 2000 Jan 1d, 12h UT1

    @return lmst   Local mean sidereal time [hours]
    """
    gmst = self.gmst(epoch=epoch)
    lmst = gmst + longitude/15.
    if lmst < 0: lmst += 24
    return lmst

  def last(self, longitude, epoch=2000):
    """
    Local Apparent Sidereal Time [hours]

    @param  self      This object
    @param  longitude Longitude of location (floating point value with East as possitive)
    @param  epoch     Julian epoch, default assumed to be 2000 Jan 1d, 12h UT1

    @return last   Local apparent sidereal time [hours]
    """
    gast = self.gast(epoch=epoch)
    last = gast + longitude/15.
    if last < 0: last += 24
    return last


def sec2hms(seconds):
  """ return seconds as (hh,mm,ss)"""
  hr = int(seconds)
  mn = int((seconds*3600 - hr*3600)/60)
  sec = seconds*3600 - hr*3600 - mn*60
  return [hr, mn, sec]

if __name__ == '__main__':

  ## get sidereal times for specified datetime
  # date='2012-07-08' time='09:44:30' lon=-80.3821638889 assuming East as positive
  d = SiderealTime(2012,7,8,9,44,30)
  print
  print 'Date', d.strftime('%Y-%m-%d %H:%M:%S')
  print 'Day of week', d.weekday()
  print 'GMST %f [hour]' % d.gmst()
  [h, m, s] = sec2hms(d.gmst())
  print 'GMST at input UT: %dh %dm %f' % (h, m, s)
  # if local longitude is provided, compute LMST
  lon=-80.3821638889
  [h, m, s] = sec2hms(abs(lon/15.))
  if lon < 0: print 'Local longitude: %dh %dm %f West' % (h, m, s)
  else: print 'Local longitude: %dh %dm %f East' % (h, m, s)
  print 'LMST %f [hour]' % d.lmst(lon)
  [h, m, s] = sec2hms(d.lmst(lon))
  print 'LMST at input UT: %dh %dm %f' % (h, m, s)
  print 'GAST %f [hour]' % d.gast()
  [h, m, s] = sec2hms(d.gast())
  print "GAST at required UT: %dh %dm %f" % (h, m, s)
  [h, m, s] = sec2hms(d.last(lon))
  print 'LAST at required UT: %dh %dm %f' % (h, m, s)

  ## get sidereal times of "now"
  # uct from number of seconds since unix epoch as seconds in UTC
  d = SiderealTime(epoch=2000)
  print
  print 'Date', d.strftime('%Y-%m-%d %H:%M:%S')
  print 'Day of week', d.weekday()
  [h, m, s] = sec2hms(d.gmst())
  print 'GMST at input UT: %dh %dm %f' % (h, m, s)
  # if local longitude is provided, compute LMST
  lon=21.41069444444
  [h, m, s] = sec2hms(abs(lon/15.))
  if lon < 0: print 'Local longitude: %dh %dm %f West' % (h, m, s)
  else: print 'Local longitude: %dh %dm %f East' % (h, m, s)
  [h, m, s] = sec2hms(d.lmst(lon))
  print 'LMST at input UT: %dh %dm %f' % (h, m, s)
  [h, m, s] = sec2hms(d.gast())
  print "GAST at required UT: %dh %dm %f" % (h, m, s)
  [h, m, s] = sec2hms(d.last(lon))
  print 'LAST at required UT: %dh %dm %f' % (h, m, s)

# -fin-


