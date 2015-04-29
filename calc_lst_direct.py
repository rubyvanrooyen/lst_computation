#!/usr/bin/python

import optparse

import datetime
import time
import numpy

# Compute LMST and LAST of a datetime object

def sec2hms(seconds):
  hr = int(seconds)
  mn = int((seconds*3600 - hr*3600)/60)
  sec = seconds*3600 - hr*3600 - mn*60
  return [hr, mn, sec]

##Astronomical Calendar (Explanatory Supplement) P 604
# Converting between Gregorian Calendar Date and Julian Date and Day number
def ymd2JD (year, month, day, hour=0.0, minute=0.0, second=0.0):
  """
  Converts a calender date (year, month, day) to a Julian date, number of days and fraction
  """
  julian_day = ((1461*(year + 4800 + ((month - 14)/12.)))/4.) \
             + ((367*(month - 2 - 12*((month - 14)/12.)))/12.) \
             - ((3*((year + 4900 + ((month - 14)/12.))/100.))/4.) \
             + day - 32074.5 # 0h based
  julian_hour = (hour*3600. + minute*60. + second)/3600./24.
  julian_date = julian_day + julian_hour
  return julian_date
def ymd2jd (year, month, day):
  """
  Converts a calender date (year, month, day) to a Julian day number
  """
  julian_day = int((1461*(year + 4800 + int((month - 14)/12.)))/4.) \
             + int((367*(month - 2 - 12*int((month - 14)/12.)))/12.) \
             - int((3*int((year + 4900 + int((month - 14)/12.))/100.))/4.) \
             + day - 32075 # 12h based
  return julian_day

##ftp://maia.usno.navy.mil/ser7/ser7.dat
# Universal Time (UT1)
def utc2ut1(UTC_timestamp):
  """
  Approximate dT = UTC-UT1 [seconds]
  """
  date_obj = datetime.datetime.utcfromtimestamp(UTC_timestamp)
  JD = ymd2jd(date_obj.year, date_obj.month, date_obj.day) # days -- starting at noon
  MJD = JD - 2400000.5 # days -- starting at midnight
  t = 2000.0 + (MJD - 51544.03) / 365.2422
  # t is the date in Besselian day fraction
  UT2_minus_UT1 = 0.022*numpy.sin(2*numpy.pi*t) - 0.012*numpy.cos(2*numpy.pi*t) - 0.006*numpy.sin(4*numpy.pi*t) + 0.007*numpy.cos(4*numpy.pi*t)
  # DUT1= (UT1-UTC) = -0.0614 - 0.00107 (MJD - 56625) - (UT2-UT1)
  DUT1 = -0.0614 - 0.00107*(MJD - 56625) - UT2_minus_UT1 # seconds
  return (UTC_timestamp+DUT1)

##Astronomical almanac B3-B10
s_day = 23.9344695917 # hours
m_day = 24.065709825  # hours
m_angle = 360./s_day  # degrees
# Computing sidereal times
def gmst(dt, epoch):
  """
  Greenwich Mean Sidereal Time [hours]
  """
  # GMST at 0h UT
  JD = ymd2JD(dt.year, 1, 0)    # beginning of the year on day 0 at 12h (noon)
  Du = (int(JD) - 0.5) - float(epoch) # at 0h (midnight)
  D=dt.timetuple().tm_yday
  H = (dt.hour*3600+dt.minute*60+(dt.second+dt.microsecond/1e6))/3600.
  ERA = 24.*((0.7790572732640 + 0.00273781191135448*(Du*m_day%m_day)*m_angle) + 0.00273781191135448*D + (1./s_day)*H) %24.
  T=((Du*24*3600+67)/3600./24. + D + 24./s_day*H)/36525. #UT
  part = (0.014506 + 4612.156534*T + 1.3915817*T*T - 0.00000044*T**3 - 0.000029956*T**4 - 3.68e-8*T**5)/3600./15.
  GMST = (ERA+part)%24
  return GMST

def eqeq(dt, epoch):
  """
  Equation of the equinoxes: Eq = GMST - GAST
  """
  # GMST at 0h UT
  JD = ymd2JD(dt.year, 1, 0) # beginning of the year on day 0 at 12h (noon)
  Du = (int(JD) - 0.5) - 2451545.0 # at 0h (midnight)
  D=dt.timetuple().tm_yday
  H = (dt.hour*3600+dt.minute*60+(dt.second+dt.microsecond/1e6))/3600.
  T=((Du*24*3600+67)/3600./24. + D + 24./s_day*H)/36525. # UT
  e = numpy.deg2rad((23.4393 - 0.0000004*T)%360)         # obliquity
  L = numpy.deg2rad((280.4665 + 36000.7698*T)%360)       # mean longitude of the Sun
  L1 = numpy.deg2rad((218.3165 + 481267.8813*T)%360)     # mean longitude of the Moon
  Om = numpy.deg2rad((125.04452 - 1934.136261*T)%360)    # longitude of the ascending node of the Moon
  delta_phi = -17.2*numpy.sin(Om) -1.32*numpy.sin(2*L) -0.23*numpy.sin(2*L1) +0.21*numpy.sin(2*Om) # nutation in longitude
  eqeq = delta_phi*numpy.cos(e)/15. # seconds
#   de = 9.2*numpy.cos(Om) + 0.57*numpy.cos(2*L) + 0.1*numpy.cos(2*L1) - 0.09*numpy.cos(2*Om)
#   eqeq = delta_phi*numpy.cos(e+de)/15. # seconds
  return eqeq



if __name__ == '__main__':
# ./calc_lst.py --date='2013-12-20' --time='11:53:59' --lon=21.41069444444 --apparent
# ./calc_lst.py --now --lon=21.41069444444 --apparent
  parser = optparse.OptionParser(usage=' \
\n  %prog --date <YYYY-MM-DD> [options] \nor \n  %prog --now \
\n\nExamples: \
\n  %prog calc_lst.py --date=\'2012-07-08\' \
\n  %prog calc_lst.py --date=\'2012-07-08\' --time=\'09:44:30\' --lon=-80.3821638889 --apparent \
\n  %prog calc_lst.py --now', \
                                 version="%prog 1.0")
  parser.add_option('--date',
                    action='store',
                    dest='date',
                    type=str,
                    default=None,
                    help='Date in format \'YYYY-MM-DD\'')
  parser.add_option('--time',
                    action='store',
                    dest='time',
                    type=str,
                    default=None,
                    help='Time (SAST) of day in format \'HH:MM:SS\'')
  parser.add_option('--lon',
                    action='store',
                    dest='lon',
                    type=float,
                    default=0.,
                    help='Local longitude \'degrees of arc\' East as positive')
  parser.add_option("--now",
                    dest='now',
                    action="store_true",
                    default=False,
                    help="Sidereal time for current date and time.")
  parser.add_option("--apparent",
                    dest='apparent',
                    action="store_true",
                    default=False,
                    help="Apparent sidereal times: GAST and LAST")
  (opts, args) = parser.parse_args()

  if not opts.now and opts.date is None: raise SystemExit(parser.print_help())

  dt = None
  if opts.now:
    # uct from number of seconds since unix epoch as seconds in UTC
    dt = datetime.datetime.utcfromtimestamp(utc2ut1(time.time()))
    print 'UTC = \t', dt
  elif opts.date is not None:
    if opts.time is not None: date_obj = datetime.datetime.strptime('%s %s' %(opts.date, opts.time), "%Y-%m-%d %H:%M:%S")
    else: date_obj = datetime.datetime.strptime(opts.date, "%Y-%m-%d")
    print "User defined date:", datetime.datetime.fromtimestamp(time.mktime(date_obj.timetuple()))
    dt = datetime.datetime.utcfromtimestamp(time.mktime(date_obj.timetuple()))

  # Assume epoch 2000 January 1d 12h UT1
  epoch = ymd2jd(2000,1,1)
  GMST = gmst(dt, epoch)
  [h, m, s] = sec2hms(GMST)
  print 'GMST at required UT: %dh %dm %f (%f)' % (h, m, s, GMST)

  # if local longitude is provided, compute LMST
  if (opts.lon != 0.0):
    [h, m, s] = sec2hms(abs(opts.lon/15.))
    if opts.lon < 0: print 'Local longitude: %dh %dm %f West' % (h, m, s)
    else: print 'Local longitude: %dh %dm %f East' % (h, m, s)

    LMST = GMST + opts.lon/15.
    if LMST < 0: LMST += 24
    [h, m, s] = sec2hms(LMST)
    print 'LMST at required UT: %dh %dm %f (%f)' % (h, m, s, LMST)

  if opts.apparent:
    eqeq = eqeq(dt, epoch)
    print "Equation of the equinoxes correction = ", eqeq

    GAST = GMST + eqeq/3600.
    [h, m, s] = sec2hms(GAST)
    print "GAST at required UT: %dh %dm %f (%f)" % (h, m, s, GAST)

    LAST = GAST + opts.lon/15.
    if LAST < 0: LAST += 24
    [h, m, s] = sec2hms(LAST)
    print 'LAST at required UT: %dh %dm %f (%f)' % (h, m, s, LAST)

# -fin-


