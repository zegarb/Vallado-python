#
# calculate algorithm 34 quantities
#

import numpy as np
import spacemath_utils as smu
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import re


rs = 696000.0
au = 149597870.0
angumb = np.arctan((rs - re) / au)
angpen = np.arctan((rs + re) / au)
reci1 = np.array([- 41221.79149309,8864.59854079,0.0])
veci1 = np.array([- 0.646416796,- 3.005940793,- 0.0])
# +50 and -80 seem to work here to get the proper angles
dtsec = 21000
reci,veci,error = obu.kepler(reci1,veci1,dtsec)
year = 2008
mon = 3
day = 16
hr = 6
min = 13
sec = 0.0
jd,jdfrac = stu.jday(year,mon,day,hr,min,sec)
rsun,rtasc,decl = obu.sun(jd + jdfrac)
np.dot(reci,rsun)
umbvert = 0.0
penvert = 0.0
umb = 'n'
pen = 'n'
# now for algorihtm 34
if np.dot(reci,rsun) < 0.0:
    ang1 = smu.angl(- rsun,reci)
    sathoriz = smu.mag(reci) * np.cos(ang1)
    satvert = smu.mag(reci) * np.sin(ang1)
    x = re / np.sin(angpen)
    penvert = np.tan(angpen) * (x + sathoriz)
    if satvert <= penvert:
        pen = 'y'
        y = re / np.sin(angumb)
        umbvert = np.tan(angumb) * (y - sathoriz)
        if satvert <= umbvert:
            umb = 'y'

print(' %11.7f  %11.4f  %11.4f  %11.4f  %11.4f U %c  P %c \n' % (ang1 * 180.0 / np.pi,sathoriz,satvert,penvert,umbvert,umb,pen))
s = 2.0 * smu.mag(reci) * ang1
print(' %11.7f  %11.7f \n' % (s,(s / smu.mag(veci)) / 60.0))
