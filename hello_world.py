import numpy as np
import datetime as datetime
import matplotlib.pyplot as plt

print "hello world"
tmp = np.arange(1, 366, 1)

year = 1996
print 'tmp =', tmp
# want to generate an array that gives the month, given year, doy
days = tmp

f = [(datetime.datetime(year, 1, 1) + datetime.timedelta(days - 1)).month for days in tmp]
print f

plt.plot(tmp, f)
plt.show()
