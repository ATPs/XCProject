# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:50:01 2016

@author: k
"""

import datetime as dt
import time as tm

tm.time()

dtnow = dt.datetime.fromtimestamp(tm.time())
dtnow

dtnow.year, dtnow.month, dtnow.day, dtnow.hour, dtnow.minute, dtnow.second # get year, month, day, etc.from a datetime

today = dt.date.today()

delta = dt.timedelta(days = 100) # create a timedelta of 100 days
delta

today - delta # the date 100 days ago



