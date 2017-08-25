# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 13:52:35 2016

@author: k
"""

import numpy
import peakutils
from peakutils.plot import plot as pplot
from matplotlib import pyplot as plt

folder = "E:\\Lab\\works\\20160212AliciaMSMS\\RAW\\"
files = ["1-1402_PGRP1_40hr.ms1","2-1402_DAP_PG.ms1","3-1402_LYS_PG.ms1","4-1402_PGRP1_DAP_PG_40hr.ms1",\
"5-1402_PGRP1_LYS_PG_40hr.ms1","6-1402_PGRP1_DAP_PG_70hr.ms1","7-1402_PGRP1_LYS_PG_70hr.ms1"]
filenamedic ={"1":files[0],"D":files[1],"L":files[2],"1D4":files[3],"1L4":files[4],"1D7":files[5],"1L7":files[6]}

