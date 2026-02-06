from datetime import datetime

import horno.image
import math

def exposuretime(header):
    return float(header["EXPTIME"])

def datamax(header):
    return 4095

def trimxslice(header):
    return slice(4,4109)

def trimyslice(header):
    return slice(0,2997)

def flatmax(header):
    return 3000