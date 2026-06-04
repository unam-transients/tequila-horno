from datetime import datetime

import horno.image
import math


def exposuretime(header):
    return float(header["EXPTIME"])


def detectortemperature(header):
    return 0.5 * (float(header["SDTTM"]) + float(header["SDTTM"]))


def datamax(header):
    return 4095

def offsetxslice(header):
    return slice(4110, 4116)

def offsetyslice(header):
    return slice(0, 2998)

def trimyslice(header):
    return slice(0, 2998)

def trimxslice(header):
    return slice(4, 4110)

def trimyslice(header):
    return slice(0, 2998)

def flatmax(header):
    return 3000

def gain(header):
    return 2.6

def readnoise(header):
    return 2.7