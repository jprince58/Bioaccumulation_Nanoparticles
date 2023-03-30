#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_experimental_data_extractor=1.0

import os
import numpy as np
import matplotlib.pyplot as plt
import docx
from docx.shared import Pt
from datetime import datetime
import matplotlib.animation as anim
from matplotlib.animation import FuncAnimation
import pandas as pd


def BCcalc(t,fit_coeff):

    to=0.5 #guess at dimensionless conversion factor
    kconv=10 #guess ay AU to particle conversion factor
    t=t*to #make dimensionless time dimensional
    out=fit_coeff[0]*t**6+fit_coeff[1]*t**5+fit_coeff[2]*t**4+fit_coeff[3]*t**3+fit_coeff[4]*t**2+fit_coeff[5]*t**1+fit_coeff[6]*t**0
    out=out/kconv #convert 
    
    return out

# %%
"""
Created on 11/2/2022

@author: joshuaprince

Purpose; Script to extract data from Tseng fits. Hmm, once I did the edited on excel this is actually a pretty simple code. Oh well. 
Version 1.0

    
"""