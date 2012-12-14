#! /usr/bin/env python

############################### Parameters for regressions #################################

# data type
type = 'trig'  # trig: sum of trig functions. 
               # step: sum of step functions

# arguments for data generation (both for trig and step)
frequencies = [10, 20, 25]
sin_amplitudes = [2, 1, 4]
cos_amplitudes = [1, -1, 1]

total_samples = 1024                                # data point to generate
time_between_samples = 1.0 / float(total_samples)   # sampling rate
ahead = 20                                          # points to predict ahead
pow_thres = 0.05                                    # threshold to filter harmonics out

# regression parameters
quantile = 0.75         # Between 0 and 1 for quantile
level = 0.95            # confidence levels
verbose = True          # print a lot or not

########################################################################################################

##################################### Python imports ###################################################

import sys
import rpy2.robjects as robjects
import math as math

rObj = robjects.r
stats = robjects.packages.importr(name = 'stats', suppress_messages = True)     # load stats package
base = robjects.packages.importr(name = 'base', suppress_messages = True)       # load stats package
qreg = robjects.packages.importr(name = 'quantreg', suppress_messages = True)   # load quantreg package

########################################################################################################

########################################## Functions ###################################################
#
# generates a data point
#
def get_value(args, time, type):
  sin_amplitude, cos_amplitude, frequency = args
  value = 0.0
  for i in range(0, len(frequency)):
    y = sin_amplitude[i] * math.sin(2.0 * math.pi * frequency[i] * time) + cos_amplitude[i] * math.cos(2.0 * math.pi * frequency[i] * time)
    if type == 'trig':
      value += y
    elif type == 'step':
        if y > 0:
            value += math.ceil(y)
    else:
        print "Invalid type: %s" %type
        quit()
  return value + rObj.rnorm(1, 0, 1)[0]

#
# generates input data
#
def set_data(args, type, total_samples, time_between_samples):
  x, y, xvalue, ymin, ymax = [], [], 0, 1000000.0, 0.0
  yhash = dict()
  yhash["ts"] = []
  for i in range(0, total_samples):
    x.append(xvalue)
    yvalue = get_value(args, x[-1], type)
    y.append(yvalue)
    yhash["ts"].append(yvalue)
    xvalue += time_between_samples
    ymin = min(ymin, yvalue)
    ymax = max(ymax, yvalue)
  return x, y, yhash, min, max
 
#
# a generic sort function for lists
# 
def sort_function(x, index, increasing):    # examples:
    if index < len(x):                      #     sort_function(x, 0, True) sorts list x by increasing in first component
        if increasing == True:              #     sort_function(x, 1, False) sorts list x by decreasing in second component
            return x[index]                 #     ...
        else:
            return -x[index]
    else:
        return x[0]

#
# increasing sort by component #1
#   
def sfunc(x):
    return sort_function(x, 1, False)
             
########################################################################################################

#
# input check and readiness
#
args = []
if len(frequencies) != len(sin_amplitudes) or len(frequencies) != len(cos_amplitudes):
  print "type: %s. Wrong parameters" %type
  quit()

args.append(sin_amplitudes)
args.append(cos_amplitudes)
args.append(frequencies)

#
# print headings
#
if verbose == True:
  print "Input Data"
  print "Sampling Interval: %f secs" %time_between_samples
  print "Total Samples %d" %total_samples

print "Frequencies: " + ", ".join(map(str, frequencies))
print "Amplitudes:: sin: " + ", ".join(map(str, sin_amplitudes)) + "  cos: " + ", ".join(map(str, cos_amplitudes))
if verbose == True:
    print "\nSpectral Analysis"

#
# generate data
#
x, y, yhash, min, max = set_data(args, type, total_samples, time_between_samples)

#
# assign input data to R objects
#
rObj.assign('xr', robjects.FloatVector(x))
rObj.assign('yr', robjects.FloatVector(y))
rObj("dt <- data.frame(xr,yr)")

#
# find the input data spectrum
#
spec = rObj("spec <- spectrum(yr, plot = FALSE)")
freq = rObj("spec$freq")    # normalized frequencies
pow = rObj("spec$spec")     # power
spec_power, spectrum = [], []
for k in range(0, len(freq)):
    spec_power.append([freq[k],pow[k]])
spec_power.sort(key = sfunc)    # sort in decreasing power
for k in range(0, len(spec_power)):
    if spec_power[k][1] > pow_thres * spec_power[0][1]:            # keep only freqs with power larger than pow_thres * max_power
        spectrum.append(spec_power[k][0] / time_between_samples)    # actual frequencies
        print "freq: %f power: %f" %(spectrum[-1], spec_power[k][1])

#       
# R assignments
#
rObj.assign('rsp', "response")
rObj.assign('prd', "prediction")
rObj.assign('se',"boot")        # boot, ker, nid, iid, rank. Only boot or iid work. needed for quantile
rObj.assign('t', "none")        # needed for quantile
rObj.assign('c',"confidence")   # needed for quantile
rObj.assign('conf', level)

#
# formula
#
model = ""
vw, vf= [], [1]
for f in spectrum:
    w = 2.0 * math.pi * f
    model += "sin(%f * xr) + cos(%f * xr) + " %(w, w)
    vw.append(w)
    vf.append('sin(%f x)' %w)
    vf.append('cos(%f x)' %w)   
model = model[:-2]  # remove the last spaces and +
fmla = robjects.Formula('yr ~ %s' % model)
rObj.assign('fmla', fmla)

#
# MSE regression
#
print '\n################# MSE REGRESSION #######################'
lmfit = rObj("lmfit <- lm(formula = fmla, data = dt)")
smry = rObj("smry <- summary(lmfit, se = se, level = conf)")
print "F statistic: %f" %smry.rx('fstatistic')[0][0]
print "R squared: %f" %smry.rx('r.squared')[0][0]
print "Adj R squared: %f" %smry.rx('adj.r.squared')[0][0]
print "Stdev: %f" %smry.rx('sigma')[0][0]
print '\n'
print "Coefficients"
print rObj("smry$coefficients")
print "Covariances"
print rObj('smry$cov.unscaled')

#
# prepare prediction array
#
xpred = []
for idx in range(total_samples, total_samples + ahead):
    xpred.append(idx)
xpredr = robjects.FloatVector(xpred)
rObj.assign('xr', xpredr)
rObj('ndt = data.frame(x = xr)')

#
# prediction
#
pred = rObj('pred = predict(lmfit, newdata = ndt, se = TRUE, interval = prd, type = rsp, level = conf)') # prediction
print "\nPredictions"
print rObj('pred$fit')

#
# QUANTILE regression
#
print '\n############################ %s QUANTILE REGRESSION #########################\n' %quantile
rObj.assign('quant', quantile)
qfit = rObj("mdl <- rq(formula = fmla, tau = quant, data = dt, model = TRUE)")
smry = rObj("smry <- summary(mdl, se = se, level = conf)")
print "Coefficients"
print rObj('smry$coefficients')
print "\nPredictions"
print rObj('pred = predict(mdl, newdata = ndt, se = se, interval = c, type = t, level = conf)') # prediction
