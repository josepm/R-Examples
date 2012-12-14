#! /usr/bin/env python

############################### Parameters for regressions #################################

points = 100        # number of data points
prediction = 2      # prediction points ahead
degree = 2          # polynomial degree
quantile = 0.75     # between 0 and 1
conf = 0.95         # conf level for regression
print "Input:: points: %d degree: %d quantile: %f conf level: %f" %(points, degree, quantile, conf)

##################################### Python imports ###################################################

import sys
import rpy2.robjects as robjects
import math as math

rObj = robjects.r
stats = robjects.packages.importr(name = 'stats', suppress_messages = True)     # load stats package
base = robjects.packages.importr(name = 'base', suppress_messages = True)       # load stats package
qreg = robjects.packages.importr(name = 'quantreg', suppress_messages = True)   # load quantreg package

########################################################################################################

############################################ Generate data ##############################################
x, y = [], []
for idx in range(0, points):
    x.append(idx)
    z = rObj.rnorm(1, 0, 10)
    y.append(idx ** (2 + degree) + z[0])  # anything here really

xr = robjects.FloatVector(x)
yr = robjects.FloatVector(y)

#          
# prediction data
#
xpred = []
if degree <= prediction:        # predict() must have at more points to predict than the poly degree (R issues)
    prediction = degree + 1
    print "\nPREDICTION: prediction increased to %d to avoid R issue" %prediction
for idx in range(points, points + prediction):
    xpred.append(idx)
xpredr = robjects.FloatVector(xpred)
newdata = base.data_frame(xr = xpredr)

#########################################################################################################

############################################## Formula ##################################################
model = 'poly(xr, %s, raw = T)'   %(str(degree)) 
fmla = robjects.Formula('yr ~ %s' % model)
env = fmla.environment
env['xr'] = xr
env['yr'] = yr

#
# MSE Regression
#
print "\n############## MSE REGRESSION ####################"
lmfit = stats.lm(fmla)
smry = base.summary(lmfit)
print "F statistic: %f" %smry.rx('fstatistic')[0][0]
print "R squared: %f" %smry.rx('r.squared')[0][0]
print "Adj R squared: %f" %smry.rx('adj.r.squared')[0][0]
print "Stdev: %f" %smry.rx('sigma')[0][0]
print '\n'
print "Coefficients"
coefs = smry.rx('coefficients')
print coefs
print 'Covariances'
covs = smry.rx('cov.unscaled')
print covs
print "\nPrediction"
pred = stats.predict(lmfit, newdata = newdata, se = True, interval = 'prediction', type = 'response', level = conf)
print pred.rx('fit')
print pred.rx('se.fit')
        
#
# quantile regression
#
print '\n######################## %s QUANTILE REGRESSION #####################' %quantile
qreg = robjects.packages.importr(name = 'quantreg', suppress_messages = True)   # load quantreg package
qfit = qreg.rq(fmla, tau = quantile, model = True)
smry = base.summary(qfit, se = 'boot', level = conf)     # se = boot or iid

print "Coefficients"
coefs = smry.rx('coefficients')
print coefs

print "\nPrediction"
# for quantile prediction must use:
#     type = none
#     se = boot or iid
#     interval = confidence
# Must use direct assignments to make it work

rObj.assign('xr', xpredr)
rObj('ndt = data.frame(x = xr)')
rObj.assign('s', 'boot')
rObj.assign('t', 'none')
rObj.assign('c', 'confidence')
rObj.assign('conf', conf)
rObj.assign('qfit', qfit)
pred = rObj('pred = predict(qfit, newdata = ndt, se = s, interval = c, type = t, level = conf)')
print pred
print "Done"
