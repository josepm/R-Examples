#!/usr/bin/ruby

# Examples: /usr/lib/ruby/gems/1.8/gems/rsruby-0.5.1.1/examples
# Examples: :~/R/i686-pc-linux-gnu-library/2.13/quantreg/demo
# http://rubyforscientificresearch.blogspot.com/search?q=RSruby

# this example shows how to pass data between Ruby and R
# Regression followed by a prediction based on the regression
# Both mse and quantile regressions are in the example
# quantile = nil does mse and quantile between 0 and 1 does quantile regression

######################### input data ###############################
imax = 20         # number of data points
iahead = 2        # points to predict past imax
degree = 2        # polynomial degree
level = 0.95      # regression confidence level                                
quantile = nil   # regression type:  nil for mse, between 0 and 1 for quantile
print = true      # true for verbose

puts "\nParameters::\n\tdata points: #{imax}\n\tprediction points: #{iahead}\n\tdegree: #{degree}\n\tconfidence: #{level}\n\tquantile: #{quantile.inspect}"

################################ Load Gems ###########################

ENV['R_HOME'] = "/usr/lib/R"  # need to define the R_HOME

# some basic error checking when loading gems
def load_gem(a_gem)
  begin
    require a_gem
  rescue LoadError
    $stderr.print "#{a_gem} missing. Exiting"
    exit
  end
end

load_gem('rubygems')
load_gem('rsruby')
r = RSRuby.instance 

#######################################################################

######################### FUNCTIONS ###################################

# this is needed to put the fit coefficients as the order is not guaranteed
def get_coefs(mdl, degree, coefs)
  kcoefs = [coefs["(Intercept)"]]                 # polynomial coefficient keys
  if degree == 1                                  # set polynomial coefficient keys
    kcoefs[1] = coefs[mdl]
  else
    (1..degree).each {|i| kcoefs[i] = coefs[mdl + i.to_s]}
  end
  return kcoefs
end

def poly_val(coefs, x)
  y, xp = 0, 1.0
  coefs.each do |c|
    y += c * xp
    xp *= x
  end
  return y
end

def print_coefs(smry, size, print)  # prints coef stats and returns their values
  vcoefs = smry["coefficients"]
  coefs, upr_coefs, lwr_coefs = [], [], []
  puts "\tCoefficients\n\tDegree\t\tEstimate\t\t\tError\t\t\t\ttValue\t\t\t\tProb >|t|" if print == true
  (0...size).each do |d|
    puts "\t#{d}\t\t#{vcoefs[d][0]}\t\t#{vcoefs[d][1]}\t\t#{vcoefs[d][2]}\t\t#{vcoefs[d][3]}" if print == true
    coefs << vcoefs[d][0]
    upr_coefs << vcoefs[d][0] + vcoefs[d][1]
    lwr_coefs << vcoefs[d][0] - vcoefs[d][1]
  end
  return coefs, upr_coefs, lwr_coefs
end

def set_data(imax, iahead, d)
  x, y = [], [] 
  (0...imax).each do |i|
    x << i
    y << 2 * rand() * (x[-1] ** d) + rand()
  end
  return x, y
end

def do_pred(degree, x, y, level, r, imax, iahead, quantile, print)
  # assign to R objects to input data
  r.assign('xr', x)
  r.assign('yr', y)
  dt = r.eval_R("dt <- data.frame(xr,yr)")
  
  ################# model ###################
  s = quantile == nil ? "\nMSE Regression" : "\nQuantile Regression"
  puts s          # type of regression 
  r.assign('rsp', "response")
  r.assign('prd', "prediction")
  r.assign('se',"boot")  # boot, ker, nid, iid, rank. Only boot or iid work. needed for quantile
  r.assign('t', "none")       # needed for quantile
  r.assign('c',"confidence")  # needed for quantil
  model = "poly(xr, #{degree}, raw = T)"    # regression model
  
  if quantile != nil      # quantile regression
    lib = "library(\'quantreg\')"
    r.eval_R("suppressPackageStartupMessages(#{lib})")          
    r.assign('quant', quantile)
    fit = r.eval_R("mdl <- rq(formula = yr ~ #{model}, tau = quant, data = dt, model = TRUE)")
    smry = r.eval_R("smry <- summary(mdl,se = se, level = #{level})")
    pred1 = r.eval_R("p1 <- predict(mdl, newdata = dt, type = t, level = #{level}, interval = c, se = se)")  # works only for se = boot, iid, 
  else                      # mse regression
    fit = r.eval_R("mdl <- lm(formula = yr ~ #{model},data = dt)") 
    smry = r.eval_R("smry <- summary(mdl,se = se, level = #{level})")
    pred1 = r.eval_R("p1 <- predict(mdl, se = TRUE, newdata = dt, type = rsp, level = #{level}, interval = prd)")
  end
  
  coefs, upr_coefs, lwr_coefs = print_coefs(smry, 1 + degree, print)  # print and extract coefs
  puts "\nF-statistic: #{smry["fstatistic"]["value"].inspect}" if quantile == nil and print == true
  
  ################### Model validation #########################
  puts "\n\nData and Model Fit" if print == true
  puts "\ty: actual data" if print == true
  puts "\tyh: model data directly computed." if print == true
  puts "\tyR, yRupr, yRlwr: model data computed by R" if print == true
  (0...imax).each do |i|
    yh = poly_val(coefs, x[i])
    if quantile == nil          # mse
      yR, yRlwr, yRupr = pred1["fit"][i]
      se = pred1["se.fit"][(1+i).to_s]
    else                        # quantile
      yR, yRlwr, yRupr = pred1[i]
      se = (yRupr - yRlwr) / 2.0
    end
    puts "\t#{i}  x: #{x[i]}  y: #{y[i]}  yh: #{yh} yR: #{yR} yRlwr: #{yRlwr} yRupr: #{yRupr} se: #{se}" if print == true
  end
  
  ################# Prediction ##################
  x = []
  iahead = degree + 1 if degree <= iahead        # predict() must have at more points to predict than the poly degree (R issues)
  (imax...imax + iahead).each {|idx| x << idx}
  r.assign('xr', x)
  ndt = r.eval_R("ndt <- data.frame(xr)") # define the prediction data range
  if quantile == nil    # mse
    pred2 = r.eval_R("p2 <- predict(mdl, newdata = ndt, se = TRUE, type = rsp, level = #{level}, interval = prd)")
  else                  # quantile
    pred2 = r.eval_R("p2 <- predict(mdl, newdata = ndt, type = t, level = #{level}, interval = c, se = se)")
  end
 
  puts "\nPredictions"
  (imax...imax + iahead).each do |idx|
    xval = idx
    if quantile == nil    # mse
      yR, yRlwr, yRupr = pred2['fit'][idx - imax]
      se = pred2["se.fit"][(1 + idx - imax).to_s]
    else    # quantile
      yR, yRlwr, yRupr = pred2[idx  - imax]
      se = (yRupr - yRlwr) / 2.0
    end
    #yh = poly_val(coefs, xval)
    puts "\t#{idx}  x: #{xval}  y: N/A  yR: #{yR} yRlwr: #{yRlwr} yRupr: #{yRupr} se: #{se}"
  end
end

########################################################################################################

############################################## MAIN ####################################################
x, y = set_data(imax, iahead, degree)
do_pred(degree, x, y, level, r, imax, iahead, quantile, print)
puts "Done"

########################################################################################################
