#!/usr/bin/ruby

# Examples: /usr/lib/ruby/gems/1.8/gems/rsruby-0.5.1.1/examples
# Examples: :~/R/i686-pc-linux-gnu-library/2.13/quantreg/demo
# http://rubyforscientificresearch.blogspot.com/search?q=RSruby
#
# this example shows how to pass data from Ruby to R
# Linear regression using trigonometric functions
# First find spectrum, then regress and finally predict
#
# Both mse and quantile regressions are in the example
# Set quantile = nil to do mse regression
# Set quantile between 0 and 1 to do quantile regression
# Set trig to regress on trigonometric functions
# Set step to regress on Hadamard type functions
#

############################### Input Values #################################

type = 'trig'   # trig: sum of trig functions. 
                # step: sum of step functions

# arguments for data generation (both for trig and step)
frequencies = [10, 20, 25]
sin_amplitudes = [2, 1, 4]
cos_amplitudes = [3, -1, 0]
total_samples = 1024                  # data point to generate
ahead = 100                           # points to predict ahead
pow_thres = 0.05                      # threshold to filter harmonics out
quantile = 0.75                        # nil for mse, between 0 and 1 for quantile
level = 0.95                          # confidence levels
print = true                          # true or false. true for lots of printing and graphs

###############################################################################

########################### Load RSRuby and gnuplot###################################
ENV['R_HOME'] = "/usr/lib/R"
ENV['RB_GNUPLOT'] = "/usr/bin/gnuplot"

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
load_gem('gnuplot') if print == true
load_gem('rsruby')
r = RSRuby.instance 


# check args for data generation
args = []
if frequencies.length != sin_amplitudes.length or frequencies.length != cos_amplitudes.length
  puts "type: #{type}. Wrong parameters"
  exit
end
args << sin_amplitudes 
args << cos_amplitudes
args << frequencies 

############ print headings ###############
if print == true
  puts "Input Data"
  puts "Total Samples: #{total_samples}"
end
puts "Frequencies: #{frequencies.inspect}"
puts "Amplitudes\n\tsin: #{sin_amplitudes.inspect}\n\tcos: #{cos_amplitudes.inspect}"

######################### FUNCTIONS ########################

def print_coefs(smry, size, freqs = nil)  # prints coef stats and returns their values
  vcoefs = smry["coefficients"]
  coefs, upr_coefs, lwr_coefs = [], [], []
  puts "\tCoefficients\n\t#\t\t\t\tEstimate\t\t\tError\t\t\t\ttValue\t\t\t\tProb >|t|"
  (0...size).each do |d|
    str = d.to_s + "\t\t"
    if freqs != nil and d > 0
      q, r = (d - 1).divmod(2)
      f = freqs[q]
      s = d % 2 == 0 ? "cos" : "sin"
      str = "#{s}( 2 * PI * #{f})"
    end
    puts "\t#{str}\t\t#{vcoefs[d][0]}\t\t#{vcoefs[d][1]}\t\t#{vcoefs[d][2]}\t\t#{vcoefs[d][3]}"
    coefs << vcoefs[d][0]
    upr_coefs << vcoefs[d][0] + vcoefs[d][1]
    lwr_coefs << vcoefs[d][0] - vcoefs[d][1]
  end
  return coefs, upr_coefs, lwr_coefs
end

# gnu plot function
def plot(vtime, vdata, start_date, end_date, app_name, min, max, xlabel)  # vdata: hash of data sets to plot
  Gnuplot.open do |gp|
      Gnuplot::Plot.new(gp) do |plot|
      plot.title "#{app_name} from #{start_date} to #{end_date}"
      plot.xrange "[#{vtime[0]}:#{vtime[-1]}]"
      plot.yrange "[#{min}:#{max}]"
      plot.xlabel xlabel
      plot.ylabel "#{app_name}"
      plot.grid "x y"
      plot.data = []
      vdata.each_key do |name|
        plot.data << Gnuplot::DataSet.new( [vtime, vdata[name]] ) { |ds|
          ds.with = "lines"    # "lines", "points", "linespoints", "dots"
          ds.title = "#{name}"
          ds.linewidth = 1
          ds.axes = "x1y1"
        }
      end
    end
  end
end

# generates a data point from input data
def get_value(args, time, type)
  sin_amplitude, cos_amplitude, frequency = args
  value = 0
  frequency.each_index do |i| 
    y = sin_amplitude[i] * Math.sin(2.0 * Math::PI * frequency[i] * time) + cos_amplitude[i] * Math.cos(2.0 * Math::PI * frequency[i] * time)
    case type
    when 'trig'
      value += y + 2 * rand()
    when 'step'
      value += y > 0 ? y.ceil : 0.0
    end
  end
  return value  
end

# generates input data
def set_data(args, type, total_samples, time_between_samples) 
  x, y, xvalue, min, max = [], [], 0, 1000000.0, 0.0
  yhash = Hash.new
  yhash["ts"] = []
  (0...total_samples).each do |i|
    x << xvalue
    yvalue = get_value(args, x[-1], type)
    y << yvalue
    yhash["ts"] << yvalue
    xvalue += time_between_samples
    min = yvalue if yvalue < min
    max = yvalue if yvalue > max
  end
  return x, y, yhash, min, max
end

# fits the data and does the prediction 
def do_pred(level, rObj, total_samples, time_between_samples, quantile, print, pow_thres, ahead, args, type)
  x, y, yhash, min, max = set_data(args, type, total_samples, time_between_samples)
  
  # assign input data to R objects to 
  rObj.assign('xr', x)
  rObj.assign('yr', y)
  dt = rObj.eval_R("dt <- data.frame(xr,yr)")
  
  # find the input data spectrum
  puts "\nSpectral Analysis"
  s = rObj.eval_R("spec <- spectrum(yr, plot = FALSE)")
  freq = rObj.eval_R("spec$freq")
  pow = rObj.eval_R("spec$spec")
  bw = rObj.eval_R("spec$bandwidth")
  spec_power, freqs = Hash.new, []
  freq.each_index {|k| spec_power[freq[k]] = pow[k]}
  sspec_power = spec_power.sort{|a, b| b[1] <=> a[1]}  # sorted by decreasing power
  sspec_power.each_index do |k| 
    if sspec_power[k][1] > pow_thres * sspec_power[0][1]
      freqs << sspec_power[k][0] / time_between_samples             # these are the actual frequencies
      puts "\tfreq: #{freqs[-1]} power: #{sspec_power[k][1]}"   # keep only freqs with power larger than pow_thres * max_power
    end
  end

  ################# model ###################
  s = quantile == nil ? "\nMSE Regression" : "\n#{quantile} Quantile Regression"
  puts s        # type of regression 
  rObj.assign('rsp', "response")
  rObj.assign('prd', "prediction")
  rObj.assign('se',"boot")  # boot, ker, nid, iid, rank. Only boot or iid work. needed for quantile
  rObj.assign('t', "none")       # needed for quantile
  rObj.assign('c',"confidence")  # needed for quantile
  
  # lm formula
  model = ""
  vw = []
  freqs.each {|f| w = 2.0 * Math::PI * f; model << "sin(#{w} * xr) + cos(#{w} * xr) + "; vw << w}
  model.slice!(-2)  # remove the last spaces and +
  
  ndt = rObj.eval_R("ndt <- data.frame(xr)")
  # regression
  if quantile != nil      # quantile regression
    lib = "library(\'quantreg\')"
    rObj.eval_R("suppressPackageStartupMessages(#{lib})")              
    rObj.assign('quant', quantile)
    fit = rObj.eval_R("mdl <- rq(formula = yr ~ #{model}, tau = quant, data = dt, model = TRUE)")
    smry = rObj.eval_R("smry <- summary(mdl, se = se, level = #{level})")
    pred1 = rObj.eval_R("p1 <- predict(mdl, newdata = ndt, type = t, level = #{level}, interval = c, se = se)")  # works only for se = boot, iid, 
  else                      # mse regression
    fit = rObj.eval_R("mdl <- lm(formula = yr ~ #{model}, data = dt)") 
    smry = rObj.eval_R("smry <- summary(mdl, se = se, level = #{level})")
    pred1 = rObj.eval_R("p1 <- predict(mdl, se = TRUE, newdata = ndt, type = rsp, level = #{level}, interval = prd)")
    coef = rObj.eval_R("c <- coef(mdl)")
  end

  coefs, upr_coefs, lwr_coefs = print_coefs(smry, 1 + 2 * freqs.length, freqs)  # print and extract coefs
  puts "\nF-statistic: #{smry["fstatistic"]["value"]}" if quantile == nil
  puts "sigma: #{smry["sigma"]}" if quantile == nil
  puts "R2: #{(100.0 * smry["r.squared"])}%" if quantile == nil
  
  # prepare data for ploting
  yhash['fit'], yhash['pred'] = [], []
  (1..total_samples).each do |i|
    if quantile == nil          # mse
      yR, yRlwr, yRupr = pred1["fit"][i - 1]
      se = pred1["se.fit"][i-1]
    else                        # quantile
      yR, yRlwr, yRupr = pred1[i - 1]
      se = (yRupr - yRlwr) / 2.0
    end
    yhash['fit'] << yR
    yhash['pred'] << yR
  end unless print == false
  
  ################# Prediction ##################
  xvalue = x[-1]
  xp = []
  (0...ahead).each do |i|
    xvalue += time_between_samples
    xp << xvalue
  end
  rObj.assign('xr', xp)
  ndt = rObj.eval_R("ndt <- data.frame(xr)") # define the prediction data range
  if quantile == nil    # mse
    pred2 = rObj.eval_R("p2 <- predict(mdl, newdata = ndt, se = TRUE, type = rsp, level = #{level}, interval = prd)")
  else                  # quantile
    pred2 = rObj.eval_R("p2 <- predict(mdl, newdata = ndt, type = t, level = #{level}, interval = c, se = se)")
  end

  len = y.length
  (len...len + ahead).each do |idx|
    xval = x[ idx - len]
    if quantile == nil    # mse
      yR, yRlwr, yRupr = pred2['fit'][idx - len]
      se = pred2["se.fit"][(1 + idx - len).to_s]
    else                  # quantile
      yR, yRlwr, yRupr = pred2[idx - len]
      se = (yRupr - yRlwr) / 2.0
    end
    yhash['pred'] << yR
    yhash['fit'] << 0.0   # padding for plot
    yhash['ts'] << 0.0    # padding for plot
  end
  plot(x + xp, yhash, "", "", "", min * 1.1, max * 1.1 , "time") if print == true     # plot the data
end

############################################## MAIN ####################################################
do_pred(level, r, total_samples, 1.0 / total_samples, quantile, print, pow_thres, ahead, args, type)   # assume total duration is 1 so that time between samples is 1/samples
puts "Done"
