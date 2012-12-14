R Examples with Ruby and Python

See http://josepferrandiz.blogspot.com/ for more details.
The r-*.{py,rb} files contain examples on how perform regressions 
using rpy2 or rsruby.
The files can be executed from the command line or using, 
> [python, ruby] file.[py,rb]
The parameters are included at the top of the file.

Trigonometric Regressions and Predictions
The r-trig-regs-pred.[py,rb] files perform least squares and/or quantile 
regressions on a trigonometric or step function derived from its frequency 
components.
Configuration
  type: defines whether the data to be fit and predicted will be a trigonometric 
        function or a step function.   
        trig: generates a trig function. 
        step: generates a step function.
        White noise is added to the data generated.
  frequencies: an array of positive values that contains the frequencies used to 
        generate the input data
  sin_amplitudes, cos_amplitudes: array that contain the amplitudes of the sin()
        and cos() components. Values can be positive or negative. Each array must 
        be the same length as the freqeuncy array.       
  total_samples: the number of data points to generate.
  ahead: number of points to predict ahead
  pow_thres: threshold to filter harmonics out. Should be 0.1 or less.
  quantile: nil for mse regression in ruby. between 0 and 1 for quantile
  level: confidence levels for the regression
  print: true or false. Set to true for lots of printing and graphs
  
Notes
1) The data points are of the form:
Sum_{k=0}^N a_k sin(2 * pi * f_k * t)+ b_k cos(2 * pi * f_k * t) + noise
The step functions are obtained by taking the ceiling value of each point.
2) There seems to be a bug in R that makes the pvalues of the parameters be 
misreported.

Polynomial Regressions and Predictions
The r-poly-regs-pred.[py,rb] files perform least squares and/or quantile 
regressions on a polynomial of arbitrart degree. 
Configuration
  imax: number of data points to generate
  iahead: number points to predict past imax
  degree: polynomial degree
  level: confidence levels for the regression                              
  quantile: nil for mse regression in ruby. between 0 and 1 for quantile
  print: true or false. Set to true for lots of printing. No graphs.

Notes:
1) There are no graphs available for the polynomial functions.
2) Due to a bug in R we must set the number of prediction points to be larger 
than the polynomial degree. The code adjusts for that and issues a warning.
