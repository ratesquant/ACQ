# ACQ - Excel Add-in for interpolation (and other things)
uses Excel-DNA, https://exceldna.codeplex.com/ (https://github.com/Excel-DNA)

# Abstract 
ACQ is lightweight .NET add-in that adds interpolation routines to Excel. It is open source and free to use for any purpose. Add-in is based on Excel-DNA, and uses .NET 4.0. 

The distribution contains: 

	1. 32-bit and 64-bit versions of add-in: ACQ32.xll, ACQ64.xll
	2. Excel spreadsheet with examples of add-in functionality: Demo.xlsx
	3. Description of numerical methods used in add-in: Primer.pdf 
	

# Installation
In order to install ACQ add-in follow the installation procedures for Custom Add-ins:
https://support.office.com/en-us/article/Add-or-remove-add-ins-0af570c4-5cf3-4fa9-9b88-403625a0b460

	1. Click the File tab, click Options, and then click the Add-Ins category 
	2. In the Manage box, click Excel Add-ins, and then click Go. The Add-Ins dialog box appears
	3. Browse to the location of the ACQ32.xll/ACQ64.xll files, and pick the xll file based on Excel bitness.
	4. Make sure your Excel security settings allow you to run Add-ins 
	5. All ACQ functions have prefix "acq_" and located in ACQ category.
	6. Let me know if you have any issues: rates.quant@gmail.com
    
# Interpolation
The following functions for interpolation are currently implemented:

	1. acq_interpolator_create - creates interpolator object (returns a handle)
	2. acq_interpolator_eval - evaluates interpolation at specified point
	3. acq_interpolation - evaluates interpolation at specified point(in-situ without constructing interpolator object)
	
Interpolaton method is specified using method argument. The argument is optional, linear interpolation is used by default. Currently implemented methods are

	0. Nearest - nearest point interpolation (also Backward and Forward)
	1. Linear - linear spline interpolation
	2. Quadratic - quadratic  spline interpolation	
	3. Cubic - natural cubic spline
	4. Hermite and HermiteQS - local cubic spline (aka Catmull-Rom spline)
	5. Steffen - monotonic cubic spline
	6. Akima - Akima spline (cubic spline with special condition for derivatives)
	7. AkimaPeriodic - Akima spline with periodic boundary conditions
	8. Multiquadrics - Interpolation based on multiquadrics radial basis function 

Bounds is optional argument that controls interpolation outside of interpolation range. When interpolating outside of range num Excel error will be returned if bounds is false, while closest point is returned when bounds is true (default).  

# 2D Interpolaton
	1. BiLinear - interpolation on rectangular grid (linear in each dimension)
	2. BiCubic - cubic interpolation on rectangular grid (natural cubic spline in each dimension)
	3. BiSteffen - Steffen interpolation on rectangular grid (Steffen in each dimension)
	4. BiAkima - Akima interpolation on rectangular grid (Akima in each dimension)
	5. BiHermite - Hermite interpolation on rectangular grid (Hermite in each dimension)
	
	
	1. acq_interpolator2d_create - creates 2D interpolator object
	2. acq_interpolator2d_eval - evaluates interpolation at specified point
	3. acq_interpolation2d - evaluates interpolation at specified point(in-situ without constructing interpolator object)

	
# Scattered Data Interpolaton (ND)
Scattered data interpolation is based on radial basis functions and currently limited to 512 interpolation nodes. The following radial basis functions are implemented (first three are the most common). Optional scale factors can be provided for each dimension

	1. Linear
	2. Cubic
	3. Multiquadrics
	4. Gaussian
	5. Thinplate,
	6. InverseQuadratic,
	7. InverseMultiquadric

the list of excel functions for scattered data interpolation:

	1. acq_interpolator_scattered_create - creates RBF interpolator
	2. acq_interpolator_scattered_eval - evaluates RBF interpolator at specified point
	3. acq_interpolator_scattered_eval_x5 - evaluates RBF interpolator at specified point, coordinates are specified individually (up to 5D)

# Linear Regression
Weighted Linear regression (with and without intercept). Compared to LINEST Excel function it allows to specify weights, does not use array formula and allows to directly compute regression estimate. 
	1. acq_regression_linear_create - creates linear regression (returns regression object)
	2. acq_regression_eval - computes regression estimate (using regression object)

# Lowess
Lowess is scatter plot smoothing based on locally-weighted linear regression (based on R code)

	1. acq_regression_lowess_create - creates lowess smoother
	2. acq_regression_lowess_eval - evaluates lowess smoother at specified point
	3. acq_regression_lowess - evaluates lowess smoother at specified point(in-situ without constructing lowess object)
	
		
# Mersenne Twister
Random number generator based on MT19937 by Takuji Nishimura and Makoto Matsumoto. 
Excel interface allows to initialize generator with single or array seed

	1. acq_random_vector(seed, size) - generate random sample using MT19937, [0, 1) 
	2. acq_random_vector_ex(seed, size) - generate random sample using array seed
	3. acq_vector_element(vector, index) - get a random number from vector (random numbers are returned as vector)
	4. acq_vector_size(vector) - get size of the vector 


# Black pricing formulas for European options
Input arguments: forward: forward price of the underlying, strike : option strike, time: time until expiration in years, rate: risk free rate, dividend: (dividend yield continuously compounded), sigma: implied volatility (annual), isCall: TRUE for call options, FALSE for puts

	1. acq_options_black_price - Compute option price  
	2. acq_options_black_vol - Compute implied volatility (to match specified price)
	3. acq_options_black_greeks - Compute option greeks: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho, Theta (analytical)
	
	
# Bachelier pricing formulas for European options
Input arguments: same as for Black options pricing formulas, but implied volatility (sigma) is normal

	1. acq_options_black_price - Compute option price  
	2. acq_options_bachelier_vol - Compute implied volatility (to match specified price) 
	3. acq_options_bachelier_greeks - Compute option greeks: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho, Theta (analytical)


# Black-Scholes pricing formulas for European options
Input arguments: same as for Black options pricing formulas, but uses spot: spot price of the underlying, instead of forward.

	1. acq_options_blackscholes_price - Compute option price  
	2. acq_options_blackscholes_vol - Compute implied volatility (to match specified price) 
	3. acq_options_blackscholes_greeks - Compute option greeks: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho, Theta, Charm, Epsilon (analytical)	


# Bjerksund and Stensland (2002) approximation for American options
	1. acq_options_bjerksund_price - Compute option price  
	2. acq_options_bjerksund_greeks - Compute option greeks: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho, Theta (numerical)	


# Binomial option pricing method for American options
Input arguments: same as for Black-Scholes with additional optional argument -  time_steps: number of times steps 

	1. acq_options_binomial_american_price - Compute option price   	
	2. acq_options_binomial_american_greeks - Compute option greeks: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho, Theta (numerical)	


# Trinomial option pricing method for American options
Input arguments: same as for Binomial option pricing 

	1. acq_options_trinomial_american_price - Compute option price  	
	2. acq_options_trinomial_american_greeks -  Compute option greeks: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho, Theta (numerical)	

# Misc Utils
Various convenience functions to compute min/max/mean ignoring non-numerical values (mean can be weighted) 
	1. acq_mean - Compute mean(x) (with optional weighted, can ignore non-numerical values)  	
	2. acq_max,acq_min, acq_absmax and acq_absmin - Compute max(x), min(x), max(abs(x)) and  min(abs(x)) (ignores non-numerical values)  
    3. acq_tostring - converts a value to string representation
    4. acq_isinteger - checks if value is integer
    5. acq_isprime - checks if number is prime
    6. acq_join	 - concatenate the elements of the range to string (with specified separator)	
	

# Shotcuts
	1. Ctrl+Shift+H - shows ACQ log window
	2. Ctrl+Shift+A - shows ACQ introspection window (for ACQ handles)
	
