# ACQ - Excel Add-in for interpolation
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

	1. acq_interpolator_create(x, y, method, bounds) - creates interpolator object (returns a handle)
	2. acq_interpolator_eval(interpolator, x) - evaluates interpolation at specified point (first argument is a handle created using acq_interpolator_create function)
	3. acq_interpolation(xi, x, y, method, bounds) - creates interpolator and computes interpolation at specified point
	
Interpolaton method is specified using method argument. The argument is optional, linear interpolation is used by default. Currently implemented methods are

	0. Nearest - nearest point interpolation (also Backward and Forward)
	1. Linear - linear spline interpolation
	2. Quadratic - quadratic  spline interpolation	
	3. Cubic - natural cubic spline
	4. Hermite - local cubic spline (aka Catmull-Rom spline)
	5. Steffen - monotonic cubic spline
	6. Akima - Akima spline (cubic spline with special condition for derivatives)
	7. AkimaPeriodic - Akima spline with periodic boundary conditions
	8. Multiquadrics - Interpolation based on multiquadrics radial basis function 

Bounds is optional argument that controls interpolation outside of interpolation range. When interpolating outside of range num Excel error will be returned if bounds is false, while closest point is returned when bounds is true (default).  

# 2D Interpolaton
	1. BiLinear - interpolation on rectangular grid (linear in each dimension)
	2. BiCubic - cubic interpolation on rectangular grid (hermite cubic spline in each dimension)
	3. BiSteffen - Steffen interpolation on rectangular grid (Steffen in each dimension)
	4. BiAkima - Akima interpolation on rectangular grid (Akima in each dimension)
	4. BiHermite - Hermite interpolation on rectangular grid (Hermite in each dimension)
	
	
# Scattered Data Interpolaton (ND)
Scattered data interpolation is based on radial basis functions and currently limited to 512 interpolation nodes. The following radial basis functions are implemented (first three are the most common). Optional scale factors can be provided for each dimension
	1. Linear
	2. Cubic
	3. Multiquadrics
	4. Gaussian
	5. Thinplate,
	6. InverseQuadratic,
	8. InverseMultiquadric

the list of excel functions for scattered data interpolation:
	1. acq_interpolator_scattered_create
	2. acq_interpolator_scattered_eval
	3. acq_interpolator_scattered_eval_x5

# Mersenne Twister
Random number generator based on MT19937 by Takuji Nishimura and Makoto Matsumoto. 
Excel interface allows to initialize generator with single or array seed

	1. acq_random_vector(seed, size) - generate random sample using MT19937, [0, 1) 
	2. acq_random_vector_ex(seed, size) - generate random sample using array seed
	3. acq_vector_element(vector, index) - get a random number from vector (random numbers are returned as vector)
	4. acq_vector_size(vector) - get size of the vector 
	
	
# Shotcuts
	1. Ctrl+Shift+H - shows ACQ log window
	2. Ctrl+Shift+A - shows ACQ introspection window (for ACQ handles)
	
