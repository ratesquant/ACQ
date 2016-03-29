# ACQ - Excel Add-in for interpolation 
uses Excel-DNA, https://exceldna.codeplex.com/

# Abstract 
ACQ is lightweight .NET Add-in that adds interpolation routines to Excel. It is open source and free to use for any purpose. This Add-in is based on Excel-DNA, and uses .NET 4.0. 

The distribution contains: 

	1. 32-bit and 64-bit versions of Add-in: ACQ32.xll, ACQ64.xll
	2. Excel spreadsheet with examples of Add-in functionality: Demo.xlsx
	3. Description of all numerical methods using in Add-in: Primer.pptx 
	

# Installation
In order to install ACQ Add-in follow the installation procedures for Custom Add-ins:
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

	0. Nearest - nearest point interpolation
	1. Linear - linear spline interpolation
	2. Quadratic - quadratic  spline interpolation	
	3. Cubic - natural cubic spline
	4. Hermite - local cubic spline (aka Catmull-Rom spline)
	5. Steffen - monotonic cubic spline
	6. Akima - Akima spline (cubic spline with special condition for derivatives)
	7. AkimaPeriodic - Akima spline with periodic boundary conditions

Bounds is optional argument that controls interpolation outside of interpolation range. When interpolating outside of range num Excel error will be returned if bounds is false, while closest point is returned when bounds is true (default).  