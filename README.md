# ACQ - Excel Add-in for interpolation 
uses Excel-DNA, https://exceldna.codeplex.com/

# Abstract 
Excel is very useful tool for working with small datasets. However, some essential functionality such as interpolation is missing. Add-ins and VBA are used to extend standard Excel functionality to turn it into truly powerful tool. All code developed as part of work remains the property of the company, and have to be implemented when changing jobs. Therefore, I decided to develop ACQ: free (for any purpose), open source implementation of frequently used excel functionality: interpolation, etc. 

The Add-in is based on Excel-DNA, and uses .NET 4.0.

The distribution contains: 

	1. 32-bit and 64-bit versions of add-in: ACQ32.xll, ACQ64.xll
	2. Excel spreadsheet with examples of add-in functionality: Demo.xlsx
	3. Description of all numerical methods using in Add-in: Primer.pptx 
	

# Installation
In order to install ACQ follow installation procedures for Custom add-ins:
https://support.office.com/en-us/article/Add-or-remove-add-ins-0af570c4-5cf3-4fa9-9b88-403625a0b460

	1. Click the File tab, click Options, and then click the Add-Ins category 
	2. In the Manage box, click Excel Add-ins, and then click Go. The Add-Ins dialog box appears
	3. Browse to the location of the ACQ32.xll/ACQ64.xll file, and pick the xll file based on bitness of your Excel version. If you select wrong file Excel will be unable to load it.
	4. Make sure your Excel security settings allow you to run Add-ins 
	5. All ACQ functions have prefix "acq_" and located in ACQ category.
	6. Let me know if you have any issues: rates.quant@gmail.com
    
# Interpolation
The following functions for interpolation are currently implemented:

	1. acq_interpolator_create(x, y, method, bounds) - create interpolator object
	2. acq_interpolator_eval(interpolator, x) - evaluate interpolation at specified point
	3. acq_interpolation(xi, x, y, method, bounds) - create interpolater and compute interpolation at specified point
	
Interpolaton method is specified using method argument. The argument is optional, linear interpolation is used by default. Currently implemented methods are
    0. Nearest - nearest point interpolation
	1. Linear - linear spline interpolation
    2. Quadratic - quadratic  spline interpolation	
	3. Cubic - natural cubic spline
	4. Hermite - local cubic spline (aka Catmull-Rom spline)
	5. Steffen - monotonic cubic spline
	6. Akima - cubic spline with Akima conditions

Bounds is optional argument that controls interpolation results outside of interpolation range. When interpolating outside of range num excel error will be returned if bounds is false, and closest point when bounds is true (default).  