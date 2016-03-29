# ACQ, Excel Add-in for interpolation 
uses Excel-DNA, https://exceldna.codeplex.com/

# Abstract 
Excel is very useful tool for working with small datasets. However, some essential functionality such as interpolation is missing. Add-ins and VBA are used to extend standard Excel functionality to turn it into truly powerful tool. All code developed as part of work remains the property of the company, and have to be implemented when changing jobs. Therefore, I decided to develop ACQ: free (for any purpose), open source implementation of frequently used excel functionality: interpolation, etc. 

The Add-in is based on Excel-DNA, and uses .NET 4.0.

The distribution contains: 

	1. 32 and 64-bit version of Add-in: ACQ32.xll, ACQ64.xll
	2. Excel spreadsheet that demonstrates Add-in functionality : Demo.xlsx
	3. Short description of all numerical methods using in Add-in: Primer.pptx 
	

# Installation 
	1. In Excel go File->Options-> Add-ins. Select "Excel Add-in" in Manage (bottom of the option screen), Click "Go"
	2. Browse to the location of the ACQ32.xll/ACQ64.xll file, and pick the xll file based on bitness of your Excel version. If you select wrong file Excel will be unable to load it.
	3. Make sure your Excel security settings allow you to run Add-ins 
	4. All ACQ functions have prefix "acq_" and located in ACQ category.
    
# Functions