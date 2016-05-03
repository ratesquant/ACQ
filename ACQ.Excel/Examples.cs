using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    public static class Examples
    {
        /*
        [ExcelFunction(Description = "acq_sum", Category = AddInInfo.Category)]
        public static double acq_sum(double a, double b)
        {
            return a + b;
        }         

        [ExcelFunction(Description = "Compute interpolated value", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_linear_interpolation(double xi, double[] x, double[] y, object bounds)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                ACQ.Math.Interpolation.InterpolationInterface interpolator = null;
                try
                {
                    interpolator = new ACQ.Math.Interpolation.LinearInterpolation(x, y); //create linear interpolator 
                    interpolator.Bounds = ExcelHelper.CheckValue(bounds, true);
                }
                catch (Exception ex)
                {
                    //LogDisplay.WriteLine("Error: " + ex.ToString());                    
                }

                if (interpolator != null)
                    return ExcelHelper.CheckNan(interpolator.Eval(xi)); //interpolate
                else
                    return ExcelError.ExcelErrorNA;
            }
        }
         */
    }
}
