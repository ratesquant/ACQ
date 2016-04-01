using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

namespace ACQ.Excel.Objects
{
    public class ExcelInterpolator2D
    {
        private static readonly string m_tag = "#Interpolator2D";
        private static readonly string m_defaultInterpolator = "Bilinear";

        [ExcelFunction(Description = "Create Interpolator object", Category = AddInInfo.Category)]
        public static object acq_interpolator2d_create(
            [ExcelArgument(Description = "Array of horizontal nodes")] double[] x1,
            [ExcelArgument(Description = "Array of vertical nodes ")]  double[] x2,
            [ExcelArgument(Description = "Table of function values ")]  double[,] y,
            [ExcelArgument(Description = "bilinear")] object method)
          {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag, new object[] { x1, x2, y, method, "acq_interpolator2d_create" },
                    (objectType, paramaters) =>
                    {
                        ACQ.Math.Interpolation.InterpolationInterface2D interpolator = construct_interpolator(x1, x2, y, method);

                        if (interpolator == null)
                            return ExcelError.ExcelErrorNull;
                        else
                            return interpolator;
                    });
            }
        }


        [ExcelFunction(Description = "Evaluate 2D interpolation at specified point", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_interpolator2d_eval(
            [ExcelArgument(Description = "Interpolator object")] string handle,
            [ExcelArgument(Description = "Interpolation point")] double x1,
            [ExcelArgument(Description = "Interpolation point")] double x2)
        {
            ACQ.Math.Interpolation.InterpolationInterface2D interpolator;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Interpolation.InterpolationInterface2D>(handle, out interpolator))
            {
                if (interpolator != null)
                {
                    return ExcelHelper.CheckNan(interpolator.Eval(x1, x2));
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Compute interpolated value (2d)", Category = AddInInfo.Category)]
        public static object acq_interpolation2d(double px1, double px2, double[] x1, double[] x2, double[,] y, object method)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                ACQ.Math.Interpolation.InterpolationInterface2D interpolator = construct_interpolator(x1, x2, y, method);

                if (interpolator != null)
                    return interpolator.Eval(px1, px2);
                else
                    return ExcelError.ExcelErrorNum;
            }
        }

        private static ACQ.Math.Interpolation.InterpolationInterface2D construct_interpolator(double[] x1, double[] x2, double[,] y, object method)
        {
            ACQ.Math.Interpolation.InterpolationInterface2D interpolator = null;
            try
            {
                string interpolation_method = ExcelHelper.Check(method, m_defaultInterpolator);
            
                interpolator = ACQ.Math.Interpolation.InterpolationFactory2D.GetInterpolator(interpolation_method, x1, x2, y);
            }
            catch (Exception ex)
            {
                LogDisplay.WriteLine("Error: " + ex.ToString());
            }
            return interpolator;
        }
    }
}
