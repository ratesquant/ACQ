using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

namespace ACQ.Excel.Objects
{
    public class ExcelInterpolator
    {
        private static readonly object m_sync = new object();
        private static readonly string m_tag = "#acqInterpolator";
        private static readonly string m_defaultInterpolator = "Linear";

        [ExcelFunction(Description = "Create Interpolator object", Category = AddInInfo.Category)]
        public static object acq_interpolator_create(
            [ExcelArgument(Description = "Array of nodes")] double[] x,
            [ExcelArgument(Description = "Array of values")]  double[] y,
            [ExcelArgument(Description = "linear, quadratic, cubic, hermite, akima, steffen etc")] object method,
            [ExcelArgument(Description = "Out of range value: false (num error), true (closest)")] object bounds)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag, new object[] { x, y, method, bounds, "acq_interpolator_create" },
                    (objectType, parameters) =>
                    {
                        ACQ.Math.Interpolation.InterpolationInterface interpolator = construct_interpolator(x, y, method, bounds);

                        if (interpolator == null)
                            return ExcelError.ExcelErrorNull;
                        else
                            return interpolator;
                    });
            }
        }
       

        [ExcelFunction(Description = "Evaluate interpolation at specified point", Category = AddInInfo.Category, IsThreadSafe = false)]
        public static object acq_interpolator_eval(
            [ExcelArgument(Description = "Interpolator object")] string handle,
            [ExcelArgument(Description = "Interpolation point")] double x)
        {
            ACQ.Math.Interpolation.InterpolationInterface interpolator;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Interpolation.InterpolationInterface>(handle, out interpolator))
            {
                if (interpolator != null)
                {
                    return ExcelHelper.CheckNan(interpolator.Eval(x));
                }
            }
            return ExcelError.ExcelErrorRef;
        }


        [ExcelFunction(Description = "Evaluate interpolation at specified point (thread safe version)", Category = AddInInfo.Category, IsThreadSafe = false)]
        public static object acq_interpolator_eval_tsafe(
            [ExcelArgument(Description = "Interpolator object")] string handle,
            [ExcelArgument(Description = "Interpolation point")] double x)
        {

            Tuple<bool, double> results = ACQ.Excel.Handles.GlobalCache.TryReadObject<ACQ.Math.Interpolation.InterpolationInterface, double, double>(handle, (interpolator, point) =>
            {
                return interpolator.Eval(point);
            }, x);

            if (results.Item1)
            {
                return ExcelHelper.CheckNan(results.Item2);
            }
            else
            {
                return ExcelError.ExcelErrorRef;

            }
        }

        [ExcelFunction(Description = "Compute interpolated value", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_interpolation(double xi, double[] x, double[] y, object method, object bounds)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                ACQ.Math.Interpolation.InterpolationInterface interpolator = construct_interpolator(x, y, method, bounds);

                if (interpolator != null)
                    return ExcelHelper.CheckNan(interpolator.Eval(xi));
                else
                    return ExcelError.ExcelErrorNA;
            }

        }

        private static ACQ.Math.Interpolation.InterpolationInterface construct_interpolator(double[] x, double[] y, object method, object bounds)
        {
            ACQ.Math.Interpolation.InterpolationInterface interpolator = null;
            try
            {
                string interpolation_method = ExcelHelper.Check(method, m_defaultInterpolator);
                bool interpolation_bounds = ExcelHelper.CheckValue(bounds, true);

                interpolator = ACQ.Math.Interpolation.InterpolationFactory.GetInterpolator(interpolation_method, x, y, interpolation_bounds);
            }
            catch (Exception ex)
            {
                lock (m_sync)
                {
                    LogDisplay.WriteLine("Error: " + ex.ToString());
                }
            }
            return interpolator;
        }
    }
}
