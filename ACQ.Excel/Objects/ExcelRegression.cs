using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

namespace ACQ.Excel.Objects
{

    public class ExcelRegression
    {
        private static readonly object m_sync = new object();
        private static readonly string m_tag = "#acqRegression";

        [ExcelFunction(Description = "Create lowess smoother, locally-weighted linear regression", Category = AddInInfo.Category)]
        public static object acq_regression_lowess_create(
            [ExcelArgument(Description = "Array of nodes")] double[] x,
            [ExcelArgument(Description = "Array of values")]  double[] y,
            [ExcelArgument(Description = "smoother span [0, 1] (optional, default 2/3)")] object span,
            [ExcelArgument(Description = "number of robust iterations (optional, default 3) ")] object nsteps,
            [ExcelArgument(Description = "regression step (optional)")] object delta)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {

                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag, new object[] { x, y, span, nsteps, delta, "acq_regression_lowess_create" },
                   (objectType, parameters) =>
                   {
                       ACQ.Math.Regression.Lowess lowess = construct_lowess(x, y, span, nsteps, delta);

                       if (lowess == null)
                           return ExcelError.ExcelErrorNull;
                       else
                           return lowess;
                   });
            }
        }

        [ExcelFunction(Description = "Evaluate lowess smoother at specified point", Category = AddInInfo.Category)]
        public static object acq_regression_lowess_eval(
            [ExcelArgument(Description = "Lowess object")] string handle, double x)
        {
            Tuple<bool, double> results = ACQ.Excel.Handles.GlobalCache.TryReadObject<ACQ.Math.Regression.Lowess, double, double>(handle, 
                (lowess, point) =>
            {
                return lowess.Eval(point);
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

        [ExcelFunction(Description = "Lowess smoother, locally-weighted linear regression", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_regression_lowess(
            [ExcelArgument(Description = "Array of nodes")] double[] x,
            [ExcelArgument(Description = "Array of values")]  double[] y, double xp,
            [ExcelArgument(Description = "smoother span [0, 1] (optional, default 2/3)")] object span,
            [ExcelArgument(Description = "number of robust iterations (optional, default 3) ")] object nsteps,
            [ExcelArgument(Description = "regression step (optional)")] object delta)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {   
                ACQ.Math.Regression.Lowess lowess = construct_lowess(x, y, span, nsteps, delta);

                return ExcelHelper.CheckNan(lowess.Eval(xp));
            }
        }

        private static ACQ.Math.Regression.Lowess construct_lowess(double[] x, double[] y, object span, object nsteps, object delta)
        {
            ACQ.Math.Regression.Lowess lowess = null;
            try
            {
                double span_input = ExcelHelper.CheckValue<double>(span, 2.0 / 3.0);
                double delta_input = ExcelHelper.CheckValue<double>(delta, 0.0); //set to zero so that lowess will pick default based on data range
                int nsteps_input = (int)ExcelHelper.CheckValue<double>(nsteps, 3);

                lowess = new Math.Regression.Lowess(x, y, span_input, nsteps_input, delta_input);
            }
            catch (Exception ex)
            {
                lock (m_sync)
                {
                    LogDisplay.WriteLine("Error: " + ex.ToString());
                }
            }
            return lowess;
        }


    }
}
