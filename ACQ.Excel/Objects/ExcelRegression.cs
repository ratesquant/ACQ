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
            [ExcelArgument(Description = "x")] double[] x,
            [ExcelArgument(Description = "y")]  double[] y,
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

        [ExcelFunction(Description = "Estimate Regression value at specified point", Category = AddInInfo.Category)]
        public static object acq_regression_estimate1d(
            [ExcelArgument(Description = "Regression object")] string handle, double x)
        {
            Tuple<bool, double> results = ACQ.Excel.Handles.GlobalCache.TryReadObject<ACQ.Math.Regression.IRegression, double, double>(handle, 
                (regression, point) =>
            {
                return regression.Estimate(point);
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

        [ExcelFunction(Description = "Estimate Regression value at specified point", Category = AddInInfo.Category,IsThreadSafe = true)]
        public static object acq_regression_estimate(
            [ExcelArgument(Description = "Regression object")] string handle, double[] x)
        {
            Tuple<bool, double> results = ACQ.Excel.Handles.GlobalCache.TryReadObject<ACQ.Math.Regression.IRegression, double, double[]>(handle,
                (regression, point) =>
                {
                    return regression.Estimate(point);
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

        [ExcelFunction(Description = "Extract regression parameters", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_regression_param(
            [ExcelArgument(Description = "Regression object")] string handle,
            [ExcelArgument(Description = "Parameter name")] string name)
        {
            Tuple<bool, double> results = ACQ.Excel.Handles.GlobalCache.TryReadObject<ACQ.Math.Regression.IRegressionParam, double, string>(handle,
                (regression, param_name) =>
                {
                    return regression.GetParam(param_name);
                }, name);

            if (results.Item1)
            {
                return ExcelHelper.CheckNan(results.Item2);
            }
            else
            {
                return ExcelError.ExcelErrorRef;
            }
        }
        /* //not needed
        [ExcelFunction(Description = "Regression summary", Category = AddInInfo.Category)]
        public static object acq_regression_summary(
            [ExcelArgument(Description = "Regression object")] string handle)
        {
            //extract regression summary
            Tuple<bool, Dictionary<string, double>> results = ACQ.Excel.Handles.GlobalCache.TryReadObject<ACQ.Math.Regression.IRegressionSummary, Dictionary<string, double>>(handle,
                (regression) =>
                {
                    return regression.Summary;
                });

            if (results.Item1)
            {                
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag_summary, new object[] { handle, "acq_regression_summary" },
                   (objectType, parameters) =>
                   {
                       return results.Item2;                    
                   });
            }
            else
            {
                return ExcelError.ExcelErrorRef;
            }
        }*/

        [ExcelFunction(Description = "Lowess smoother, locally-weighted linear regression", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_regression_lowess(
            [ExcelArgument(Description = "x")] double[] x,
            [ExcelArgument(Description = "y")]  double[] y, double xp,
            [ExcelArgument(Description = "smoother span [0, 1] (optional, default 2/3)")] object span,
            [ExcelArgument(Description = "number of robust iterations (optional, default 3) ")] object nsteps,
            [ExcelArgument(Description = "regression step (optional)")] object delta)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {   
                ACQ.Math.Regression.Lowess lowess = construct_lowess(x, y, span, nsteps, delta);

                return ExcelHelper.CheckNan(lowess.Estimate(xp));
            }
        }


        [ExcelFunction(Description = "Create Least Angle Regression (LARS)", Category = AddInInfo.Category)]
        public static object acq_regression_lars_create(
            [ExcelArgument(Description = "x")] double[,] x,
            [ExcelArgument(Description = "y")]  double[] y)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {

                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag, new object[] { x, y, "acq_regression_lars_create" },
                   (objectType, parameters) =>
                   {
                       ACQ.Math.Regression.Lars lars = new Math.Regression.Lars(x, y);

                       if (lars == null)
                           return ExcelError.ExcelErrorNull;
                       else
                           return lars;
                   });
            }
        }

        [ExcelFunction(Description = "Create Linear Regression", Category = AddInInfo.Category)]
        public static object acq_regression_linear_create(
            [ExcelArgument(Description = "x")] double[,] x,
            [ExcelArgument(Description = "y")]  double[] y,
            [ExcelArgument(Description = "Intercept, optional(default=true)")]  object intercept,
            [ExcelArgument(Description = "Weights, optional(default = none)")]  object weights)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                bool include_intercept = ExcelHelper.CheckValue(intercept, true);
                double[] reg_weights = ExcelHelper.CheckArray<double>(weights);

                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag, new object[] { x, y, reg_weights, include_intercept, "acq_regression_linear_create" },
                   (objectType, parameters) =>
                   {
                       ACQ.Math.Regression.LinearRegression regression = null;

                       try
                       {
                           regression = new Math.Regression.LinearRegression(x, y, reg_weights, include_intercept);
                       }
                       catch (Exception ex)
                       {
                           LogDisplay.WriteLine("Error: " + ex.ToString());      
                       }

                       if (regression == null)
                           return ExcelError.ExcelErrorNull;
                       else
                           return regression;
                   });
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
