using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;


namespace ACQ.Excel
{
    public static class StatUtils
    {
        /// <summary>
        /// TODO: ignore empty elements, right now they are zeros
        /// </summary>
        /// <param name="x"></param>
        /// <param name="w"></param>
        /// <param name="ignore_na"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute mean (optionaly weighted)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_mean(
            [ExcelArgument(Description = "Array")] object x,
            [ExcelArgument(Description = "Array of Weights [optional]")] object w,
            [ExcelArgument(Description = "Ignore NA flag [optional] (default: false)")] object ignore_na)
        {
            double[] x_input = ExcelHelper.CheckArray<double>(x);
            double[] w_input = ExcelHelper.CheckArray<double>(w);

            bool na_rm = ExcelHelper.CheckValue<bool>(ignore_na, false);            

            double result = ACQ.Math.Stats.Utils.WeightedMean(x_input, w_input, na_rm);

            return ExcelHelper.CheckNan(result);
        }

        /// <summary>
        /// Compute Max
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute max(x) (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_max(
            [ExcelArgument(Description = "Array")] object x)
        {               
            double[] array = ExcelHelper.CheckArray<double>(x);

            double max = ACQ.Math.Stats.Utils.Max(array);

            return ExcelHelper.CheckNan(max);
        }

        /// <summary>
        /// Compute Max
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute max(abs(x)) (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_absmax(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x);

            double absmax = ACQ.Math.Stats.Utils.AbsMax(array);

            return ExcelHelper.CheckNan(absmax);
        }

        /// <summary>
        /// Compute Min
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute min(x) (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_min(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x);

            double min = ACQ.Math.Stats.Utils.Min(array);

            return ExcelHelper.CheckNan(min);
        }

        /// <summary>
        /// Compute Min
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute min(abs(x)) (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_absmin(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x);

            double min = ACQ.Math.Stats.Utils.AbsMin(array);

            return ExcelHelper.CheckNan(min);
        }

        /// <summary>
        /// Compute Min
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute sum(x) (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_sum(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x);

            double sum = ACQ.Math.Stats.Utils.Sum(array);

            return ExcelHelper.CheckNan(sum);
        }

        /// <summary>
        /// Compute Min
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Total number of numeric elements in x (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_count(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x); //it is a bit wasteful to create a new array, just to count number of numeric values, but it keeps code cleaner 

            return array.Length;
        }

        /// <summary>
        /// Compute standard deviation
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute standard deviation of x (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_stdev(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x);

            double stdev = ACQ.Math.Stats.Utils.Std(array);

            return ExcelHelper.CheckNan(stdev);
        }


        /// <summary>
        /// Compute standard deviation
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [ExcelFunction(Description = "Compute sum of squares of x (ignores non-numeric values)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_sumofsquares(
            [ExcelArgument(Description = "Array")] object x)
        {
            double[] array = ExcelHelper.CheckArray<double>(x);

            double sum_sqr = ACQ.Math.Stats.Utils.SumOfSquares(array);

            return ExcelHelper.CheckNan(sum_sqr);
        }

        [ExcelFunction(Description = "Binomial proportion confidence interval (Wilson - method)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_binomtest(
            [ExcelArgument(Description = "Vector of number of successes in the binomial experiment")] int x,
            [ExcelArgument(Description = "Vector of number of independent trials in the binomial experiment.")] int n,
            [ExcelArgument(Description = "Confidence level (e.g. 0.95)")] double conf_level,
            [ExcelArgument(Description = "Lower (TRUE) or Upper bound (FALSE)")] bool isLower)
        {
            
            var result = ACQ.Math.Stats.BinomTest.ConfidenceInterval(x, n, conf_level);

            return ExcelHelper.CheckNan(isLower ? result.Item1 : result.Item2);
        }
    }
}
