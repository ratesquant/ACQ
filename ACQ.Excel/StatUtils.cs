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
            double max = Double.NegativeInfinity;

            object[,] array = x as object[,];

            if (array != null)
            {
                int n = array.GetLength(0);
                int m = array.GetLength(1);
                
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    { 
                        object item = array[i, j];
                        if (item is Double)
                        {
                            double temp = (Double)item;
                            if (temp > max)
                            {
                                max = temp;
                            }
                        }                 
                    }
                }
            }
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
            double max = Double.NegativeInfinity;

            object[,] array = x as object[,];

            if (array != null)
            {
                int n = array.GetLength(0);
                int m = array.GetLength(1);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        object item = array[i, j];
                        if (item is Double)
                        {
                            double temp = System.Math.Abs( (Double)item);
                            if (temp > max)
                            {
                                max = temp;
                            }
                        }
                    }
                }
            }
            return ExcelHelper.CheckNan(max);
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
            double min = Double.PositiveInfinity;

            object[,] array = x as object[,];

            if (array != null)
            {
                int n = array.GetLength(0);
                int m = array.GetLength(1);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        object item = array[i, j];
                        if (item is Double)
                        {
                            double temp = (Double)item;
                            if (temp < min)
                            {
                                min = temp;
                            }
                        }
                    }
                }
            }
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
            double min = Double.PositiveInfinity;

            object[,] array = x as object[,];

            if (array != null)
            {
                int n = array.GetLength(0);
                int m = array.GetLength(1);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        object item = array[i, j];
                        if (item is Double)
                        {
                            double temp = System.Math.Abs((Double)item);
                            if (temp < min)
                            {
                                min = temp;
                            }
                        }
                    }
                }
            }
            return ExcelHelper.CheckNan(min);
        }
    }
}
