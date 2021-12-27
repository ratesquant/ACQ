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
        public static object acq_mean(object x, object w, object ignore_na)
        {
            double[] x_input = ExcelHelper.CheckArray<double>(x);
            double[] w_input = ExcelHelper.CheckArray<double>(w);

            bool na_rm = ExcelHelper.CheckValue<bool>(ignore_na, false);            

            double result = ACQ.Math.Stats.Utils.WeightedMean(x_input, w_input, na_rm);

            return ExcelHelper.CheckNan(result);
        }
       
    }
}
