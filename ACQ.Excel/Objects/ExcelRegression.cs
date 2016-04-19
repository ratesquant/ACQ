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

        [ExcelFunction(Description = "Lowess Regression", Category = AddInInfo.Category)]
        public static object acq_regression_lowess(
            [ExcelArgument(Description = "Array of nodes")] double[] x,
            [ExcelArgument(Description = "Array of values")]  double[] y, double xp, object span, object nsteps, object delta)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                double span_input = ExcelHelper.CheckValue<double>(span, 2.0 / 3.0);
                double delta_input = ExcelHelper.CheckValue<double>(delta, 0.01);
                int nsteps_input = ExcelHelper.CheckValue<int>(nsteps, 3);

                ACQ.Math.Regression.Lowess lowess = new Math.Regression.Lowess(x, y, span_input, nsteps_input, delta_input);

                return lowess.Eval(xp);
            }
        }
    }
}
