using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    public static class MathUtils
    {
        [ExcelFunction(Description = "Compute first derivative: finite difference central 3pt", Category = AddInInfo.Category)]
        public static object acq_diff1_c3pt(double[] x, double[] y)
        {
            object result = ExcelError.ExcelErrorNA;

            if (x == null || y == null || x.Length != y.Length || x.Length != 3)
            {
                result = ExcelError.ExcelErrorNA;            
            }
            else
            {
                if (x[2] > x[1] && x[1] > x[0])
                {
                    double dx = x[2] - x[0];
                    double dx0 = x[1] - x[0];
                    double dx1 = x[2] - x[1];

                    double dy0 = (y[1] - y[0]) / dx0;
                    double dy1 = (y[2] - y[1]) / dx1;

                    double a = dx0 / dx;
                    double b = 1.0 - a;

                    result = dy0 * b + dy1 * a;
                }
            }
            return result;
        }

        [ExcelFunction(Description = "Compute second derivative: finite difference central 3pt", Category = AddInInfo.Category)]
        public static object acq_diff2_c3pt(double[] x, double[] y)
        {
            object result = ExcelError.ExcelErrorNA;

            if (x == null || y == null || x.Length != y.Length || x.Length != 3)
            {
                result = ExcelError.ExcelErrorNA;
            }
            else
            {
                if (x[2] > x[1] && x[1] > x[0])
                {
                    double dx = x[2] - x[0];
                    double dx0 = x[1] - x[0];
                    double dx1 = x[2] - x[1];

                    double dy0 = (y[1] - y[0]) / dx0;
                    double dy1 = (y[2] - y[1]) / dx1;

                    result = (dy1 - dy0) / dx;
                }
            }
            return result;
        }
    }
}
