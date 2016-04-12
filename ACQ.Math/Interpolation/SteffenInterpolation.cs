using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class SteffenInterpolation : HermiteInterpolation
    {
        public SteffenInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, bounds)
        {
        }

        protected override void compute_coefficients(double[] x, double[] y, out double[] a)
        {
            int n = x.Length;

            a = new double[n];


            for (int i = 1; i < n - 1; i++)
            {
                double dx = x[i + 1] - x[i - 1];
                double dx0 = x[i] - x[i - 1];
                double dx1 = x[i + 1] - x[i];

                double dy0 = (y[i] - y[i - 1]) / dx0;
                double dy1 = (y[i + 1] - y[i]) / dx1;

                double p = (dy0 * dx1 + dy1 * dx0) / dx;

                //steffen conditions that limit derivatives

                if (dy0 * dy0 <= 0.0)
                {
                    a[i] = 0.0;
                }
                else if (System.Math.Abs(p) > 2 * System.Math.Abs(dy0) || System.Math.Abs(p) > 2 * System.Math.Abs(dy1))
                {
                    a[i] = 2.0 * Utils.Sign(dy0) * System.Math.Min(System.Math.Abs(dy0), System.Math.Abs(dy1));

                }else
                {
                    a[i] = p;
                }
            }

            a[0] = (y[0] - y[1]) / (x[0] - x[1]);
            a[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

        }

    }
}
