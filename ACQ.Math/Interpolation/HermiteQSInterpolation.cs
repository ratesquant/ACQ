using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Hermite Interpolation using quadratic finite difference slope
    /// </summary>
    public class HermiteQSInterpolation : HermiteInterpolation
    {
        public HermiteQSInterpolation(double[] x, double[] y, bool bounds = true)
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

                a[i] = (dy0 * dx1 + dy1 * dx0) / dx; // fit parabola to the three points and compute derivative of the parabola
            }

            a[0] = (y[1] - y[0]) / (x[1] - x[0]);
            a[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
        }

    }
}
