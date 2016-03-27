using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Hermite Cubic Spline Interpolation
    /// </summary>
    public class HermiteInterpolation : InterpolationBase
    {
        private readonly double[] m_c;

        public HermiteInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, bounds)
        {
            compute_coefficients(m_x, m_y, out m_c);
        }

        public override double Eval(double x)
        {
            double value;

            int index = FindIndex(x, out value);

            if (index > 0)
            {
                double x0 = m_x[index - 1];
                double x1 = m_x[index];
                double y0 = m_y[index - 1];
                double y1 = m_y[index];
                double dx = x1 - x0;
                double a = (x1 - x) / dx;
                double b = 1.0 - a;

                double h1, h2, h3, h4;

                hermite_basis(b, out h1, out h2, out h3, out h4);

                value = h1 * y0 + h2 * y1 + dx * (m_c[index - 1] * h3 + m_c[index] * h4);
            }

            return value;
        }

        protected void hermite_basis(double u, out double h1, out double h2, out double h3, out double h4)
        {
            double v = 1.0 - u;
            double z = u * v;
            h4 = -u * z; // h11 = -u^2*(1-u)
            h3 = z + h4; // h10 = u*(1-u)^2
            h1 = v + h3 + h4; //h00 = (1+2t)(1-t)^2
            h2 = 1.0 - h1; // h01 = t^2*(3-2t)    
        }

        private static void compute_coefficients(double[] x, double[] y, out double[] a)
        {
            int n = x.Length;

            a = new double[n];

            if (n == 2)
            {
                a[0] = a[1] = (y[0] - y[1]) / (x[0] - x[1]);
            }
            else
            {

                for (int i = 1; i < n - 1; i++)
                {
                    a[i] = 0.5 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) + (y[i] - y[i - 1]) / (x[i] - x[i - 1]));
                }

                a[0] = 1.5 * (y[0] - y[1]) / (x[0] - x[1]) - 0.5 * a[1];
                a[n - 1] = 1.5 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) - 0.5 * a[n - 2];
            }

        }
    }
}
