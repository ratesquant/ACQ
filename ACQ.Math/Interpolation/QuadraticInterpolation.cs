using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Quadratic Interpolation
    /// </summary>
    public class QuadraticInterpolation : InterpolationBase
    {
        private readonly double[] m_u;

        public QuadraticInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, bounds)
        {
            compute_coefficients(m_x, m_y, out m_u);
        }

        public override double Eval(double x)
        {
            double value;

            int index = FindIndex(x, out value);

            if (index > 0)
            {
                double x0, x1, y0, y1;

                Bracket(index, out x0, out x1, out y0, out y1);

                double dx = x1 - x0;
                double b = (x - x0) / dx;
                double q1, q2, q3;

                quadratic_basis(b, out q1, out q2, out q3);

                value = q1 * y0 + q3 * y1 + q2 * m_u[index - 1];
            }

            return value;
        }

        private void compute_coefficients(double[] x, double[] y, out double[] u)
        {
            int n = x.Length;

            u = new double[n - 1];

            double dy = y[1] - y[0];

            for (int i = 0; i < n - 1; i++)
            {
                u[i] = 0.25 * (dy + 3 * y[i] + y[i + 1]);
                dy = (y[i] + 3 * y[i + 1] - 4 * u[i]);
            }
 
        }

        private void quadratic_basis(double x, out double q1, out double q2, out double q3)
        {
            double a = 2.0 * x * x;
            q1 = 1 - 3 * x + a;
            q2 =  4 * x - 2 * a;
            q3 = -x + a;
        }
    }
}
