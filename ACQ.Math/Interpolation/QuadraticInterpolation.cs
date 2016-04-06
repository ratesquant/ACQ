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

            if (n == 2)
            {
                //linear interpolation, in case we only have two points
                u[0] = 0.5 * (y[1] + y[0]);
            }
            else if (n == 3)
            {
                //fit parabola using 3 points (i.e. second order Lagrange polynomial) 
                for (int i = 0; i < 2; i++)
                {
                    double m = 0.5 * (x[i+1] + x[i]);
                    double c1 = ((m - x[1]) / (x[0] - x[1])) * ((m - x[2]) / (x[0] - x[2]));
                    double c2 = ((m - x[0]) / (x[1] - x[0])) * ((m - x[2]) / (x[1] - x[2]));
                    double c3 = ((m - x[1]) / (x[2] - x[0])) * ((m - x[0]) / (x[2] - x[1]));

                    u[i] = y[0] * c1 + y[1] * c2 + y[2] * c3;
                }
            }
            else
            {
                double dx1 = x[1] - x[0];
                double dx2 = x[2] - x[1];
                double dx3 = x[2] - x[0];

                double c1 = -(2 * dx1 + dx2) / (dx1 * dx3);
                double c2 = (dx1 + dx2) / (dx1 * dx2);
                double c3 = -dx1 / (dx2 * dx3) ;
 
                double dy = y[0] * c1 + y[1] * c2 + y[2] * c3;

                for (int i = 0; i < n - 1; i++)
                {
                    double dx = x[i + 1] - x[i];
                    u[i] = 0.25 * (dy * dx + 3 * y[i] + y[i + 1]);
                    dy = (y[i] + 3 * y[i + 1] - 4 * u[i])/dx;
                }
            }
 
        }

        private void quadratic_basis(double x, out double q1, out double q2, out double q3)
        {
            q1 = (1 - 2 * x) * (1 - x);
            q2 =  4 * x * (1 - x);
            q3 = (2 * x - 1) * x;
        }
    }
}
