using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{

    /// <summary>
    /// The Akima interpolation is a continuously differentiable sub-spline interpolation. 
    /// Uses hermit(cubic) basis functions. 
    /// Hiroshi Akima, Journal of the ACM, Vol. 17, No. 4, October 1970, pages 589-602
    /// </summary>
    public class AkimaInterpolation : InterpolationBase
    {
        private readonly double[] m_b, m_c, m_d;
        private readonly bool m_periodic;

        public AkimaInterpolation(double[] x, double[] y, bool bounds = true)
            : this(x, y, bounds, false)
        {
        }
        public AkimaInterpolation(double[] x, double[] y, bool bounds = true, bool periodic = false)
            : base(x, y, bounds)
        {
            m_periodic = periodic;

            compute_coefficients(m_x, m_y, m_periodic, out m_b, out m_c, out m_d);
        }

        public override double Eval(double x)
        {
            double value;

            int index = FindIndex(x, out value);

            if(index > 0 )
            {
                index = index - 1;

                double dx = x - m_x[index];

                value = m_y[index] + dx * (m_b[index] + dx * (m_c[index] + dx * m_d[index]));
            }

            return value;
        }

        private static void compute_coefficients(double[] x, double[] y, bool periodicBoundary, out double[] b, out double[] c, out double[] d)
        {
            int n = x.Length;

            b = new double[n - 1];
            c = new double[n - 1];
            d = new double[n - 1];

            double[] s = new double[n + 3]; //slope of the function

            for (int i = 0; i < n - 1; i++)
            {
                s[i + 2] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            }

            if (periodicBoundary)
            {
                // periodic boundary conditions
                s[0] = s[n - 1];
                s[1] = s[n];
                s[n + 1] = s[2];
                s[n + 2] = s[3];
            }
            else
            {
                // non-periodic boundary conditions
                s[0] = 3.0 * s[2] - 2.0 * s[3];
                s[1] = 2.0 * s[2] - s[3];
                s[n + 1] = 2.0 * s[n] - s[n - 1];
                s[n + 2] = 3.0 * s[n] - 2.0 * s[n - 1];
            }

            //compute coefficients 
            for (int i = 0; i < n - 1; i++)
            {
                double ne = System.Math.Abs(s[i + 3] - s[i + 2]) + System.Math.Abs(s[i + 1] - s[i]);

                if (ne == 0.0)
                {
                    b[i] = s[i + 2];
                    c[i] = 0.0;
                    d[i] = 0.0;
                }
                else
                {
                    double dx = x[i + 1] - x[i];
                    double ne_next = System.Math.Abs(s[i + 4] - s[i + 3]) + System.Math.Abs(s[i + 2] - s[i + 1]);
                    double alpha_i = System.Math.Abs(s[i + 1] - s[i]) / ne;
                    double alpha_ip1;
                    double tl_ip1;
                    if (ne_next == 0.0)
                    {
                        tl_ip1 = s[i + 2];
                    }
                    else
                    {
                        alpha_ip1 = System.Math.Abs(s[i + 2] - s[i + 1]) / ne_next;
                        tl_ip1 = (1.0 - alpha_ip1) * s[i + 2] + alpha_ip1 * s[i + 3];
                    }
                    b[i] = (1.0 - alpha_i) * s[i + 1] + alpha_i * s[i + 2];
                    c[i] = (3.0 * s[i + 2] - 2.0 * b[i] - tl_ip1) / dx;
                    d[i] = (b[i] + tl_ip1 - 2.0 * s[i + 2]) / (dx * dx);
                }
            }
        }
    }

    public class AkimaPeriodicInterpolation : AkimaInterpolation
    {
        public AkimaPeriodicInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, bounds, true)
        { }
    }
}
