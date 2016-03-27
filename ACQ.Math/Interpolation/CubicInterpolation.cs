using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Cubic Spline Interpolation (C2) with continuous first and second derivatives
    /// (natural boundary conditions)
    /// </summary>
    public class CubicInterpolation : InterpolationBase
    {
        private readonly double[] m_c;
        private readonly bool m_periodic;

        public CubicInterpolation(double[] x, double[] y, bool bounds = true)
            : this(x, y, bounds, false)
        {
        }
        public CubicInterpolation(double[] x, double[] y, bool bounds = true, bool periodic = false)
            : base(x, y, bounds)
        {
            m_periodic = periodic;
            compute_coefficients(m_x, m_y, m_periodic, out m_c);
        }

        public override double Eval(double x)
        {  
            double value;

            int index = FindIndex(x, out value);

            if(index > 0 )
            {
                double x0 = m_x[index - 1];
                double x1 = m_x[index];
                double y0 = m_y[index - 1];
                double y1 = m_y[index];
                double dx = x1 - x0;
                double dy = y1 - y0;

                double t = (x - x0) / dx;
                double a = m_c[index - 1] * dx - dy;
                double b = -m_c[index] * dx + dy;

                value = (1 - t) * y0 + t * y1 + t * (1 - t) * (a * (1 - t) + b * t);
            }

            return value;
        }

        //see [Engeln-Mullges + Uhlig, p. 254]
        private static void compute_coefficients(double[] x, double[] y, bool periodic, out double[] c)
        {
            int n = x.Length;

            c = new double[n];

            //there are n-2 equations, that can be written as symmetric tridiagonal matrix 
            //2 equations at the boundary are given by natural conditions 
            // dx[i] = x[i+1] - x[i]
            // k[0] + 0.5 * k[1] = 1.5 * dy[0] /dx[0]
            // k[i-1] * dx[i] + k[i] * 2 * (dx[i] + dx[i-1]) + k[i+1] * dx[i-1] = 3 * (dy[i-1]*dx[i]/dx[x-1] + dy[i]*dx[i-1]/dx[i])
            // k[n-2] + 0.5 * k[n-1] = 1.5 * dy[n-2]/dx[n-2]

            double[] g = new double[n]; //diag elements

            c[0] = 0.5;
            g[0] = 1.5 * (y[1] - y[0])/(x[1] - x[0]);

            //do gauss elimination (TODO: choletsky)
            //1. eliminate lower off diagonal
            //2. scale diagonal elements to 1
            //3. keep track of upper off diagonal elements in c
            //4. rhs store in array: g
            for (int i = 1; i < n - 1; i++)
            {
                double dx0 = x[i] - x[i - 1];
                double dx1 = x[i + 1] - x[i];
                double dy0 = y[i] - y[i - 1];
                double dy1 = y[i + 1] - y[i];

                double di = 2 * (dx0 + dx1);
                double ei = c[i-1];
                double ec = dx1;
                double gi = 3 * (dy0 * dx1 / dx0 + dy1 * dx0 / dx1);

                double scale = 1.0/(di - ei * ec);
                g[i] = (gi - g[i - 1] * ec) * scale;
                c[i] = dx0 * scale; 
            }


            //back substitution
            double dxn = x[n - 1] - x[n - 2];
            double dyn = y[n - 1] - y[n - 2];
            double dn = 1.0 - 0.5 * c[n - 2];
            double gn = (1.5 * dyn / dxn - g[n - 2] * 0.5);

            c[n - 1] = gn / dn;

            for (int i = n - 2; i >= 0 ; i--)
            {
                c[i] = g[i] - c[i + 1] * c[i];
            }
        }
    }

    public class CubicPeriodicInterpolation : CubicInterpolation
    {
        public CubicPeriodicInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, bounds, false)
        { }
    }
}
