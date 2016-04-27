using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Exponential Tension Spline Interpoaltion. 
    /// P. Rentrop, "An Algorithm for the Computation of the Exponential Tension Spline", Numerische Mathematik 35(1):81-93 · February 1980    
    /// </summary>
    public class ExpTensionInterpolation : InterpolationBase
    {
        private readonly double[] m_p;  //names of variables are the same as in article (except for prefix m_) 
        private readonly double[] m_hp;
        private readonly double[] m_ph;
        private readonly double[] m_d;

        public ExpTensionInterpolation(double[] x, double[] y, double tension = 1.0)
            : base(x, y)
        {
            m_p = new double[x.Length - 1]; // tension paramaters are for intervals between nodes, therefore there are n-1 of them 

            ACQ.Math.Utils.FillArray(m_p, tension); //init all tensions

            compute_coefficients(x, y, m_p, out m_d, out m_hp, out m_ph);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="p">tension paramaters, size = n - 1 </param>
        public ExpTensionInterpolation(double[] x, double[] y, double[] tension)
            : base(x, y)
        {
            //if user specified longer array of tensions, dont throw an exception   
            if (tension == null || tension.Length < x.Length - 1)  
            {
                throw new ArgumentNullException("tension");
            }
            
            m_p = new double[x.Length - 1];

            for(int i = 0; i<m_p.Length; i++)
            {
                m_p[i] = tension[i];
            }

            compute_coefficients(x, y, m_p, out m_d, out m_hp, out m_ph);
        }

        public override double Eval(double x)
        {
            double value;

            int index = FindIndex(x, out value);

            if (index > 0)
            {
                int i1 = index - 1;
                double x0 = m_x[i1];
                double x1 = m_x[index];
                
                double dx = x1 - x0;                
                double t0 = (x - x0) / dx;
                double t1 = 1.0 - t0;

                if (m_hp[index - 1] > 0.5)
                {
                    double e0 = System.Math.Exp(-t0 * m_hp[i1]);
                    double e1 = System.Math.Exp(-t1 * m_hp[i1]);
                    double c = 1.0 - m_ph[i1] * m_ph[i1];

                    value = m_y[index] * t0 + m_y[i1] * t1 + (m_d[index] * (e1 * (1.0 - e0 * e0) / c - t0) + m_d[i1] * (e0 * (1.0 - e1 * e1) / c - t1)) / (m_p[i1] * m_p[i1]);
                }
                else
                {
                    double e0 = t0 * m_hp[i1];
                    double e1 = t1 * m_hp[i1];
                    double c = dx * dx / (1.0 + m_hp[i1] * m_hp[i1] * m_ph[i1]);

                    value = t0 * (m_y[index] + m_d[index] * c * (t0 * t0 * phi(e0 * e0) - m_ph[i1] )) + t1 * (m_y[i1] + m_d[i1] * c * (t1 * t1 * phi(e1 * e1) - m_ph[i1]));
                }
            }

            return value;
        }

        private static void compute_coefficients(double[] x, double[] y, double[] p, out double[] d, out double[] hp, out double[] ph)
        {
            //converted from procedure exsplcoeff, Rentrop article
            int n = x.Length - 1;

            d = new double[n + 1]; //[0, n]
            hp = new double[n];    //[0, n - 1]
            ph = new double[n];    //[0, n - 1]

            double[] dq = new double[n]; //[1, n - 1]
            double[] q = new double[n + 1]; //[0, n]
            double[] r = new double[n + 1]; //[0, n]


            double u = y[0];
            double c, c1, c2, w;

            for (int i = 1; i <= n; i++)
            {
                int i1 = i - 1;
                double h = x[i] - x[i1]; // should be > 0 
                double v = y[i];
                hp[i1] = System.Math.Abs(h * p[i1]);
                d[i] = (v - u) / h;
                u = v;

                if (hp[i1] > 0.5)
                {
                    ph[i1] = System.Math.Exp(-hp[i1]);
                    c = ph[i1] * ph[i1];
                    c1 = 1.0 - c;
                    c2 = c1 / hp[i1];
                    c1 = c1 * hp[i1] / h;
                    q[i] = (1.0 - c2 + c) / c1;
                    r[i] = (c2 - 2.0 * ph[i1]) / c1;
                }
                else
                {
                    c = hp[i1] * hp[i1];
                    ph[i1] = phi(c);
                    w = h / (1.0 + c * ph[i1]);
                    c = c * 0.25;
                    c2 = 1.0 + c * ph[i1];
                    q[i] = (0.5 * c2 * c2 - ph[i1]) * w;
                    r[i] = ph[i1] * w; 
                }
            }

            d[0] = u = 0.0;
            for (int i = 1; i < n; i++)
            {
                q[i] = q[i] + q[i + 1] - u * r[i];
                dq[i] = d[i + 1] - d[i];
                d[i] = dq[i] - u * d[i - 1];
                u = r[i + 1] / q[i];
            }

            d[n] = 0.0;
            for (int i = n - 1; i >= 1; i--)
            {
                d[i] = (d[i] - r[i + 1] * d[i + 1]) / q[i];
            }
        }
        /// <summary>
        /// (sinh(x) - x) / x^3, for x [0, 0.5]
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        private static double phi(double a)
        {
            return ((0.27713991169e-5 * a + 0.19840927713e-3) * a + 0.83333336379e-2)*a + 0.16666666666;
        }
    }
}
