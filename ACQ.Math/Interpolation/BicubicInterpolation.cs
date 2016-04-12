using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class BicubicInterpolation : InterpolationBase2D
    {
        public BicubicInterpolation(double[] x1, double[] x2, double[,] y)
            : base(x1, x2, y, false)
        {
        }
        public BicubicInterpolation(double[] x1, double[] x2, double[,] y, bool copyData)
            : base(x1, x2, y, copyData)
        {
        }

        public override double Eval(double x1, double x2)
        {
            double value = Double.NaN;

            int index1, index2;


            FindIndex(x1, x2, out index1, out index2);

            if (index1 > 0 && index2 > 0)
            {
                //      x10       x11  
                //x20   y10  y1e  y11
                //           val
                //x21   y20  y2e  y21 

                int n1 = m_x1.Length;
                //interpolation along x1, between index2-1 and index2
                double x10 = m_x1[index1 - 1];
                double x11 = m_x1[index1];
                double x20 = m_x2[index2 - 1];
                double x21 = m_x2[index2];
                double y10 = m_y[index2 - 1, index1 - 1];
                double y11 = m_y[index2 - 1, index1];
                double y20 = m_y[index2, index1 - 1];
                double y21 = m_y[index2, index1];

                //x1 derivatives
                double y10dx1 = dydx1(m_y, m_x1, index1 - 1, index2 - 1);
                double y11dx1 = dydx1(m_y, m_x1, index1    , index2 - 1);
                double y20dx1 = dydx1(m_y, m_x1, index1 - 1, index2);
                double y21dx1 = dydx1(m_y, m_x1, index1,     index2);

                //x2 derivatives
                double y10dx2 = dydx2(m_y, m_x2, index1 - 1, index2 - 1);
                double y11dx2 = dydx2(m_y, m_x2, index1,     index2 - 1);
                double y20dx2 = dydx2(m_y, m_x2, index1 - 1, index2);
                double y21dx2 = dydx2(m_y, m_x2, index1,     index2);

                //TODO: cross derivatives
                double y10dx1dx2 = 0;
                double y11dx1dx2 = 0;
                double y20dx1dx2 = 0;
                double y21dx1dx2 = 0;

                double b1 = (x1 - x10) / (x11 - x10);
                double b2 = (x2 - x20) / (x21 - x20);

                double h11, h21, h31, h41, h12, h22, h32, h42;
                HermiteInterpolation.hermite_basis(b1, out h11, out h21, out h31, out h41);
                HermiteInterpolation.hermite_basis(b2, out h12, out h22, out h32, out h42);

                //h1 * y0 + h2 * y1 + dx * (m_c[index - 1] * h3 + m_c[index] * h4);
                value  = h11 * h12 * y10 + h21 * h12 * y11 + h12 * h11 * y20 + h22 * h21 * y20; // values at the nodes 
                //TODO: add derivatives and cross derivatives
                
            }

            return value;
        }

        /// <summary>
        /// compute first derivatives using finite difference (central in the middle, forward at the boundary)
        /// </summary>
        /// <param name="y"></param>
        /// <param name="x"></param>
        /// <param name="index1"></param>
        /// <param name="index2"></param>
        /// <returns></returns>
        private static double dydx1(double[,] y, double[] x, int index1, int index2)
        {
            double dy = 0;

            if (index1 == 0)
            {
                dy = (y[index2, index1 + 1] - y[index2, index1]) / (x[index1 + 1] - x[index1]);
            }
            else if (index1 == x.Length - 1)
            {
                dy = (y[index2, index1] - y[index2, index1 - 1]) / (x[index1] - x[index1 - 1]);
            }
            else 
            {
                double dx  = x[index1 + 1] - x[index1 - 1];
                double dx0 = x[index1]     - x[index1 - 1];
                double dx1 = x[index1 + 1] - x[index1];

                double dy0 = (y[index2, index1] - y[index2, index1 - 1])/dx0;
                double dy1 = (y[index2, index1 + 1] - y[index2, index1])/dx1;

                dy = (dy0 * dx1 + dy1 * dx0) / dx;
             }

            return dy;
        }

        private static double dydx2(double[,] y, double[] x, int index1, int index2)
        {
            //can be combined with dyx1, but separate it for now for clarity 
            double dy = 0;

            if (index2 == 0)
            {
                dy = (y[index2 + 1, index1] - y[index2, index1]) / (x[index2 + 1] - x[index2]);
            }
            else if (index2 == x.Length - 1)
            {
                dy = (y[index2, index1] - y[index2 - 1, index1]) / (x[index2] - x[index2 - 1]);
            }
            else
            {
                double dx = x[index2 + 1] - x[index2 - 1];
                double dx0 = x[index2] - x[index2 - 1];
                double dx1 = x[index2 + 1] - x[index2];

                double dy0 = (y[index2, index1] - y[index2 - 1, index1]) / dx0;
                double dy1 = (y[index2 + 1, index1] - y[index2, index1]) / dx1;

                dy = (dy0 * dx1 + dy1 * dx0) / dx;
            }

            return dy;
        }
    }
}
