using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class BiLinearInterpolation : InterpolationBase2D
    {
        public BiLinearInterpolation(double[] x1, double[] x2, double[,] y)
            : base(x1, x2, y, true)
        {
        }
        public BiLinearInterpolation(double[] x1, double[] x2, double[,] y, bool copyData)
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

                //linear interpolation along x1, for index2-1 and index2
                double x10 = m_x1[index1 - 1];
                double x11 = m_x1[index1];
                double x20 = m_x2[index2 - 1];
                double x21 = m_x2[index2];
                double y10 = m_y[index2 - 1, index1 - 1];
                double y11 = m_y[index2 - 1, index1];
                double y20 = m_y[index2, index1 - 1];
                double y21 = m_y[index2, index1];

                double b1 = (x1 - x10) / (x11 - x10);
                double b2 = (x2 - x20) / (x21 - x20);

                double y1e = y10 + b1 * (y11 - y10);
                double y2e = y20 + b1 * (y21 - y20);

                value = y1e + b2 * (y2e - y1e);
            }

            return value;
        }
    }
}
