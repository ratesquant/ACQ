using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    class BackwardInterpolation : InterpolationBase
    {
        public BackwardInterpolation(double[] x, double[] y)
            : base(x, y)
        { }

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

                value = x < x1 ? y0 : y1;
            }

            return value;
        }
    }
}
