using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    class ForwardInterpolation : InterpolationBase
    {
        public ForwardInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, bounds)
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

                value = y1;
            }

            return value;
        }
    }
}
