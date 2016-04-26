using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Linear Interpolation
    /// </summary>
    public class LinearInterpolation : InterpolationBase
    {
        public LinearInterpolation(double[] x, double[] y)
            : base(x, y)
        { }

        public override double Eval(double x)
        {  
            double value;

            int index = FindIndex(x, out value);

            if(index > 0 )
            {
                double x0, x1, y0, y1;

                Bracket(index, out x0, out x1, out y0, out y1);
                
                double dx = x1 - x0;
                double b = (x - x0) / dx;

                value = y0 + b * (y1 - y0);
            }

            return value;
        }
    }
}
