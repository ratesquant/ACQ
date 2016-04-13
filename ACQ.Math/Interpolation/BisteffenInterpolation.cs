using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class BisteffenInterpolation : BiInterpolation<SteffenInterpolation>
    {
        public BisteffenInterpolation(double[] x1, double[] x2, double[,] y)
            : base(x1, x2, y, true)
        {
        }
        public BisteffenInterpolation(double[] x1, double[] x2, double[,] y, bool copyData)
            : base(x1, x2, y, copyData)
        {
        }

    }
}
