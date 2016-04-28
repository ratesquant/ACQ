using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class SteffenSoftInterpolation : SteffenInterpolation
    {        
        public SteffenSoftInterpolation(double[] x, double[] y)
            : base(x, y)
        {
        }

        protected override double beta()
        {
            return 3.0;//2.0 is default beta of Steffen method, use 3.0 for softer results             
        }
    }
}
