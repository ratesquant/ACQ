using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    class MultiquadricsInterpolation : RbfInterpolation1D
    {
        public MultiquadricsInterpolation(double[] x, double[] y, bool bounds = true)
            : base(x, y, enRadialBasisFunction.Multiquadrics)
        { }
    }
}
