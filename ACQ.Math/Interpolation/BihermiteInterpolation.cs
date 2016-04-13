using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// Hermite cubic cpline interpolation on 2D rectangular grid
    /// this is sometimes called bicubic interpolation 
    /// </summary>
    public class BiHermiteInterpolation : BiInterpolation<HermiteInterpolation>
    {
        public BiHermiteInterpolation(double[] x1, double[] x2, double[,] y)
            : base(x1, x2, y, true)
        {
        }
        public BiHermiteInterpolation(double[] x1, double[] x2, double[,] y, bool copyData)
            : base(x1, x2, y, copyData)
        {
        }
    }
}
