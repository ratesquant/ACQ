using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ACQ.Math;

namespace ACQ.Test
{
    class Program
    {
        static void Main(string[] args)
        {
            TestInterpolation();
        }

        public static void TestInterpolation()
        {
            ACQ.Math.Interpolation.InterpolationFactory factory = new Math.Interpolation.InterpolationFactory();

            double[] x = new double[]{1,2,3};
            double[] y = new double[]{1,2,3};
            double[,] f = new double[3, 3];

            ACQ.Math.Interpolation.InterpolationFactory.GetInterpolator("Linear", x, y, true);

            ACQ.Math.Interpolation.InterpolationFactory2D.GetInterpolator("Bilinear", x, y, f);
 
        }
    }
}
