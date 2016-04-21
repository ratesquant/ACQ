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
            TestSort();
            //TestInterpolation();
        }

        public static void TestSort()
        {
            System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();

            int n = 1024;

            double[] x = new double[n];

            ACQ.Math.Random.RandomBase rng = new ACQ.Math.Random.MersenneTwister();

            for (int i = 0; i < x.Length; i++)
            {
                x[i] = rng.NextDouble();
            }

            double[] xs = (double[])x.Clone();

            Array.Sort<double>(xs);

            for (int i = 0; i < n; i++)
            {
                int m = (int)System.Math.Floor(rng.NextDouble() * n);

                timer.Restart();

                ACQ.Math.Sort.Psort(x, m);
                Console.WriteLine("x[{0}] = {1}, {2}, psorted in {3} ", m, x[m], x[m] - xs[m], timer.ElapsedMilliseconds);
            }
 
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
