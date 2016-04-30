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
            TestLinearAlgebra();
            //TestSort();
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
        public static void TestLinearAlgebra()
        {  
            double eps = ACQ.Math.Const.epsilon;
            
            for (int n = 3; n <= 128; n++) 
            {
                ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateMagic(n);

                double trace = M.Trace(); // trace   = diagonal sum, should be the magic sum, (n^3 + n)/2.

                ACQ.Math.Linalg.EigenvalueDecomposition E = new ACQ.Math.Linalg.EigenvalueDecomposition(0.5 * (M+(M.Transpose())));
                double[] d = E.GetRealEigenvalues();    //maximum eigenvalue of (A + A')/2, should equal trace             

                //int r = M.rank();   //should equal n if n is odd, be less than n if n is even
                //double c = M.cond();

                ACQ.Math.Linalg.LuDecomposition LU = new ACQ.Math.Linalg.LuDecomposition(M);
                ACQ.Math.Linalg.Matrix L = LU.GetL();
                ACQ.Math.Linalg.Matrix U = LU.GetU();
                int[] p = LU.GetPivotIndex();
                ACQ.Math.Linalg.Matrix R = L * U - M.Submatrix(p, 0, n - 1);
                double lu_res = R.Norm1() / (n * eps);

                ACQ.Math.Linalg.QrDecomposition QR = new ACQ.Math.Linalg.QrDecomposition(M);
                ACQ.Math.Linalg.Matrix Q = QR.GetQ();
                R = QR.GetR();
                R = Q * R - M;
                double qr_res = R.Norm1() / (n * eps);

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                U = SV.GetU();
                ACQ.Math.Linalg.Matrix S = SV.GetS();
                ACQ.Math.Linalg.Matrix V = SV.GetV();                
                R = U * S * V.Transpose() - M;
                double sv_res = R.Norm1() / (n * eps);
                double rank = SV.Rank();
                double cond = SV.Cond();


                Console.WriteLine("n:{0}, trace:{1}, ev:{2}, lu:{3}, qr:{4}, sv:{5}", n, trace, d[d.Length - 1], lu_res, qr_res, sv_res); 
            } 
        }
        public static void TestInterpolation()
        {
            ACQ.Math.Interpolation.InterpolationFactory factory = new Math.Interpolation.InterpolationFactory();

            double[] x = new double[]{1,2,3};
            double[] y = new double[]{1,2,3};
            double[,] f = new double[3, 3];

            ACQ.Math.Interpolation.InterpolationFactory.GetInterpolator("Linear", x, y);

            ACQ.Math.Interpolation.InterpolationFactory2D.GetInterpolator("Bilinear", x, y, f);
 
        }
    }
}
