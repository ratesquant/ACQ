using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

using ACQ.Math;

namespace ACQ.Test
{
    class Program
    {
        static void Main(string[] args)
        {
            TestInterpolation2();
            //TestThreadSafety();
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

        public static void TestThreadSafety()
        {
            ThreadPool.SetMaxThreads(16, 16);

            for (int i = 0; i < 32*1024; i++)
            {
                Console.WriteLine("Adding {0}", i);
                ThreadPool.QueueUserWorkItem(addin_function, i); //simulate excel calling our add-in function                                
            }

            //System.Threading.Thread thread = new System.Threading.Thread(ExecuteInForeground);
            //thread.Start(10);                       
            do
            {
                System.Threading.Thread.Sleep(500);                
            } while (true);         
        }

        private static void addin_function(object a)
        {
            Console.WriteLine("Thread {0}: {1}", Thread.CurrentThread.ManagedThreadId, a);

            //here we capture variable "a" and then compare it with "b" which is passed later
            System.Threading.Thread thread = new System.Threading.Thread(
                (object b) =>
                {
                    Console.WriteLine("Thread {0}: {1}, a:{2} b:{3}, {4}",
                        Thread.CurrentThread.ManagedThreadId,
                        Thread.CurrentThread.ThreadState, a, b, (int)a - (int)b);
                    if (a != b)
                    {
                        throw new Exception("Different values detected");
                    }                    
                });
            System.Threading.Thread.Sleep(20); //sleep to ensure that other threads can enter before we can finish
            thread.Start(a);            
        }

        public static void TestLinearAlgebra()
        {  
            double eps = ACQ.Math.Const.epsilon;
            
            //test magic matrices
            for (int n = 3; n <= 32; n++) 
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
                double lu_res = R.Norm1();

                ACQ.Math.Linalg.QrDecomposition QR = new ACQ.Math.Linalg.QrDecomposition(M);
                ACQ.Math.Linalg.Matrix Q = QR.GetQ();
                R = QR.GetR();
                R = Q * R - M;
                double qr_res = R.Norm1();

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                U = SV.GetU();
                ACQ.Math.Linalg.Matrix S = SV.GetS();
                ACQ.Math.Linalg.Matrix V = SV.GetV();                
                R = U * S * V.Transpose() - M;
                double sv_res = R.Norm1();
                double rank = SV.Rank();
                double cond = SV.Cond();

                Console.WriteLine("n:{0,3} trace:{1:E4} ev:{2:E4} lu:{3:E4} qr:{4:E4} sv:{5:E4}", n, trace, d[d.Length - 1], lu_res, qr_res, sv_res); 
            }

            //test random matrices
            ACQ.Math.Random.IRandomGenerator rng = new ACQ.Math.Random.MersenneTwister();

            for (int n = 1; n <= 8; n++)
            {
                ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateRandom(n, n, rng);

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                ACQ.Math.Linalg.Matrix U = SV.GetU();
                ACQ.Math.Linalg.Matrix S = SV.GetS();
                ACQ.Math.Linalg.Matrix V = SV.GetV();
                ACQ.Math.Linalg.Matrix R = U * S * V.Transpose() - M;
                ACQ.Math.Linalg.Matrix pinv = SV.PseudoInverse();                
                double pinv_res1 = (M * pinv * M - M).Norm1();
                double pinv_res2 = (pinv * M * pinv - pinv).Norm1();                
                double sv_res = R.Norm1();
                double rank = SV.Rank();
                double cond = SV.Cond();

                Console.WriteLine("n:{0,3} rank:{1:E4} cond:{2:E4} sv:{3:E4} pinv:{4:E4} {5:E4}", n, rank, cond, sv_res, pinv_res1, pinv_res2); 
            }

            //non square matrixes, m > n
            for (int n = 1; n <= 32; n++)
            {
                for (int m = n; m <= 32; m++)
                {
                    ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateRandom(m, n, rng);

                    double err_scale = 1.0 / ( System.Math.Max(m, n) * eps );
                    ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                    ACQ.Math.Linalg.Matrix U = SV.GetU();
                    ACQ.Math.Linalg.Matrix S = SV.GetS();
                    ACQ.Math.Linalg.Matrix V = SV.GetV();
                    ACQ.Math.Linalg.Matrix R = U * S * V.Transpose() - M;                    
                    ACQ.Math.Linalg.Matrix pinv = SV.PseudoInverse();                    
                    double pinv_res1 = (M * pinv * M - M).Norm1() * err_scale;
                    double pinv_res2 = (pinv * M * pinv - pinv).Norm1() * err_scale;
                    double sv_res = R.Norm1() * err_scale;
                    double rank = SV.Rank();
                    double cond = SV.Cond();

                    Console.WriteLine("{0} x {1}, rank:{2:F4}, cond:{3:F4}, sv:{4:F4}, pinv:{5:F4}, {6:F4}", m, n, rank, cond, sv_res, pinv_res1, pinv_res2);
                }
            }

            //linear systems: square with full rank
            for (int n = 1; n <= 32; n++)
            {
                int m = 10;
                ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateRandom(n, n, rng);
                ACQ.Math.Linalg.Matrix B = new Math.Linalg.Matrix(n, m);

                for (int i = 0; i < M.Rows; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < M.Columns; j++)
                    {
                        sum += M[i, j];
                    }

                    for (int k = 0; k < B.Columns; k++)
                    {
                        B[i, k] = sum * (k + 1);                        
                    }
                }

                double err_scale = 1.0 / (n * eps);

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                ACQ.Math.Linalg.Matrix X1 = SV.Solve(B);                                
                ACQ.Math.Linalg.Matrix X2 = SV.PseudoInverse() * B;
                double sv_res1 = (M * X1 - B).Norm1() * err_scale;
                double sv_res2 = (M * X2 - B).Norm1() * err_scale;

                ACQ.Math.Linalg.LuDecomposition LU = new ACQ.Math.Linalg.LuDecomposition(M);
                X1 = LU.Solve(B);                
                double lu_res1 = (M * X1 - B).Norm1() * err_scale;

                ACQ.Math.Linalg.QrDecomposition QR = new ACQ.Math.Linalg.QrDecomposition(M);                
                X1 = QR.Solve(B);
                double qr_res1 = (M * X1 - B).Norm1() * err_scale;

                Console.WriteLine("{0}, svd:{1:F4}, sv:{2:F4} lu:{3:F4} qr:{4:F4}", n, sv_res1, sv_res2, lu_res1, qr_res1);       
            }

            //linear systems: symmetric with full rank
            for (int n = 1; n <= 32; n++)
            {
                int m = 10;
                ACQ.Math.Linalg.Matrix A = ACQ.Math.Linalg.Matrix.CreateRandom(n, n, rng);
                ACQ.Math.Linalg.Matrix M = A * A.Transpose();
                ACQ.Math.Linalg.Matrix B = new Math.Linalg.Matrix(n, m);

                for (int i = 0; i < M.Rows; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < M.Columns; j++)
                    {
                        sum += M[i, j];
                    }

                    for (int k = 0; k < B.Columns; k++)
                    {
                        B[i, k] = sum * (k + 1);
                    }
                }

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                ACQ.Math.Linalg.Matrix X1 = SV.Solve(B);
                ACQ.Math.Linalg.Matrix X2 = SV.PseudoInverse() * B;
                double sv_res1 = (M * X1 - B).Norm1();
                double sv_res2 = (M * X2 - B).Norm1();

                ACQ.Math.Linalg.LuDecomposition LU = new ACQ.Math.Linalg.LuDecomposition(M);
                X1 = LU.Solve(B);
                double lu_res1 = (M * X1 - B).Norm1();

                ACQ.Math.Linalg.QrDecomposition QR = new ACQ.Math.Linalg.QrDecomposition(M);
                X1 = QR.Solve(B);
                double qr_res1 = (M * X1 - B).Norm1();

                ACQ.Math.Linalg.CholeskyDecomposition CD = new ACQ.Math.Linalg.CholeskyDecomposition(M);
                X1 = CD.Solve(B);
                double cd_res1 = (M * X1 - B).Norm1();

                Console.WriteLine("{0,3}, sv:{1:E4} pinv:{2:E4} lu:{3:E4} qr:{4:E4}, cd:{5:E4}", n, sv_res1, sv_res2, lu_res1, qr_res1, cd_res1);
            }

            ACQ.Math.Random.IRandomGenerator rnorm = new ACQ.Math.Random.BoxMuller(rng);

            //linear systems: overdetermined, full rank (least squares solution)
            for (int n = 1; n <= 32; n++)
            {
                ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateRandom(2 * n, n, rnorm);
                ACQ.Math.Linalg.Matrix B = ACQ.Math.Linalg.Matrix.CreateRandom(2 * n, 1, rnorm);

                for (int i = 0; i < M.Rows; i++)
                {
                    double sum = 0;
                    for (int k = 0; k < M.Columns; k++)
                    {
                        sum += M[i, k];
                    }
                    B[i, 0] = sum + 0.1 * B[i, 0];
                }

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                ACQ.Math.Linalg.Matrix Xsv = SV.Solve(B);
                double sv_res1 = (M * Xsv - B).Norm1();
                
                ACQ.Math.Linalg.QrDecomposition QR = new ACQ.Math.Linalg.QrDecomposition(M);
                ACQ.Math.Linalg.Matrix Xqr = QR.Solve(B);
                double qr_res1 = (M * Xqr - B).Norm1();

                double res = (Xsv - Xqr).Norm1();

                Console.WriteLine("{0,3}, sv:{1:E4} qr:{2:E4} res:{3:E4}", n, sv_res1, qr_res1, res);
            }

            //linear systems: overdetermined but not full rank (least squares solution), 
            //QR does not work in this case,
            for (int n = 1; n <= 32; n++)
            {
                ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateRandom(2 * n, n + 1, rnorm);
                ACQ.Math.Linalg.Matrix B = ACQ.Math.Linalg.Matrix.CreateRandom(2 * n, 1, rnorm);

                for (int i = 0; i < M.Rows; i++)
                {
                    M[i, n] = M[i, n - 1];//dublicate column
                    double sum = 0;
                    for (int k = 0; k < M.Columns; k++)
                    {
                        sum += M[i, k];
                    }
                    B[i, 0] = sum + 0.1 * B[i, 0];
                }

                //solution shold be approximatly all ones, last two values should be the same
                //condition number is large (infinite actually, sometimes finite because of truncation errors) 
                //rank = M.Columns - 1 (i.e. rank = n)
                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M);
                ACQ.Math.Linalg.Matrix Xsv = SV.Solve(B); 
                double sv_res1 = (M * Xsv - B).Norm1();

                ACQ.Math.Linalg.QrDecomposition QR = new ACQ.Math.Linalg.QrDecomposition(M);
                double qr_res1 = Double.NaN;
                double res = Double.NaN;
                if (QR.IsFullRank)
                {
                    //occasionally QR thinks that matrix is full rank, 
                    //but it is poorly conditioned so solution is wrong 
                    ACQ.Math.Linalg.Matrix Xqr = QR.Solve(B);
                    qr_res1 = (M * Xqr - B).Norm1();
                    res = (Xsv - Xqr).Norm1();
                }

                Console.WriteLine("{0,3}, sv:{1:E4} qr:{2:E4} res:{3:E4}", n, sv_res1, qr_res1, res);
            }

            //linear systems: underdetermined (least norm solution)
            for (int n = 1; n <= 32; n++)
            {
                ACQ.Math.Linalg.Matrix M = ACQ.Math.Linalg.Matrix.CreateRandom(n, 2 * n + 1, rnorm);
                ACQ.Math.Linalg.Matrix B = ACQ.Math.Linalg.Matrix.CreateRandom(n, 1, rnorm);

                for (int i = 0; i < M.Rows; i++)
                {
                    M[i, 2 * n] = M[i, 2 * n - 1];//dublicate column
                }

                ACQ.Math.Linalg.SvDecomposition SV = new ACQ.Math.Linalg.SvDecomposition(M.Transpose());
                ACQ.Math.Linalg.Matrix Xsv = SV.Solve(B, true); //last two values should be the same
                double sv_res1 = (M * Xsv - B).Norm1(); //this should be zero 
                double x_norm = Xsv.Norm1();

                Console.WriteLine("{0,3}, sv:{1:E4} xnorm:{2:E4}", n, sv_res1, x_norm);
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

        public static void TestInterpolation2()
        {            
            double[] x = new double[] { 0, 1 };
            double[] y = new double[] { 0, 1 };            

            ACQ.Math.Interpolation.InterpolationInterface interpolator = new ACQ.Math.Interpolation.LinearInterpolation(x, y);

            double xi = 0.5;
            double yi = interpolator.Eval(xi);
        }
    }
}
