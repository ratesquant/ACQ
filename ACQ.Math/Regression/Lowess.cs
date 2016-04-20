using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Regression
{
    /*
    *  R : A Computer Langage for Statistical Data Analysis
    *  Copyright (C) 1996  Robert Gentleman and Ross Ihaka
    *  Copyright (C) 1999-2016 The R Core Team
    *
    *  This program is free software; you can redistribute it and/or modify
    *  it under the terms of the GNU General Public License as published by
    *  the Free Software Foundation; either version 2 of the License, or
    *  (at your option) any later version.
    *
    *  This program is distributed in the hope that it will be useful,
    *  but WITHOUT ANY WARRANTY; without even the implied warranty of
    *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    *  GNU General Public License for more details.
    *
    *  You should have received a copy of the GNU General Public License
    *  along with this program; if not, a copy is available at
    *  https://www.R-project.org/Licenses/
    */
    /// <summary>
    /// Lowess regression, adapted from R source code,http://svn.r-project.org/R/trunk/src/library/stats/src/lowess.c
    /// it attempts to faithfully reporoduce R code without any "improvments" 
    /// </summary>
    public class Lowess
    {
        ACQ.Math.Interpolation.LinearInterpolation m_interpolator;

        /// <summary>
        /// Construct Lowess
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="nsteps"> number of robustifying iterations which should be performed. (default = 3) </param>
        /// <param name="span">smoother span [0, 1], default = 2/3. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness</param>
        /// <param name="delta">Defaults to 1/100th of the range of x</param>
        public Lowess(double[] x, double[] y, double span = 2.0/3.0, int nsteps = 3, double delta = 0.0)
        {
            if (x == null || y == null)
                throw new ArgumentNullException("Lowess input arrays can not be null");

            if (x.Length != y.Length)
                throw new ArgumentException("Lowess input arrays x and y should have the same length");

            if (x.Length < 2)
                throw new ArgumentException("Lowess input arrays should have at least 2 nodes");

            if (nsteps < 0)
                throw new ArgumentException("Lowess nsteps must be >= 0");

            //the code below actually will work for any value of span (because it checks later), but R checks for negative span, so we as well for consistency 
            if (span < 0)
                throw new ArgumentException("Lowess span must be > 0");
            
            int n = x.Length;

            //check that x is sorted
            bool ordered = true;
            for (int i = 0; i < n - 1; i++)
            {
                if (x[i + 1] < x[i]) //same values are allowed
                {
                    ordered = false;
                    break;
                }
            }

            double[] x_in = x;
            double[] y_in = y;

            if (!ordered)
            {
                //make copy and sort if provided array is not ordered
                x_in = (double[])x.Clone();
                y_in = (double[])y.Clone();

                Array.Sort<double, double>(x_in, y_in);
            }

            if (delta <= 0.0)
            {
                delta = 0.01 * (x_in[n - 1] - x_in[0]); //same default as in R
            }

            double[] ys = new double[n];
           
            clowess(x_in, y_in, ys, span, nsteps, delta); //this is main r - function 

            //check that all values of x_in are uniques
            bool unique = true;
            for (int i = 0; i < n - 1; i++)
            {
                if (x[i + 1] == x[i]) //same values are allowed
                {
                    unique = false;
                    break;
                }
            }

            if (unique)
            {
                m_interpolator = new ACQ.Math.Interpolation.LinearInterpolation(x_in, ys);
            }
            else 
            {
                List<double> xu = new List<double>(n);
                List<double> yu = new List<double>(n);

                xu.Add(x_in[0]);
                yu.Add(ys[0]);

                double prev_x = x_in[0];

                for (int i = 1; i < n; i++)
                {
                    if(x_in[i] > prev_x)
                    {
                        xu.Add(x_in[i]);
                        yu.Add(ys[i]);
                        prev_x = x_in[i];
                    }
                }
                m_interpolator = new ACQ.Math.Interpolation.LinearInterpolation(xu.ToArray(), yu.ToArray());
            }
        }

        public double Eval(double xp)
        {
            double result = m_interpolator.Eval(xp);
            return result;
        }

        #region R routines
        private static bool lowest(double[] x, double[] y, double xs, out double ys, int nleft, int nright, double[] w, bool userw, double[] rw)
        {
            int n = x.Length;
            int nrt, j;
            double a, b, c, h, h1, h9, r, range;
            bool success = false;

            range = x[n-1] - x[0];
            h = System.Math.Max(xs - x[nleft], x[nright] - xs);
            h9 = 0.999 * h;
            h1 = 0.001 * h;

            // sum of weights 

            a = 0.0;
            j = nleft;
            while (j < n)
            {
                // compute weights 
                // (pick up all ties on right) 
                w[j] = 0.0;
                r = System.Math.Abs(x[j] - xs);
                if (r <= h9)
                {
                    if (r <= h1)
                        w[j] = 1.0;
                    else
                        w[j] = Utils.Cube(1.0 - Utils.Cube(r / h));
                    if (userw)
                    {
                        w[j] *= rw[j];
                    }
                    a += w[j];
                }
                else if (x[j] > xs)
                    break;
                j = j + 1;
            }

            // rightmost pt (may be greater
            // than nright because of ties)

            nrt = j - 1;
            if (a <= 0.0)
            {
                success = false;
            }
            else
            {
                success = true;
            }

            // weighted least squares 
            // make sum of w[j] == 1 

            for (j = nleft; j <= nrt; j++)
                w[j] /= a;
            if (h > 0.0)
            {
                a = 0.0;

                //  use linear fit 
                // weighted center of x values

                for (j = nleft; j <= nrt; j++)
                {
                    a += w[j] * x[j];
                }
                b = xs - a;
                c = 0.0;
                for (j = nleft; j <= nrt; j++)
                {
                    c += w[j] * Utils.Sqr(x[j] - a);
                }
                if (System.Math.Sqrt(c) > 0.001 * range)
                {
                    b /= c;

                    // points are spread out
                    // enough to compute slope 

                    for (j = nleft; j <= nrt; j++)
                    {
                        w[j] *= (b * (x[j] - a) + 1.0);
                    }
                }
            }
            ys = 0.0;
            for (j = nleft; j <= nrt; j++)
            {
                ys += w[j] * y[j];
            }

            return success;
        }

        private static void clowess(double[] x, double[] y, double[] ys, double f, int nsteps, double delta)
        {
            int i, j, last, m1, m2, nleft, nright, ns;
            double alpha, c1, c9, cmad, cut, d1, d2, denom, r;
            int n = x.Length; // n >= 2

            double[] rw = new double[n];
            double[] res = new double[n];

            // at least two, at most n points
            ns = System.Math.Max(2, System.Math.Min(n, (int)(f * n + 1.0e-7)));


            // robustness iterations
            for (int iter = 0; iter <= nsteps; iter++)
            {
                nleft = 0;
                nright = ns - 1;
                last = -1;	// index of prev estimated point 
                i = 0;		// index of current point

                for (; ; )
                {
                    if (nright < n - 1)
                    {
                        // move nleft,  nright to right 
                        // if radius decreases 

                        d1 = x[i] - x[nleft];
                        d2 = x[nright + 1] - x[i];

                        if (d1 > d2)
                        {
                            nleft++;
                            nright++;
                            continue;
                        }
                    }

                    // fitted value at x[i]
                    double y_temp;
                    bool success = lowest(x, y, x[i], out y_temp, nleft, nright, res, iter > 0, rw);

                    if (success)
                    {
                        ys[i] = y_temp;
                    }
                    else
                    {
                        ys[i] = y[i];
                    }

                    // all weights zero 
                    // copy over value (all rw==0) 

                    if (last < i - 1)
                    {
                        denom = x[i] - x[last];

                        // skipped points -- interpolate
                        for (j = last + 1; j < i; j++)
                        {
                            alpha = (x[j] - x[last]) / denom;
                            ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last];
                        }
                    }

                    // last point actually estimated 
                    last = i;

                    // x coord of close points 
                    cut = x[last] + delta;
                    for (i = last + 1; i < n; i++)
                    {
                        if (x[i] > cut)
                            break;
                        if (x[i] == x[last])
                        {
                            ys[i] = ys[last];
                            last = i;
                        }
                    }
                    i = System.Math.Max(last + 1, i - 1);
                    if (last >= n - 1)
                    {
                        break;
                    }
                }//loop over all points

                // residuals
                for (i = 0; i < n; i++)
                {
                    res[i] = y[i] - ys[i];
                }

                // overall scale estimate
                double sc = 0.0;
                for (i = 0; i < n; i++)
                {
                    sc += System.Math.Abs(res[i]);
                }
                sc /= n;

                // compute robustness weights, except last time 
                if (iter < nsteps)
                {
                    for (i = 0; i < n; i++)
                    {
                        rw[i] = System.Math.Abs(res[i]);
                    }

                    //Compute cmad = 6 * median(rw[], n) 
                    m1 = n / 2;
                    Array.Sort<double>(rw); //R does partial sort here to figure out median, do full a sort here for clarity
                    if (n % 2 == 0)
                    {
                        m2 = n - m1 - 1;
                        cmad = 3.0 * (rw[m1] + rw[m2]);
                    }
                    else
                    { // n odd 
                        cmad = 6.0 * rw[m1];
                    }

                    if (cmad < 1e-7 * sc) // effectively zero 
                    {
                        break;
                    }

                    c9 = 0.999 * cmad;
                    c1 = 0.001 * cmad;

                    for (i = 0; i < n; i++)
                    {
                        r = System.Math.Abs(res[i]);
                        if (r <= c1)
                            rw[i] = 1.0;
                        else if (r <= c9)
                            rw[i] = Utils.Sqr(1.0 - Utils.Sqr(r / cmad));
                        else
                            rw[i] = 0.0;
                    }
                }
            }
        }
        #endregion

    }
}
