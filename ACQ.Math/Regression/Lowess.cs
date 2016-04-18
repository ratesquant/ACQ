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
    /// Lowess regression, adapted from R source code, 
    /// </summary>
    public class Lowess
    {
        double[] m_x;
        double[] m_y;
        double[] m_rw;

        /// <summary>
        /// Construct Loess
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="nsteps"> number of robustifying iterations which should be performed. (default = 3) </param>
        /// <param name="span">smoother span. default(2/3) This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness</param>
        public Lowess(double[] x, double[] y, int nsteps = 3, double span = 2.0/3.0)
        {
            double delta = 0.001;
            m_x = (double[])x.Clone();
            m_y = (double[])y.Clone();

            double[] ys = new double[x.Length];
            double[] rw = new double[x.Length];
            double[] res = new double[x.Length];

            clowess(m_x, m_y, span, nsteps, delta, ys, rw, res);
 
        }
        private static void lowest(double[] x, double[] y, int n, double xs, out double ys, int nleft, int nright, double[] w, bool userw, double[] rw, out bool ok)
        {
            int nrt, j;
            double a, b, c, h, h1, h9, r, range;

            //x--;
            //y--;
            //w--;
            //rw--;

            range = x[n] - x[1];
            h = System.Math.Max(xs - x[nleft], x[nright] - xs);
            h9 = 0.999 * h;
            h1 = 0.001 * h;

            // sum of weights 

            a = 0.0;
            j = nleft;
            while (j <= n)
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
                ok = false;
            }
            else
            {
                ok = true;
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
        }
        private static void rPsort(double[] x, int n, int k)
        {
            rPsort2(x, 0, n - 1, k);
        }
        private static void rPsort2(double[] x, int lo, int hi, int k)
        {
            double v, w;
            bool nalast = true;
            int L, R, i, j;

            for (L = lo, R = hi; L < R; )
            {
                v = x[k];
                for (i = L, j = R; i <= j; )
                {
                    while (rcmp(x[i], v, nalast) < 0) i++;
                    while (rcmp(v, x[j], nalast) < 0) j--;
                    if (i <= j) { w = x[i]; x[i++] = x[j]; x[j--] = w; }
                }
                if (j < k) L = i;
                if (k < i) R = j;
            }
        }
        private static int rcmp(double x, double y, bool nalast)
        {
            bool nax = Double.IsNaN(x);
            bool nay = Double.IsNaN(y);
            if (nax && nay) return 0;
            if (nax) return nalast ? 1 : -1;
            if (nay) return nalast ? -1 : 1;
            if (x < y) return -1;
            if (x > y) return 1;
            return 0;
        }
        private static void clowess(double[] x, double[] y, double span, int nsteps, double delta, double[] ys, double[] rw, double[] res)
        {
            int i, iter, j, last, m1, m2, nleft, nright, ns;
            bool ok;
            double alpha, c1, c9, cmad, cut, d1, d2, denom, r, sc;
            int n = x.Length;

            if (n < 2)
            {
                ys[0] = y[0];
                return;
            }


            // at least two, at most n points
            ns = System.Math.Max(2, System.Math.Min(n, (int)(span * n + 1e-7)));


            // robustness iterations

            iter = 1;
            while (iter <= nsteps + 1)
            {
                nleft = 1;
                nright = ns;
                last = 0;	/* index of prev estimated point */
                i = 1;		/* index of current point */

                for (; ; )
                {
                    if (nright < n)
                    {
                        // move nleft,  nright to right 
                        // if radius decreases 

                        d1 = x[i] - x[nleft];
                        d2 = x[nright + 1] - x[i];

                        // if d1 <= d2 with 
                        // x[nright+1] == x[nright], 
                        // lowest fixes 

                        if (d1 > d2)
                        {

                            // radius will not 
                            // decrease by 
                            // move right 

                            nleft++;
                            nright++;
                            continue;
                        }
                    }

                    // fitted value at x[i]
                    double y_temp;
                    lowest(x, y, n, x[i], out y_temp, nleft, nright, res, iter > 1, rw, out ok);
                    if (ok)
                    {
                        ys[i] = y_temp;
                    }else
                    {
                        ys[i] = y[i];
                    }

                    // all weights zero 
                    // copy over value (all rw==0) 

                    if (last < i - 1)
                    {
                        denom = x[i] - x[last];

                        // skipped points -- interpolate 
                        // non-zero - proof? 

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
                    for (i = last + 1; i <= n; i++)
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
                    if (last >= n)
                    {
                        break;
                    }
                }
                // residuals
                for (i = 0; i < n; i++)
                {
                    res[i] = y[i + 1] - ys[i + 1];
                }

                // overall scale estimate */
                sc = 0.0;
                for (i = 0; i < n; i++)
                {
                    sc += System.Math.Abs(res[i]);
                }
                sc /= n;

                // compute robustness weights 
                // except last time 

                if (iter > nsteps)
                    break;

                for (i = 0; i < n; i++)
                {
                    rw[i] = System.Math.Abs(res[i]);
                }

                /* Compute   cmad := 6 * median(rw[], n)  ---- */
                m1 = n / 2;
                // partial sort, for m1 & m2 
                rPsort(rw, n, m1);
                if (n % 2 == 0)
                {
                    m2 = n - m1 - 1;
                    rPsort(rw, n, m2);
                    cmad = 3.0 * (rw[m1] + rw[m2]);
                }
                else
                { // n odd 
                    cmad = 6.0 * rw[m1];
                }

                if (cmad < 1e-7 * sc) // effectively zero 
                    break;
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
                iter++;
            }
        }
    }
}
