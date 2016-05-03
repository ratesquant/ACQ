using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math
{
    public class Stats
    {
        public static double SumOfSquares(double[] x)
        {
            double mean;
            double ssd = SumOfSquaredDev(x, out mean);
            return ssd + mean * mean * x.Length;
        }

        public static double SumOfSquaredDev(double[] x)
        {
            double mean;
            return SumOfSquaredDev(x, out mean);
        }

        public static double SumOfSquaredDev(double[] x, out double mean)
        {
            double tss;

            if (x != null && x.Length > 0)
            {
                int n = x.Length;

                double sum = 0.0;

                for (int i = 0; i < n; i++)
                {
                    sum += x[i];
                }

                mean = sum / n;
                double sum2 = 0.0;
                double sum3 = 0.0;

                for (int i = 0; i < n; i++)
                {
                    double dx = (x[i] - mean);
                    sum2 += dx * dx;
                    sum3 += dx;
                }

                tss = sum2 - sum3 * sum3 / n;
            }
            else
            {
                tss = mean = Double.NaN;                
            }
            return tss;
        }

        public static double SumOfSquaredDev(double[] x, double[] w)
        {
            double mean;
            return SumOfSquaredDev(x, w, out mean);
        }

        public static double SumOfSquaredDev(double[] x, double[] w, out double mean)
        {
            double tss;

            if (w == null)
            {
                return SumOfSquaredDev(x, out mean); //unweighted version
            }

            if (x != null && x.Length > 0 && w.Length == x.Length)
            {
                int n = x.Length;

                double sum = 0.0;
                double sw = 0.0;

                for (int i = 0; i < n; i++)
                {
                    sum += x[i] * w[i];
                    sw += w[i];
                }

                mean = sum / sw;
                double sum2 = 0.0;
                double sum3 = 0.0;

                for (int i = 0; i < n; i++)
                {
                    double dx = (x[i] - mean);
                    sum2 += w[i] * dx * dx;
                    sum3 += w[i] * dx;
                }

                tss = sum2 - sum3 * sum3 / sw;
            }
            else
            {
                tss = mean = Double.NaN;
            }
            return tss;
        }  
    }
}
