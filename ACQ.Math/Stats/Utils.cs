using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;

namespace ACQ.Math.Stats
{
    public static class Utils
    {
       /// <summary>
       /// Mean with compensated summation, 
       /// Kahan summation algorithm minimizes the error when adding a sequence of finite precision floating point numbers.
       /// </summary>
       /// <param name="x"></param>
       /// <returns></returns>
        public static double KahanSum(double[] x)
        {
            double dSum = x[0];
            double c = 0.0;
            double y, t;

            for (int i = 1; i < x.Length; i++)
            {
                y = x[i] - c;  
                t = dSum + y;      //Sum is big, y small, so low-order digits of y are lost.
                c = (t-dSum) - y;  //(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
                dSum = t;
            }
            return dSum; 
        } 
        /*
        public static double KahanMean(double[] x)
        {
            return KahanSum(x) / x.Length;
        }
       
        public static double Sum(double[] x)
        {
            double dSum = x.Length == 0 ? Double.NaN : x[0];

            for (int i = 1; i < x.Length; i++)
                dSum += x[i];

            return dSum;
        }
        public static double Mean(double[] x)
        {
            return Sum(x) / x.Length;
        }*/

        public static double Max(double[] x)
        {
            if (x == null || x.Length == 0)
                return Double.NaN;

            double max = x[0];

            for (int i = 1; i < x.Length; i++)
            {
                if (x[i] > max)
                {
                    max = x[i];
                }
            }

            return max;
        }

        /// <summary>
        /// Standard deviation of a sample(unbiased), normalized by n-1 where n is the length of sample. 
        /// </summary>
        /// <param name="x">sample</param>
        /// <returns></returns>
        public static double Std(double[] x)
        {
            return Sqrt(Var(x));
        }

        /// <summary>
        /// Varience of a sample (unbiased), normalized by n-1 where n is the length of sample. 
        /// </summary>
        /// <param name="x">sample</param>
        /// <returns></returns>
        public static double Var(double[] x)
        {
            int n = x.Length;

            if (n == 1) return 0;

            double dVar, ep, mean, s;
            
            mean = Mean(x);

            dVar = 0;    ep = 0;
            for (int i = 0; i < n; i++)
            {
                s = x[i] - mean;
                ep += s;
                dVar += s * s;
            }
            dVar = (dVar - ep * ep / n) / (n - 1);

            return dVar;
        }

        /// <summary>
        /// Computes mean and variance by using roundoff stable algorithm
        /// </summary>        
        public static void MeanAndVar(double[] x, out double mean, out double var)
        {
            int count = 0;
            double alpha = x[0];
            double beta = 0;
            double temp;

            for (int i = 1; i < x.Length; i++)
            {
                count = i + 1;
                temp = (x[i] - alpha);
                alpha += temp / count;
                beta += i * temp * temp / count;
            }
            mean = alpha;
            var = beta / (count - 1);
        }


        /// <summary>
        /// Sample skewness (biased), E[(x-mean)^3]/std^3
        /// </summary>
        /// <param name="x">sample</param>
        /// <returns></returns>
        public static double Skewness(double[] x)
        {
            int n = x.Length;

            if (n == 1) return Double.NaN;

            double dSkew, ep, mean, s, t, dVar;
            
            mean = Mean(x);
            dSkew = 0; dVar = 0;  ep = 0;
            for (int i = 0; i < n; i++)
            {
                s = x[i] - mean;
                t = s * s;
                ep += s;                
                dVar += t;
                dSkew += t * s;
            }
            dVar = (dVar - ep * ep / n) / n;
            if (dVar == 0)
                dSkew = Double.NaN;
            else
                dSkew /= (n*dVar*Sqrt(dVar));

            return dSkew;
        }

        /// <summary>
        /// Sample kurtosis (biased), E[(x-mean)^4]/std^4
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Kurtosis(double[] x)
        {
            int n = x.Length;

            if (n == 1) return Double.NaN;

            double dKurt, ep, mean, s, t, dVar;

            mean = Mean(x);

            dKurt = 0; dVar = 0; ep = 0;
            for (int i = 0; i < n; i++)
            {
                s = x[i] - mean;
                t = s * s;
                ep += s;
                dVar += t;
                dKurt += t * t;
            }
            dVar = (dVar - ep * ep / n) / (n);

            if (dVar == 0)
                dKurt = Double.NaN;
            else
                dKurt /= (n * dVar * dVar);

            return dKurt;
        }

        /// <summary>
        /// Central moment of all orders(biased), E[(x-mean)^order]
        /// </summary>
        /// <param name="x"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public static double Moment(double[] x, int order)
        {
            int n = x.Length;            

            double dMoment, mean;

            mean = Mean(x);

            dMoment = 0;
            for (int i = 0; i < n; i++)
            {                
                dMoment += Pow(x[i] - mean, order);
            }


            return dMoment/n;            
        }


        public static double KahanBabushkaNeumaierSum(double[] x)
        {
            if (x == null || x.Length == 0)
                return Double.NaN;

            double sum = 0.0;
            double c = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                double xi = x[i];
                double t = sum + xi;

                if (System.Math.Abs(sum) >= System.Math.Abs(xi))
                    c += (sum - t) + xi;
                else
                    c += (xi - t) + sum;
                sum = t;
            }
            return sum + c;
        }

        /// <summary>
        /// Uses Kahan summation
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Sum(double[] x, bool ignore_na = false)
        {
            if (x == null || x.Length == 0)
                return Double.NaN;

            double sum = 0.0;
            double c = 0.0;

            if (ignore_na)
            {
                for (int i = 0; i < x.Length; i++)
                {
                    if (!Double.IsNaN(x[i]))
                    {
                        double y = x[i] - c;
                        double t = sum + y;

                        c = (t - sum) - y;
                        sum = t;
                    }
                }
            }
            else
            {
                for (int i = 0; i < x.Length; i++)
                {
                    double y = x[i] - c;
                    double t = sum + y;

                    c = (t - sum) - y;
                    sum = t;
                }
            }
            return sum;
        }

        /// <summary>
        /// Uses Kahan summation
        /// </summary>
        /// <param name="x"></param>
        /// <param name="ignore_na"></param>
        /// <returns></returns>
        public static double Mean(double[] x, bool ignore_na = false)
        {
            if (x == null || x.Length == 0)
                return Double.NaN;

            double mean;
            double sum = 0.0;
            double c = 0.0;

            if (ignore_na)
            {
                int count = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    if (!Double.IsNaN(x[i]))
                    {
                        double y = x[i] - c;
                        double t = sum + y;

                        c = (t - sum) - y;
                        sum = t;
                        count++;
                    }
                }
                mean = sum / count;
            }
            else
            {
                mean = Sum(x) / x.Length;
            }
            return mean;
        }

        public static double WeightedMean(double[] x, double[] w, bool ignore_na = false)
        {
            if (x == null || x.Length == 0)
                return Double.NaN;

            double result = 0;

            if (w != null) //weighted sum
            {
                if (x.Length != w.Length)
                    return Double.NaN;

                if (ignore_na)
                {
                    double sum = 0, sum_w = 0;
                    double c = 0, c_w = 0;
                    for (int i = 0; i < x.Length; i++)
                    {
                        if (!Double.IsNaN(x[i]) && !Double.IsNaN(w[i]))
                        {
                            double y = x[i] * w[i] - c;
                            double t = sum + y;

                            c = (t - sum) - y;
                            sum = t;

                            y = w[i] - c_w;
                            t = sum_w + y;

                            c = (t - sum_w) - y;
                            sum_w = t;
                        }
                    }

                    result = sum / sum_w;
                }
                else
                {
                    double sum = 0;
                    double c = 0;

                    for (int i = 0; i < x.Length; i++)
                    {
                        double y = x[i] * w[i] - c;
                        double t = sum + y;

                        c = (t - sum) - y;
                        sum = t;
                    }

                    result = sum / Sum(w, ignore_na);
                }
            }
            else
            {
                result = Mean(x, ignore_na); //unweighted mean
            }
            return result;
        }

        public static double SumOfSquares(double[] x)
        {
            double mean;
            double ssd = SumOfSquaredDev(x, out mean);
            return ssd + mean * mean * x.Length;
        }

        public static double SumOfSquares(double[] x, double[] w)
        {
            double mean;
            double ssd = SumOfSquaredDev(x, w, out mean);

            if (w == null)
                ssd = ssd + mean * mean * x.Length;
            else
                ssd = ssd + mean * mean * Sum(w);

            return ssd;
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
