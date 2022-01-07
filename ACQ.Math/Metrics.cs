using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Math
{
    /// <summary>
    /// An implementation of evaluation metrics in R that are commonly
    /// used in supervised machine learning.It implements metrics for
    /// regression, time series, binary classification, classification,and information retrieval problems.
    /// </summary>
    public class Metrics
    {
        private static bool IsValidInput(double[] actual, double[] predicted)
        {
            if (actual == null || predicted == null)
                return false;

            if (actual.Length == 0 || actual.Length != predicted.Length)
                return false;

            return true;
        }
        /// <summary>
        ///  Mean Absolute Percent Error between two numeric vectors
        /// </summary>
        /// <param name="actual"></param>
        /// <param name="predicted"></param>
        /// <returns></returns>
        public static double MAPE(double[] actual, double[] predicted)
        {
            double result = Double.NaN;

            if (IsValidInput(actual, predicted))
            {
                result = 0;
                for (int i = 0; i < actual.Length; i++)
                {
                    result += Abs(actual[i] - predicted[i]) / actual[i];
                }
                result = result / actual.Length;
            }

            return result;
        }

        #region Binomial Metrics
        /// <summary>
        /// Computes AUC, actual should contain only 0 and 1
        /// The function only checks if actual is not equal to zero, assumes that all other values indicate that even happened (i.e. 1)
        /// </summary>
        /// <param name="actual"></param>
        /// <param name="predicted"></param>
        /// <returns></returns>
        public static double AUC(double[] actual, double[] predicted)
        {
            double result = Double.NaN;

            if (IsValidInput(actual, predicted))
            {
                int[] index = new int[actual.Length];
                int[] rank = new int[actual.Length];

                for (int i = 0; i < actual.Length; i++)
                {
                    index[i] = i;
                }
                
                Array.Sort(predicted, index);

                //ties are not handled
                for (int i = 0; i < index.Length; i++)
                {
                    rank[index[i]] = i;
                }

                int n_pos = 0;
                double sum = 0;
                for (int i = 0; i < actual.Length; i++)
                {
                    if (actual[i] != 0d)
                    {
                        n_pos++;
                        sum += rank[i];
                    }
                }
                int n_neg = actual.Length - n_pos;

                if (n_neg == 0)
                    result = 1.0;
                else
                    result = (sum/n_pos - 0.5*(n_pos - 1)) / n_neg;
            }

            return result;
        }

        /// <summary>
        /// Gini Coefficient is 2*AUC – 1, and its purpose is to normalize the AUC so that a random classifier scores 0, and a perfect classifier scores 1
        /// </summary>
        /// <param name="actual"></param>
        /// <param name="predicted"></param>
        /// <returns></returns>
        public static double Gini(double[] actual, double[] predicted)
        {
            double auc = AUC(actual, predicted);
            return 2 * auc - 1d;
        }
        
        public static double LogLoss(double[] actual, double[] predicted)
        {
            double result = Double.NaN;

            if (IsValidInput(actual, predicted))
            {
                result = 0;
                for (int i = 0; i < actual.Length; i++)
                {
                    result += actual[i] * Log(predicted[i]) + (1 - actual[i]) * Log(1 - predicted[i]);
                }
                result = -result / actual.Length;
            }

            return result;
        }

        /// <summary>
        /// Kolmogorov-Smirnov Statistics: [0, 1]
        /// </summary>
        /// <param name="actual"></param>
        /// <param name="predicted"></param>
        /// <returns></returns>
        public static double KS(double[] actual, double[] predicted)
        {
            double result = Double.NaN;

            if (IsValidInput(actual, predicted))
            {
                Array.Sort(predicted, actual);

                int m0 = 0;
                int m1 = 0;
                int max_diff = 0;
                int m0_max = 0, m1_max = 0;

                for (int i = 0; i < predicted.Length; i++)
                {
                    if (actual[i] != 0d)
                    {
                        m1++;
                    }
                    else
                    {
                        m0++;
                    }

                    int diff = Abs(m1 - m0);

                    if (diff > max_diff)
                    {
                        m0_max = m0;
                        m1_max = m1;
                        max_diff = diff;
                    }
                }

                if (m1 != 0 && m0 != 0)
                {
                    result = Abs( (double)m1_max / m1 - (double)m0_max / m0);
                }
            }

            return result;
        }
        #endregion 
    }
}
