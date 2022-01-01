using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Stats
{
    /// <summary>
    /// Binomial confidence intervals
    /// https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    /// </summary>
    public class BinomTest
    {
        /// <summary>
        /// Confidence interval for the binomial probability.
        /// </summary>
        /// <param name="x">Vector of number of successes in the binomial experiment</param>
        /// <param name="n">Vector of number of independent trials in the binomial experiment.</param>
        /// <param name="conf_level"></param>
        /// <returns></returns>
        public static Tuple<double, double> ConfidenceInterval(int x, int n, double conf_level = 0.95)
        {
            double p = (double) x / n;
            double alpha = 1 - conf_level;
            double z = ACQ.Math.Special.InverseNormalCdf(1 - 0.5 * alpha);
            double z2 = z * z;

            double p1 = p + 0.5 * z2 / n;
            double p2 = z * System.Math.Sqrt((p * (1 - p) + 0.25 * z2 / n) / n);
            double p3 = 1 + z2 / n;
            double lcl = (p1 - p2) / p3;
            double ucl = (p1 + p2) / p3;

            return new Tuple<double, double>(lcl, ucl);
        }
    }
}
