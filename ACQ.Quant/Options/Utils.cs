using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Quant.Options
{
    public enum enOptionGreeks
    {
        Price,
        /// <summary>
        /// dP/dF
        /// </summary>
        Delta,
        /// <summary>
        /// d(dP/dF)/dF
        /// </summary>
        Gamma,
        /// <summary>
        /// dP/ds
        /// </summary>
        Vega,
        /// <summary>
        /// d(dP/ds)/ds
        /// </summary>
        Vomma,
        Vanna,
        /// <summary>
        /// dP/dr
        /// </summary>
        Rho,
        /// <summary>
        /// -dP/dt
        /// </summary>
        Theta,
        /// <summary>
        /// -d(dP/dF)dt
        /// </summary>
        Charm,
        /// <summary>
        /// dP/d(dividend)
        /// </summary>
        Epsilon
    }

    public class Utils
    {
        public static double ImpliedVol(Func<double, double> opt_price, double target_price)
        {
            double min_sigma = 0;
            double max_sigma = 1;
            const double sigma_limit = 1e6;
            
            //Func<double, double> opt_price = x => Price(forward, strike, time, rate, x, isCall);

            if (opt_price(min_sigma) > target_price)
                return Double.NaN;

            //find right limit for volatility. 
            while (opt_price(max_sigma) < target_price && max_sigma < sigma_limit)
            {
                max_sigma *= 2;
            }

            if (opt_price(max_sigma) < target_price)
                return Double.NaN;

            var solver = new ACQ.Math.Roots.Brent();

            ACQ.Math.Roots.IterationResults results = solver.Solve(delegate (double x)
            {
                return opt_price(x) - target_price;
            }, min_sigma, max_sigma);

            return results.Root;
        }
    }

}
