using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Quant.Options
{

    /// <summary>
    /// Option greeks (derivatives): here U - underlying, s - volatility, r - risk free rate,
    /// </summary>
    public enum enOptionGreeks
    {
        Price,
        /// <summary>
        /// dP/dU
        /// </summary>
        Delta,
        /// <summary>
        /// d(dP/dU)/dU = d(Delta)/dU
        /// </summary>
        Gamma,
        /// <summary>
        /// dP/ds
        /// </summary>
        Vega,
        /// <summary>
        /// Volga, d(dP/ds)/ds = d(Vega)/ds
        /// </summary>
        Vomma,
        /// <summary>
        /// d(dP/ds)/dU = d(Vega)/dU
        /// </summary>
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
        /// -d(dP/dU)dt = d(Delta) / dt
        /// </summary>
        Charm,
        /// <summary>
        /// dP/d(dividend)
        /// </summary>
        Epsilon
    }

    public class Utils
    {
        public delegate double OptionPriceDelegate(double underlying, double strike, double time, double rate, double dividend, double sigma);
        public static double ImpliedVol(Func<double, double> opt_price, double target_price)
        {
            double min_sigma = 1e-15;
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

        
        public static double NumericalGreeks(OptionPriceDelegate price, enOptionGreeks greek, double underlying, double strike, double time, double rate, double dividend, double sigma)
        {
            const double dx = 1e-3; //shifts to compute numeric greeks
            
            double value = Double.NaN;
            double q = dividend;
            double S = underlying;
            double K = strike;
            double t = time;
            double r = rate;

            switch (greek)
            {
                case enOptionGreeks.Price:
                    value = price(S, K, time, r, q, sigma);
                    break;
                case enOptionGreeks.Delta:
                    value = (price(S + dx, K, t, r, q, sigma) - price(S - dx, K, t, r, q, sigma)) / (2 * dx);
                    break;
                case enOptionGreeks.Gamma:
                    value = (price(S + dx, K, t, r, q, sigma) + price(S - dx, K, t, r, q, sigma) - 2*price(S, K, t, r, q, sigma)) / (dx * dx);
                    break;
                case enOptionGreeks.Vega:
                    double sigma_up = sigma + dx;
                    double sigma_dn = System.Math.Max(1e-12, sigma - dx);
                    value = (price(S, K, t, r, q, sigma_up) - price(S, K, t, r, q, sigma_dn)) / (sigma_up - sigma_dn);
                    break;
                case enOptionGreeks.Vomma:
                    value = (price(S, K, t, r, q, sigma + dx) + price(S, K, t, r, q, sigma - dx) - 2 * price(S, K, t, r, q, sigma)) / (dx * dx); ;
                    break;
                case enOptionGreeks.Rho:
                    value = (price(S, K, t, r + dx, q, sigma) - price(S, K, t, r - dx, q, sigma)) / (2 * dx);
                    break;
                case enOptionGreeks.Theta:
                    double dt = System.Math.Min(dx, t);
                    value = (price(S, K, t - dt, r, q, sigma) - price(S, K, t, r, q, sigma)) / (dt);
                    break;
                case enOptionGreeks.Vanna:                    
                    value = (price(S+dx, K, t, r, q, sigma+dx) - price(S+dx, K, t, r, q, sigma-dx) - price(S - dx, K, t, r, q, sigma + dx) + price(S - dx, K, t, r, q, sigma - dx)) / (4 * dx * dx);
                    break;

            }
            return value;
        }
    }

}
