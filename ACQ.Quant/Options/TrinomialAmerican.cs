using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{
    /// <summary>
    /// Implements original Cox, Ross, & Rubinstein (CRR) method for pricing american options on binomial tree (log-normal)
    /// https://en.wikipedia.org/wiki/Binomial_options_pricing_model
    /// </summary>
    public class TrinomialAmerican
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="dividend"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall, int time_steps = 1024)
        {
            double r = rate;
            double q = dividend;
            double S = spot;
            double K = strike;

            int n = time_steps; //number of time steps
            double dt = time / n;
            double df = Exp(-r * dt);
            double up = Exp(sigma * Sqrt(2 * dt));
            double up2 = up * up; //up/dn ratio
            double dn = 1d / up;
            double k_up = Sqrt(up);
            double k_dn = 1d/ k_up;
            double carry = Exp(0.5 * (r - q) * dt);

            double p_up = (k_up * carry - 1d) / (up - 1d);
            double p_dn = (up - k_up * carry) / (up - 1d);
            p_up = df * p_up * p_up;
            p_dn = df * p_dn * p_dn;
            double p_md = df - (p_up + p_dn);

            if (p_up < 0 || p_dn < 0 || p_md < 0)
            {
                //number of steps needs to be increased to keep dt below the following threshold,
                //we are not going to do this inside the function - since it is computationaly costly, 
                //dt < 2*sigma^2 /(r - q)^2

                return Double.NaN;
            }


            double[] v = new double[2 * n + 1]; //option values
            double[] p = new double[2 * n + 1]; //underlying asset prices

            p[0] = S * Pow(up, -n);
            for (int i = 1; i <= 2 * n; i++)
            {
                p[i] = up * p[i - 1]; //compute terminal distribution of prices, (there might be some error accumulation)
            }
            // options values at expiration
            for (int i = 0; i <= 2 * n; i++)
            {
                v[i] = isCall ? Max(0d, p[i] - K) : Max(0d, K - p[i]);
            }

            for (int j = 2 * n - 1; j >= n; j--)
            {
                double v_dn = v[2 * n - j - 1]; //dn price
                for (int i = 2 * n - j; i <= j; i++)
                {                    
                    double v_i = p_up * v[i + 1] + p_dn * v_dn + p_md * v[i];
                    
                    v_dn = v[i]; //save mid-price, it will be used as down price on next iteration

                    v[i] = isCall ? Max(v_i, p[i] - K) : Max(v_i, K - p[i]);
                }
            }

            return v[n];
        }

        public static double Greeks(enOptionGreeks greek, double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall, int time_steps)
        {
            double value = Double.NaN;

            if (greek == enOptionGreeks.Price)
            {
                value = Price(spot, strike, time, rate, dividend, sigma, isCall, time_steps);
            }
            else
            {
                Utils.OptionPriceDelegate price_function = delegate (double S, double K, double t, double r, double q, double v) {
                    return TrinomialAmerican.Price(S, K, t, r, q, v, isCall, time_steps);
                };

                value = Utils.NumericalGreeks(price_function, greek, spot, strike, time, rate, dividend, sigma);
            }

            return value;
        }
    }
}
