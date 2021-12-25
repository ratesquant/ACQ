using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{

    /// <summary>
    /// Implements original Cox, Ross, & Rubinstein (CRR) method for pricing options on binomial tree (log-normal)
    /// https://en.wikipedia.org/wiki/Binomial_options_pricing_model
    /// </summary>
    public class Binomial
    {
        public delegate double OptionPayoff(double price, double strike);
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
        public static double Price(double spot, double strike, double time, double rate, double dividend, double sigma, OptionPayoff payoff, bool isAmerican = true, int time_steps = 1024)
        {
            double r = rate;
            double q = dividend;
            double S = spot;
            double K = strike;

            int n = time_steps; //number of time steps
            double dt = time / n;
            double df = Exp(-r * dt);
            double up = Exp(sigma * Sqrt(dt));
            double up2 = up * up; //up/dn ratio
            double dn = 1d / up;
            double p_up = df * (up * Exp((r - q) * dt) - 1d) / (up2 - 1d);
            //double p_up = (up * Exp(- q* dt) - df) / (up2 - 1d);
            double p_dn = df - p_up;

            if (p_up < 0 || p_dn < 0)
            {
                //number of steps needs to be increased to keep dt below the following threshold,
                //we are not going to do this inside the function - since it is computationaly costly, 
                //dt < sigma^2 /(r - q)^2

                return Double.NaN;
            }         


            double[] v = new double[n + 1]; //option values
            double[] p = new double[n + 1]; //underlying asset prices

            p[0] = S * Pow(up, -n);
            for (int i = 1; i <= n; i++)
            {
                p[i] = up2 * p[i - 1]; //compute terminal distribution of prices, (there might be some error accumulation)
            }
            // options values at expiration
            for (int i = 0; i <= n; i++)
            {
                v[i] = payoff(p[i], K);
            }

            for (int j = n - 1; j >= 0; j--)
            {
                for (int i = 0; i <= j; i++)
                {
                    v[i] = p_up * v[i + 1] + p_dn * v[i];
                    p[i] = dn * p[i + 1];

                    if (isAmerican)
                    {
                        v[i] = Max(v[i], payoff(p[i], K));
                    }
                }
            }

            return v[0];
        }

        public static double Greeks(enOptionGreeks greek, double spot, double strike, double time, double rate, double dividend, double sigma, OptionPayoff payoff, bool isAmerican, int time_steps)
        {
            double value = Double.NaN;

            if (greek == enOptionGreeks.Price)
            {
                value = Price(spot, strike, time, rate, dividend, sigma, payoff, isAmerican, time_steps);
            }
            else
            {
                Utils.OptionPriceDelegate price_function = delegate (double S, double K, double t, double r, double q, double v) {
                    return Binomial.Price(S, K, t, r, q, v, payoff, isAmerican, time_steps);
                };

                value = Utils.NumericalGreeks(price_function, greek, spot, strike, time, rate, dividend, sigma);
            }

            return value;
        }
    }
}
