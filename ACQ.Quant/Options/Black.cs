using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{
    public enum enOptionGreeks
    {
        Price,
        Delta,
        Gamma,
        Vega,
        Vomma,
        Vanna,
        Rho,
        Theta
    }
    /// <summary>
    /// Black option greeks, https://en.wikipedia.org/wiki/Greeks_(finance)
    /// Reference implementation
    /// </summary>
    public class Black
    {
        public static double Greeks(enOptionGreeks greek, double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double value = Double.NaN;
            
            //TODO: redo using reflection
            switch(greek)
            {
                case enOptionGreeks.Price:
                    value = Price(forward, strike, time, rate, sigma, isCall);
                    break;
                case enOptionGreeks.Delta:
                    value = Delta(forward, strike, time, rate, sigma, isCall);
                    break;
                case enOptionGreeks.Gamma:
                    value = Gamma(forward, strike, time, rate, sigma);
                    break;
                case enOptionGreeks.Vega:
                    value = Vega(forward, strike, time, rate, sigma);
                    break;
                case enOptionGreeks.Vomma:
                    value = Vomma(forward, strike, time, rate, sigma);
                    break;
                case enOptionGreeks.Vanna:
                    value = Vanna(forward, strike, time, rate, sigma);
                    break;
                case enOptionGreeks.Rho:
                    value = Rho(forward, strike, time, rate, sigma, isCall);
                    break;
                case enOptionGreeks.Theta:
                    value = Theta(forward, strike, time, rate, sigma, isCall);
                    break;
            }
            return value;
        }
        public static double Price(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double v = sigma * Sqrt(t);

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);

            double price = Double.NaN;

            if (isCall)
            {
                price = df * (F * Math.Special.NormalCdf(d1) - K * Math.Special.NormalCdf(d2));
            }
            else
            {
                price = -df * (F * Math.Special.NormalCdf(-d1) - K * Math.Special.NormalCdf(-d2));
            }
            return price;
        }

        public static double ImpliedVol(double forward, double strike, double time, double rate, double price, bool isCall)
        {
            double min_sigma = 0;
            double max_sigma = 1;
            const double sigma_limit = 1e6;

            Func<double, double> opt_price = x => Price(forward, strike, time, rate, x, isCall);

            if (opt_price(min_sigma) > price)
                return Double.NaN;

            //find right limit for volatility. 
            while (opt_price(max_sigma) < price && max_sigma < sigma_limit)
            {
                max_sigma *= 2;
            }

            if (opt_price(max_sigma) < price)
                return Double.NaN;

            var solver = new ACQ.Math.Roots.Brent();

            ACQ.Math.Roots.IterationResults results = solver.Solve(delegate (double x)
            {
                return opt_price(x) - price;
            }, min_sigma, max_sigma);

            return results.Root;
        }

        /// <summary>
        /// Delta - first oder derivative with respect to forward dP/dF 
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Delta(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double v = sigma * Sqrt(t);

            double d1 = (Log(F / K) + 0.5 * v * v) / v;            
            double df = Exp(-r * t);

            double delta = Double.NaN;

            if (isCall)
            {
                delta = df * Math.Special.NormalCdf(d1);
            }
            else
            {
                delta = -df * Math.Special.NormalCdf(-d1);
            }
            return delta;
        }

        /// <summary>
        /// Gamma  - second oder derivative with respect to forward d2P/dF^2, same for call and put 
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Gamma(double forward, double strike, double time, double rate, double sigma)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double v = sigma * Sqrt(t);

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double df = Exp(-r * t);

            double gamma = df * Math.Special.NormalPdf(d1) / (F * v);
            
            return gamma;
        }

        /// <summary>
        /// Vega  - derivative with respect to volatility (sigma) dP/dsigma, same for call and put 
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Vega(double forward, double strike, double time, double rate, double sigma)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double df = Exp(-r * t);

            double vega = df * F * Math.Special.NormalPdf(d1) * sqrt_t;

            return vega;
        }

        /// <summary>
        /// Theta  - derivative with respect to time
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Theta(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);

            double theta = Double.NaN;

            if (isCall)
            {
                theta = df * (-0.5 * F * Math.Special.NormalPdf(d1) * sigma/ sqrt_t - r * K * Math.Special.NormalCdf(d2) + r * F * Math.Special.NormalCdf(d1));
            }
            else
            {
                theta = df * (-0.5 * F * Math.Special.NormalPdf(d1) * sigma / sqrt_t + r * K * Math.Special.NormalCdf(-d2) - r * F * Math.Special.NormalCdf(-d1));
            }
            return theta;
        }

        public static double Rho(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);

            double rho = Double.NaN;

            if (isCall)
            {
                rho = -t * df * (F * Math.Special.NormalCdf(d1) - K * Math.Special.NormalCdf(d2));
            }
            else
            {
                rho = -t * df * (-F * Math.Special.NormalCdf(-d1) + K * Math.Special.NormalCdf(-d2));
            }
            return rho;
        }

        /// <summary>
        /// Vanna - second order derivative with respect to forward and volatility
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Vanna(double forward, double strike, double time, double rate, double sigma)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);

            double vanna = - df * Math.Special.NormalPdf(d1) * d2 / sigma;
            
            return vanna;
        }

        /// <summary>
        /// Vomma (Volga) - second order derivative with respect to the volatility
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Vomma(double forward, double strike, double time, double rate, double sigma)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(F / K) + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);

            double vomma =  df * F * Math.Special.NormalPdf(d1) * sqrt_t * d2 * d1 / sigma;

            return vomma;
        }
    }
}
