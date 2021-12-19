using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{
    /// <summary>
    /// Bachelier option pricing (same as Black but with Normal distribution instead of Log-Normal)
    /// Kazuhiro Iwasawa, "Analytic Formula for the European Normal Black Scholes Formula". New York University. 2 December 2001.    
    /// </summary>
    public class Bachelier
    {

        public static double Greeks(enOptionGreeks greek, double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double value = Double.NaN;
            
            switch (greek)
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
                case enOptionGreeks.Theta:
                    value = Theta(forward, strike, time, rate, sigma, isCall);
                    break;
                case enOptionGreeks.Rho:
                    value = Rho(forward, strike, time, rate, sigma, isCall);
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

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double price;

            if (isCall)
            {
                price = df * ((F - K) * Math.Special.NormalCdf(d) + v * Math.Special.NormalPdf(d));
            }
            else
            {
                price = df *((K - F) * Math.Special.NormalCdf(-d) + v * Math.Special.NormalPdf(-d));
            }
            return price;
        }

        public static double ImpliedVol(double forward, double strike, double time, double rate, double option_price, bool isCall)
        {
            Func<double, double> opt_price = x => Price(forward, strike, time, rate, x, isCall);

            double implied_vol = Utils.ImpliedVol(opt_price, option_price);

            return implied_vol;
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

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double delta;

            if (isCall)
            {
                delta = df * Math.Special.NormalCdf(d);
            }
            else
            {
                delta = -df * Math.Special.NormalCdf(-d);
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

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double gamma = df * Math.Special.NormalPdf(d) / v;

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

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double vega = df * Math.Special.NormalPdf(d) * sqrt_t;

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

            double d = (F - K) / v;            
            double df = Exp(-r * t);

            double theta;

            if (isCall)
            {
                theta = r * df *((F-K)* Math.Special.NormalCdf(d) + v * Math.Special.NormalPdf(d)) - 0.5 * df * Math.Special.NormalPdf(d) * sigma / sqrt_t;
            }
            else
            {
                theta = r * df * ((K - F) * Math.Special.NormalCdf(-d) + v * Math.Special.NormalPdf(d)) - 0.5 * df * Math.Special.NormalPdf(d) * sigma / sqrt_t;
            }
            return theta;
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

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double vomma = df * sqrt_t * Math.Special.NormalPdf(d) * d * d / sigma;

            return vomma;
        }

        /// <summary>
        /// Rho  - derivative with respect to the risk free rate
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Rho(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double K = strike;
            double F = forward;
            double t = time;
            double r = rate;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double rho;

            if (isCall)
            {
                rho = -t * df * ((F - K) * Math.Special.NormalCdf(d) + v * Math.Special.NormalPdf(d));
            }
            else
            {
                rho = -t * df * ((K - F) * Math.Special.NormalCdf(-d) + v * Math.Special.NormalPdf(-d));
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

            double d = (F - K) / v;
            double df = Exp(-r * t);

            double vanna = -df * d * Math.Special.NormalPdf(d) / sigma;

            return vanna;
        }

    }
}
