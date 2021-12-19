using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{
    /// <summary>
    /// Black-Scholes option greeks, https://en.wikipedia.org/wiki/Greeks_(finance)
    /// Reference implementation
    /// </summary>
    public class BlackScholes
    {
        
        public static double Greeks(enOptionGreeks greek, double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double value = Double.NaN;

            //TODO: redo using reflection
            switch (greek)
            {
                case enOptionGreeks.Price:
                    value = Price(spot, strike, time, rate, dividend, sigma, isCall);
                    break;
                case enOptionGreeks.Delta:
                    value = Delta(spot, strike, time, rate, dividend, sigma, isCall);
                    break;
                case enOptionGreeks.Gamma:
                    value = Gamma(spot, strike, time, rate, dividend, sigma);
                    break;
                case enOptionGreeks.Vega:
                    value = Vega(spot, strike, time, rate, dividend, sigma);
                    break;
                case enOptionGreeks.Vomma:
                    value = Vomma(spot, strike, time, rate, dividend, sigma);
                    break;
                case enOptionGreeks.Vanna:
                    value = Vanna(spot, strike, time, rate, dividend, sigma);
                    break;
                case enOptionGreeks.Rho:
                    value = Rho(spot, strike, time, rate, dividend, sigma, isCall);
                    break;
                case enOptionGreeks.Theta:
                    value = Theta(spot, strike, time, rate, dividend, sigma, isCall);
                    break;
                case enOptionGreeks.Charm:
                    value = Charm(spot, strike, time, rate, dividend, sigma, isCall);
                    break;
                case enOptionGreeks.Epsilon:
                    value = Epsilon(spot, strike, time, rate, dividend, sigma, isCall);
                    break;
            }
            return value;
        }
        
        public static double Price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double v = sigma * Sqrt(t);

            double d1 = (Log(S / K) + (r - q)  * t + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);
            double dfq = Exp(-q * t);

            double price;

            if (isCall)
            {
                price = S * dfq * Math.Special.NormalCdf(d1) - K * df * Math.Special.NormalCdf(d2);
            }
            else
            {
                price = -S * dfq * Math.Special.NormalCdf(-d1) + K * df * Math.Special.NormalCdf(-d2);
            }
            return price;
        }

        public static double ImpliedVol(double spot, double strike, double time, double rate, double dividend, double option_price, bool isCall)
        {
            Func<double, double> opt_price = x => Price(spot, strike, time, rate, dividend, x, isCall);

            double implied_vol = Utils.ImpliedVol(opt_price, option_price);

            return implied_vol;
        }

        public static double Delta(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double v = sigma * Sqrt(t);

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double dfq = Exp(-q * t);

            double delta;

            if (isCall)
            {
                delta = dfq * Math.Special.NormalCdf(d1);
            }
            else
            {
                delta = -dfq * Math.Special.NormalCdf(-d1);
            }
            return delta;
        }

        /// <summary>
        /// Gamma  - second oder derivative with respect to forward d2P/dS^2, same for call and put 
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Gamma(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double v = sigma * Sqrt(t);

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;            
            double dfq = Exp(-q * t);

            double gamma = dfq * Math.Special.NormalPdf(d1) / (S * v);

            return gamma;
        }

        /// <summary>
        /// Vega  - derivative with respect to volatility (sigma) dP/dsigma, same for call and put 
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Vega(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double dfq = Exp(-q * t);

            double vega = dfq * S * Math.Special.NormalPdf(d1) * sqrt_t;

            return vega;
        }

        /// <summary>
        /// Vomma (Volga) - second order derivative with respect to the volatility
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="dividend"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Vomma(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double d2 = d1 - v;            
            double dfq = Exp(-q * t);

            double vomma = dfq * S * Math.Special.NormalPdf(d1) * sqrt_t * d2 * d1 / sigma;

            return vomma;
        }

        /// <summary>
        /// Vanna - second order derivative with respect to spot and volatility
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public static double Vanna(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double d2 = d1 - v;            
            double dfq = Exp(-q * t);

            double vanna = -dfq * Math.Special.NormalPdf(d1) * d2 / sigma;

            return vanna;
        }


        /// <summary>
        /// Theta  - derivative with respect to time
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Theta(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);
            double dfq = Exp(-q * t);

            double theta;

            if (isCall)
            {
                theta = -0.5 * dfq * S * Math.Special.NormalPdf(d1) * sigma / sqrt_t - r * K * df * Math.Special.NormalCdf(d2) + q * S * dfq * Math.Special.NormalCdf(d1);
            }
            else
            {
                theta = -0.5 * dfq * S * Math.Special.NormalPdf(d1) * sigma / sqrt_t + r * K * df * Math.Special.NormalCdf(-d2) - q * S * dfq * Math.Special.NormalCdf(-d1);
            }
            return theta;
        }

        public static double Rho(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double df = Exp(-r * t);
            double dfq = Exp(-q * t);

            double rho;

            if (isCall)
            {
                rho = t * df * K * Math.Special.NormalCdf(d2);
            }
            else
            {
                rho = -t * df * K * Math.Special.NormalCdf(-d2);
            }
            return rho;
        }

        /// <summary>
        /// Second order derivative with respect to underlying and time
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="dividend"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Charm(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double d2 = d1 - v;
            double dfq = Exp(-q * t);

            double charm;

            if (isCall)
            {
                charm = q * dfq * Math.Special.NormalCdf(d1) - 0.5 * dfq * Math.Special.NormalPdf(d1) * (2 * (r - q) * t - d2 * v) / (t * v);
            }
            else
            {
                charm = -q * dfq * Math.Special.NormalCdf(-d1) - 0.5*dfq * Math.Special.NormalPdf(d1) * (2 * (r - q) * t - d2 * v) / (t * v);
            }
            return charm;
        }

        /// <summary>
        /// Derivative with respect to the dividend yield
        /// </summary>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="dividend"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Epsilon(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double q = dividend;
            double sqrt_t = Sqrt(t);
            double v = sigma * sqrt_t;

            double d1 = (Log(S / K) + (r - q) * t + 0.5 * v * v) / v;
            double dfq = Exp(-q * t);

            double epsilon;

            if (isCall)
            {
                epsilon = -S * t * dfq * Math.Special.NormalCdf(d1);
            }
            else
            {
                epsilon = S * t * dfq * Math.Special.NormalCdf(-d1);
            }
            return epsilon;
        }

    }
}
