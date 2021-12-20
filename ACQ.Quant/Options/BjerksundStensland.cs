using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{
    /// <summary>
    /// Bjerksund and Stensland (2002) Approximation for American options
    /// Adopted from "The Complete Guide to Option Pricing Formulas", Espen Gaarder Haug 
    /// "Numerical Methods versus Bjerksund and Stensland Approximations for American Options Pricing" (2014) by Marasovic Branka, Aljinovic Zdravka, Poklepovic Tea
    /// </summary>
    public class BjerksundStensland
    {
        /// <summary>
        /// Compute finite difference greeks 
        /// </summary>
        /// <param name="greek"></param>
        /// <param name="spot"></param>
        /// <param name="strike"></param>
        /// <param name="time"></param>
        /// <param name="rate"></param>
        /// <param name="dividend"></param>
        /// <param name="sigma"></param>
        /// <param name="isCall"></param>
        /// <returns></returns>
        public static double Greeks(enOptionGreeks greek, double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double value = Double.NaN;

            if (greek == enOptionGreeks.Price)
            {
                value = Price(spot, strike, time, rate, dividend, sigma, isCall);
            } else
            { 
                Utils.OptionPriceDelegate price_function = delegate (double S, double K, double t, double r, double q, double v) {
                    return BjerksundStensland.Price(S, K, t, r, q, v, isCall); };

                value = Utils.NumericalGreeks(price_function, greek, spot, strike, time, rate, dividend, sigma);
            }

            return value;
        }
        public static double Price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double price;
            if (isCall)
                price = PriceCall(spot, strike, time, rate, dividend, sigma);
            else
                price = PriceCall(strike,spot, time, rate- dividend, -dividend, sigma);
            return price;
        }


        private static double PriceCall(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double K = strike;
            double S = spot;
            double t = time;
            double r = rate;
            double b = dividend;
            double v = sigma;
            double v2 = v * v;
            

            double price;

            double BInfinity, B0;
            double ht1, ht2, I1, I2;
            double alfa1, alfa2, Beta, t1;


            t1 = 0.5 * (Sqrt(5) - 1) * t;

            if (b >= r)// Never optimal to exercise before maturity
                price = ACQ.Quant.Options.BlackScholes.Price(S, K, t, r, b, sigma, true);
            else
            {
                Beta = (0.5 - b / v2) + Sqrt((b / v2 - 0.5) * (b / v2 - 0.5) + 2 * r / v2);
                BInfinity = Beta / (Beta - 1) * K;
                B0 = Max(K, r / (r - b) * K);

                ht1 = -(b * t1 + 2 * v * Sqrt(t1)) * K * K / ((BInfinity - B0) * B0);
                ht2 = -(b * t + 2 * v * Sqrt(t)) * K * K / ((BInfinity - B0) * B0);
                I1 = B0 + (BInfinity - B0) * (1 - Exp(ht1));
                I2 = B0 + (BInfinity - B0) * (1 - Exp(ht2));
                alfa1 = (I1 - K) * Pow(I1, -Beta);
                alfa2 = (I2 - K) * Pow(I2, -Beta);

                if (S >= I2)
                    price = S - K;
                else
                {
                    price = alfa2 * Pow(S, Beta) - 
                        alfa2 * phi(S, t1, Beta, I2, I2, r, b, v) + 
                        phi(S, t1, 1, I2, I2, r, b, v) - 
                        phi(S, t1, 1, I1, I2, r, b, v) - 
                        K * phi(S, t1, 0, I2, I2, r, b, v) + 
                        K * phi(S, t1, 0, I1, I2, r, b, v) + 
                        alfa1 * phi(S, t1, Beta, I1, I2, r, b, v) -
                        alfa1 * ksi(S, t, Beta, I1, I2, I1, t1, r, b, v) + 
                        ksi(S, t, 1, I1, I2, I1, t1, r, b, v) - 
                        ksi(S, t, 1, K, I2, I1, t1, r, b, v) - 
                        K * ksi(S, t, 0, I1, I2, I1, t1, r, b, v) + 
                        K * ksi(S, t, 0, K, I2, I1, t1, r, b, v);
                }
            }
          
            return price;
        }

        private static double phi(double S, double t, double gamma, double h, double i, double r, double b, double v)
        {
            double lambda, kappa, d;
            double v2 = v * v;
            double v_sqrt_t = v * Sqrt(t);
            double phi;

            lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1) * v2) * t;
            d = -(Log(S / h) + (b + (gamma - 0.5) * v2) * t) / v_sqrt_t;
            kappa = 2 * b / v2 + 2 * gamma - 1;
            phi = Exp(lambda) * Pow(S, gamma) * (Math.Special.NormalCdf(d) - Pow(i / S, kappa) * Math.Special.NormalCdf(d - 2 * Log(i / S) / v_sqrt_t));

            return phi;
        }

        private static double ksi(double S, double t2, double gamma, double h, double I2, double I1, double t1, double r, double b, double v)
        {
            double e1, e2, e3, e4;
            double f1, f2, f3, f4;
            double rho, kappa, lambda;
            double v_sqrt_t1 = v * Sqrt(t1);
            double v_sqrt_t2 = v * Sqrt(t2);
            double v2 = v * v;

            e1 = (Log(S / I1) + (b + (gamma - 0.5) * v2) * t1) / v_sqrt_t1;
            e2 = (Log(I2* I2 / (S * I1)) + (b + (gamma - 0.5) * v2) * t1) / v_sqrt_t1;
            e3 = (Log(S / I1) - (b + (gamma - 0.5) * v2) * t1) / v_sqrt_t1;
            e4 = (Log(I2 * I2 / (S * I1)) - (b + (gamma - 0.5) * v2) * t1) / v_sqrt_t1;


            f1 = (Log(S / h) + (b + (gamma - 0.5) * v2) * t2) / v_sqrt_t2;
            f2 = (Log(I2 * I2 / (S * h)) + (b + (gamma - 0.5) * v2) * t2) / v_sqrt_t2;
            f3 = (Log(I1 * I1 / (S * h)) + (b + (gamma - 0.5) * v2) * t2) / v_sqrt_t2;
            f4 = (Log(S * I1 * I1 / (h * I2 * I2)) + (b + (gamma - 0.5) * v2) * t2) / v_sqrt_t2;

            rho = Sqrt(t1 / t2);
            lambda = -r + gamma * b + 0.5 * gamma * (gamma - 1) * v2;
            kappa = 2 * b / (v2) + (2 * gamma - 1);

            double ksi = Exp(lambda * t2) * Pow(S, gamma) * (Math.Special.CBND(-e1, -f1, rho) - Pow(I2 / S, kappa) * Math.Special.CBND(-e2, -f2, rho) -
                Pow(I1 / S, kappa) * Math.Special.CBND(-e3, -f3, -rho) + Pow(I1 / I2, kappa) * Math.Special.CBND(-e4, -f4, -rho));
         
            return ksi;
        }


        public static double ImpliedVol(double forward, double strike, double time, double rate, double dividend, double option_price, bool isCall)
        {
            Func<double, double> opt_price = x => Price(forward, strike, time, rate, dividend, x, isCall);

            double implied_vol = Utils.ImpliedVol(opt_price, option_price);

            return implied_vol;
        }
    }
}
