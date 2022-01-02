using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using static System.Math;

namespace ACQ.Quant.Options
{
    public class SABRParams    
    {
        double m_alpha; //initial vol
        double m_beta; //CEV exponent
        double m_rho; //correlation 
        double m_nu; //vol of vol

        public SABRParams(double alpha, double beta, double rho, double nu, bool check_params = false)
        {
            if (check_params)
            {
                if (!(alpha > 0.0))
                {
                    throw new ArgumentException("alpha > 0.0", "alpha");
                }

                if (!(beta >= 0.0 && beta <= 1.0))
                {
                    throw new ArgumentException("beta >= 0.0 && beta <= 1.0", "beta");
                }

                if (!(rho * rho < 1.0))
                {
                    throw new ArgumentException("-1 < rho < 1.0", "rho");
                }

                if (!(nu >= 0.0))
                {
                    throw new ArgumentException("nu >= 0.0", "nu");
                }
            }

            m_alpha = alpha;
            m_beta = beta;
            m_rho = rho;
            m_nu = nu;            
        }

        public bool isValid
        {
            get 
            {
                if (!(m_alpha > 0.0))
                {
                    return false;//throw new ArgumentException("alpha > 0.0", "alpha");
                }

                if (!(m_beta >= 0.0 && m_beta <= 1.0))
                {
                    return false; //throw new ArgumentException("beta >= 0.0 && beta <= 1.0", "beta");
                }

                if (!(m_rho * m_rho < 1.0))
                {
                    return false; //throw new ArgumentException("-1 < rho < 1.0", "rho");
                }

                if (!(m_nu >= 0.0))
                {
                    return false; //throw new ArgumentException("nu >= 0.0", "nu");
                }
                return true;
            }
            
        }

        public double Alpha
        {
            get
            {
                return m_alpha;
            }
        }
        public double Beta
        {
            get
            {
                return m_beta;
            }
        }
        public double Rho
        {
            get
            {
                return m_rho;
            }
        }
        public double Nu
        {
            get
            {
                return m_nu;
            }
        }
    }

    /// <summary>
    /// Hagan approximation for classical SABR model: "Managing Smile Risk" by PATRICK S. HAGAN, DEEP KUMAR, ANDREW S. LESNIEWSKI, AND DIANA E. WOODWARD
    /// </summary>
    public class HaganSABR
    {
        public static double Price(double forward, double strike, double time, double rate, SABRParams sabr_params, bool isCall, bool useNormalApprox = false)
        {
            double price = Double.NaN;

            if (sabr_params.isValid)
            {
                if (useNormalApprox)
                {
                    double normal_vol = HaganSABR.NormalVol(forward, strike, time, sabr_params);
                    price = Bachelier.Price(forward, strike, time, rate, normal_vol, isCall);
                }
                else
                {
                    double black_vol = HaganSABR.BlackVol(forward, strike, time, sabr_params);
                    price = Black.Price(forward, strike, time, rate, black_vol, isCall);
                }
            }

            return price;
        }

        public static double BlackVol(double forward, double strike, double time, SABRParams sabr_params)
        {
            if (!sabr_params.isValid)
            {
                return Double.NaN;
            }

            double alpha = sabr_params.Alpha;
            double beta = sabr_params.Beta;
            double rho = sabr_params.Rho;
            double nu = sabr_params.Nu;

            double K = strike;
            double F = forward;
            double t = time;            
            
            double A = Pow(F * K, 1.0 - beta);
            double sqrt_a = Sqrt(A);
            double log_r;

            if (Abs(F - K) > 1e-12)
            {
                log_r = Log(F / K);
            }
            else
            {
                double epsilon = (F - K) / K;
                log_r = epsilon - 0.5 * epsilon * epsilon;
            }
            double z = (nu / alpha) * sqrt_a * log_r;
            double B = 1.0 - 2.0 * rho * z + z * z;
            double C = (1.0 - beta) * (1.0 - beta) * log_r * log_r;
            double tmp = (Sqrt(B) + z - rho) / (1.0 - rho);
            double xx = Log(tmp);
            double D = sqrt_a * (1.0 + C / 24.0 + C * C / 1920.0);
            double d = 1.0 + t * ((1.0 - beta) * (1.0 - beta) * alpha * alpha / (24.0 * A) + 0.25 * rho * beta * nu * alpha / sqrt_a + (2.0 - 3.0 * rho * rho) * (nu * nu / 24.0));

            double coef;

            if (Abs(z * z) > 2.22045e-15)
            {
                coef = z / xx;
            }
            else
            {
                coef = 1.0 - 0.5 * rho * z - (3.0 * rho * rho - 2.0) * z * z / 12.0;
             }
            return ((alpha / D) * coef * d);
        }

        public static double NormalVol(double forward, double strike, double time, SABRParams sabr_params)
        {
            if (!sabr_params.isValid)
            {
                return Double.NaN;
            }

            double alpha = sabr_params.Alpha;
            double beta = sabr_params.Beta;
            double rho = sabr_params.Rho;
            double nu = sabr_params.Nu;

            double K = strike;
            double F = forward;
            double t = time;

            double A = Pow(F * K, 1.0 - beta);
            double sqrtA = Sqrt(A);
            double log_r;

            if (Abs(F - K) > 1e-12)
            {
                log_r = Log(F / K);
            }
            else
            {
                double epsilon = (F - K) / K;
                log_r = epsilon - 0.5 * epsilon * epsilon;
            }

            double z = (nu / alpha) * sqrtA * log_r;
            double B = 1.0 - 2.0 * rho * z + z * z;
            double C = (1.0 - beta) * (1.0 - beta) * log_r * log_r;
            double D = log_r * log_r;
            double tmp = (Sqrt(B) + z - rho) / (1.0 - rho);
            double xx = Log(tmp);
            double E_1 = (1.0 + D / 24.0 + D * D / 1920.0);
            double E_2 = (1.0 + C / 24.0 + C * C / 1920.0);
            double E = E_1 / E_2;
            double d = 1.0 + t * (-beta * (2 - beta) * alpha * alpha / (24.0 * A) +
                                0.25 * rho * beta * nu * alpha / sqrtA +
                                (2.0 - 3.0 * rho * rho) * (nu * nu / 24.0));
            double coef;

            if (Abs(z * z) > 2.22045e-15)
            {
                coef = z / xx;
            }
            else
            {
                coef = 1.0 - 0.5 * rho * z - (3.0 * rho * rho - 2.0) * z * z / 12.0;
            }

            return (alpha * Pow(F * K, beta / 2.0) * E * coef * d);
        }
    }
}
