//#define EXPORT_GREEKS_FUNCTIONS
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

using ACQ.Quant;

namespace ACQ.Excel.Finance
{
    public static class Options
    {
        #region Black Option Greeks

        [ExcelFunction(Description = "Black option price", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_price(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double price = ACQ.Quant.Options.Black.Price(forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(price);
        }

        [ExcelFunction(Description = "Black option implied vol", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_vol(double forward, double strike, double time, double rate, double option_price, bool isCall)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                double vol = ACQ.Quant.Options.Black.ImpliedVol(forward, strike, time, rate, option_price, isCall);
                return ExcelHelper.CheckNan(vol);
            }
        }

        [ExcelFunction(Description = "Black option greeks", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_greeks(
            [ExcelArgument(Description = "Name of the option greek: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho,Theta")] string greek_name, 
            double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            var greek = ExcelHelper.CheckEnum<ACQ.Quant.Options.enOptionGreeks>(greek_name, ACQ.Quant.Options.enOptionGreeks.Price);

            double value = ACQ.Quant.Options.Black.Greeks(greek, forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(value);
        }
        
        #if EXPORT_GREEKS_FUNCTIONS
        [ExcelFunction(Description = "Black option delta (dP/dF)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_delta(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double delta = ACQ.Quant.Options.Black.Delta(forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(delta);
        }

        [ExcelFunction(Description = "Black option gamma d(dP/dF)/dF", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_gamma(double forward, double strike, double time, double rate, double sigma)
        {
            double gamma = ACQ.Quant.Options.Black.Gamma(forward, strike, time, rate, sigma);
            return ExcelHelper.CheckNan(gamma);
        }

        [ExcelFunction(Description = "Black option vega dP/dsigma", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_vega(double forward, double strike, double time, double rate, double sigma)
        {
            double vega = ACQ.Quant.Options.Black.Vega(forward, strike, time, rate, sigma);
            return ExcelHelper.CheckNan(vega);
        }
        [ExcelFunction(Description = "Black option vanna d(dP/dsigma)/dF", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_vanna(double forward, double strike, double time, double rate, double sigma)
        {
            double vanna = ACQ.Quant.Options.Black.Vanna(forward, strike, time, rate, sigma);
            return ExcelHelper.CheckNan(vanna);
        }
        [ExcelFunction(Description = "Black option vomma d(dP/dsigma)/dsigma", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_vomma(double forward, double strike, double time, double rate, double sigma)
        {
            double vamma = ACQ.Quant.Options.Black.Vomma(forward, strike, time, rate, sigma);
            return ExcelHelper.CheckNan(vamma);
        }

        [ExcelFunction(Description = "Black option theta  dP/dt", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_theta(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double theta = ACQ.Quant.Options.Black.Theta(forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(theta);            
        }

        [ExcelFunction(Description = "Black option rho  dP/dr", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_black_rho(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double rho = ACQ.Quant.Options.Black.Rho(forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(rho);
        }
        #endif
        #endregion
        
        #region Black-Scholes Option Greeks        

        [ExcelFunction(Description = "Black-Scholes option price", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double price = ACQ.Quant.Options.BlackScholes.Price(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(price);
        }

        [ExcelFunction(Description = "Black-Scholes option implied vol", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_vol(double spot, double strike, double time, double rate, double dividend, double option_price, bool isCall)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                double vol = ACQ.Quant.Options.BlackScholes.ImpliedVol(spot, strike, time, rate, dividend, option_price, isCall);
                return ExcelHelper.CheckNan(vol);
            }
        }

        [ExcelFunction(Description = "Black-Scholes option greeks", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_greeks(
        [ExcelArgument(Description = "Name of the option greek: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho,Theta")] string greek_name,
        double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            var greek = ExcelHelper.CheckEnum<ACQ.Quant.Options.enOptionGreeks>(greek_name, ACQ.Quant.Options.enOptionGreeks.Price);

            double value = ACQ.Quant.Options.BlackScholes.Greeks(greek, spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(value);
        }
        #if EXPORT_GREEKS_FUNCTIONS
        [ExcelFunction(Description = "Black-Scholes option delta (dP/dS)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_delta(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double delta = ACQ.Quant.Options.BlackScholes.Delta(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(delta);
        }

        [ExcelFunction(Description = "Black-Scholes option gamma d(dP/dS)/dS", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_gamma(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double gamma = ACQ.Quant.Options.BlackScholes.Gamma(spot, strike, time, rate, dividend, sigma);
            return ExcelHelper.CheckNan(gamma);
        }

        [ExcelFunction(Description = "Black-Scholes option vega dP/dsigma", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_vega(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double vega = ACQ.Quant.Options.BlackScholes.Vega(spot, strike, time, rate, dividend, sigma);
            return ExcelHelper.CheckNan(vega);
        }

        [ExcelFunction(Description = "Black-Scholes option Vomma d(dP/dsigma)/dsigma", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_vomma(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double vega = ACQ.Quant.Options.BlackScholes.Vomma(spot, strike, time, rate, dividend, sigma);
            return ExcelHelper.CheckNan(vega);
        }

        [ExcelFunction(Description = "Black-Scholes option charm d(dP/dS)/dt", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_charm(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double charm = ACQ.Quant.Options.BlackScholes.Charm(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(charm);
        }

        [ExcelFunction(Description = "Black-Scholes option epsilon dP/d(dividend)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_epsilon(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double epsilon = ACQ.Quant.Options.BlackScholes.Epsilon(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(epsilon);
        }

        [ExcelFunction(Description = "Black-Scholes option rho  dP/dr", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_rho(double spot, double strike, double time, double rate, double dividend,double sigma, bool isCall)
        {
            double rho = ACQ.Quant.Options.BlackScholes.Rho(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(rho);
        }

        [ExcelFunction(Description = "Black-Scholes option theta  dP/dt", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_theta(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double theta = ACQ.Quant.Options.BlackScholes.Theta(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(theta);
        }

        [ExcelFunction(Description = "Black-Scholes option vanna d(dP/dsigma)/dS", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_blackscholes_vanna(double spot, double strike, double time, double rate, double dividend, double sigma)
        {
            double vanna = ACQ.Quant.Options.BlackScholes.Vanna(spot, strike, time, rate, dividend, sigma);
            return ExcelHelper.CheckNan(vanna);
        }
        #endif
        
        #endregion
        
        #region Bachelier Option Greeks

        [ExcelFunction(Description = "Bachelier option price", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_bachelier_price(double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            double price = ACQ.Quant.Options.Bachelier.Price(forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(price);
        }

        [ExcelFunction(Description = "Bachelier option implied vol", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_bachelier_vol(double forward, double strike, double time, double rate, double option_price, bool isCall)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                double vol = ACQ.Quant.Options.Bachelier.ImpliedVol(forward, strike, time, rate, option_price, isCall);
                return ExcelHelper.CheckNan(vol);
            }
        }

        [ExcelFunction(Description = "Bachelier option greeks", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_bachelier_greeks(
            [ExcelArgument(Description = "Name of the option greek: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho,Theta")] string greek_name,
            double forward, double strike, double time, double rate, double sigma, bool isCall)
        {
            var greek = ExcelHelper.CheckEnum<ACQ.Quant.Options.enOptionGreeks>(greek_name, ACQ.Quant.Options.enOptionGreeks.Price);

            double value = ACQ.Quant.Options.Bachelier.Greeks(greek, forward, strike, time, rate, sigma, isCall);
            return ExcelHelper.CheckNan(value);
        }
        
        #endregion
        
        #region Bjerksund and Stensland (2002) Approximation for American options
        [ExcelFunction(Description = "Bjerksund and Stensland (2002) Approximation for American options", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_bjerksund_price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            double price = ACQ.Quant.Options.BjerksundStensland.Price(spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(price);
        }

        [ExcelFunction(Description = "Numerical Greeks using Bjerksund and Stensland (2002) Approximation for American options", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_bjerksund_greeks(
           [ExcelArgument(Description = "Name of the option greek: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho,Theta")] string greek_name,
           double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall)
        {
            var greek = ExcelHelper.CheckEnum<ACQ.Quant.Options.enOptionGreeks>(greek_name, ACQ.Quant.Options.enOptionGreeks.Price);

            double value = ACQ.Quant.Options.BjerksundStensland.Greeks(greek, spot, strike, time, rate, dividend, sigma, isCall);
            return ExcelHelper.CheckNan(value);
        }
        #endregion

        #region Binomial Price for American options
        [ExcelFunction(Description = "Binomial Price for American options", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_binomial_american_price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall, object time_steps)
        {
            //System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
            //timer.Restart();

            int n_steps = (int)ExcelHelper.CheckValue<double>(time_steps, 500);
            double price = ACQ.Quant.Options.BinomialAmerican.Price(spot, strike, time, rate, dividend, sigma, isCall, n_steps);

            //LogDisplay.WriteLine("Binomial American Price: {0}, {1}", n_steps, timer.Elapsed.TotalMilliseconds);            

            return ExcelHelper.CheckNan(price);
        }

        [ExcelFunction(Description = "Numerical Greeks using Binomial Price for American options", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_binomial_american_greeks(
            [ExcelArgument(Description = "Name of the option greek: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho,Theta")] string greek_name,
            double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall, object time_steps)
        {
            var greek = ExcelHelper.CheckEnum<ACQ.Quant.Options.enOptionGreeks>(greek_name, ACQ.Quant.Options.enOptionGreeks.Price);
            int n_steps = ExcelHelper.CheckValue<int>(time_steps, 500);

            double value = ACQ.Quant.Options.BinomialAmerican.Greeks(greek, spot, strike, time, rate, dividend, sigma, isCall, n_steps);
            return ExcelHelper.CheckNan(value);
        }
        #endregion

        #region Trinomial Price for American options
        [ExcelFunction(Description = "Trinomial Price for American options", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_trinomial_american_price(double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall, object time_steps)
        {
            int n_steps = (int)ExcelHelper.CheckValue<double>(time_steps, 500);
            double price = ACQ.Quant.Options.TrinomialAmerican.Price(spot, strike, time, rate, dividend, sigma, isCall, n_steps);
            
            return ExcelHelper.CheckNan(price);
        }

        [ExcelFunction(Description = "Numerical Greeks using Trinomial Price for American options", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_options_trinomial_american_greeks(
            [ExcelArgument(Description = "Name of the option greek: Price, Delta, Gamma, Vega, Vomma, Vanna, Rho,Theta")] string greek_name,
            double spot, double strike, double time, double rate, double dividend, double sigma, bool isCall, object time_steps)
        {
            var greek = ExcelHelper.CheckEnum<ACQ.Quant.Options.enOptionGreeks>(greek_name, ACQ.Quant.Options.enOptionGreeks.Price);
            int n_steps = ExcelHelper.CheckValue<int>(time_steps, 500);

            double value = ACQ.Quant.Options.TrinomialAmerican.Greeks(greek, spot, strike, time, rate, dividend, sigma, isCall, n_steps);
            return ExcelHelper.CheckNan(value);
        }
        #endregion
                
         //System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
         //timer.Restart();
         //LogDisplay.WriteLine("Elapsed: {0}, {1}", n_steps, timer.Elapsed.TotalMilliseconds);                    
    }
}

 