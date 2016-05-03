using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Regression
{
    public interface IRegression
    {
        double Estimate(double[] x);
    }

    enum enRegressionResults
    {
 
    }

    public abstract class RegressionBase
    {
    }
}
