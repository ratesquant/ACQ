using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Regression
{
    public interface IRegression
    {
        double Eval(params double[] x);
    }

    public interface IRegressionParam
    {
        double GetParam(string name);
    }

    public interface IRegressionSummary
    {
        Dictionary<string, double> Summary { get; }
    }
}
