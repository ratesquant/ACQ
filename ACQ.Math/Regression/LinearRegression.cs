using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Regression
{
    public class LinearRegression : IRegression, IRegressionSummary, IRegressionParam
    {
        readonly bool m_intercept;
        readonly string[] m_names;
        readonly bool m_weighted;

        //paramaters
        double[] m_coeffs;
        double[] m_coeffs_stderr;
        
        
        int m_observations;

        Dictionary<string, double> m_summary;   

        /// <summary>
        /// 1D constructor
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="w"></param>
        /// <param name="intercept"></param>
        public LinearRegression(double[] x, double[] y, double[] w, bool intercept)
        {
            if (x == null || y == null || x.Length != y.Length)            
            {
                throw new ArgumentException("LinearRegression: check input arguments (x and y can't be null and must have the same size)");
            }

            m_weighted = (w != null);

            if (m_weighted && w.Length != y.Length)            
            {
                throw new ArgumentException("LinearRegression: array of weights should have the same length as x and y");
            }

            m_intercept = intercept;

            //init variable names
            m_names = new string[1 + (m_intercept ? 1 : 0)];
            m_names[0] = "c1";

            if (m_intercept)
            {
                m_names[1] = "(intercept)";
            }

            Math.Linalg.Matrix A = null;

            int n = x.Length;

            //check input for NaN
            for (int i = 0; i < x.Length; i++)
            { 
                if(Double.IsNaN(x[i]) || Double.IsNaN(y[i]))
                    throw new ArgumentException("LinearRegression: there should not be NaN values in x or y");

                if(m_weighted && Double.IsNaN(w[i]))
                    throw new ArgumentException("LinearRegression: weights vector should not have NaN values");
            }

            if (m_intercept)
            {               
                int m = 1;

                A = new Linalg.Matrix(n, m + 1);

                for (int i = 0; i < n; i++)
                {
                    A[i, 0] = x[i];
                    A[i, 1] = 1.0; //intercept
                }
            }
            else
            {
                A = new Linalg.Matrix(x);
            }  

            compute_coefficients(A, y, w);
        }

        public LinearRegression(double[,] x, double[] y, double[] w, bool intercept)
        {
            if (x == null || y == null || x.GetLength(0) != y.Length)
            {
                throw new ArgumentException("LinearRegression: check input arguments (x and y can't be null and length of y should be the same as number of x rows)");
            }

            m_weighted = (w != null);

            if (m_weighted && w.Length != y.Length)
            {
                throw new ArgumentException("LinearRegression: array of weights should have the same length as y");
            }

            m_intercept = intercept;

            //init variable names
            int n = x.GetLength(0);
            int m = x.GetLength(1);            
            
            m_names = new string[m + (m_intercept ? 1 : 0)];

            for (int j = 0; j < m; j++)
            {
                m_names[j] = String.Format("c{0}", j + 1);
            }

            if (m_intercept)
            {
                m_names[m] = "(intercept)";
            }

            //check input for NaN
            for (int i = 0; i < n; i++)
            {
                if (Double.IsNaN(y[i]))
                    throw new ArgumentException("LinearRegression: there should not be NaN values in y");

                if (m_weighted && Double.IsNaN(w[i]))
                    throw new ArgumentException("LinearRegression: weights vector should not have NaN values");

                for (int j = 0; j < m; j++)
                {
                    if (Double.IsNaN(x[i, j]))
                        throw new ArgumentException("LinearRegression: there should not be NaN values in x");

                }
            }


            compute_regression(x, y, w);
        }        

        public double Estimate(params double[] x)
        {
            double result = Double.NaN;

            if(x != null)
            {
                if (x.Length + (m_intercept ? 1 : 0) == m_coeffs.Length)
                {
                    double sum = m_intercept ? m_coeffs[m_coeffs.Length - 1] : 0.0;

                    for (int i = 0; i < x.Length; i++)
                    {
                        sum += m_coeffs[i] * x[i];
                    }
                    result = sum;                    
                }
            }
            return result;
        }

        public Dictionary<string, double> Summary
        {
            get
            {
                return m_summary;
            }
        }

        public double GetParam(string name)
        {
            double value ;

            if (!m_summary.TryGetValue(name, out value))
            {
                value = Double.NaN;
            }
            return value;
        }

        private void compute_regression(double[,] x, double[] y, double[] w)
        {
            int m = x.GetLength(1);

            Math.Linalg.Matrix A = null;

            int n = x.GetLength(0);

            if (m_intercept)
            {
                A = new Linalg.Matrix(n, m + 1);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < m; j++)
                    {
                        A[i, j] = x[i, j];
                    }
                    A[i, m] = 1.0;
                }
            }
            else
            {
                A = new Linalg.Matrix(x);
            }

            compute_coefficients(A, y, w);
        }

        private void compute_coefficients(Linalg.Matrix A, double[] y, double[] w)
        {
            Linalg.Matrix b = new Linalg.Matrix(y, false);

            double[] w_sqrt = null;
            if (m_weighted)
            {
                w_sqrt = new double[b.Rows];
                for (int i = 0; i < b.Rows; i++)
                {
                    w_sqrt[i] = System.Math.Sqrt(w[i]);

                    b[i, 0] = b[i, 0] * w_sqrt[i]; //here we assume that specified weights are sigmas (i.e. not sigma square)

                    for (int j = 0; j < A.Columns; j++)
                    {
                       A[i, j] = A[i, j] * w_sqrt[i];
                    }
                }
            }

            Linalg.Matrix coeffs;
            Linalg.Matrix cov;
            Linalg.SvDecomposition sv;
            if (A.Rows > A.Columns)
            {
                sv = new Linalg.SvDecomposition(A);
                coeffs = sv.Solve(b, false);
            }
            else
            {
                sv = new Linalg.SvDecomposition(A.Transpose());
                coeffs = sv.Solve(b, true);
            }
            cov = sv.Cov();

            Linalg.Matrix fit = A * coeffs;
            Linalg.Matrix residuals = fit - b;

            double aic_weight_correction = 0;
            double weight_sum = fit.Rows;
            if (m_weighted)
            {
                for (int i = 0; i < fit.Rows; i++)
                {
                    fit[i, 0] = fit[i, 0] / w_sqrt[i];
                    residuals[i, 0] = fit[i, 0] - b[i, 0] / w_sqrt[i];

                    aic_weight_correction += System.Math.Log(w[i]);
                }
                weight_sum = Math.Stats.Utils.Sum(w);
            }

            double mss = Math.Stats.Utils.SumOfSquaredDev(fit.RowPackedData(), w);
            double rss = Math.Stats.Utils.SumOfSquares(residuals.RowPackedData(), w); //without intercept residuals dont add up to zero            
            double tss = Math.Stats.Utils.SumOfSquaredDev(y, w);


            m_observations = A.Rows;
            m_coeffs = coeffs.RowPackedData();

            int p = m_coeffs.Length;
            int n = m_observations;
            int dof = n - p;
            int dof_intercept = m_intercept ? 1 : 0;

            m_coeffs_stderr = new double[p]; //coefficient standard errors

            double stderr = System.Math.Sqrt( rss / dof ); //residual sum of squares

            for (int i = 0; i < p; i++)
            {
                m_coeffs_stderr[i] = stderr * System.Math.Sqrt(cov[i, i]);
            }

            double r2 = mss / (mss + rss);
            double r2_adj = 1.0 - (n - dof_intercept) * (1.0 - r2) / dof;
            double fvalue = (mss / (p - dof_intercept)) / (rss / dof);
            double aic_bic = n + n * System.Math.Log(2 * System.Math.PI) + n * System.Math.Log(rss / n) -aic_weight_correction;
            double aic = aic_bic + 2 * (m_coeffs.Length + 1);
            double bic = aic_bic + System.Math.Log(m_observations) * (m_coeffs.Length + 1);

            //generate summary
            m_summary = new Dictionary<string, double>();

            m_summary["n"] = m_observations;
            m_summary["p"] = m_coeffs.Length;
            m_summary["dof"] = dof;
            m_summary["stderr"] = stderr;
            m_summary["rss"] = rss;
            m_summary["tss"] = tss;
            m_summary["mss"] = mss;
            m_summary["r2"] = r2;
            m_summary["r2_adj"] = r2_adj;
            m_summary["fvalue"] = fvalue;
            m_summary["aic"] = aic;
            m_summary["bic"] = bic;

            for (int i = 0; i < p; i++)
            {
                double coef = m_coeffs[i];
                double coef_std = m_coeffs_stderr[i];
                double tvalue = coef / coef_std;
                double pvalue = Math.Special.incbet(0.5 * dof, 0.5, dof / (dof + tvalue * tvalue));  

                m_summary[m_names[i]] = coef;
                m_summary[m_names[i] + "_se"] = coef_std;
                m_summary[m_names[i] + "_tvalue"] = tvalue;
                m_summary[m_names[i] + "_pvalue"] = pvalue;
            }

            m_summary["intercept"] = m_intercept ? 1.0 : 0.0;
        }
    }
}
