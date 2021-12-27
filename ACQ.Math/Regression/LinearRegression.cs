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
                throw new ArgumentException("LinearRegression: check input arguments");
            }

            bool use_weights = (w != null);

            if (use_weights && w.Length != y.Length)            
            {
                throw new ArgumentException("LinearRegression: array of weights should have the same length");
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

            if (m_intercept)
            {               
                int m = 1;

                A = new Linalg.Matrix(n, m + 1);

                for (int i = 0; i < n; i++)
                {
                    A[i, 0] = x[0];
                    A[i, 1] = 1.0; //intercept
                }
            }
            else
            {
                A = new Linalg.Matrix(x, false);
            }  

            compute_coefficients(A, y, w);
        }

        public LinearRegression(double[,] x, double[] y, double[] w, bool intercept)
        {
            if (x == null || y == null || x.GetLength(0) != y.Length)
            {
                throw new ArgumentException("LinearRegression: check input arguments");
            }
            
            if ((w != null) && w.Length != y.Length)
            {
                throw new ArgumentException("LinearRegression: array of weights should have the same length");
            }

            m_intercept = intercept;

            //init variable names
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
                A = new Linalg.Matrix(x, false);
            }

            compute_coefficients(A, y, w);
        }

        private void compute_coefficients(Linalg.Matrix A, double[] y, double[] w)
        {
            Linalg.Matrix b = new Linalg.Matrix(y, false);

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


            double mss = m_intercept ? Math.Stats.Utils.SumOfSquaredDev(fit.RowPackedData(), w) : Math.Stats.Utils.SumOfSquares(fit.RowPackedData());
            double tss = m_intercept ? Math.Stats.Utils.SumOfSquaredDev(y, w) : Math.Stats.Utils.SumOfSquares(y);
            double rss = Math.Stats.Utils.SumOfSquares(residuals.RowPackedData()); //without intercept residuals dont add up to zero

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
            double aic = n + n * System.Math.Log(2 * System.Math.PI) + n * System.Math.Log(rss / n) + 2 * (m_coeffs.Length + 1); // sum(log(w))

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
