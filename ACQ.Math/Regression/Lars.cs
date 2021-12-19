using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Regression
{
    /// <summary>
    /// Least Angle Regression (with optional LASSO solution) 
    /// 
    /// Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). Least angle regression. The Annals of statistics, 32(2), 407-499.
    /// </summary>
    public class Lars
    {
        List<double[]> m_beta; //coefficients 
        List<double> m_gamma; //step sizes 
 
        /// <summary>
        /// First Lars model, (includes intercept) 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        public Lars(double[,] x, double[] y)
        {
            bool lasso = false;            

            int n = x.GetLength(0);  //number of observations 
            int p = x.GetLength(1); //number of variables 

            int maxvars = System.Math.Min(n - 1, p); //maximum number of variables 
            int maxit = 8 * maxvars;

            m_beta = new List<double[]>(); // this are coefficients on each step 
            m_gamma = new List<double>();
            
            SortedSet<int> c_set = new SortedSet<int>(); //candidate set
            SortedSet<int> a_set = new SortedSet<int>(); //active set

            //initialize candidate set with all available variables
            for (int i = 0; i < p; i++)
            {
                c_set.Add(i);
            }

            //compute gramm matrix (gram = x' * x) 
            Linalg.Matrix full_gram = new Linalg.Matrix(p, p);             

            for (int i = 0; i < p; i++)
            {
                for (int j = 0; j < p; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < n; k++)
                    {
                        sum += x[k, i] * x[k, j];
                    }
                    full_gram[i, j] = sum;
                }
            }

            //bool stop_flag = false;

            double[] mu = new double[n]; //lars regression vector
            double[] c = new double[p]; //correlations 

            for (int it = 0; it < maxit; it++)
            {
                //compute correlations 
                for (int j = 0; j < p; j++)
                {
                    double sum = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        sum += x[i, j] * (y[i] - mu[i]);
                    }
                    c[j] = sum;
                }

                //find abs max corr from candidate set
                double max_abs_c = 0.0;
                int max_abs_c_index = -1;

                foreach (int i in c_set)
                {
                    double abs_c = System.Math.Abs(c[i]);

                    if (abs_c > max_abs_c)
                    {
                        max_abs_c = abs_c;
                        max_abs_c_index = i;
                    }
                }

                //exit if there is no correlation
                if (max_abs_c < ACQ.Math.Const.epsilon)
                {
                    break;
                }

                a_set.Add(max_abs_c_index);
                c_set.Remove(max_abs_c_index);

                int vars = a_set.Count;

                double[] s = new double[vars];

                foreach (int i in a_set)
                {
                    s[i] = ACQ.Math.Utils.Sign(c[i]);
                }

                //compute partical Gram matrix, Gram = X(active_columns)' * X(active_columns)
                int[] active_indices = a_set.ToArray();
                Linalg.Matrix gram = full_gram.Submatrix(active_indices, active_indices);
                Linalg.CholeskyDecomposition gram_chol = new Linalg.CholeskyDecomposition(gram);
                Linalg.Matrix inv_gram = gram_chol.Solve(s);

                //compute coefficients of equiangular vector
                double[] w = new double[vars];

                double norm = 0.0;

                for (int i = 0; i < vars; i++)
                {
                    w[i] = s[i] * inv_gram[i, 0];
                    norm += w[i];
                }

                double scale = 1.0 / System.Math.Sqrt(norm);

                for (int i = 0; i < vars; i++)
                {
                    w[i] = scale * w[i];
                }
                //compute equiangular vector
                double[] u = new double[n];

                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < vars; j++)
                    {
                        sum += x[i, active_indices[j]] * w[j];
                    }
                    u[i] = sum;
                }

                double gamma = max_abs_c / scale; // set gamma to the largest value (i.e. use regular least squares)    

                //correlation (angle) between equiangular vector and all remaining variables
                foreach (int i in c_set)
                {
                    double angle = 0.0;
                    for (int j = 0; j < n; j++)
                    {
                        angle += x[j, i] * u[j];
                    }

                    double t1 = (max_abs_c - c[i]) / (scale - angle);
                    double t2 = (max_abs_c + c[i]) / (scale + angle);

                    if (t1 > 0)
                    {
                        gamma = System.Math.Min(t1, gamma);
                    }

                    if (t2 > 0)
                    {
                        gamma = System.Math.Min(t2, gamma);
                    }
                }

                //LASSO code here
                if (lasso)
                {
 
                }

                //update coefficients 
                double[] beta = new double[p];

                if (m_beta.Count > 0)
                {
                    double[] pev_beta = m_beta[m_beta.Count - 1];

                    for (int i = 0; i < vars; i++)
                    {
                        int index = active_indices[i];
                        beta[index] = pev_beta[index];
                    }
                }

                for (int i = 0; i < vars; i++)
                {
                    beta[active_indices[i]] += gamma * w[i];
                }

                m_beta.Add(beta);
                m_gamma.Add(gamma);


                //update lars vector
                for (int i = 0; i < n; i++)
                {
                    mu[i] += gamma * u[i];
                }
            }
        }
    }
}
