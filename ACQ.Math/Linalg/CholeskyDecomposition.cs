using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Linalg
{
    /// <summary>
    ///		Cholesky Decomposition of a symmetric, positive definite matrix. (Adapted from JAMA)
    ///	</summary>
    /// <remarks>
    ///		For a symmetric, positive definite matrix A, the Cholesky decomposition is a
    ///		lower triangular matrix L so that >A = L * L'.
    ///		If the matrix is not symmetric or positive definite, the constructor returns a partial 
    ///		decomposition and sets two internal variables that can be queried using the
    ///		 PositiveDefinite property.
    ///		
    /// TODO: consider adding tolerance to check for positivity, and pivoting 
    ///	</remarks>
    public class CholeskyDecomposition
    {
        private Matrix m_l;
        private bool m_positiveDefinite; // positive definite (only symmetric matrix can be positive definite)
        //it is possible to extend definition of "positive definite" to include some non-symmetric real matrices, but we dont need it here 
        
        /// <summary>
        /// Construct a Cholesky Decomposition
        /// </summary>
        /// <param name="A"></param>
        /// <param name="checkSymmetry">check if matrix is symmetric, if not use elements below diagonal and assume that it is symmetric</param>
        public CholeskyDecomposition(Matrix A, bool checkSymmetry = false)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }

            if (!A.IsSquare)
            {
                throw new ArgumentException("Matrix is not square.", "value");
            }

            int n = A.Rows;
            m_l = new Matrix(n, n);

            double[,] a = A.Data;
            double[,] l = m_l.Data;

            m_positiveDefinite = (A.Columns == n); //only symmetric can be positive definite
            
            for (int j = 0; j < n; j++)
            {
                double d = 0.0;
                for (int k = 0; k < j; k++)
                {   
                    double s = 0.0;
                    for (int i = 0; i < k; i++)
                    {
                        s += l[k, i] * l[j, i];
                    }
                    l[j, k] = s = (a[j, k] - s) / l[k, k];
                    d = d + s * s;

                    if (checkSymmetry) //Symmetry check, does not work very well with floating numbers because of truncation errors  
                    {
                        m_positiveDefinite = m_positiveDefinite && (a[k, j] == a[j, k]);
                    }
                }

                d = a[j, j] - d;

                //Positive-definite: d > 0,  Positive-semidefinite d >= 0                
                m_positiveDefinite = m_positiveDefinite & (d > 0.0); 

                l[j,j] = System.Math.Sqrt(System.Math.Max(d, 0.0));

                for (int k = j + 1; k < n; k++)
                {
                    l[j, k] = 0.0;
                }
            }
        }

        /// <summary>Returns true if the matrix is symmetric and positive definite.</summary>
        public bool IsPositiveDefinite
        {
            get
            {
                return m_positiveDefinite;
            }
        }

      
        /// <summary>Returns the left triangular factor L so that A = L * L'.</summary>
        public Matrix GetL()
        {
            return m_l.Clone();            
        }

        /// <summary>Solves a set of equation systems of type A * X = B.</summary>
        /// <param name="B">Right hand side matrix with as many rows as A and any number of columns.</param>
        /// <returns>Matrix X so that L * L' * X = B.</returns>
        public Matrix Solve(Matrix B)
        {
            if (B == null)
            {
                throw new ArgumentNullException("B");
            }

            if (B.Rows != m_l.Rows)
            {
                throw new ArgumentException("Matrix dimensions do not match.");
            }

            if (!IsPositiveDefinite)
            {
                throw new InvalidOperationException("Matrix is not positive definite.");
            }

            int n = m_l.Rows;
            int nx = B.Columns;

            Matrix X = (Matrix)B.Clone();
            double[,] x = X.Data;
            double[,] l = m_l.Data;

            // Solve L*Y = B;
            for (int k = 0; k < n; k++)
            {
                for (int j = 0; j < nx; j++)
                {
                    for (int i = 0; i < k; i++)
                    {
                        x[k, j] -= x[i, j] * l[k, i];
                    }
                    x[k, j] /= l[k, k];
                }
            }

            // Solve L'*X = Y;
            for (int k = n - 1; k >= 0; k--)
            {
                for (int j = 0; j < nx; j++)
                {
                    for (int i = k + 1; i < n; i++)
                    {
                        x[k, j] -= x[i, j] * l[i, k];
                    }
                    x[k, j] /= l[k, k];
                }
            }

            return X;
        }
    }

}
