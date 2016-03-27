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
    ///		Symmetric and PositiveDefinite properties.
    ///	</remarks>
    public class CholeskyDecomposition
    {
        private Matrix m_l;
        private bool m_symPosDef; //symmetric and positive definite

        /// <summary>Construct a Cholesky Decomposition.</summary>
        public CholeskyDecomposition(Matrix A)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }

            if (!A.Square)
            {
                throw new ArgumentException("Matrix is not square.", "value");
            }

            int n = A.Rows;
            m_l = new Matrix(n, n);

            double[,] a = A.Data;
            double[,] l = m_l.Data;

            m_symPosDef = (A.Columns == n);

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

                    m_symPosDef = m_symPosDef && (a[k, j] == a[j, k]);
                }

                d = a[j, j] - d;

                m_symPosDef = m_symPosDef & (d > 0.0);
                l[j,j] = System.Math.Sqrt(System.Math.Max(d, 0.0));
                for (int k = j + 1; k < n; k++)
                {
                    l[j, k] = 0.0;
                }
            }
        }

        /// <summary>Returns true if the matrix is symmetric.</summary>
        public bool IsSymmetricAndPositiveDefinite
        {
            get
            {
                return m_symPosDef;
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

            if (!IsSymmetricAndPositiveDefinite)
            {
                throw new InvalidOperationException("Matrix is not symmetric and positive definite.");
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
