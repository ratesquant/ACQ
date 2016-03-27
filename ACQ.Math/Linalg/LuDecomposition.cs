using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Linalg
{
    /// <summary>
    ///   LU decomposition of a rectangular matrix. (Adapted from JAMA)
    /// </summary>
    /// <remarks>
    ///   For an m x n matrix A with m >= n, the LU decomposition is an m x n
    ///   unit lower triangular matrix L, an n x n upper triangular matrix U,
    ///   and a permutation vector piv of length m so that A(piv)=L*U.
    ///   If m < n, then L is m x m and U is m x n.
    ///   The LU decompostion with pivoting always exists, even if the matrix is
    ///   singular, so the constructor will never fail.  The primary use of the
    ///   LU decomposition is in the solution of square systems of simultaneous
    ///   linear equations. This will fail if NonSingular returns false.
    /// </remarks>
    public class LuDecomposition
    {
        private Matrix m_lu;
        private int m_psign;
        private int[] m_pivot;

        /// <summary>Construct a LU decomposition.</summary>	
        public LuDecomposition(Matrix A)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }
            //"left-looking", dot-product, Crout/Doolittle algorithm.

            this.m_lu = (Matrix)A.Clone();
            double[,] lu = m_lu.Data;
            int rows = A.Rows;
            int columns = A.Columns;
            m_pivot = new int[rows];
            for (int i = 0; i < rows; i++)
            {
                m_pivot[i] = i;
            }

            m_psign = 1;
            double[] LUcolj = new double[rows];

            // Outer loop.
            for (int j = 0; j < columns; j++)
            {
                // Make a copy of the j-th column to localize references.
                for (int i = 0; i < rows; i++)
                {
                    LUcolj[i] = lu[i, j];
                }

                // Apply previous transformations.
                for (int i = 0; i < rows; i++)
                {
                    // Most of the time is spent in the following dot product.
                    int kmax = System.Math.Min(i, j);
                    double s = 0.0;
                    for (int k = 0; k < kmax; k++)
                    {
                        s += lu[i, k] * LUcolj[k];
                    }
                    lu[i, j] = LUcolj[i] -= s;
                }

                // Find pivot and exchange if necessary.
                int p = j;
                for (int i = j + 1; i < rows; i++)
                {
                    if (System.Math.Abs(LUcolj[i]) > System.Math.Abs(LUcolj[p]))
                    {
                        p = i;
                    }
                }

                if (p != j)
                {
                    for (int k = 0; k < columns; k++)
                    {
                        double t = lu[p, k];
                        lu[p, k] = lu[j, k];
                        lu[j, k] = t;
                    }

                    int v = m_pivot[p];
                    m_pivot[p] = m_pivot[j];
                    m_pivot[j] = v;

                    m_psign = -m_psign;
                }

                // Compute multipliers.

                if (j < rows & lu[j, j] != 0.0)
                {
                    for (int i = j + 1; i < rows; i++)
                    {
                        lu[i, j] /= lu[j, j];
                    }
                }
            }
        }

        /// <summary>Returns if the matrix is non-singular.</summary>
        public bool IsNonSingular
        {
            get
            {
                for (int j = 0; j < m_lu.Columns; j++)
                {
                    if (m_lu[j, j] == 0)
                    {
                        return false;
                    }
                }
                return true;
            }
        }

        public bool IsSingular
        {
            get
            {
                return !IsNonSingular;
            }
        }

        /// <summary>Returns the determinant of the matrix.</summary>
        public double Determinant
        {
            get
            {
                if (m_lu.Rows != m_lu.Columns)
                {
                    throw new ArgumentException("Matrix must be square.");
                }
                double det = (double)m_psign;
                for (int j = 0; j < m_lu.Columns; j++)
                {
                    det *= m_lu[j, j];
                }
                return det;
            }
        }

        /// <summary>Returns the lower triangular factor L with A=LU.</summary>
        public Matrix GetL()
        {
            int rows = m_lu.Rows;
            int columns = m_lu.Columns;
            Matrix L = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (i > j)
                        L[i, j] = m_lu[i, j];
                    else if (i == j)
                        L[i, j] = 1.0;
                    else
                        L[i, j] = 0.0;
                }
            }
            return L;
        }

        /// <summary>Return upper triangular factor U</summary>
        public Matrix GetU()
        {
            int rows = m_lu.Rows;
            int columns = m_lu.Columns;
            Matrix U = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (i <= j)
                        U[i, j] = m_lu[i, j];
                    else
                        U[i, j] = 0.0;
                }
            }
            return U;
           
        }

        /// <summary>Returns the pivot permuation vector.</summary>
        public double[] GetPivot()
        {
            int rows = m_lu.Rows;

            double[] p = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                p[i] = (double)m_pivot[i];
            }

            return p;
        }

        /// <summary>Solves a set of equation systems of type A * X = B.</summary>
        /// <param name="B">Right hand side matrix with as many rows as A and any number of columns.</param>
        /// <returns>Matrix X so that L * U * X = B</returns>
        public Matrix Solve(Matrix B)
        {
            if (B == null)
            {
                throw new ArgumentNullException("B");
            }

            if (B.Rows != m_lu.Rows)
            {
                throw new ArgumentException("Invalid matrix dimensions.", "value");
            }

            if (IsSingular)
            {
                throw new InvalidOperationException("Matrix is singular");
            }

            // Copy right hand side with pivoting
            int count = B.Columns;
            Matrix X = B.Submatrix(m_pivot, 0, count - 1);

            int rows = m_lu.Rows;
            int columns = m_lu.Columns;
            double[,] lu = m_lu.Data;

            // Solve L*Y = B(piv,:)
            for (int k = 0; k < columns; k++)
            {
                for (int i = k + 1; i < columns; i++)
                {
                    for (int j = 0; j < count; j++)
                    {
                        X[i, j] -= X[k, j] * lu[i, k];
                    }
                }
            }

            // Solve U*X = Y;
            for (int k = columns - 1; k >= 0; k--)
            {
                for (int j = 0; j < count; j++)
                {
                    X[k, j] /= lu[k, k];
                }

                for (int i = 0; i < k; i++)
                {
                    for (int j = 0; j < count; j++)
                    {
                        X[i, j] -= X[k, j] * lu[i, k];
                    }
                }
            }

            return X;
        }
    }
}
