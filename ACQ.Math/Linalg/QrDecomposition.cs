using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Linalg
{
    /// <summary>
    ///	  QR decomposition for a rectangular matrix. (Adapted from JAMA)
    /// </summary>
    /// <remarks>
    ///   For an m x n matrix A with m >= n, the QR decomposition is an m x n
    ///   orthogonal matrix Q and an n x n upper triangular 
    ///   matrix R so that A = Q * R.
    ///   The QR decompostion always exists, even if the matrix does not have
    ///   full rank, so the constructor will never fail. The primary use of the
    ///   QR decomposition is in the least squares solution of nonsquare systems
    ///   of simultaneous linear equations.
    ///   This will fail if FullRank returns false.
    /// </remarks>
    public class QrDecomposition
    {
        private Matrix m_qr;
        private double[] m_rdiag;

        /// <summary>Construct a QR decomposition.</summary>	
        public QrDecomposition(Matrix A)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }

            m_qr = (Matrix)A.Clone();
            double[,] qr = m_qr.Data;
            int m = A.Rows;
            int n = A.Columns;
            m_rdiag = new double[n];

            for (int k = 0; k < n; k++)
            {
                // Compute 2-norm of k-th column without under/overflow.
                double nrm = 0;
                for (int i = k; i < m; i++)
                {
                    nrm = Utils.Hypotenuse(nrm, qr[i, k]);
                }

                if (nrm != 0.0)
                {
                    // Form k-th Householder vector.
                    if (qr[k, k] < 0)
                    {
                        nrm = -nrm;
                    }

                    for (int i = k; i < m; i++)
                    {
                        qr[i, k] /= nrm;
                    }

                    qr[k, k] += 1.0;

                    // Apply transformation to remaining columns.
                    for (int j = k + 1; j < n; j++)
                    {
                        double s = 0.0;

                        for (int i = k; i < m; i++)
                        {
                            s += qr[i, k] * qr[i, j];
                        }

                        s = -s / qr[k, k];

                        for (int i = k; i < m; i++)
                        {
                            qr[i, j] += s * qr[i, k];
                        }
                    }
                }

                m_rdiag[k] = -nrm;
            }
        }

        /// <summary>Least squares solution of A * X = B</summary>
        /// <param name="B">Right-hand-side matrix with as many rows as A and any number of columns.</param>
        /// <returns>A matrix that minimized the two norm of Q * R * X - B.</returns>
        public Matrix Solve(Matrix B)
        {
            if (B == null)
            {
                throw new ArgumentNullException("value");
            }

            if (B.Rows != m_qr.Rows)
            {
                throw new ArgumentException("Matrix row dimensions must agree.");
            }

            if (!IsFullRank)
            {
                throw new InvalidOperationException("Matrix is rank deficient.");
            }

            // Copy right hand side
            int count = B.Columns;
            Matrix X = B.Clone();
            int m = m_qr.Rows;
            int n = m_qr.Columns;
            double[,] qr = m_qr.Data;

            // Compute Y = transpose(Q)*B
            for (int k = 0; k < n; k++)
            {
                for (int j = 0; j < count; j++)
                {
                    double s = 0.0;

                    for (int i = k; i < m; i++)
                    {
                        s += qr[i, k] * X[i, j];
                    }

                    s = -s / qr[k, k];

                    for (int i = k; i < m; i++)
                    {
                        X[i, j] += s * qr[i, k];
                    }
                }
            }

            // Solve R*X = Y;
            for (int k = n - 1; k >= 0; k--)
            {
                for (int j = 0; j < count; j++)
                {
                    X[k, j] /= m_rdiag[k];
                }

                for (int i = 0; i < k; i++)
                {
                    for (int j = 0; j < count; j++)
                    {
                        X[i, j] -= X[k, j] * qr[i, k];
                    }
                }
            }

            return X.Submatrix(0, n - 1, 0, count - 1);
        }

        /// <summary>Is the matrix full rank?</summary>
        public bool IsFullRank
        {
            get
            {
                int columns = m_qr.Columns;
                for (int i = 0; i < columns; i++)
                {
                    if (m_rdiag[i] == 0)
                    {
                        return false;
                    }
                }

                return true;
            }
        }

        /// <summary>Returns R</summary>
        public Matrix GetH()
        {
            int m = m_qr.Rows;
            int n = m_qr.Columns;
            Matrix H = new Matrix(m, n);
            double[,] h = H.Data;
            double[,] qr = m_qr.Data;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i >= j)
                    {
                        h[i, j] = qr[i, j];
                    }
                    else
                    {
                        h[i, j] = 0.0;
                    }
                }
            }
            return H;            
        }
        /// <summary>Returns R</summary>
        public Matrix GetR()
        {
            int n = this.m_qr.Columns;
            Matrix R = new Matrix(n, n);
            double[,] r = R.Data;
            double[,] qr = m_qr.Data;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i < j)
                    {
                        r[i, j] = qr[i, j];
                    }
                    else if (i == j)
                    {
                        r[i, j] = m_rdiag[i];
                    }
                    else
                    {
                        r[i, j] = 0.0;
                    }
                }
            }
            return R;            
        }

        /// <summary>Returns the orthogonal factor <c>Q</c>.</summary>
        public Matrix GetQ()
        {
            int m = m_qr.Rows;
            int n = m_qr.Columns;
            Matrix Q = new Matrix(m, n);
            double[,] q = Q.Data;
            double[,] qr = m_qr.Data;
            for (int k = n - 1; k >= 0; k--)
            {
                for (int i = 0; i < m; i++)
                {
                    q[i, k] = 0.0;
                }

                q[k, k] = 1.0;
                for (int j = k; j < n; j++)
                {
                    if (qr[k, k] != 0)
                    {
                        double s = 0.0;

                        for (int i = k; i < m; i++)
                        {
                            s += qr[i, k] * q[i, j];
                        }

                        s = -s / qr[k, k];

                        for (int i = k; i < m; i++)
                        {
                            q[i, j] += s * qr[i, k];
                        }
                    }
                }
            }

            return Q;
        }
    }
}
