using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Linalg
{
    /// <summary>
    /// Singular Value Decomposition. (Adapted from JAMA and LAPACK)
    ///  
    /// For an m-by-n matrix A with m >= n, the singular value decomposition is
    /// an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
    /// an n-by-n orthogonal matrix V so that A = U*S*V'.
    /// 
    /// The singular values, sigma[k] = S[k][k], are ordered so that
    /// sigma[0] >= sigma[1] >= ... >= sigma[n-1].
    /// 
    /// The singular value decompostion always exists, so the constructor will
    /// never fail (prodived that m >= n). The matrix condition number and the effective numerical
    /// rank can be computed from this decomposition.
    /// 
    ///             A    =    U    *     S      *      Vt
    /// general   m x n     m x m      m x n         n x n
    /// m>=n      m x n     m x n      n x n         n x n (this is what is being computed in this class )
    /// compact   m x n     m x k      k x k         k x n , k = min(m, n) 
    /// 
    /// in case when m less than n, transpose the matrix before computing SVD
    /// </summary>
    public class SvDecomposition
    {
        private Matrix m_U;
        private Matrix m_V;
        private double[] m_s;
        
        public SvDecomposition(Matrix A)
        {            
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }
            
            // Initialize.
            int m = A.Rows;
            int n = A.Columns;  
            
            if (m < n)
            {
                throw new ArgumentException("SVD only works for m >= n, (transpose input matrix)");
            }

            double[,] a = (double[,])A.Data.Clone();
            
            m_s = new double[n]; //n
            m_U = new Matrix(m, n); // m x n
            m_V = new Matrix(n, n);
            double[] e = new double [n];
            double[] work = new double [m];
            bool wantu = true;
            bool wantv = true;
            
            double[,] u = m_U.Data;
            double[,] v = m_V.Data;
            
            // Reduce A to bidiagonal form, storing the diagonal elements
            // in s and the super-diagonal elements in e.
            
            int nct = System.Math.Min(m - 1,n);
            int nrt = System.Math.Max(0, System.Math.Min(n - 2, m));
            for (int k = 0; k < System.Math.Max(nct, nrt); k++) 
            {
                if (k < nct) 
                {
                    // Compute the transformation for the k-th column and
                    // place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    m_s[k] = 0;
                    for (int i = k; i < m; i++) 
                    {
                        m_s[k] = Math.Utils.Hypotenuse(m_s[k],a[i, k]);
                    }
                    if (m_s[k] != 0.0) 
                    {
                        if (a[k, k] < 0.0) 
                        {
                            m_s[k] = -m_s[k];
                        }
                        for (int i = k; i < m; i++) 
                        {
                            a[i, k] /= m_s[k];               
                        }
                        a[k, k] += 1.0;
                    }
                    m_s[k] = -m_s[k];
                }

                for (int j = k+1; j < n; j++) 
                {
                    if ((k < nct) & (m_s[k] != 0.0))  
                    {
                        // Apply the transformation
                        double t = 0;
                        for (int i = k; i < m; i++) 
                        {
                            t += a[i, k] * a[i, j];
                        }
                        t = -t / a[k, k];
                        for (int i = k; i < m; i++) 
                        {
                            a[i, j] += t * a[i, k];
                        }
                    }
                    // Place the k-th row of A into e for the
                    // subsequent calculation of the row transformation.

                    e[j] = a[k, j];
                }
                if (wantu & (k < nct)) {

                    // Place the transformation in U for subsequent back
                    // multiplication.

                    for (int i = k; i < m; i++) 
                    {
                        u[i, k] = a[i, k];
                    }
                }
                if (k < nrt) 
                {
                    // Compute the k-th row transformation and place the
                    // k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0;
                    for (int i = k+1; i < n; i++)
                    {
                        e[k] = Math.Utils.Hypotenuse(e[k],e[i]);
                    }
                    if (e[k] != 0.0) 
                    {
                        if (e[k+1] < 0.0) 
                        {
                            e[k] = -e[k];
                        }
                        for (int i = k+1; i < n; i++) 
                        {
                            e[i] /= e[k];
                        }
                        e[k+1] += 1.0;
                    }
                    e[k] = -e[k];
                    if ((k+1 < m) & (e[k] != 0.0)) 
                    {
                        // Apply the transformation.

                        for (int i = k+1; i < m; i++) 
                        {
                            work[i] = 0.0;
                        }
                        for (int j = k+1; j < n; j++) 
                        {
                            for (int i = k+1; i < m; i++) 
                            {
                                work[i] += e[j]*a[i, j];
                            }
                        }
                        for (int j = k+1; j < n; j++) 
                        {
                            double t = -e[j]/e[k+1];
                            for (int i = k+1; i < m; i++) 
                            {
                                a[i, j] += t*work[i];
                            }
                        }
                    }
                    if (wantv) 
                    {
                        // Place the transformation in V for subsequent
                        // back multiplication.
                        for (int i = k+1; i < n; i++)                 
                        {
                            v[i, k] = e[i];
                        }
                    }
                }
            }
            // Set up the final bidiagonal matrix or order p.
            int p = System.Math.Min(n,m+1);
            if (nct < n) 
            {
                m_s[nct] = a[nct, nct];
            }
            if (m < p) 
            {
                m_s[p-1] = 0.0;
            }
            if (nrt+1 < p) 
            {
                e[nrt] = a[nrt, p-1];
            }
            e[p-1] = 0.0;
            
            // If required, generate U.
            if (wantu) 
            {
                for (int j = nct; j < n; j++) 
                {
                    for (int i = 0; i < m; i++) 
                    {
                        u[i, j] = 0.0;
                    }
                    u[j, j] = 1.0;
                }
                for (int k = nct-1; k >= 0; k--) 
                {
                    if (m_s[k] != 0.0) 
                    {
                        for (int j = k+1; j < n; j++) 
                        {
                            double t = 0;
                            for (int i = k; i < m; i++) 
                            {
                                t += u[i, k] * u[i, j];
                            }
                            t = -t/u[k, k];
                            for (int i = k; i < m; i++) 
                            {
                                u[i, j] += t * u[i, k];
                            }
                        }
                        for (int i = k; i < m; i++ ) 
                        {
                            u[i, k] = -u[i, k];
                        }
                        u[k, k] = 1.0 + u[k, k];
                        for (int i = 0; i < k-1; i++) 
                        {
                            u[i, k] = 0.0;
                        }
                    } else 
                    {
                        for (int i = 0; i < m; i++) 
                        {
                            u[i, k] = 0.0;
                        }
                        u[k, k] = 1.0;
                    }
                }
            }
            
            // If required, generate V.
            if (wantv) 
            {
                for (int k = n-1; k >= 0; k--) 
                {
                    if ((k < nrt) & (e[k] != 0.0)) {
                        for (int j = k+1; j < n; j++) {
                            double t = 0;
                            for (int i = k+1; i < n; i++) {
                                t += v[i, k] * v[i, j];
                            }
                            t = -t/v[k+1, k];
                            for (int i = k+1; i < n; i++) {
                                v[i, j] += t*v[i, k];
                            }
                        }
                    }
                    for (int i = 0; i < n; i++) 
                    {
                        v[i, k] = 0.0;
                    }
                    v[k, k] = 1.0;
                }
            }
            // Main iteration loop for the singular values.
            int pp = p-1;
            int iter = 0;
            double eps = Math.Const.epsilon;
            double tiny = Math.Const.tiny;

            while (p > 0) 
            {
                int k,kase;
                // Here is where a test for too many iterations would go.
                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.

                // kase = 1     if s(p) and e[k-1] are negligible and k<p
                // kase = 2     if s(k) is negligible and k<p
                // kase = 3     if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step).
                // kase = 4     if e(p-1) is negligible (convergence).
                
                for (k = p - 2; k >= -1; k--) 
                {
                    if (k == -1) 
                    {
                        break;
                    }
                    if (System.Math.Abs(e[k]) <= tiny + eps*(System.Math.Abs(m_s[k]) + System.Math.Abs(m_s[k+1]))) 
                    {
                        e[k] = 0.0;
                        break;
                    }
                }
                if (k == p - 2) 
                {
                    kase = 4;
                } 
                else 
                {
                    int ks;
                    for (ks = p-1; ks >= k; ks--) 
                    {
                        if (ks == k) 
                        {
                            break;
                        }
                        double t = (ks != p ? System.Math.Abs(e[ks]) : 0.0) +  (ks != k+1 ? System.Math.Abs(e[ks-1]) : 0.0);

                        if (System.Math.Abs(m_s[ks]) <= tiny + eps*t)  
                        {
                            m_s[ks] = 0.0;
                            break;
                        }
                    }
                    if (ks == k) 
                    {
                        kase = 3;
                    } else if (ks == p-1) 
                    {
                        kase = 1;
                    } else 
                    {
                        kase = 2;
                        k = ks;
                    }
                }
                k++;

                // Perform the task indicated by kase.
                switch (kase) 
                {
                    // Deflate negligible s(p).
                    case 1: 
                        {
                            double f = e[p-2];
                            e[p-2] = 0.0;
                            for (int j = p-2; j >= k; j--) 
                            {
                                double t = Utils.Hypotenuse(m_s[j], f);
                                double cs = m_s[j]/t;
                                double sn = f/t;
                                m_s[j] = t;
                                if (j != k) 
                                {
                                    f = -sn*e[j-1];
                                    e[j-1] = cs*e[j-1];
                                }
                                if (wantv) 
                                {
                                    for (int i = 0; i < n; i++) 
                                    {
                                        t = cs * v[i, j] + sn * v[i, p-1];
                                        v[i, p-1] = -sn * v[i, j] + cs * v[i, p-1];
                                        v[i, j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Split at negligible s(k).
                    case 2: 
                        {
                            double f = e[k-1];
                            e[k-1] = 0.0;
                            for (int j = k; j < p; j++) 
                            {
                                double t = Utils.Hypotenuse(m_s[j], f);
                                double cs = m_s[j]/t;
                                double sn = f/t;
                                m_s[j] = t;
                                f = -sn*e[j];
                                e[j] = cs*e[j];
                                if (wantu) 
                                {
                                    for (int i = 0; i < m; i++) 
                                    {
                                        t = cs*u[i, j] + sn*u[i, k-1];
                                        u[i, k-1] = -sn*u[i, j] + cs*u[i, k-1];
                                        u[i, j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Perform one qr step.
                    case 3: 
                        {
                            // Calculate the shift.
                            double scale = Math.Utils.Max(System.Math.Abs(m_s[p-1]),
                                System.Math.Abs(m_s[p-2]),System.Math.Abs(e[p-2]), 
                                System.Math.Abs(m_s[k]), System.Math.Abs(e[k]));
                            double sp = m_s[p-1]/scale;
                            double spm1 = m_s[p-2]/scale;
                            double epm1 = e[p-2]/scale;
                            double sk = m_s[k]/scale;
                            double ek = e[k]/scale;
                            double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
                            double c = (sp*epm1)*(sp*epm1);
                            double shift = 0.0;
                            if ((b != 0.0) | (c != 0.0)) 
                            {
                                shift = System.Math.Sqrt(b*b + c);
                                if (b < 0.0) 
                                {
                                    shift = -shift;
                                }
                                shift = c/(b + shift);
                            }
                            double f = (sk + sp)*(sk - sp) + shift;
                            double g = sk*ek;
                            
                            // Chase zeros.
                            for (int j = k; j < p-1; j++) 
                            {
                                double t = Utils.Hypotenuse(f, g);
                                double cs = f/t;
                                double sn = g/t;
                                if (j != k) 
                                {
                                    e[j-1] = t;
                                }
                                f = cs*m_s[j] + sn*e[j];
                                e[j] = cs*e[j] - sn*m_s[j];
                                g = sn*m_s[j+1];
                                m_s[j+1] = cs*m_s[j+1];
                                if (wantv) 
                                {
                                    for (int i = 0; i < n; i++) 
                                    {
                                        t = cs*v[i, j] + sn * v[i, j+1];
                                        v[i, j+1] = -sn * v[i, j] + cs * v[i, j+1];
                                        v[i, j] = t;
                                    }
                                }
                                t = Utils.Hypotenuse(f, g);
                                cs = f/t;
                                sn = g/t;
                                m_s[j] = t;
                                f = cs*e[j] + sn*m_s[j+1];
                                m_s[j+1] = -sn*e[j] + cs*m_s[j+1];
                                g = sn*e[j+1];
                                e[j+1] = cs*e[j+1];
                                if (wantu && (j < m-1)) 
                                {
                                    for (int i = 0; i < m; i++) 
                                    {
                                        t = cs * u[i, j] + sn*u[i, j+1];
                                        u[i, j+1] = -sn * u[i, j] + cs * u[i, j+1];
                                        u[i, j] = t;
                                    }
                                }
                            }
                            e[p-2] = f;
                            iter = iter + 1;
                        }
                        break;

                    // Convergence.
                    case 4: 
                        {
                            // Make the singular values positive.   
                            if (m_s[k] <= 0.0) 
                            {
                                m_s[k] = (m_s[k] < 0.0 ? -m_s[k] : 0.0);
                                if (wantv) 
                                {
                                    for (int i = 0; i <= pp; i++) 
                                    {
                                        v[i, k] = -v[i, k];
                                    }
                                }
                            }
                            // Order the singular values.
   
                            while (k < pp) 
                            {
                                if (m_s[k] >= m_s[k+1]) 
                                {
                                    break;
                                }
                                double t = m_s[k];
                                m_s[k] = m_s[k+1];
                                m_s[k+1] = t;
                                if (wantv && (k < n-1)) 
                                {
                                    for (int i = 0; i < n; i++) 
                                    {
                                        t = v[i, k+1]; 
                                        v[i, k+1] = v[i, k]; 
                                        v[i, k] = t;
                                    }
                                }
                                if (wantu && (k < m-1)) 
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = u[i, k+1]; 
                                        u[i, k+1] = u[i, k]; 
                                        u[i, k] = t;
                                    }
                                }
                                k++;
                            }
                            iter = 0;
                            p--;
                        }
                        break;
                }
            }
        }


        /// <summary>
        /// Return the left singular vectors
        /// </summary>
        /// <returns></returns>
        public Matrix GetU()
        {
            return m_U.Clone();
        }

        /// <summary>
        /// Return the right singular vectors
        /// </summary>
        /// <returns></returns>
        public Matrix GetV()
        {
            return m_V.Clone();
        }

        public Matrix GetS()
        {
            Matrix S = new Matrix(m_U.Columns, m_V.Rows, 0.0);

            for (int i = 0; i < Math.Utils.Min(S.Rows, S.Columns, m_s.Length); i++)
            {
                S[i, i] = m_s[i];
            }
            return S;
        }

        /// <summary>
        /// Return the one-dimensional array of singular values
        /// </summary>
        /// <returns></returns>
        public double[] GetSingularValues()
        {
            return (double[])m_s.Clone();
        }

        public double Norm2()
        {
            return m_s[0];
        }

        /// <summary>
        /// Two norm condition number, max(S)/min(S)
        /// </summary>
        /// <returns></returns>
        public double Cond()
        {
            int m = m_U.Rows;
            int n = m_V.Columns;
            return m_s[0] / m_s[Math.Utils.Min(m, n, m_s.Length) - 1]; //all other elements if presents are zero
        }

        /// <summary>
        /// Effective numerical matrix rank
        /// </summary>
        /// <returns></returns>
        public int Rank()
        {
            int m = m_U.Rows;
            int n = m_V.Columns;

            double tol = System.Math.Max(m, n) * m_s[0] * Math.Const.epsilon;
            int r = 0;
            for (int i = 0; i < m_s.Length; i++)
            {
                if (m_s[i] > tol)
                {
                    r++;
                }
            }
            return r;
        }

        
        /// <summary>
        /// Returns the pseudoinverse of a matrix, such that
        /// X = PseudoInverse(A) produces a matrix X of the same dimensions as A' so that A*X*A = A, X*A*X = X.
        /// A = U * S * V', pinv(A) = V * S^-1 * U'
        /// </summary>
        public Matrix PseudoInverse()
        {
            int m = m_U.Rows;
            int n = m_V.Columns;            

            Matrix res = new Matrix(n, m);

            double[] scale = new double[n];

            double tol = System.Math.Max(m, n) * Math.Utils.Max(m_s) * Math.Const.epsilon;

            for (int k = 0; k < n; k++)
            {
                scale[k] = m_s[k] > tol ? 1.0 / m_s[k] : 0.0;
            }            
            
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {                    
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                    {                        
                        sum += m_V[i, k] * scale[k] * m_U[j, k];
                    }
                    res[i, j] = sum;
                }
            }

            return res;
        }

        public double[] Solve(double[] b, bool transposed = false)
        {
            Matrix x = Solve(new Matrix(b, false));

            return x.RowPackedData();
        }

        /// <summary>
        /// Solve Ax = B, given SVD decomposition of A, since decomposition currently implemented for A with rows >= cols
        /// linear system with A rows &lt cols can be solved using SVN decomposition of transposed A i.e. M = A' 
        /// </summary>
        /// <param name="rhs"></param>
        /// <param name="transposed">true is decomposition was contructed using transpose of A</param>
        /// <returns></returns>
        public Matrix Solve(Matrix B, bool transposed = false)
        {
            int m = m_U.Rows;
            int n = m_V.Columns;

            if ((B.Rows != m && !transposed) || (B.Rows != n && transposed))
            {
                throw new ArgumentException("B has wrong number of rows");
            }

            double tol = Math.Const.epsilon * System.Math.Max(m, n) * m_s[0];

            double[] scales = new double[n];

            for (int i = 0; i < n; i++) //m_s should have size n
            {
                scales[i] = ((System.Math.Abs(m_s[i]) > tol) ? 1.0 / m_s[i] : 0.0);
            }            
            
            if (transposed)
            {
                Matrix VTB = new Matrix(n, B.Columns);

                for (int j = 0; j < B.Columns; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        double sum = 0;
                        for (int k = 0; k < n; k++)
                        {
                            sum += m_V[k, i] * B[k, j];
                        }
                        VTB[i, j] = sum * scales[i];
                    }
                }
                return m_U * VTB; 
            }
            else
            {
                Matrix UTB = new Matrix(n, B.Columns);

                for (int j = 0; j < B.Columns; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        double sum = 0;
                        for (int k = 0; k < m; k++)
                        {
                            sum += m_U[k, i] * B[k, j];
                        }
                        UTB[i, j] = sum * scales[i];
                    }
                }
                return m_V * UTB; 
            }
        }

        public Matrix Cov()
        {
            int m = m_U.Rows;
            int n = m_V.Columns;

            double tol = Math.Const.epsilon * System.Math.Max(m, n) * m_s[0];

            double[] scales = new double[n];

            for (int i = 0; i < n; i++) //m_s should have size n
            {
                scales[i] = ((System.Math.Abs(m_s[i]) > tol) ? 1.0 / (m_s[i] * m_s[i]) : 0.0);
            }
            
            Matrix cov = new Matrix(n, n);
            double[,] v = m_V.Data;

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                    {
                        sum += v[i, k] * scales[k] * v[j, k];
                    }
                    cov[i, j] = sum;
                }
            }
            return cov;
        }
    }
}
