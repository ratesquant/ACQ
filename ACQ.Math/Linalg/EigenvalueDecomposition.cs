using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Linalg
{
    /// <summary>
    /// Determines the eigenvalues and eigenvectors of a real square matrix. (Adapted from JAMA)
    /// </summary>
    /// <remarks>
    /// If <c>A</c> is symmetric, then <c>A = V * D * V'</c> and <c>A = V * V'</c>
    /// where the eigenvalue matrix <c>D</c> is diagonal and the eigenvector matrix <c>V</c> is orthogonal.
    /// If <c>A</c> is not symmetric, the eigenvalue matrix <c>D</c> is block diagonal
    /// with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    /// <c>lambda+i*mu</c>, in 2-by-2 blocks, <c>[lambda, mu; -mu, lambda]</c>.
    /// The columns of <c>V</c> represent the eigenvectors in the sense that <c>A * V = V * D</c>.
    /// The matrix V may be badly conditioned, or even singular, so the validity of the equation
    /// <c>A=V*D*inverse(V)</c> depends upon the condition of <c>V</c>.
    /// </remarks>
    public class EigenvalueDecomposition
    {
        private readonly int n;           	// matrix dimension
        private readonly double[] m_d, m_e; 		// storage of eigenvalues.
        private readonly Matrix m_v; 			// storage of eigenvectors.
        private readonly Matrix m_h;  			// storage of nonsymmetric Hessenberg form.
        private readonly double[] m_ort;    	// storage for nonsymmetric algorithm.
        private readonly bool m_symmetric;

        /// <summary>Construct an eigenvalue decomposition.</summary>
        public EigenvalueDecomposition(Matrix A)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }

            if (A.Rows != A.Columns)
            {
                throw new ArgumentException("Matrix is not a square matrix.", "value");
            }

           

            n = A.Columns;
            m_v = new Matrix(n, n);
            m_d = new double[n];
            m_e = new double[n];

            double[,] a = A.Data;
            double[,] v = m_v.Data;

            // Check for symmetry.
            m_symmetric = true;
            for (int j = 0; (j < n) & m_symmetric; j++)
            {
                for (int i = 0; (i < n) & m_symmetric; i++)
                {
                    m_symmetric = (a[i, j] == a[j, i]);
                }
            }

            if (m_symmetric)
            {
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        v[i, j] = a[i, j];
                    }
                }
                tred2(); // Tridiagonalize.
                tql2();// Diagonalize.
            }
            else
            {
                m_h = new Matrix(n, n);
                m_ort = new double[n];

                for (int j = 0; j < n; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        m_h[i, j] = a[i, j];
                    }
                }
                orthes();// Reduce to Hessenberg form.
                hqr2();// Reduce Hessenberg to real Schur form.
            }
        }

        private void tred2()
        {
            // Symmetric Householder reduction to tridiagonal form.
            // This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            double[,] v = m_v.Data;
            for (int j = 0; j < n; j++)
            {
                m_d[j] = v[n - 1, j];
            }

            // Householder reduction to tridiagonal form.
            for (int i = n - 1; i > 0; i--)
            {
                // Scale to avoid under/overflow.
                double scale = 0.0;
                double h = 0.0;
                for (int k = 0; k < i; k++)
                    scale = scale + System.Math.Abs(m_d[k]);

                if (scale == 0.0)
                {
                    m_e[i] = m_d[i - 1];
                    for (int j = 0; j < i; j++)
                    {
                        m_d[j] = v[i - 1, j];
                        v[i, j] = 0.0;
                        v[j, i] = 0.0;
                    }
                }
                else
                {
                    // Generate Householder vector.
                    for (int k = 0; k < i; k++)
                    {
                        m_d[k] /= scale;
                        h += m_d[k] * m_d[k];
                    }

                    double f = m_d[i - 1];
                    double g = System.Math.Sqrt(h);
                    if (f > 0)
                    {
                        g = -g;
                    }

                    m_e[i] = scale * g;
                    h = h - f * g;
                    m_d[i - 1] = f - g;
                    for (int j = 0; j < i; j++)
                        m_e[j] = 0.0;

                    // Apply similarity transformation to remaining columns.
                    for (int j = 0; j < i; j++)
                    {
                        f = m_d[j];
                        v[j, i] = f;
                        g = m_e[j] + v[j, j] * f;
                        for (int k = j + 1; k <= i - 1; k++)
                        {
                            g += v[k, j] * m_d[k];
                            m_e[k] += v[k, j] * f;
                        }
                        m_e[j] = g;
                    }

                    f = 0.0;
                    for (int j = 0; j < i; j++)
                    {
                        m_e[j] /= h;
                        f += m_e[j] * m_d[j];
                    }

                    double hh = f / (h + h);
                    for (int j = 0; j < i; j++)
                        m_e[j] -= hh * m_d[j];

                    for (int j = 0; j < i; j++)
                    {
                        f = m_d[j];
                        g = m_e[j];
                        for (int k = j; k <= i - 1; k++)
                        {
                            v[k, j] -= (f * m_e[k] + g * m_d[k]);
                        }

                        m_d[j] = v[i - 1, j];
                        v[i, j] = 0.0;
                    }
                }
                m_d[i] = h;
            }

            // Accumulate transformations.
            for (int i = 0; i < n - 1; i++)
            {
                v[n-1, i] = v[i, i];
                v[i, i] = 1.0;
                double h = m_d[i + 1];
                if (h != 0.0)
                {
                    for (int k = 0; k <= i; k++)
                        m_d[k] = v[k, i + 1] / h;

                    for (int j = 0; j <= i; j++)
                    {
                        double g = 0.0;
                        for (int k = 0; k <= i; k++)
                            g += v[k, i + 1] * v[k, j];
                        for (int k = 0; k <= i; k++)
                            v[k, j] -= g * m_d[k];
                    }
                }

                for (int k = 0; k <= i; k++)
                    v[k, i + 1] = 0.0;
            }

            for (int j = 0; j < n; j++)
            {
                m_d[j] = v[n - 1, j];
                v[n - 1, j] = 0.0;
            }

            v[n - 1, n - 1] = 1.0;
            m_e[0] = 0.0;
        }

        private void tql2()
        {
            double[,] v = m_v.Data;

            // Symmetric tridiagonal QL algorithm.
            // This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            for (int i = 1; i < n; i++)
                m_e[i - 1] = m_e[i];

            m_e[n - 1] = 0.0;

            double f = 0.0;
            double tst1 = 0.0;
            double eps = System.Math.Pow(2.0, -52.0);

            for (int l = 0; l < n; l++)
            {
                // Find small subdiagonal element.
                tst1 = System.Math.Max(tst1, System.Math.Abs(m_d[l]) + System.Math.Abs(m_e[l]));
                int m = l;
                while (m < n)
                {
                    if (System.Math.Abs(m_e[m]) <= eps * tst1)
                        break;
                    m++;
                }

                // If m == l, d[l] is an eigenvalue, otherwise, iterate.
                if (m > l)
                {
                    int iter = 0;
                    do
                    {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        double g = m_d[l];
                        double p = (m_d[l + 1] - g) / (2.0 * m_e[l]);
                        double r = Utils.Hypotenuse(p, 1.0);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        m_d[l] = m_e[l] / (p + r);
                        m_d[l + 1] = m_e[l] * (p + r);
                        double dl1 = m_d[l + 1];
                        double h = g - m_d[l];
                        for (int i = l + 2; i < n; i++)
                        {
                            m_d[i] -= h;
                        }

                        f = f + h;

                        // Implicit QL transformation.
                        p = m_d[m];
                        double c = 1.0;
                        double c2 = c;
                        double c3 = c;
                        double el1 = m_e[l + 1];
                        double s = 0.0;
                        double s2 = 0.0;
                        for (int i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * m_e[i];
                            h = c * p;
                            r = Utils.Hypotenuse(p, m_e[i]);
                            m_e[i + 1] = s * r;
                            s = m_e[i] / r;
                            c = p / r;
                            p = c * m_d[i] - s * g;
                            m_d[i + 1] = h + s * (c * g + s * m_d[i]);

                            // Accumulate transformation.
                            for (int k = 0; k < n; k++)
                            {
                                h = v[k, i + 1];
                                v[k, i + 1] = s * v[k, i] + c * h;
                                v[k, i] = c * v[k, i] - s * h;
                            }
                        }

                        p = (-s) * s2 * c3 * el1 * m_e[l] / dl1;
                        m_e[l] = s * p;
                        m_d[l] = c * p;

                        // Check for convergence.
                    }
                    while (System.Math.Abs(m_e[l]) > eps * tst1);
                }
                m_d[l] = m_d[l] + f;
                m_e[l] = 0.0;
            }

            // Sort eigenvalues and corresponding vectors.
            for (int i = 0; i < n - 1; i++)
            {
                int k = i;
                double p = m_d[i];
                for (int j = i + 1; j < n; j++)
                {
                    if (m_d[j] < p)
                    {
                        k = j;
                        p = m_d[j];
                    }
                }

                if (k != i)
                {
                    m_d[k] = m_d[i];
                    m_d[i] = p;
                    for (int j = 0; j < n; j++)
                    {
                        p = v[j, i];
                        v[j, i] = v[j, k];
                        v[j, k] = p;
                    }
                }
            }
        }

        private void orthes()
        {
            // Nonsymmetric reduction to Hessenberg form.
            // This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
            int low = 0;
            int high = n - 1;

            for (int m = low + 1; m <= high - 1; m++)
            {
                // Scale column.

                double scale = 0.0;
                for (int i = m; i <= high; i++)
                    scale = scale + System.Math.Abs(m_h[i, m - 1]);

                if (scale != 0.0)
                {
                    // Compute Householder transformation.
                    double h = 0.0;
                    for (int i = high; i >= m; i--)
                    {
                        m_ort[i] = m_h[i, m - 1] / scale;
                        h += m_ort[i] * m_ort[i];
                    }

                    double g = System.Math.Sqrt(h);
                    if (m_ort[m] > 0) g = -g;

                    h = h - m_ort[m] * g;
                    m_ort[m] = m_ort[m] - g;

                    // Apply Householder similarity transformation
                    // H = (I - u * u' / h) * H * (I - u * u') / h)
                    for (int j = m; j < n; j++)
                    {
                        double f = 0.0;
                        for (int i = high; i >= m; i--)
                            f += m_ort[i] * m_h[i, j];

                        f = f / h;
                        for (int i = m; i <= high; i++)
                            m_h[i, j] -= f * m_ort[i];
                    }

                    for (int i = 0; i <= high; i++)
                    {
                        double f = 0.0;
                        for (int j = high; j >= m; j--)
                            f += m_ort[j] * m_h[i, j];

                        f = f / h;
                        for (int j = m; j <= high; j++)
                            m_h[i, j] -= f * m_ort[j];
                    }

                    m_ort[m] = scale * m_ort[m];
                    m_h[m, m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    m_v[i, j] = (i == j ? 1.0 : 0.0);

            for (int m = high - 1; m >= low + 1; m--)
            {
                if (m_h[m, m - 1] != 0.0)
                {
                    for (int i = m + 1; i <= high; i++)
                        m_ort[i] = m_h[i, m - 1];

                    for (int j = m; j <= high; j++)
                    {
                        double g = 0.0;
                        for (int i = m; i <= high; i++)
                            g += m_ort[i] * m_v[i, j];

                        // Double division avoids possible underflow.
                        g = (g / m_ort[m]) / m_h[m, m - 1];
                        for (int i = m; i <= high; i++)
                            m_v[i, j] += g * m_ort[i];
                    }
                }
            }
        }

        private void cdiv(double xr, double xi, double yr, double yi, out double cdivr, out double cdivi)
        {
            // Complex scalar division.
            double r;
            double d;
            if (System.Math.Abs(yr) > System.Math.Abs(yi))
            {
                r = yi / yr;
                d = yr + r * yi;
                cdivr = (xr + r * xi) / d;
                cdivi = (xi - r * xr) / d;
            }
            else
            {
                r = yr / yi;
                d = yi + r * yr;
                cdivr = (r * xr + xi) / d;
                cdivi = (r * xi - xr) / d;
            }
        }

        private void hqr2()
        {
            // Nonsymmetric reduction from Hessenberg to real Schur form.   
            // This is derived from the Algol procedure hqr2, by Martin and Wilkinson, Handbook for Auto. Comp.,
            // Vol.ii-Linear Algebra, and the corresponding  Fortran subroutine in EISPACK.
            int nn = this.n;
            int n = nn - 1;
            int low = 0;
            int high = nn - 1;
            double eps = System.Math.Pow(2.0, -52.0);
            double exshift = 0.0;
            double p = 0;
            double q = 0;
            double r = 0;
            double s = 0;
            double z = 0;
            double t;
            double w;
            double x;
            double y;

            // Store roots isolated by balanc and compute matrix norm
            double norm = 0.0;
            for (int i = 0; i < nn; i++)
            {
                if (i < low | i > high)
                {
                    m_d[i] = m_h[i, i];
                    m_e[i] = 0.0;
                }

                for (int j = System.Math.Max(i - 1, 0); j < nn; j++)
                    norm = norm + System.Math.Abs(m_h[i, j]);
            }

            // Outer loop over eigenvalue index
            int iter = 0;
            while (n >= low)
            {
                // Look for single small sub-diagonal element
                int l = n;
                while (l > low)
                {
                    s = System.Math.Abs(m_h[l - 1, l - 1]) + System.Math.Abs(m_h[l, l]);
                    if (s == 0.0) s = norm;
                    if (System.Math.Abs(m_h[l, l - 1]) < eps * s)
                        break;

                    l--;
                }

                // Check for convergence
                if (l == n)
                {
                    // One root found
                    m_h[n, n] = m_h[n, n] + exshift;
                    m_d[n] = m_h[n, n];
                    m_e[n] = 0.0;
                    n--;
                    iter = 0;
                }
                else if (l == n - 1)
                {
                    // Two roots found
                    w = m_h[n, n - 1] * m_h[n - 1, n];
                    p = (m_h[n - 1, n - 1] - m_h[n, n]) / 2.0;
                    q = p * p + w;
                    z = System.Math.Sqrt(System.Math.Abs(q));
                    m_h[n, n] = m_h[n, n] + exshift;
                    m_h[n - 1, n - 1] = m_h[n - 1, n - 1] + exshift;
                    x = m_h[n, n];

                    if (q >= 0)
                    {
                        // Real pair
                        z = (p >= 0) ? (p + z) : (p - z);
                        m_d[n - 1] = x + z;
                        m_d[n] = m_d[n - 1];
                        if (z != 0.0)
                            m_d[n] = x - w / z;
                        m_e[n - 1] = 0.0;
                        m_e[n] = 0.0;
                        x = m_h[n, n - 1];
                        s = System.Math.Abs(x) + System.Math.Abs(z);
                        p = x / s;
                        q = z / s;
                        r = System.Math.Sqrt(p * p + q * q);
                        p = p / r;
                        q = q / r;

                        // Row modification
                        for (int j = n - 1; j < nn; j++)
                        {
                            z = m_h[n - 1, j];
                            m_h[n - 1, j] = q * z + p * m_h[n, j];
                            m_h[n, j] = q * m_h[n, j] - p * z;
                        }

                        // Column modification
                        for (int i = 0; i <= n; i++)
                        {
                            z = m_h[i, n - 1];
                            m_h[i, n - 1] = q * z + p * m_h[i, n];
                            m_h[i, n] = q * m_h[i, n] - p * z;
                        }

                        // Accumulate transformations
                        for (int i = low; i <= high; i++)
                        {
                            z = m_v[i, n - 1];
                            m_v[i, n - 1] = q * z + p * m_v[i, n];
                            m_v[i, n] = q * m_v[i, n] - p * z;
                        }
                    }
                    else
                    {
                        // Complex pair
                        m_d[n - 1] = x + p;
                        m_d[n] = x + p;
                        m_e[n - 1] = z;
                        m_e[n] = -z;
                    }

                    n = n - 2;
                    iter = 0;
                }
                else
                {
                    // No convergence yet	 

                    // Form shift
                    x = m_h[n, n];
                    y = 0.0;
                    w = 0.0;
                    if (l < n)
                    {
                        y = m_h[n - 1, n - 1];
                        w = m_h[n, n - 1] * m_h[n - 1, n];
                    }

                    // Wilkinson's original ad hoc shift
                    if (iter == 10)
                    {
                        exshift += x;
                        for (int i = low; i <= n; i++)
                            m_h[i, i] -= x;

                        s = System.Math.Abs(m_h[n, n - 1]) + System.Math.Abs(m_h[n - 1, n - 2]);
                        x = y = 0.75 * s;
                        w = -0.4375 * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30)
                    {
                        s = (y - x) / 2.0;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = System.Math.Sqrt(s);
                            if (y < x) s = -s;
                            s = x - w / ((y - x) / 2.0 + s);
                            for (int i = low; i <= n; i++)
                                m_h[i, i] -= s;
                            exshift += s;
                            x = y = w = 0.964;
                        }
                    }

                    iter = iter + 1;

                    // Look for two consecutive small sub-diagonal elements
                    int m = n - 2;
                    while (m >= l)
                    {
                        z = m_h[m, m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / m_h[m + 1, m] + m_h[m, m + 1];
                        q = m_h[m + 1, m + 1] - z - r - s;
                        r = m_h[m + 2, m + 1];
                        s = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                        p = p / s;
                        q = q / s;
                        r = r / s;
                        if (m == l)
                            break;
                        if (System.Math.Abs(m_h[m, m - 1]) * (System.Math.Abs(q) + System.Math.Abs(r)) < eps * (System.Math.Abs(p) * (System.Math.Abs(m_h[m - 1, m - 1]) + System.Math.Abs(z) + System.Math.Abs(m_h[m + 1, m + 1]))))
                            break;
                        m--;
                    }

                    for (int i = m + 2; i <= n; i++)
                    {
                        m_h[i, i - 2] = 0.0;
                        if (i > m + 2)
                            m_h[i, i - 3] = 0.0;
                    }

                    // Double QR step involving rows l:n and columns m:n
                    for (int k = m; k <= n - 1; k++)
                    {
                        bool notlast = (k != n - 1);
                        if (k != m)
                        {
                            p = m_h[k, k - 1];
                            q = m_h[k + 1, k - 1];
                            r = (notlast ? m_h[k + 2, k - 1] : 0.0);
                            x = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                            if (x != 0.0)
                            {
                                p = p / x;
                                q = q / x;
                                r = r / x;
                            }
                        }

                        if (x == 0.0) break;

                        s = System.Math.Sqrt(p * p + q * q + r * r);
                        if (p < 0) s = -s;

                        if (s != 0)
                        {
                            if (k != m)
                                m_h[k, k - 1] = -s * x;
                            else
                                if (l != m)
                                    m_h[k, k - 1] = -m_h[k, k - 1];

                            p = p + s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q = q / p;
                            r = r / p;

                            // Row modification
                            for (int j = k; j < nn; j++)
                            {
                                p = m_h[k, j] + q * m_h[k + 1, j];
                                if (notlast)
                                {
                                    p = p + r * m_h[k + 2, j];
                                    m_h[k + 2, j] = m_h[k + 2, j] - p * z;
                                }

                                m_h[k, j] = m_h[k, j] - p * x;
                                m_h[k + 1, j] = m_h[k + 1, j] - p * y;
                            }

                            // Column modification
                            for (int i = 0; i <= System.Math.Min(n, k + 3); i++)
                            {
                                p = x * m_h[i, k] + y * m_h[i, k + 1];
                                if (notlast)
                                {
                                    p = p + z * m_h[i, k + 2];
                                    m_h[i, k + 2] = m_h[i, k + 2] - p * r;
                                }

                                m_h[i, k] = m_h[i, k] - p;
                                m_h[i, k + 1] = m_h[i, k + 1] - p * q;
                            }

                            // Accumulate transformations
                            for (int i = low; i <= high; i++)
                            {
                                p = x * m_v[i, k] + y * m_v[i, k + 1];
                                if (notlast)
                                {
                                    p = p + z * m_v[i, k + 2];
                                    m_v[i, k + 2] = m_v[i, k + 2] - p * r;
                                }

                                m_v[i, k] = m_v[i, k] - p;
                                m_v[i, k + 1] = m_v[i, k + 1] - p * q;
                            }
                        }
                    }
                }
            }

            // Backsubstitute to find vectors of upper triangular form
            if (norm == 0.0)
            {
                return;
            }

            for (n = nn - 1; n >= 0; n--)
            {
                p = m_d[n];
                q = m_e[n];

                // Real vector
                if (q == 0)
                {
                    int l = n;
                    m_h[n, n] = 1.0;
                    for (int i = n - 1; i >= 0; i--)
                    {
                        w = m_h[i, i] - p;
                        r = 0.0;
                        for (int j = l; j <= n; j++)
                            r = r + m_h[i, j] * m_h[j, n];

                        if (m_e[i] < 0.0)
                        {
                            z = w;
                            s = r;
                        }
                        else
                        {
                            l = i;
                            if (m_e[i] == 0.0)
                            {
                                m_h[i, n] = (w != 0.0) ? (-r / w) : (-r / (eps * norm));
                            }
                            else
                            {
                                // Solve real equations
                                x = m_h[i, i + 1];
                                y = m_h[i + 1, i];
                                q = (m_d[i] - p) * (m_d[i] - p) + m_e[i] * m_e[i];
                                t = (x * s - z * r) / q;
                                m_h[i, n] = t;
                                m_h[i + 1, n] = (System.Math.Abs(x) > System.Math.Abs(z)) ? ((-r - w * t) / x) : ((-s - y * t) / z);
                            }

                            // Overflow control
                            t = System.Math.Abs(m_h[i, n]);
                            if ((eps * t) * t > 1)
                                for (int j = i; j <= n; j++)
                                    m_h[j, n] = m_h[j, n] / t;
                        }
                    }
                }
                else if (q < 0)
                {
                    double cdivr, cdivi;

                    // Complex vector
                    int l = n - 1;

                    // Last vector component imaginary so matrix is triangular
                    if (System.Math.Abs(m_h[n, n - 1]) > System.Math.Abs(m_h[n - 1, n]))
                    {
                        m_h[n - 1, n - 1] = q / m_h[n, n - 1];
                        m_h[n - 1, n] = -(m_h[n, n] - p) / m_h[n, n - 1];
                    }
                    else
                    {
                        cdiv(0.0, -m_h[n - 1, n], m_h[n - 1, n - 1] - p, q, out cdivr, out cdivi);
                        m_h[n - 1, n - 1] = cdivr;
                        m_h[n - 1, n] = cdivi;
                    }

                    m_h[n, n - 1] = 0.0;
                    m_h[n, n] = 1.0;
                    for (int i = n - 2; i >= 0; i--)
                    {
                        double ra, sa, vr, vi;
                        ra = 0.0;
                        sa = 0.0;
                        for (int j = l; j <= n; j++)
                        {
                            ra = ra + m_h[i, j] * m_h[j, n - 1];
                            sa = sa + m_h[i, j] * m_h[j, n];
                        }

                        w = m_h[i, i] - p;

                        if (m_e[i] < 0.0)
                        {
                            z = w;
                            r = ra;
                            s = sa;
                        }
                        else
                        {
                            l = i;
                            if (m_e[i] == 0)
                            {
                                cdiv(-ra, -sa, w, q, out cdivr, out cdivi);
                                m_h[i, n - 1] = cdivr;
                                m_h[i, n] = cdivi;
                            }
                            else
                            {
                                // Solve complex equations
                                x = m_h[i, i + 1];
                                y = m_h[i + 1, i];
                                vr = (m_d[i] - p) * (m_d[i] - p) + m_e[i] * m_e[i] - q * q;
                                vi = (m_d[i] - p) * 2.0 * q;
                                if (vr == 0.0 & vi == 0.0)
                                    vr = eps * norm * (System.Math.Abs(w) + System.Math.Abs(q) + System.Math.Abs(x) + System.Math.Abs(y) + System.Math.Abs(z));
                                cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, out cdivr, out cdivi);
                                m_h[i, n - 1] = cdivr;
                                m_h[i, n] = cdivi;
                                if (System.Math.Abs(x) > (System.Math.Abs(z) + System.Math.Abs(q)))
                                {
                                    m_h[i + 1, n - 1] = (-ra - w * m_h[i, n - 1] + q * m_h[i, n]) / x;
                                    m_h[i + 1, n] = (-sa - w * m_h[i, n] - q * m_h[i, n - 1]) / x;
                                }
                                else
                                {
                                    cdiv(-r - y * m_h[i, n - 1], -s - y * m_h[i, n], z, q, out cdivr, out cdivi);
                                    m_h[i + 1, n - 1] = cdivr;
                                    m_h[i + 1, n] = cdivi;
                                }
                            }

                            // Overflow control
                            t = System.Math.Max(System.Math.Abs(m_h[i, n - 1]), System.Math.Abs(m_h[i, n]));
                            if ((eps * t) * t > 1)
                                for (int j = i; j <= n; j++)
                                {
                                    m_h[j, n - 1] = m_h[j, n - 1] / t;
                                    m_h[j, n] = m_h[j, n] / t;
                                }
                        }
                    }
                }
            }

            // Vectors of isolated roots
            for (int i = 0; i < nn; i++)
                if (i < low | i > high)
                    for (int j = i; j < nn; j++)
                        m_v[i, j] = m_h[i, j];

            // Back transformation to get eigenvectors of original matrix
            for (int j = nn - 1; j >= low; j--)
                for (int i = low; i <= high; i++)
                {
                    z = 0.0;
                    for (int k = low; k <= System.Math.Min(j, high); k++)
                        z = z + m_v[i, k] * m_h[k, j];
                    m_v[i, j] = z;
                }
        }

        /// <summary>Returns the real parts of the eigenvalues.</summary>
        public double[] GetRealEigenvalues()
        {
            return (double[])m_d.Clone();
        }

        /// <summary>Returns the imaginary parts of the eigenvalues.</summary>	
        public double[] GetImagEigenvalues()
        {
            return (double[])m_e.Clone();            
        }

        /// <summary>Returns the eigenvector matrix.</summary>
        public Matrix GetV()
        {
            return m_v.Clone();            
        }

        /// <summary>Returns the block diagonal eigenvalue matrix.</summary>
        public Matrix GetD()
        {   
            Matrix D = new Matrix(n, n);
            double[,] d = D.Data;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    d[i, j] = 0.0;
                }

                d[i, i] = m_d[i];
                if (m_e[i] > 0)
                {
                    d[i, i + 1] = m_e[i];
                }
                else if (m_e[i] < 0)
                {
                    d[i, i - 1] = m_e[i];
                }
            }

            return D;            
        }
    }
}
