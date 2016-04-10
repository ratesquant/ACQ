using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Globalization;

namespace ACQ.Math.Linalg
{
    /// <summary>Matrix 
    /// clarity over performance
    /// </summary>
    public class Matrix
    {
        private double[,] m_data;

        #region Constructors
        /// <summary>Constructs an empty matrix of the given size.</summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        public Matrix(int rows, int columns)
        {
            m_data = new double[rows, columns];     
        }

        public Matrix(Matrix A) : this(A.Rows, A.Columns)
        {
            m_data = (double[,])A.m_data.Clone();
        }

        /// <summary>Constructs a matrix of the given size and assigns a given value to all elements.</summary>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        /// <param name="value">Value to assign to all elements.</param>
        public Matrix(int rows, int columns, double value) : this(rows, columns)
        {
            for (int i = 0; i < rows ; i++)
            {
                for (int j = 0; j <  cols; j++)
                {
                    m_data[i, j] = value;
                }
            }
        }

        /// <summary>Constructs a matrix from the given array.</summary>
        /// <param name="value">The array the matrix gets constructed from.</param>        
        public Matrix(double[][] data)
        {
            if (data == null || data[0] == null)
            {
                throw new ArgumentNullException("data");
            }

            m_data = new double[data.Length, data[0].Length];

            for (int i = 0; i < rows; i++)
            {
                if (data[i].Length != cols)
                    throw new ArgumentException("Argument out of range.");

                for (int j = 0; j < cols; j++)
                {
                    m_data[i, j] = data[i][j];
                }                
            }          
        }

        public Matrix(double[,] data)
        {
            m_data = (double[,])data.Clone();
        }

        public Matrix(double[] data, bool isRow = true)
            : this(isRow ? 1 : data.Length, isRow ? data.Length : 1)
        {
            if (isRow)
            {
                for (int i = 0; i < data.Length; i++)
                {
                    m_data[0, i] = data[i];
                }
            }
            else 
            {
                for (int i = 0; i < data.Length; i++)
                {
                    m_data[i, 0] = data[i];
                }
            }
        }
        #endregion

        private int rows
        {
            get
            {
                return m_data.GetLength(0);
            }
        }

        private int cols
        {
            get
            {
                return m_data.GetLength(1);
            }
        }



        /// <summary>Determines weather two instances are equal.</summary>
        public override bool Equals(object obj)
        {
            return Equals(this, (Matrix)obj);
        }

        /// <summary>Determines weather two instances are equal.</summary>
        public static bool Equals(Matrix A, Matrix B)
        {
            if (((object)A) == ((object)B))
            {
                return true;
            }

            if ((((object)A) == null) || (((object)B) == null))
            {
                return false;
            }

            if ((A.Rows != B.Rows) || (A.Columns != B.Columns))
            {
                return false;
            }

            for (int i = 0; i < A.Rows; i++)
            {
                for (int j = 0; j < A.Columns; j++)
                {
                    if (A[i, j] != B[i, j])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        public override int GetHashCode()
        {
            return m_data.GetHashCode();
        }

        internal double[,] Data
        {
            get
            {
                return m_data;
            }
        }

        /// <summary>Returns the number of columns.</summary>
        public int Rows
        {
            get
            {
                return rows;
            }
        }

        /// <summary>Returns the number of columns.</summary>
        public int Columns
        {
            get
            {
                return cols;
            }
        }

        /// <summary>Return true if the matrix is a square matrix.</summary>
        public bool Square
        {
            get
            {
                return (rows == cols);
            }
        }

        /// <summary>Returns true if the matrix is symmetric.</summary>
        public bool Symmetric
        {
            get
            {
                if (Square)
                {
                    for (int i = 0; i < rows; i++)
                    {
                        for (int j = 0; j <= i; j++)
                        {
                            if (m_data[i, j] != m_data[i, j])
                            {
                                return false;
                            }
                        }
                    }

                    return true;
                }

                return false;
            }
        }

        /// <summary>Access the value at the given location.</summary>
        public double this[int row, int column]
        {
            set
            {
                m_data[row, column] = value;
            }

            get
            {
                return m_data[row, column];
            }
        }

        public bool CheckIndexes(int row, int column)
        {
            return row >= 0 && column >= 0 && row < rows && column < cols;
        }
       

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="startRow">Start row index</param>
        /// <param name="endRow">End row index</param>
        /// <param name="startColumn">Start column index</param>
        /// <param name="endColumn">End column index</param>
        public Matrix Submatrix(int startRow, int endRow, int startColumn, int endColumn)
        {
            if ((startRow > endRow) || (startColumn > endColumn) || 
                (startRow < 0) || (startRow >= rows) || 
                (endRow < 0) || (endRow >= rows) || (startColumn < 0) || 
                (startColumn >= cols) || (endColumn < 0) || (endColumn >= cols))
            {
                throw new ArgumentException("Argument out of range.");
            }

            Matrix X = new Matrix(endRow - startRow + 1, endColumn - startColumn + 1);
            double[,] x = X.Data;
            for (int i = startRow; i <= endRow; i++)
            {
                for (int j = startColumn; j <= endColumn; j++)
                {
                    x[i - startRow, j - startColumn ] = m_data[i, j];
                }
            }

            return X;
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="rowIndexes">Array of row indices</param>
        /// <param name="columnIndexes">Array of column indices</param>
        public Matrix Submatrix(int[] rowIndexes, int[] columnIndexes)
        {
            Matrix X = new Matrix(rowIndexes.Length, columnIndexes.Length);
            double[,] x = X.Data;
            for (int i = 0; i < rowIndexes.Length; i++)
            {
                for (int j = 0; j < columnIndexes.Length; j++)
                {
                    if ((rowIndexes[i] < 0) || (rowIndexes[i] >= rows) || 
                        (columnIndexes[j] < 0) || (columnIndexes[j] >= cols))
                    {
                        throw new ArgumentException("Argument out of range.");
                    }

                    x[i, j] = m_data[rowIndexes[i], columnIndexes[j]];
                }
            }

            return X;
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="i0">Start row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="c">Array of row indices</param>
        public Matrix Submatrix(int i0, int i1, int[] c)
        {
            if ((i0 > i1) || (i0 < 0) || (i0 >= rows) || (i1 < 0) || (i1 >= rows))
            {
                throw new ArgumentException("Argument out of range.");
            }

            Matrix X = new Matrix(i1 - i0 + 1, c.Length);
            double[,] x = X.Data;
            for (int i = i0; i <= i1; i++)
            {
                for (int j = 0; j < c.Length; j++)
                {
                    if ((c[j] < 0) || (c[j] >= cols))
                    {
                        throw new ArgumentException("Argument out of range.");
                    }

                    x[i - i0, j] = m_data[i, c[j]];
                }
            }

            return X;
        }

        /// <summary>Returns a sub matrix extracted from the current matrix.</summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        public Matrix Submatrix(int[] r, int j0, int j1)
        {
            if ((j0 > j1) || (j0 < 0) || (j0 >= cols) || (j1 < 0) || (j1 >= cols))
            {
                throw new ArgumentException("Argument out of range.");
            }

            Matrix X = new Matrix(r.Length, j1 - j0 + 1);
            double[,] x = X.Data;
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = j0; j <= j1; j++)
                {
                    if ((r[i] < 0) || (r[i] >= rows))
                    {
                        throw new ArgumentException("Argument out of range.");
                    }

                    x[i, j - j0] = m_data[r[i], j];
                }
            }

            return X;
        }

        /// <summary>Creates a copy of the matrix.</summary>
        public Matrix Clone()
        {            
            return new Matrix(this);
        }

        /// <summary>Returns the transposed matrix.</summary>
        public Matrix Transpose()
        {
            Matrix X = new Matrix(cols, rows);
            double[,] x = X.Data;
           
            for (var j = 0; j < cols; j++)
            {
                for (var i = 0; i < rows; i++)
                {
                    x[j, i] = m_data[i, j];
                }
            }

            return X;
        }

        /// <summary>Returns the One Norm for the matrix.</summary>
        /// <value>The maximum column sum.</value>
        public double Norm1
        {
            get
            {
                double f = 0;
                for (int j = 0; j < cols; j++)
                {
                    double s = 0;
                    for (int i = 0; i < rows; i++)
                    {
                        s += System.Math.Abs(m_data[i, j]);
                    }

                    f = System.Math.Max(f, s);
                }
                return f;
            }
        }

        /// <summary>Returns the Infinity Norm for the matrix.</summary>
        /// <value>The maximum row sum.</value>
        public double InfinityNorm
        {
            get
            {
                double f = 0;
                for (int i = 0; i < rows; i++)
                {
                    double s = 0;
                    for (int j = 0; j < cols; j++)
                    {
                        s += System.Math.Abs(m_data[i, j]);
                    }
                    f = System.Math.Max(f, s);
                }
                return f;
            }
        }

        /// <summary>Returns the Frobenius Norm for the matrix.</summary>
        /// <value>The square root of sum of squares of all elements.</value>
        public double FrobeniusNorm
        {
            get
            {
                double f = 0;
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        f = Utils.Hypotenuse(f, m_data[i, j]);
                    }
                }

                return f;
            }
        }

        /// <summary>Unary minus.</summary>
        public static Matrix Negate(Matrix A)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }

            int rows = A.Rows;
            int columns = A.Columns;
            double[,] data = A.Data;

            Matrix X = new Matrix(rows, columns);
            double[,] x = X.Data;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j <  columns; j++)
                {
                    x[i, j] = -data[i, j];
                }
            }

            return X;
        }

        /// <summary>Unary minus.</summary>
        public static Matrix operator -(Matrix A)
        {
            if (A == null)
            {
                throw new ArgumentNullException("value");
            }

            return Negate(A);
        }

        /// <summary>Matrix equality.</summary>
        public static bool operator ==(Matrix A, Matrix B)
        {
            return Equals(A, B);
        }

        /// <summary>Matrix inequality.</summary>
        public static bool operator !=(Matrix A, Matrix B)
        {
            return !Equals(A, B);
        }
      

        /// <summary>Matrix addition.</summary>
        public static Matrix Add(Matrix A, Matrix B)
        {
            if (A == null || B == null)
            {
                throw new ArgumentNullException();
            }

            int rows = A.Rows;
            int columns = A.Columns;
            double[,] a = A.Data;

            if ((rows != B.Rows) || (columns != B.Columns))
            {
                throw new ArgumentException("Matrix dimension do not match.");
            }

            Matrix X = new Matrix(rows, columns);
            double[,] x = X.Data;
            double[,] b = B.Data;
            
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j <  columns; j++)
                {
                    x[i, j] = a[i, j] + b[i, j];
                }
            }
            return X;
        }

        /// <summary>Matrix addition.</summary>
        public static Matrix operator +(Matrix A, Matrix B)
        {
            return Add(A, B);
        }

        /// <summary>Matrix subtraction.</summary>
        public static Matrix Subtract(Matrix A, Matrix B)
        {
            if (A == null || B == null)
            {
                throw new ArgumentNullException();
            }

            int rows = A.Rows;
            int columns = A.Columns;
            double[,] a = A.Data;

            if ((rows != B.Rows) || (columns != B.Columns))
            {
                throw new ArgumentException("Matrix dimension do not match.");
            }

            Matrix X = new Matrix(rows, columns);
            double[,] x = X.Data;
            double[,] b = B.Data;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j <  columns; j++)
                {
                    x[i, j] = a[i, j] - b[i, j];
                }
            }
            return X;
        }

        /// <summary>Matrix subtraction.</summary>
        public static Matrix operator -(Matrix A, Matrix B)
        {
            return Subtract(A, B);
        }

        /// <summary>Matrix-scalar multiplication.</summary>
        public static Matrix Multiply(Matrix A, double c)
        {
            if (A == null)
            {
                throw new ArgumentNullException("A");
            }

            int rows = A.Rows;
            int columns = A.Columns;
            double[,] a = A.Data;

            Matrix X = new Matrix(rows, columns);
            double[,] x = X.Data;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j <  columns; j++)
                {
                    x[i, j] = a[i, j] * c;
                }
            }

            return X;
        }

        /// <summary>Matrix-scalar multiplication.</summary>
        public static Matrix operator *(Matrix A, double c)
        {
            return Multiply(A, c);
        }

        /// <summary>Matrix-scalar multiplication.</summary>
        public static Matrix operator *(double c, Matrix A)
        {
            return Multiply(A, c);
        }

        /// <summary>Matrix-matrix multiplication.</summary>
        public static Matrix Multiply(Matrix A, Matrix B)
        {
            if (A == null || B == null)
            {
                throw new ArgumentNullException();
            }

            int rows = A.Rows;
            double[,] a = A.Data;

            if (B.Rows != A.Columns)
            {
                throw new ArgumentException("Matrix dimensions are not valid.");
            }

            int columns = B.Columns;
            int r_rows = B.Rows;
            Matrix X = new Matrix(rows, columns);
            double[,] x = X.Data;
            double[,] b = B.Data;

            int size = A.Columns;           

            for (int m1 = 0; m1 < rows; m1++)
            {
                for (var n1 = 0; n1 < columns; n1++)
                {
                    double sum = 0;
                    for (var k1 = 0; k1 < r_rows; k1++)
                    {
                        sum += a[m1, k1] * b[k1, n1];
                    }
                    x[m1, n1] += sum;
                }
            }
            return X;
        }

        /// <summary>Matrix-matrix multiplication.</summary>
        public static Matrix operator *(Matrix A, Matrix B)
        {
            return Multiply(A, B);
        }

        /// <summary>Returns the LHS solution vetor if the matrix is square or the least squares solution otherwise.</summary>
        public Matrix Solve(Matrix B)
        {
            return (rows == cols) ? new LuDecomposition(this).Solve(B) : new QrDecomposition(this).Solve(B);
        }

        /// <summary>Inverse of the matrix if matrix is square, pseudoinverse otherwise.</summary>
        public Matrix GetInverse()
        {
            return Solve(Matrix.CreateDiagonal(rows, rows, 1.0));            
        }

        /// <summary>Determinant if matrix is square.</summary>
        public double GetDeterminant()
        {
            return new LuDecomposition(this).GetDeterminant();            
        }

        /// <summary>Returns the trace of the matrix.</summary>
        /// <returns>Sum of the diagonal elements.</returns>
        public double Trace
        {
            get
            {
                double trace = 0;
                for (int i = 0; i < System.Math.Min(rows, cols); i++)
                {
                    trace += m_data[i, i];
                }
                return trace;
            }
        }

        /// <summary>Returns a diagonal matrix of the given size.</summary>
        public static Matrix CreateDiagonal(int rows, int columns, double value)
        {
            Matrix X = new Matrix(rows, columns);
            double[,] x = X.Data;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    x[i, j] = ((i == j) ? value : 0.0);
                }
            }
            return X;
        }

        public static Matrix CreateDiagonal(double[] values)
        {
            Matrix res = new Matrix(values.Length, values.Length);

            for (int i = 0; i < values.Length; i++)
            {
                res[i, i] = values[i];
            }
            return res;
        }

        /// <summary>Returns the matrix in a textual form.</summary>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(rows).Append(" x ").Append(cols).Append(" : ");
            sb.AppendLine();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    sb.Append(m_data[i, j]).Append(", ");
                }
                sb.AppendLine();
            }

            return sb.ToString();            
        }
        /*
        public static void SVD(Matrix Mat_, out double[,] S_, out double[,] U_, out double[,] V_)
        {
            int Rows, Cols;
            int m, MP, n, NP;
            double[] w;
            double[,] A, v;

            m = Mat_.Rows + 1;
            n = Mat_.Columns + 1;

            if (m < n)
            {
                m = n;
                MP = NP = n;
            }
            else if (m > n)
            {
                n = m;
                NP = MP = m;
            }
            else
            {
                MP = m;
                NP = n;
            }

            A = new double[m + 1, n + 1];

            for (int row = 1; row <= Mat_.Rows + 1; row++)
                for (int col = 1; col <= Mat_.Columns + 1; col++)
                { A[row, col] = Mat_[row - 1, col - 1]; }

            const int NMAX = 100;
            v = new double[NP + 1, NP + 1];
            w = new double[NP + 1];

            int k, l, nm;
            int flag, i, its, j, jj;

            double[,] U_temp, S_temp, V_temp;
            double anorm, c, f, g, h, s, scale, x, y, z;
            double[] rv1 = new double[NMAX];

            l = 0;
            nm = 0;
            g = 0.0;
            scale = 0.0;
            anorm = 0.0;

            for (i = 1; i <= n; i++)
            {
                l = i + 1;
                rv1[i] = scale * g;
                g = s = scale = 0.0;
                if (i <= m)
                {
                    for (k = i; k <= m; k++) scale += Math.Abs(A[k, i]);
                    if (scale != 0)
                    {
                        for (k = i; k <= m; k++)
                        {
                            A[k, i] /= scale;
                            s += A[k, i] * A[k, i];
                        }
                        f = A[i, i];
                        g = -Utils.Sign(Math.Sqrt(s), f);
                        h = f * g - s;
                        A[i, i] = f - g;
                        if (i != n)
                        {
                            for (j = l; j <= n; j++)
                            {
                                for (s = 0, k = i; k <= m; k++) s += A[k, i] * A[k, j];
                                f = s / h;
                                for (k = i; k <= m; k++) A[k, j] += f * A[k, i];
                            }
                        }
                        for (k = i; k <= m; k++) A[k, i] *= scale;
                    }
                }
                w[i] = scale * g;
                g = s = scale = 0.0;
                if (i <= m && i != n)
                {
                    for (k = l; k <= n; k++) scale += Math.Abs(A[i, k]);
                    if (scale != 0)
                    {
                        for (k = l; k <= n; k++)
                        {
                            A[i, k] /= scale;
                            s += A[i, k] * A[i, k];
                        }
                        f = A[i, l];
                        g = -Utils.Sign(Math.Sqrt(s), f);
                        h = f * g - s;
                        A[i, l] = f - g;
                        for (k = l; k <= n; k++) rv1[k] = A[i, k] / h;
                        if (i != m)
                        {
                            for (j = l; j <= m; j++)
                            {
                                for (s = 0.0, k = l; k <= n; k++) s += A[j, k] * A[i, k];
                                for (k = l; k <= n; k++) A[j, k] += s * rv1[k];
                            }
                        }
                        for (k = l; k <= n; k++) A[i, k] *= scale;
                    }
                }
                anorm = Math.Max(anorm, (Math.Abs(w[i]) + Math.Abs(rv1[i])));
            }
            for (i = n; i >= 1; i--)
            {
                if (i < n)
                {
                    if (g != 0)
                    {
                        for (j = l; j <= n; j++)
                            v[j, i] = (A[i, j] / A[i, l]) / g;
                        for (j = l; j <= n; j++)
                        {
                            for (s = 0.0, k = l; k <= n; k++) s += A[i, k] * v[k, j];
                            for (k = l; k <= n; k++) v[k, j] += s * v[k, i];
                        }
                    }
                    for (j = l; j <= n; j++) v[i, j] = v[j, i] = 0.0;
                }
                v[i, i] = 1.0;
                g = rv1[i];
                l = i;
            }
            for (i = n; i >= 1; i--)
            {
                l = i + 1;
                g = w[i];
                if (i < n)
                    for (j = l; j <= n; j++) A[i, j] = 0.0;
                if (g != 0)
                {
                    g = 1.0 / g;
                    if (i != n)
                    {
                        for (j = l; j <= n; j++)
                        {
                            for (s = 0.0, k = l; k <= m; k++) s += A[k, i] * A[k, j];
                            f = (s / A[i, i]) * g;
                            for (k = i; k <= m; k++) A[k, j] += f * A[k, i];
                        }
                    }
                    for (j = i; j <= m; j++) A[j, i] *= g;
                }
                else
                {
                    for (j = i; j <= m; j++) A[j, i] = 0.0;
                }
                ++A[i, i];
            }
            for (k = n; k >= 1; k--)
            {
                for (its = 1; its <= 30; its++)
                {
                    flag = 1;
                    for (l = k; l >= 1; l--)
                    {
                        nm = l - 1;
                        if (Math.Abs(rv1[l]) + anorm == anorm)
                        {
                            flag = 0;
                            break;
                        }
                        if (Math.Abs(w[nm]) + anorm == anorm) break;
                    }
                    if (flag != 0)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i <= k; i++)
                        {
                            f = s * rv1[i];
                            if (Math.Abs(f) + anorm != anorm)
                            {
                                g = w[i];
                                h = Utils.Hypotenuse(f, g);
                                w[i] = h;
                                h = 1.0 / h;
                                c = g * h;
                                s = (-f * h);
                                for (j = 1; j <= m; j++)
                                {
                                    y = A[j, nm];
                                    z = A[j, i];
                                    A[j, nm] = y * c + z * s;
                                    A[j, i] = z * c - y * s;
                                }
                            }
                        }
                    }
                    z = w[k];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            w[k] = -z;
                            for (j = 1; j <= n; j++) v[j, k] = (-v[j, k]);
                        }
                        break;
                    }
                    if (its == 30) Console.WriteLine("No convergence in 30 SVDCMP iterations");
                    x = w[l];
                    nm = k - 1;
                    y = w[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = Utils.Hypotenuse(f, 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + Utils.Sign(g, f))) - h)) / x;
                    c = s = 1.0;
                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = Utils.Hypotenuse(f, h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y = y * c;
                        for (jj = 1; jj <= n; jj++)
                        {
                            x = v[jj, j];
                            z = v[jj, i];
                            v[jj, j] = x * c + z * s;
                            v[jj, i] = z * c - x * s;
                        }
                        z = Utils.Hypotenuse(f, h);
                        w[j] = z;
                        if (z != 0)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }
                        f = (c * g) + (s * y);
                        x = (c * y) - (s * g);
                        for (jj = 1; jj <= m; jj++)
                        {
                            y = A[jj, j];
                            z = A[jj, i];
                            A[jj, j] = y * c + z * s;
                            A[jj, i] = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0;
                    rv1[k] = f;
                    w[k] = x;
                }
            }

            S_temp = new double[NP, NP];
            V_temp = new double[NP, NP];
            U_temp = new double[MP, NP];

            for (i = 1; i <= NP; i++) S_temp[i - 1, i - 1] = w[i];

            S_ = S_temp;

            for (i = 1; i <= NP; i++)
                for (j = 1; j <= NP; j++) V_temp[i - 1, j - 1] = v[i, j];

            V_ = V_temp;

            for (i = 1; i <= MP; i++)
                for (j = 1; j <= NP; j++) U_temp[i - 1, j - 1] = A[i, j];

            U_ = U_temp;
        }

        /// <summary>
        /// Returns the pseudoinverse of a matrix, such that
        /// X = PINV(A) produces a matrix 'X' of the same dimensions
        /// as A' so that A*X*A = A, X*A*X = X.
        /// In case of an error the error is raised as an exception.
        /// </summary>
        /// <param name="Mat">An array whose pseudoinverse is to be found</param>
        /// <returns>The pseudoinverse of the array as an array</returns>
        public static Matrix PINV(Matrix Mat)
        {
            Matrix Matrix = new Matrix(Mat);
            int i, m, n, j, l;
            double[,] S;
            double EPS, MulAdd, Tol, Largest_Item = 0;

            double[,] svdU, svdS, svdV;

            m = Mat.Rows;
            n = Mat.Columns;

            SVD(Mat, out svdS, out svdU, out svdV);

            EPS = 2.2204E-16;
            m++;
            n++;

            Matrix Part_I = new Matrix(m, n);
            Matrix Part_II = new Matrix(m, n);
            S = new Double[n, n];

            MulAdd = 0;
            for (i = 0; i <= svdS.GetUpperBound(0); i++)
            {
                if (i == 0) Largest_Item = svdS[0, 0];
                if (Largest_Item < svdS[i, i]) Largest_Item = svdS[i, i];
            }

            if (m > n) Tol = m * Largest_Item * EPS;
            else Tol = n * Largest_Item * EPS;

            for (i = 0; i < n; i++) S[i, i] = svdS[i, i];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    for (l = 0; l < n; l++)
                    {
                        if (S[l, j] > Tol) MulAdd += svdU[i, l] * (1 / S[l, j]);
                    }
                    Part_I[i, j] = MulAdd;
                    MulAdd = 0;
                }
            }

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    for (l = 0; l < n; l++)
                    {
                        MulAdd += Part_I[i, l] * svdV[j, l];
                    }
                    Part_II[i, j] = MulAdd;
                    MulAdd = 0;
                }
            }
            return Part_II.Transpose();
        }*/
    }
        
  
    /*
    /// <summary>
    /// 	Singular Value Decomposition for a rectangular matrix.
    /// </summary>
    /// <remarks>
    ///	  For an m-by-n matrix <c>A</c> with <c>m >= n</c>, the singular value decomposition is
    ///   an m-by-n orthogonal matrix <c>U</c>, an n-by-n diagonal matrix <c>S</c>, and
    ///   an n-by-n orthogonal matrix <c>V</c> so that <c>A = U * S * V'</c>.
    ///   The singular values, <c>sigma[k] = S[k,k]</c>, are ordered so that
    ///   <c>sigma[0] >= sigma[1] >= ... >= sigma[n-1]</c>.
    ///   The singular value decompostion always exists, so the constructor will
    ///   never fail. The matrix condition number and the effective numerical
    ///   rank can be computed from this decomposition.
    /// </remarks>
    public class SingularValueDecomposition
    {
        private double[][] m_u;
        private double[][] m_v;
        private double[] m_s; // singular values
        private int m;
        private int n;

        private double m_eps = Math.Pow(2.0, -52.0);

        /// <summary>Construct singular value decomposition.</summary>
        public SingularValueDecomposition(Matrix value)
        {
            if (value == null)
            {
                throw new ArgumentNullException("value");
            }

            Matrix copy = (Matrix)value.Clone();
            double[] a = copy.Array;
            m = value.Rows;
            n = value.Columns;
            int nu = Math.Min(m, n);
            m_s = new double[Math.Min(m + 1, n)];
            m_u = Kamlex.MathLib.Utils.JaggedArray<double>(m, nu);
            m_v = Kamlex.MathLib.Utils.JaggedArray<double>(n, n);
            double[][] u = m_u;
            double[][] v = m_v;
            double[] e = new double[n];
            double[] work = new double[m];
            bool wantu = true;
            bool wantv = true;

            // Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e.
            int nct = Math.Min(m - 1, n);
            int nrt = Math.Max(0, Math.Min(n - 2, m));
            for (int k = 0; k < Math.Max(nct, nrt); k++)
            {
                if (k < nct)
                {
                    // Compute the transformation for the k-th column and place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    m_s[k] = 0;
                    for (int i = k; i < m; i++)
                    {
                        m_s[k] = Hypotenuse(m_s[k], a[k*m + i]);
                    }

                    if (m_s[k] != 0.0)
                    {
                        if (a[k*m+k] < 0.0)
                        {
                            m_s[k] = -m_s[k];
                        }

                        for (int i = k; i < m; i++)
                        {
                            a[k*m + i] /= m_s[k];
                        }

                        a[k*m + k] += 1.0;
                    }

                    m_s[k] = -m_s[k];
                }

                for (int j = k + 1; j < n; j++)
                {
                    if ((k < nct) & (m_s[k] != 0.0))
                    {
                        // Apply the transformation.
                        double t = 0;
                        for (int i = k; i < m; i++)
                            t += a[k*m +i] * a[j*m+i];
                        t = -t / a[k *m +k];
                        for (int i = k; i < m; i++)
                            a[j*m +i] += t * a[k*m + i];
                    }

                    // Place the k-th row of A into e for the subsequent calculation of the row transformation.
                    e[j] = a[j*m +k];
                }

                if (wantu & (k < nct))
                {
                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (int i = k; i < m; i++)
                        u[i][k] = a[k*m+i];
                }

                if (k < nrt)
                {
                    // Compute the k-th row transformation and place the k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0;
                    for (int i = k + 1; i < n; i++)
                    {
                        e[k] = Hypotenuse(e[k], e[i]);
                    }

                    if (e[k] != 0.0)
                    {
                        if (e[k + 1] < 0.0)
                            e[k] = -e[k];

                        for (int i = k + 1; i < n; i++)
                            e[i] /= e[k];

                        e[k + 1] += 1.0;
                    }

                    e[k] = -e[k];
                    if ((k + 1 < m) & (e[k] != 0.0))
                    {
                        // Apply the transformation.
                        for (int i = k + 1; i < m; i++)
                            work[i] = 0.0;

                        for (int j = k + 1; j < n; j++)
                            for (int i = k + 1; i < m; i++)
                                work[i] += e[j] * a[j*m +i];

                        for (int j = k + 1; j < n; j++)
                        {
                            double t = -e[j] / e[k + 1];
                            for (int i = k + 1; i < m; i++)
                                a[j*m+i] += t * work[i];
                        }
                    }

                    if (wantv)
                    {
                        // Place the transformation in V for subsequent back multiplication.
                        for (int i = k + 1; i < n; i++)
                            v[i][k] = e[i];
                    }
                }
            }

            // Set up the final bidiagonal matrix or order p.
            int p = Math.Min(n, m + 1);
            if (nct < n) m_s[nct] = a[nct * m + nct];
            if (m < p) m_s[p - 1] = 0.0;
            if (nrt + 1 < p) e[nrt] = a[(p - 1) * m + nrt];
            e[p - 1] = 0.0;

            // If required, generate U.
            if (wantu)
            {
                for (int j = nct; j < nu; j++)
                {
                    for (int i = 0; i < m; i++)
                        u[i][j] = 0.0;
                    u[j][j] = 1.0;
                }

                for (int k = nct - 1; k >= 0; k--)
                {
                    if (m_s[k] != 0.0)
                    {
                        for (int j = k + 1; j < nu; j++)
                        {
                            double t = 0;
                            for (int i = k; i < m; i++)
                                t += u[i][k] * u[i][j];

                            t = -t / u[k][k];
                            for (int i = k; i < m; i++)
                                u[i][j] += t * u[i][k];
                        }

                        for (int i = k; i < m; i++)
                            u[i][k] = -u[i][k];

                        u[k][k] = 1.0 + u[k][k];
                        for (int i = 0; i < k - 1; i++)
                            u[i][k] = 0.0;
                    }
                    else
                    {
                        for (int i = 0; i < m; i++)
                            u[i][k] = 0.0;
                        u[k][k] = 1.0;
                    }
                }
            }

            // If required, generate V.
            if (wantv)
            {
                for (int k = n - 1; k >= 0; k--)
                {
                    if ((k < nrt) & (e[k] != 0.0))
                    {
                        for (int j = k + 1; j < nu; j++)
                        {
                            double t = 0;
                            for (int i = k + 1; i < n; i++)
                                t += v[i][k] * v[i][j];

                            t = -t / v[k + 1][k];
                            for (int i = k + 1; i < n; i++)
                                v[i][j] += t * v[i][k];
                        }
                    }

                    for (int i = 0; i < n; i++)
                        v[i][k] = 0.0;
                    v[k][k] = 1.0;
                }
            }

            // Main iteration loop for the singular values.
            int pp = p - 1;
            int iter = 0;
            double eps = Math.Pow(2.0, -52.0);
            double tiny = Math.Pow(2.0, -966.0);
            while (p > 0)
            {
                int k, kase;

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
                        break;

                    if (Math.Abs(e[k]) <= tiny + eps * (Math.Abs(m_s[k]) + Math.Abs(m_s[k + 1])))
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
                    for (ks = p - 1; ks >= k; ks--)
                    {
                        if (ks == k)
                            break;

                        double t = (ks != p ? Math.Abs(e[ks]) : 0.0) + (ks != k + 1 ? Math.Abs(e[ks - 1]) : 0.0);
                        if (Math.Abs(m_s[ks]) <= tiny + eps * t)
                        {
                            m_s[ks] = 0.0;
                            break;
                        }
                    }

                    if (ks == k)
                        kase = 3;
                    else if (ks == p - 1)
                        kase = 1;
                    else
                    {
                        kase = 2;
                        k = ks;
                    }
                }

                k++;

                //if (iter > 10000)
                //    break;
                    //kase = 4;

                // Perform the task indicated by kase.
                switch (kase)
                {
                    // Deflate negligible s(p).
                    case 1:
                        {
                            double f = e[p - 2];
                            e[p - 2] = 0.0;
                            for (int j = p - 2; j >= k; j--)
                            {
                                double t = Hypotenuse(m_s[j], f);
                                double cs = m_s[j] / t;
                                double sn = f / t;
                                m_s[j] = t;
                                if (j != k)
                                {
                                    f = -sn * e[j - 1];
                                    e[j - 1] = cs * e[j - 1];
                                }

                                if (wantv)
                                {
                                    for (int i = 0; i < n; i++)
                                    {
                                        t = cs * v[i][j] + sn * v[i][p - 1];
                                        v[i][p - 1] = -sn * v[i][j] + cs * v[i][p - 1];
                                        v[i][j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Split at negligible s(k).
                    case 2:
                        {
                            double f = e[k - 1];
                            e[k - 1] = 0.0;
                            for (int j = k; j < p; j++)
                            {
                                double t = Hypotenuse(m_s[j], f);
                                double cs = m_s[j] / t;
                                double sn = f / t;
                                m_s[j] = t;
                                f = -sn * e[j];
                                e[j] = cs * e[j];
                                if (wantu)
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = cs * u[i][j] + sn * u[i][k - 1];
                                        u[i][k - 1] = -sn * u[i][j] + cs * u[i][k - 1];
                                        u[i][j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Perform one qr step.
                    case 3:
                        {
                            // Calculate the shift.
                            double scale = Utils.Max(Math.Abs(m_s[p - 1]), Math.Abs(m_s[p - 2]), Math.Abs(e[p - 2]), Math.Abs(m_s[k]), Math.Abs(e[k]));
                            double sp = m_s[p - 1] / scale;
                            double spm1 = m_s[p - 2] / scale;
                            double epm1 = e[p - 2] / scale;
                            double sk = m_s[k] / scale;
                            double ek = e[k] / scale;
                            double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                            double c = (sp * epm1) * (sp * epm1);
                            double shift = 0.0;
                            if ((b != 0.0) | (c != 0.0))
                            {
                                shift = Math.Sqrt(b * b + c);
                                if (b < 0.0)
                                    shift = -shift;
                                shift = c / (b + shift);
                            }

                            double f = (sk + sp) * (sk - sp) + shift;
                            double g = sk * ek;

                            // Chase zeros.
                            for (int j = k; j < p - 1; j++)
                            {
                                double t = Hypotenuse(f, g);
                                double cs = f / t;
                                double sn = g / t;
                                if (j != k)
                                    e[j - 1] = t;
                                f = cs * m_s[j] + sn * e[j];
                                e[j] = cs * e[j] - sn * m_s[j];
                                g = sn * m_s[j + 1];
                                m_s[j + 1] = cs * m_s[j + 1];
                                if (wantv)
                                {
                                    for (int i = 0; i < n; i++)
                                    {
                                        t = cs * v[i][j] + sn * v[i][j + 1];
                                        v[i][j + 1] = -sn * v[i][j] + cs * v[i][j + 1];
                                        v[i][j] = t;
                                    }
                                }

                                t = Hypotenuse(f, g);
                                cs = f / t;
                                sn = g / t;
                                m_s[j] = t;
                                f = cs * e[j] + sn * m_s[j + 1];
                                m_s[j + 1] = -sn * e[j] + cs * m_s[j + 1];
                                g = sn * e[j + 1];
                                e[j + 1] = cs * e[j + 1];
                                if (wantu && (j < m - 1))
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = cs * u[i][j] + sn * u[i][j + 1];
                                        u[i][j + 1] = -sn * u[i][j] + cs * u[i][j + 1];
                                        u[i][j] = t;
                                    }
                                }
                            }

                            e[p - 2] = f;
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
                                    for (int i = 0; i <= pp; i++)
                                        v[i][k] = -v[i][k];
                            }

                            // Order the singular values.
                            while (k < pp)
                            {
                                if (m_s[k] >= m_s[k + 1])
                                    break;

                                double t = m_s[k];
                                m_s[k] = m_s[k + 1];
                                m_s[k + 1] = t;
                                if (wantv && (k < n - 1))
                                    for (int i = 0; i < n; i++)
                                    {
                                        t = v[i][k + 1];
                                        v[i][k + 1] = v[i][k];
                                        v[i][k] = t;
                                    }

                                if (wantu && (k < m - 1))
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = u[i][k + 1];
                                        u[i][k + 1] = u[i][k];
                                        u[i][k] = t;
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

        /// <summary>Returns the condition number <c>max(S) / min(S)</c>.</summary>
        public double Condition
        {
            get
            {
                return m_s[0] / m_s[Math.Min(m, n) - 1];
            }
        }

        /// <summary>Returns the Two norm.</summary>
        public double Norm2
        {
            get
            {
                return m_s[0];
            }
        }

        /// <summary>Returns the effective numerical matrix rank.</summary>
        /// <value>Number of non-negligible singular values.</value>
        public int Rank
        {
            get
            {               
                double tol = Math.Max(m, n) * m_s[0] * m_eps;
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
        }

        public Matrix U
        {
            get 
            {
                return new Matrix(m_u);
            }
        }

        public Matrix V
        {
            get
            {
                return new Matrix(m_v);
            }
        }

        /// <summary>Return the one-dimensional array of singular values.</summary>		
        public double[] SingularValues
        {
            get
            {
                return this.m_s;
            }
        }

        public Matrix Inverse()
        {
            double tol = m_eps * Math.Max(m, n) * m_s[0];

            int vrows = m_v.Length;
            int vcols = m_v[0].Length;

            double[,] X = new double[vrows, m_s.Length];

            //X =  C * S^-1
            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < vcols; j++)
                {
                    if (Math.Abs(m_s[j]) > tol)
                        X[i, j] = m_v[i][j] / m_s[j];
                }
            }

            //Y = X*U'
            int urows = m_u.Length;
            int ucols = m_u[0].Length;

            Matrix res = new Matrix(vrows, urows);
            
            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    double sum = 0;

                    for (int k = 0; k < ucols; k++)
                        sum += X[i, k] * m_u[j][k];
                    res[i, j] = sum;
                }
            }
            return res;
        }

        /// <summary>
        /// Solve Ax = b, given SVD decomposition of A, since decomposition currently implemented for A with rows >= cols
        /// linear system with A rows &lt cols can be solved using SVN decomposition of transposed A i.e. M = A' 
        /// </summary>
        /// <param name="rhs"></param>
        /// <param name="transposed">true is decomposition was contructed using transpose of A</param>
        /// <returns></returns>
        public double[] Solve(double[] rhs, bool transposed)
        {
            Matrix Y = new Matrix(rhs, false);

            if(rhs.Length != m && rhs.Length != n)
                throw new ArgumentException("rhs has wrong size");

            System.Diagnostics.Debug.Assert(m_s.Length == n);
            
            double tol = m_eps * Math.Max(m, n) * m_s[0];
            
            double[] svals = new double[n];
            double[][] VL = Kamlex.MathLib.Utils.JaggedArray<double>(n, n);

            for (int i = 0; i < n; i++)
               svals[i] = ((Math.Abs(m_s[i]) > tol) ? 1.0 / m_s[i] : 0.0);
            
            //compute V*S^-1
            for (int i = 0; i < n; i++)
            {                
                for (int j = 0; j < n; j++)
                {                    
                    VL[i][j] = m_v[i][j] * svals[j]; 
                }
            }

            int urows = m_u.Length;
            int vrows = m_v.Length;

            double[][] VLU = Kamlex.MathLib.Utils.JaggedArray<double>(vrows, urows);

            //compute V*S^-1*U^-1
            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    double sum = 0;

                    for (int k = 0; k < n; k++)
                        sum += VL[i][k] * m_u[j][k];
                    VLU[i][j] = sum;
                }
            }

            if (transposed)
            {
                double[] res = new double[m];

                for (int i = 0; i < m; i++)
                {
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                        sum += rhs[k] * VLU[k][i];
                    res[i] = sum;
                }

                return res;
            }
            else
            {
                double[] res = new double[n];

                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int k = 0; k < m; k++)
                        sum += VLU[i][k] * rhs[k];
                    res[i] = sum;
                }
                return res;
            }
        }




        private static double Hypotenuse(double a, double b)
        {
            if (Math.Abs(a) > Math.Abs(b))
            {
                double r = b / a;
                return Math.Abs(a) * Math.Sqrt(1 + r * r);
            }

            if (b != 0)
            {
                double r = a / b;
                return Math.Abs(b) * Math.Sqrt(1 + r * r);
            }

            return 0.0;
        }
    }*/
}
