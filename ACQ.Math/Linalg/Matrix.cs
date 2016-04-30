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
                {
                    throw new ArgumentException("Argument out of range.");
                }

                for (int j = 0; j < cols; j++)
                {
                    m_data[i, j] = data[i][j];
                }                
            }          
        }

        public Matrix(double[,] data, bool copy = true)
        {
            if (copy)
            {
                m_data = (double[,])data.Clone();
            }else
            {
                m_data = data;
            }
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


        public static Matrix CreateIdentity(int n)
        {
            Matrix d = new Matrix(n, n, 0.0);

            for (int i = 0; i < n; i++)
            {
                d[i, i] = 1.0;
            }
            return d;
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

        /// <summary>
        /// Create random n x m matrix 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="m"></param>
        /// <param name="rng"></param>
        /// <returns></returns>
        public static Matrix CreateRandom(int n, int m, ACQ.Math.Random.IRandomGenerator rng)
        {
            Matrix random = new Matrix(n, m);
            double[,] d = random.m_data;

            for (int i = 0; i < random.Rows; i++)
            {
                for (int j = 0; j < random.Columns; j++)
                {
                    d[i, j] = rng.NextDouble();
                }
            }
            return random;
        }
        /// <summary>
        /// Create magic square matrix
        /// </summary>
        /// <param name="n">size of the matrix</param>
        /// <returns></returns>
        public static Matrix CreateMagic(int n)
        {
            Matrix magic = new Matrix(n, n);
            double[,] M = magic.m_data;
            double t;

            if ((n % 2) == 1)   // Odd order
            {
                int a = (n+1)/2;
                int b = (n+1);
                for (int j = 0; j < n; j++) 
                {
                    for (int i = 0; i < n; i++) 
                    {
                        M[i, j] = n*((i+j+a) % n) + ((i+2*j+b) % n) + 1;
                    }
                }                
            }
            else if ((n % 4) == 0) // Doubly Even Order
            {
                for (int j = 0; j < n; j++) 
                {
                    for (int i = 0; i < n; i++) 
                    {
                        if (((i+1)/2)%2 == ((j+1)/2)%2) 
                        {
                            M[i, j] = n*n-n*i-j;
                        } else 
                        {
                            M[i, j] = n*i+j+1;
                        }
                    }
                }
            }
            else  // Singly Even Order
            {
                int p = n/2;
                int k = (n-2)/4;
                Matrix A = CreateMagic(p);

                for (int j = 0; j < p; j++) 
                {
                    for (int i = 0; i < p; i++) 
                    {
                        double aij = A[i, j];
                        M[i, j] = aij;
                        M[i, j+p] = aij + 2*p*p;
                        M[i+p, j] = aij + 3*p*p;
                        M[i+p, j+p] = aij + p*p;
                    }
                }
                for (int i = 0; i < p; i++) 
                {
                    for (int j = 0; j < k; j++) 
                    {
                        t = M[i, j]; 
                        M[i, j] = M[i+p, j];
                        M[i+p, j] = t;
                    }
                    for (int j = n-k+1; j < n; j++) 
                    {
                        t = M[i, j]; 
                        M[i, j] = M[i+p, j];
                        M[i+p, j] = t;
                    }
                }
                t = M[k, 0]; 
                M[k, 0] = M[k+p, 0]; 
                M[k+p, 0] = t;
                t = M[k, k];
                M[k, k] = M[k+p, k];
                M[k+p, k] = t;
            }
            return magic;
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
        public bool IsSquare
        {
            get
            {
                return (rows == cols);
            }
        }

        /// <summary>Returns true if the matrix is symmetric.</summary>
        public bool IsSymmetric
        {
            get
            {
                if (IsSquare)
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

        public double[] RowPackedData()
        {
            double[] d = new double[rows * cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    d[j + i * cols] = m_data[i, j];
                }
            }
            return d;
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
        public double Norm1()
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

        /// <summary>Returns the Infinity Norm for the matrix.</summary>
        /// <value>The maximum row sum.</value>
        public double InfinityNorm()
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

        /// <summary>Returns the Frobenius Norm for the matrix.</summary>
        /// <value>The square root of sum of squares of all elements.</value>
        public double FrobeniusNorm()
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
        public double Trace()
        {
            double trace = 0;
            for (int i = 0; i < System.Math.Min(rows, cols); i++)
            {
                trace += m_data[i, i];
            }
            return trace;
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
    }
}
