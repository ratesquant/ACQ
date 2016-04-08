using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Linalg
{
    public class Vector
    {
        private double[] m_data;

        #region Constructors
        public Vector(int size)
        {
            m_data = new double[size];     
        }

        public Vector(int size, double value) : this(size)
        {
            Utils.FillArray(m_data, value);
        }

        public Vector(Vector a) : this(a.Size)
        {
            a.m_data.CopyTo(m_data, 0);
        }
        public Vector(double[] data, bool copyData = false)            
        {
            if (copyData)
            {
                m_data = (double[])data.Clone();
            }
            else
            {
                m_data = data;
            }
        }
        #endregion

        public int Size
        {
            get
            {
                return m_data.Length;
            }
        }

        /// <summary>Determines weather two instances are equal.</summary>
        public override bool Equals(object obj)
        {
            return Equals(this, (Vector)obj);
        }

        /// <summary>Determines weather two instances are equal.</summary>
        public static bool Equals(Vector a, Vector b)
        {
            if (((object)a) == ((object)b))
            {
                return true;
            }

            if ((((object)a) == null) || (((object)b) == null))
            {
                return false;
            }

            if ((a.Size != b.Size))
            {
                return false;
            }

            for (int i = 0; i < a.Size; i++)
            {
                if (a.m_data[i] != b.m_data[i])
                {
                    return false;
                }
            }

            return true;
        }

        public override int GetHashCode()
        {
            return m_data.GetHashCode();
        }

        internal double[] Data
        {
            get
            {
                return m_data;
            }
        }
        public double this[int i]
        {
            set
            {
                m_data[i] = value;
            }

            get
            {
                return m_data[i];
            }
        }

        public Vector Clone()
        {
            return new Vector(this);
        }

        /// <summary>Unary minus.</summary>
        public static Vector Negate(Vector a)
        {
            if (a == null)
            {
                throw new ArgumentNullException("a");
            }

            int size = a.Size;
            double[] data = a.m_data;

            Vector X = new Vector(size);
            double[] x = X.m_data;
            for (int i = 0; i < size; i++)
            {
                x[i] = -data[i];
            }

            return X;
        }

        /// <summary>Unary minus.</summary>
        public static Vector operator -(Vector a)
        {
            if (a == null)
            {
                throw new ArgumentNullException("a");
            }

            return Negate(a);
        }

        /// <summary>Vector equality.</summary>
        public static bool operator ==(Vector a, Vector b)
        {
            return Equals(a, b);
        }

        /// <summary>Vector inequality.</summary>
        public static bool operator !=(Vector a, Vector b)
        {
            return !Equals(a, b);
        }

        /// <summary>Vector-scalar multiplication.</summary>
        public static Vector Multiply(Vector a, double c)
        {
            if (a == null)
            {
                throw new ArgumentNullException("a");
            }

            int size = a.Size;
            double[] data = a.Data;

            Vector X = new Vector(size);
            double[] x = X.Data;

            for (int i = 0; i < size; i++)
            {
                x[i] = a[i] * c;
            }

            return X;
        }

        /// <summary>Vector-scalar addition.</summary>
        public static Vector Add(Vector a, double c)
        {
            if (a == null)
            {
                throw new ArgumentNullException("a");
            }

            int size = a.Size;
            double[] data = a.Data;

            Vector X = new Vector(size);
            double[] x = X.Data;

            for (int i = 0; i < size; i++)
            {
                x[i] = a[i] + c;
            }

            return X;
        }

        /// <summary>Vector-scalar multiplication.</summary>
        public static Vector operator *(Vector a, double c)
        {
            return Multiply(a, c);
        }

        /// <summary>Vector-scalar multiplication.</summary>
        public static Vector operator *(double c, Vector a)
        {
            return Multiply(a, c);
        }

        /// <summary>Vector-scalar addition.</summary>
        public static Vector operator +(Vector a, double c)
        {
            return Add(a, c);
        }

        /// <summary>Vector-scalar addition.</summary>
        public static Vector operator +(double c, Vector a)
        {
            return Add(a, c);
        }

        /// <summary>Vector dot product</summary>
        public static double DotProduct(Vector a, Vector b)
        {
            if (a == null || b == null)
            {
                throw new ArgumentNullException();
            }
            
            if (a.Size != b.Size)
            {
                throw new ArgumentException("Vectors should have the same size");
            }

            int size = a.Size;
            double[] ad = a.Data;
            double[] bd = b.Data;           

            double sum = 0;
            for (int i = 0; i < size; i++)
            {
                sum += ad[i] * bd[i];
            }
            return sum;
        }

        /// <summary>
        /// Vector dot product
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static double operator *(Vector A, Vector B)
        {
            return DotProduct(A, B);
        }
    }
}
