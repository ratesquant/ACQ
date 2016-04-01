using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public interface InterpolationInterface2D
    {
        double Eval(double x1, double x2);
    }

    /// <summary>
    /// Abstract Base class for 1D interpolation 
    /// </summary>
    public abstract class InterpolationBase2D : InterpolationInterface2D
    {
        #region Members
        protected readonly double[] m_x1; //horizontal nodes (along x-axis)
        protected readonly double[] m_x2; //vertical nodes (along y-axis)
        protected readonly double[,] m_y; //grid of values, number of rows should be the same as size of x2
        #endregion Members

        #region Constructors
        /// <summary>
        /// Creates Interpolation
        /// </summary>
        /// <param name="x">interpolation nodes - should be sorted i.e. x[i+1] > x[i]</param>
        /// <param name="y">function values at interpolation nodes i.e. y[i] = f(x[i])</param>        
        /// <param name="bCopyData"> Make copy of input arrays if true, otherwise save the reference</param>
        public InterpolationBase2D(double[] x1, double[] x2, double[,] y, bool bCopyData = false)
        {
            if (x1 == null || x2 == null || y == null)
                throw new ArgumentNullException("interpolation arrays can not be null");

            if (x1.Length != y.GetLength(1) || x2.Length != y.GetLength(0))
                throw new ArgumentException("interpolation array x1, x2 and y have incompatible sizes");

            if (x1.Length < 1 || x2.Length < 1)
                throw new ArgumentException("interpolation array should have at least 2 nodes");

            //check that data is ordered   
            if (!isSorted(x1) || !isSorted(x2))
            {
                throw new ArgumentException("interpolation nodes should be ordered");                
            }


            //Data is not copied to save memory
            if (bCopyData)
            {
                m_x1 = (double[])x1.Clone();
                m_x2 = (double[])x2.Clone();
                m_y  = (double[,])y.Clone();
            }
            else
            {
                m_x1 = x1;
                m_x2 = x2;
                m_y = y;
            }
        }
        #endregion Constructors

        private static bool isSorted(double[] x)
        {
            for (int i = 0; i < x.Length - 1; i++)
            {
                if (x[i + 1] <= x[i])
                {
                    return false;
                }
            }
            return true;
        }

        #region Abstruct Methods

        public abstract double Eval(double x1, double x2);

        #endregion

        #region Public Methods
        /// <summary>
        /// Computes interpolated value at the node
        /// </summary>
        /// <param name="x">node coordinate</param>
        /// <returns></returns>
        public double this[double x1, double x2]
        {
            get
            {
                return Eval(x1, x2);
            }
        }

        /// <summary>
        /// Size of interpolation data set
        /// </summary>
        public int Size1
        {
            get
            {
                return m_x1.Length;
            }
        }
        public int Size2
        {
            get
            {
                return m_x2.Length;
            }
        }
        #endregion Public Methods

        #region Protected Methods


        /// <summary>
        /// Function returns the index i such that x[i-1] < x <= x[i]         
        /// </summary>
        /// <param name="xi"></param>
        /// <returns></returns>
        protected void FindIndex(double x1, double x2, out int index1, out int index2)
        {
            index1 = FindIndex(m_x1, x1);
            index2 = FindIndex(m_x2, x2);
        }

        /// <summary>
        /// Function returns the index i such that a[i-1] < x <= a[i]         
        /// </summary>
        /// <param name="xi"></param>
        /// <returns></returns>
        private static int FindIndex(double[] a, double x)
        {
            int index;

            if (x < a[0] || x > a[a.Length - 1])
            {
                index = 0;
            }
            else
            {
                index = Array.BinarySearch<double>(a, x);

                if (index == 0) // x = a[0]
                {
                    index = 1; //we can always do this because we check in constructor that there are at least two nodes
                }
                else if (index < 0)
                {
                    index = ~index; //binary search returns compliment of the node larger than a value
                }
            }

            return index;
        }

        #endregion
    }
}
