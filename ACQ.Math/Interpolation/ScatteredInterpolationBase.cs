using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public interface ScatteredInterpolationInterface
    {
        double Eval(double[] x);
    }

    /// <summary>
    /// Abstract Base class for scattered data interpolation 
    /// </summary>
    public abstract class ScatteredInterpolationBase : ScatteredInterpolationInterface
    {
        #region Members
        protected readonly double[,] m_x; //scattered nodes - one node per row, columns are coordinates 
        protected readonly double[] m_y; //function values
        #endregion Members

        #region Constructors
        /// <summary>
        /// Creates Interpolation
        /// </summary>
        public ScatteredInterpolationBase(double[,] x, double[] y, bool bCopyData = true)
        {
            if (x == null || y == null)
            {
                throw new ArgumentNullException("interpolation arrays can not be null");
            }

            if (y.Length != x.GetLength(0))
            {
                throw new ArgumentException("interpolation array x and y should have incompatible sizes");
            }

            //Data is not copied to save memory
            if (bCopyData)
            {
                m_x = (double[,])x.Clone();
                m_y = (double[])y.Clone();
            }
            else
            {
                m_x = x;
                m_y = y;
            }
        }
        #endregion Constructors

        #region Abstruct Methods

        public abstract double Eval(double[] x);

        #endregion

        #region Public Methods
        /// <summary>
        /// Computes interpolated value at the node
        /// </summary>
        /// <param name="x">node coordinate</param>
        /// <returns></returns>
        public double this[double[] x]
        {
            get
            {
                return Eval(x);
            }
        }

        /// <summary>
        /// Size of interpolation data set
        /// </summary>
        public int Size
        {
            get
            {
                return m_y.Length;
            }
        }
        public int Dim
        {
            get
            {
                return m_x.GetLength(1);
            }
        }
        #endregion Public Methods

    }
}
