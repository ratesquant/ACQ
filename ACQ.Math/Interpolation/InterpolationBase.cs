using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public interface InterpolationInterface
    {
        double Eval(double x);
        bool Bounds { get; set; }
    }
    /// <summary>
    /// Abstract Base class for 1D interpolation 
    /// </summary>
    public abstract class InterpolationBase : InterpolationInterface
    {
        #region Members
        protected readonly double[] m_x;
        protected readonly double[] m_y;
        protected bool m_bounds = true; //this is not needed in constructor and can be switched on or off
        #endregion Members

        #region Constructors
        /// <summary>
        /// Creates Interpolation
        /// </summary>
        /// <param name="x">interpolation nodes - should be sorted i.e. x[i+1] > x[i]</param>
        /// <param name="y">function values at interpolation nodes i.e. y[i] = f(x[i])</param>        
        /// <param name="bCopyData"> Make copy of input arrays if true, otherwise save the reference</param>
        public InterpolationBase(double[] x, double[] y, bool bCopyData = true)
        {
            if (x == null || y == null)
                throw new ArgumentNullException("interpolation arrays can not be null");

            if (x.Length != y.Length)
                throw new ArgumentException("interpolation x and y arrays should have the same length");

            if (x.Length < 2)
            {
                throw new ArgumentException("interpolation array should have at least 2 nodes");
            }

            //check that data is ordered           
            for (int i = 0; i < x.Length - 1; i++)
            {
                if (x[i + 1] <= x[i])
                {
                    throw new ArgumentException("interpolation nodes should be ordered");
                }
            }


            //Data is not copied to save memory
            if (bCopyData)
            {
                m_x = (double[])x.Clone();
                m_y = (double[])y.Clone();
            }
            else
            {
                m_x = x;
                m_y = y;
            }            
        }

        public bool Bounds
        {
            set
            {
                m_bounds = value;
            }
            get
            {
                return m_bounds;
            }
        }
        #endregion Constructors

        #region Abstruct Methods

        public abstract double Eval(double x);

        #endregion

        #region Public Methods
        /// <summary>
        /// Computes interpolated value at the node
        /// </summary>
        /// <param name="x">node coordinate</param>
        /// <returns></returns>
        public double this[double x]
        {
            get
            {
                return Eval(x);
            }
        }

        /// <summary>
        /// Evaluates interpolation at x
        /// </summary>
        /// <param name="xi">interpolation nodes</param>
        /// <returns>interpolated function values</returns>
        public virtual double[] Eval(double[] x)
        {
            if (x == null)
            {
                throw new ArgumentException("Interpolation array can not be null");
            }

            double[] y = new double[x.Length];

            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Eval(x[i]);
            }

            return y;
        }

        public virtual void Eval(double[] x, double[] y)
        {
            if (x == null || y == null)
                throw new ArgumentException("Interpolation arrays can not be null");

            if (x.Length != y.Length)
                throw new ArgumentException("interpolation arrays should have the same length");

            for (int i = 0; i < x.Length; i++)
            {
                y[i] = Eval(x[i]);
            }
        }

        public double CheckError()
        {
            double error = 0;

            for (int i = 0; i < m_x.Length; i++)
            {
                error += Utils.Sqr(Eval(m_x[i]) - m_y[i]);
            }
            return System.Math.Sqrt(error);
        }

        /// <summary>
        /// Size of interpolation data set
        /// </summary>
        public int Size
        {
            get
            {
                return m_x.Length;
            }
        }

        #endregion Public Methods

        #region Protected Methods

     
        /// <summary>
        /// Function returns the index i such that x[i-1] < x <= x[i]         
        /// </summary>
        /// <param name="xi"></param>
        /// <returns></returns>
        protected int FindIndex(double x)
        {
            int index = Array.BinarySearch<double>(m_x, x);           

            return (index < 0) ? ~index : index;
        }

        /// <summary>
        /// Finds an index of interpolation region
        /// </summary>
        /// <param name="x"></param>
        /// <param name="value"></param>
        /// <returns>index of interpolation region [1, n-1], index = 0 when x outside of range</returns>
        protected int FindIndex(double x, out double value)
        {
            int index = 0; 

            value = Double.NaN;

            if (x < m_x[0]) //check left boundary
            {
                if (m_bounds)
                {
                    value = m_y[0];
                }
            }
            else if (x > m_x[m_x.Length - 1]) //check right boundary
            {
                if (m_bounds)
                {
                    value = m_y[m_y.Length - 1];
                }
            }
            else
            {
                index = Array.BinarySearch<double>(m_x, x);

                if (index == 0) // x = m_x[0]
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

        protected void Bracket(int index1, out double x0, out double x1, out double y0, out double y1)
        {
            int index0 = index1 - 1;
            
            x0 = m_x[index0];
            x1 = m_x[index1];
            y0 = m_y[index0];
            y1 = m_y[index1];
        }

        #endregion

        #region Static Methods
        /*
        public static double[] Compute(double[] x, double[] y, double[] xi, InterpolationMethod method)
        {
            Interp1DBase interp = null;

            switch (method)
            {
                case InterpolationMethod.Linear:
                    interp = new Interp1DLinear(x, y);
                    break;
                case InterpolationMethod.Akima:
                    interp = new Interp1DAkima(x, y);
                    break;
                case InterpolationMethod.AkimaPeriodic:
                    interp = new Interp1DAkima(x, y, true);
                    break;
                case InterpolationMethod.Nearest:
                    interp = new Interp1DNearest(x, y);
                    break;
                case InterpolationMethod.CubicSpline:
                    interp = new Interp1DCubicSpline(x, y);
                    break;
                case InterpolationMethod.CubicHermiteSpline:
                    interp = new Interp1DCubicHermiteSpline(x, y);
                    break;
                default:
                    interp = new Interp1DLinear(x, y);
                    break;
            }

            return interp.Compute(xi);
        }
        */
        #endregion
    }
}
