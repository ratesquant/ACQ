using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// not an efficient implementation of RBF interpolaton, sutable for small number of points (less than 512)
    /// </summary>
    public class RbfInterpolation : ScatteredInterpolationBase
    {
        private const int m_sizeLimit = 512;

        double[] m_a;
        readonly double[] m_invscale; //one over scale

        RadialBasisFunction m_basisFunction;

        public RbfInterpolation(double[,] x, double[] y, enRadialBasisFunction basisFunction)
            : this(x, y, basisFunction, null, 0.0)
        {
        }
        public RbfInterpolation(double[,] x, double[] y, enRadialBasisFunction basisFunction, double[] scale, double smooth)
            : base(x, y, true)
        {
            int n = x.GetLength(0);
            int dim = x.GetLength(1);

            if (n > m_sizeLimit)
            {
                throw new ArgumentOutOfRangeException(String.Format("RbfInterpolation does not support more than {0} nodes", m_sizeLimit));
            }

            m_basisFunction = RadialBasisFunction.CreateRadialBasisFunction(basisFunction);

            m_invscale = new double[dim];

            //scale is optional argument, use scale of one if missing
            if (scale == null)
            {
                for (int i = 0; i < m_invscale.Length; i++)
                {
                    m_invscale[i] = 1.0;
                }
            }
            else
            {
                for (int i = 0; i < m_invscale.Length; i++)
                {
                    if (i < scale.Length)
                        m_invscale[i] = 1.0 / scale[i];
                    else
                        m_invscale[i] = 1.0;
                }
            }


            ACQ.Math.Linalg.Matrix A = new ACQ.Math.Linalg.Matrix(n + dim + 1, n + dim + 1);
            ACQ.Math.Linalg.Matrix rhs = new ACQ.Math.Linalg.Matrix(n + dim + 1, 1);

            double[] dx = new double[dim];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    for (int k = 0; k < dim; k++)
                    {
                        dx[k] = (x[i, k] - x[j, k])*m_invscale[k];
                    }
                    double r = norm2(dx);
                    double func = m_basisFunction[r];
                    A[i, j] = func;
                    A[j, i] = func;
                }
                A[i, i] = A[i, i] - smooth;
                rhs[i, 0] = y[i];
            }

            //Linear Polynomial part
            for (int i = 0; i < n; i++)
            {
                A[i, n] = 1;
                A[n, i] = 1;

                for (int j = 0; j < dim; j++)
                {
                    A[i, n + 1 + j] = x[i, j] * m_invscale[j];
                    A[n + 1 + j, i] = x[i, j] * m_invscale[j];
                }
            }

            ACQ.Math.Linalg.Matrix res = A.Solve(rhs);

            m_a = new double[n + dim + 1];

            for (int i = 0; i < res.Rows; i++)
            {
                m_a[i] = res[i, 0];
            }
        }

        public override double Eval(double[] x)
        {
            double res = 0;
            int n = m_x.GetLength(0);
            int dim = m_x.GetLength(1);

            if (x == null || x.Length != dim)
            {
                return Double.NaN; // dont throw exception when it is natural to return nan
            }

            double[] dx = new double[dim];

            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < dim; k++)
                {
                    dx[k] = (m_x[i, k] - x[k]) * m_invscale[k];
                }
                res += m_a[i] * m_basisFunction[norm2(dx)];
            }
            res += m_a[n];

            for (int k = 0; k < dim; k++)
            {
                res += m_a[n + 1 + k] * x[k] * m_invscale[k];
            }

            return res;
        }

        protected double norm2(double[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; i++)
            {
                sum += x[i] * x[i];
            }
            return System.Math.Sqrt(sum);

        }
    }
}
