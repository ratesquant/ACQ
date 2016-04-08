using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace ACQ.Math.Interpolation
{

    /// <summary>
    /// not an efficient implementation of RBF interpolaton, sutable for small number of points (less than 512)
    /// </summary>
    public class RbfInterpolation1D : InterpolationBase
    {
        private const int m_sizeLimit = 512;

        double[] m_a;
        readonly double m_scale; //one over basis scale squared

        RadialBasisFunction m_basisFunction;

        public RbfInterpolation1D(double[] x, double[] y, enRadialBasisFunction basisFunction)
            : this(x, y, basisFunction, 1.0, 0.0)
        {
        }
        public RbfInterpolation1D(double[] x, double[] y, enRadialBasisFunction basisFunction, double scale, double smooth)
            : base(x, y, false)
        {
            m_basisFunction = RadialBasisFunction.CreateRadialBasisFunction(basisFunction);
            m_scale = scale;

            int n = x.Length;

            if (n > m_sizeLimit)
            {
                throw new ArgumentOutOfRangeException(String.Format("RbfInterpolation1D does not support more than {0} nodes", m_sizeLimit));
            }

            ACQ.Math.Linalg.Matrix A = new ACQ.Math.Linalg.Matrix(n + 2, n + 2);
            ACQ.Math.Linalg.Matrix rhs = new ACQ.Math.Linalg.Matrix(n + 2, 1);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double r = System.Math.Abs((x[i] - x[j])/ m_scale);
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
                A[i, n + 1] = x[i] / m_scale;

                A[n,   i] = 1;
                A[n + 1, i] = x[i] / m_scale;
            }

            ACQ.Math.Linalg.Matrix res = A.Solve(rhs);

            m_a = new double[n+2];

            for (int i = 0; i < res.Rows; i++)
            {
                m_a[i] = res[i, 0];
            }
        }  
        
        public override double Eval(double x)
        {   
            double res = 0;
            int n = m_x.Length;

            for (int i = 0; i < n; i++)
            {
                res += m_a[i] * m_basisFunction[System.Math.Abs((x - m_x[i])/m_scale)];
            }
            res += m_a[n];
            res += m_a[n + 1] * x / m_scale;

            return res;            
        }
    }
}
