using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    class RbfInterpolation
    {
    }

    /*
     * 
     * public class Interp1DRBF : Interp1DBase
    {
        public enum enRadialBasisFunction
        {
            Linear,
            Cubic,
            Multiquadrics,
            Gaussian,
            Thinplate,
            InverseQuadratic,                        
            InverseMultiquadric
        }

        delegate double BasisFunction(double distance);

        BasisFunction m_basisFunction;

        double[] m_a;
        readonly double m_basisScale; //one over basis scale squared

        public Interp1DRBF(double[] x, double[] y, enRadialBasisFunction basisFunction)
            : this(x, y, basisFunction, 1.0, 0)
        {
        }
        public Interp1DRBF(double[] x, double[] y, enRadialBasisFunction basisFunction, double basisScale, double smooth)
            : base(x, y)
        {
            m_basisFunction = GetBasisFunction(basisFunction);
            m_basisScale = 1 / (basisScale * basisScale);

            int n = x.Length;

            MathLib.Matrix A = new Matrix(n + 2, n + 2);
            MathLib.Matrix rhs = new Matrix(n + 2, 1);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double r = Math.Abs(x[i] - x[j]);
                    double func = m_basisFunction(r);
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
                A[i, n+1] = x[i];

                A[n,   i] = 1;
                A[n+1, i] = x[i];
            }

            Matrix res = A.Solve(rhs);

            m_a = new double[n+2];

            for (int i = 0; i < res.Rows; i++)
            {
                m_a[i] = res[i, 0];
            }
        }  
        
        public override double this[double xi]
        {
            get
            {
                double res = 0;
                int n = m_x.Length;

                for (int i = 0; i < n; i++)
                {
                    res += m_a[i] * m_basisFunction(Math.Abs(xi - this.m_x[i]));
                }
                res += m_a[n];
                res += m_a[n + 1] * xi;

                return res;
            }
        }
        #region Basis Functions
        BasisFunction GetBasisFunction(enRadialBasisFunction basisFunction)
        {
            switch (basisFunction)
            {
                case enRadialBasisFunction.Linear:
                    return LinearBasis;
                case enRadialBasisFunction.InverseQuadratic:
                    return InverseQuadraticBasis;
                case enRadialBasisFunction.Cubic:
                    return CubicBasis;
                case enRadialBasisFunction.Gaussian:
                    return GaussianBasis;
                case enRadialBasisFunction.Multiquadrics:
                    return MultiquadricsBasis;
                case enRadialBasisFunction.Thinplate:
                    return ThinplatesBasis;
                case enRadialBasisFunction.InverseMultiquadric:
                    return InverseMultiquadricBasis;
            }
            return null;
        }

        double LinearBasis(double r)
        {
            return r;
        }        

        double CubicBasis(double r)
        {
            return r * r * r;
        }

        double GaussianBasis(double r)
        {
            return Math.Exp(-r * r * m_basisScale);
        }

        double MultiquadricsBasis(double r)
        {
            return Math.Sqrt(1 + r * r * m_basisScale);
        }
        double InverseQuadraticBasis(double r)
        {
            return 1/(1 + r * r * m_basisScale);
        }
        double InverseMultiquadricBasis(double r)
        {
            return 1.0 / MultiquadricsBasis(r);
        }

        double ThinplatesBasis(double r)
        {
            if (r < 1e-12)
                return 0;
            else
                return r * r * Math.Log(r); 
        }

        #endregion

      
    }
     * */
}
