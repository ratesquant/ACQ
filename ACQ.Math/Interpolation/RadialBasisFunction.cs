using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
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

    public class RadialBasisFunction
    {
        public delegate double BasisFunction(double distance);

        BasisFunction m_basis;

        public RadialBasisFunction(BasisFunction basis)
        {
            m_basis = basis;
        }

        public double this[double d]
        {
            get
            {
                return m_basis(d);
            }
        }

        public static RadialBasisFunction CreateRadialBasisFunction(enRadialBasisFunction basisFunction)
        {
            RadialBasisFunction rbf = new RadialBasisFunction(GetBasisFunction(basisFunction));

            return rbf;
        }


        #region Basis Functions
        private static BasisFunction GetBasisFunction(enRadialBasisFunction basisFunction)
        {
            //TODO: use reflection 
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

        static double LinearBasis(double r)
        {
            return r;
        }

        static double CubicBasis(double r)
        {
            return r * r * r;
        }

        static double GaussianBasis(double r)
        {
            return System.Math.Exp(-r * r);
        }

        static double MultiquadricsBasis(double r)
        {
            return System.Math.Sqrt(1 + r * r);
        }
        static double InverseQuadraticBasis(double r)
        {
            return 1 / (1 + r * r);
        }
        static double InverseMultiquadricBasis(double r)
        {
            return 1.0 / MultiquadricsBasis(r);
        }

        static double ThinplatesBasis(double r)
        {
            if (r < 1e-12)
                return 0;
            else
                return r * r * System.Math.Log(r);
        }

        #endregion

      
    }
}
