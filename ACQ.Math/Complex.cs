using System;

namespace ACQ.Math
{
	/// <summary>
	/// Struct Representing A Complex Number
	/// </summary>
	public struct Complex 
	{
		#region Private Fields
		private readonly double m_real;
		private readonly double m_imag;
		#endregion
		#region Constructors
	
        public Complex(double real, double imag) 
		{ 
			m_real = real; 
			m_imag = imag; 
		}

		public Complex(double[] array)
		{
			if ((array == null) || (array.Length != 2))
				throw new ArgumentException("Input array is not a complex number.", "array");

			m_real = array[0];
			m_imag = array[1];
		}
		#endregion
		#region Polar Properties and Constructor
		/// <summary>
		/// Static Constructor method for constructing a complex number from polar coordinates
		/// </summary>
		/// <param name="mod">Modulus of the complex number</param>
		/// <param name="arg">Theta in the Argand representation of a complex number</param>
		/// <returns>A complex number</returns>
		public static Complex ComplexPolar(double mod, double arg)
		{
            return new Complex(mod * System.Math.Cos(arg), mod * System.Math.Sin(arg));
		}

		/// <summary>
		/// Returns the Argand angle (or theta in the Euler equation) of this complex number
		/// </summary>
		public double Arg
		{
			get
			{
				// this check is required to capture the right phase of the wave
				if (this.real >0)
                    return System.Math.Atan(this.imag / this.real);
				else
                    return System.Math.Atan(this.imag / this.real) + System.Math.PI;
			}
		}
		/// <summary>
		/// Returns the modulus of this complex number
		/// </summary>
		public double R
		{
			get
			{
				return this.Mod;
			}
		}
		
		/// <summary>
		/// Returns the polar coordinates of this complex number as an array (r,theta) where c = r*exp(theta)
		/// </summary>
		public double[] ToPolar
		{
			get
			{
				return new double[] {this.Mod, this.Arg};
			}
		}

		#endregion
		#region Public Properties
		/// <summary>
		/// The Real Part
		/// </summary>
		public double real
		{
			get { return m_real; }			
		}

		/// <summary>
		/// The Imaginary Part
		/// </summary>
		public double imag
		{
			get { return m_imag; }
		}

		/// <summary>
		/// Returns the modulus of this complex number
		/// Mod = Sqrt(real ^ 2 + imag ^ 2)
		/// </summary>
		public double Mod
		{
			get
			{
                return System.Math.Sqrt(m_real * m_real + m_imag * m_imag);
			}
		}
	
		/// <summary>
		/// Returns the complex conjugate of this complex number
		/// </summary>
		public Complex Conjugate
		{
			get
			{
				return new Complex(this.real, -this.imag);
			}
		}

		/// <summary>
		/// Gets the e^c
		/// </summary>
		public Complex Exp
		{
			get
			{
                double e = System.Math.Exp(this.real);
                return new Complex(e * System.Math.Cos(this.imag), e * System.Math.Sin(this.imag));
			}
		}

		/// <summary>
		/// Get the log(c)
		/// </summary>
		public Complex Log
		{
			get
			{
                return new Complex(System.Math.Log(this.Mod), this.Arg);
			}
		}

		/// <summary>
		/// Get The Sqrt(c)
		/// </summary>
		public Complex Sqrt
		{
			get
			{
				return this.Pow(0.5);
			}
		}
		#endregion
		#region Casting operators
		/// <summary>
		/// Cast to real (double)
		/// </summary>
		/// <param name="complex"></param>
		/// <returns>The real part of the complex number</returns>
		public static explicit operator double(Complex complex)
		{
			return complex.real;
		}

		/// <summary>
		/// Cast to a double[] (real, imag)
		/// </summary>
		/// <param name="complex"></param>
		/// <returns></returns>
		public static explicit operator double[](Complex complex)
		{
			return new double[] {complex.real, complex.imag};
		}

		/// <summary>
		/// Construct using a real number
		/// </summary>
		/// <param name="real">Real number</param>
		/// <returns>A complex number</returns>
		public static implicit operator Complex(double real)
		{
			return new Complex(real,0d);
		}
		#endregion	
		#region Basic unary operators
		/// <summary>
		/// Returns the Complex Number
		/// </summary>
		/// <param name="a"></param>
		/// <returns></returns>
		public static Complex operator + (Complex a)
		{
			return a;
		}

		/// <summary>
		/// Subtracts two complex numbers
		/// </summary>
		/// <param name="a"></param>
		/// <returns></returns>
		public static Complex operator - (Complex a)
		{
			return new Complex(-a.real, -a.imag);
		}
		
		/// <summary>
		/// Returns the conjugate of this complex number
		/// </summary>
		/// <param name="rhs"></param>
		/// <returns></returns>
		public static Complex operator ~ (Complex rhs)
		{
			return rhs.Conjugate;
		}
		#endregion
		#region  Basic binary operators for addition, subtraction, multiplication, and division.
		/// <summary>
		/// There are some things in life that don't require commenting
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <returns></returns>
		public static Complex operator + (Complex a, Complex b)
		{
			return new Complex(a.real + b.real, a.imag + b.imag);
		}
		/// <summary>
		/// There are some things in life that don't require commenting
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <returns></returns>
		public static Complex operator - (Complex a, Complex b)
		{
			return new Complex(a.real - b.real, a.imag - b.imag);
		}
		/// <summary>
		/// There are some things in life that don't require commenting
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <returns></returns>
		public static Complex operator * (Complex a, Complex b)
		{
			return new Complex(a.real*b.real - a.imag*b.imag, a.real*b.imag + a.imag*b.real);
		}
		/// <summary>
		/// There are some things in life that don't require commenting
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <returns></returns>
		public static Complex operator / (Complex a, Complex b)
		{
			double r1 = a.real, r2 = b.real, i1 = a.imag, i2 = b.imag;
	
			/* Complex division using a trick to avoid possible overflows, underflows, and loss of precision */

            if (System.Math.Abs(r2) >= System.Math.Abs(i2))
			{
				double d = i2/r2;
				Complex c = new Complex( r1 + i1*d , i1 - r1*d);
				c = c * (1/ (r2 + i2*d));
				return c;
			}
			else
			{
				double d = r2/i2;
				Complex c = new Complex( i1 + r1*d , i1*d - r1 );
				return c *= 1/(i2 + r2*d);
			}
		}

		#endregion
		#region Complex Functions
		/// <summary>
		/// Returns complex number raised to the power k
		/// </summary>
		/// <param name="x"></param>
		/// <param name="k"></param>
		/// <returns></returns>
		public static Complex Pow(Complex x, double k)
		{
			double r = System.Math.Pow(x.R, k);
			double arg = x.Arg;

            return new Complex(r * System.Math.Cos(arg * k), r * System.Math.Sin(arg * k));
		}
		
		/// <summary>
		/// Returns this complex number to the power of k
		/// </summary>
		/// <param name="k">The power to raise the complex number</param>
		/// <returns>This complex number raised to the power of k</returns>
		public Complex Pow(double k)
		{
			return Complex.Pow(this,k);
		}
		#endregion
		#region ToString overrides
		public static Complex Parse(string complex_number)
		{
			/* check that the number is in one of the following formats:

				a + ib
				a + bi
				a + i*b
				a + b*i
				all of the above with j instead of i
				all of the above with brackets round them (e.g. (a + bj))
			
				r,arg (polar format, with and without brackets)
			*/
			throw new NotImplementedException("Not yet implemented!");
		}
	
		public override string ToString() 
		{
			if (imag == 0)
				return this.real.ToString();
			else
				return String.Format("({0}+{1}i)", real, imag);
		}

		public string ToStringPolar() 
		{
			return String.Format("Mod({0}) Arg({1})", this.Mod.ToString(), this.Arg.ToString());
		}

		public string ToStringScalar() 
		{
			return String.Format("{0}+{1}j", real, imag);
		}
        public static bool operator ==(Complex c1, Complex c2)
        {
            return (c1.m_real == c2.m_real) && (c1.m_imag == c2.m_imag);
        }

        public static bool operator !=(Complex c1, Complex c2)
        {
            return !(c1==c2);
        }

        public override int GetHashCode()
        {
            return (m_real.GetHashCode() ^ m_imag.GetHashCode());
        }

        public override bool Equals(object o)
        {
            return (o.GetType() == this.GetType()) ? (this == (Complex)o) : false;
        }
		#endregion
	}
}