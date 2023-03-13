using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math
{
    #region cephes
    /*
        
        This file is the original README from that page, which also states the license.

           Some software in this archive may be from the book _Methods and
        Programs for Mathematical Functions_ (Prentice-Hall or Simon & Schuster
        International, 1989) or from the Cephes Mathematical Library, a
        commercial product. In either event, it is copyrighted by the author.
        What you see here may be used freely but it comes with no support or
        guarantee.

           The two known misprints in the book are repaired here in the
        source listings for the gamma function and the incomplete beta
        integral.


           Stephen L. Moshier
           moshier@na-net.ornl.gov
     * */
    #endregion
    /// <summary>
    /// Special functions, adapted from Cephes (subset)
    /// https://github.com/jeremybarnes/cephes
    /// http://www.netlib.org/cephes/
    /// function names are the same as Cephes
    /// </summary>
    public class Special
    {
        #region const
        private const double TWOPI = 6.283185307179586476925286766559;
        private const double MACHEP = 1.11022302462515654042E-16;
        private const double MAXLOG = 7.09782712893383996732E2;
        private const double MINLOG = -7.451332191019412076235E2;
        private const double MAXGAM = 171.624376956302725;
        private const double SQTPI = 2.50662827463100050242E0;
        private const double SQRTH = 7.07106781186547524401E-1;
        private const double SQRT2 = 1.41421356237309504880;
        private const double LOGPI = 1.14472988584940017414;
        private const double LS2PI = 0.91893853320467274178;// log( sqrt( 2*pi ) )
        private const double SQT2PI_INV = 0.39894228040143267793994605993438; //1/sqrt(2*pi)

        private const double MAXLGM = 2.556348e305;

        private const double big = 4.503599627370496e15;
        private const double biginv =  2.22044604925031308085e-16;

        private const double s2pi = 2.50662827463100050242E0; //sqrt(2pi) 

        #endregion 

        #region polycoeffs
        private static readonly double[] m_gamma_P = {
          1.60119522476751861407E-4,
          1.19135147006586384913E-3,
          1.04213797561761569935E-2,
          4.76367800457137231464E-2,
          2.07448227648435975150E-1,
          4.94214826801497100753E-1,
          9.99999999999999996796E-1	};
        private static readonly double[] m_gamma_Q = {
        -2.31581873324120129819E-5,
         5.39605580493303397842E-4,
        -4.45641913851797240494E-3,
         1.18139785222060435552E-2,
         3.58236398605498653373E-2,
        -2.34591795718243348568E-1,
         7.14304917030273074085E-2,
         1.00000000000000000320E0};
        private static readonly double[] m_gamma_STIR = {
        7.87311395793093628397E-4,
        -2.29549961613378126380E-4,
        -2.68132617805781232825E-3,
        3.47222221605458667310E-3,
        8.33333333333482257126E-2};
        private static readonly double[] m_gamma_A = {
        8.11614167470508450300E-4,
        -5.95061904284301438324E-4,
        7.93650340457716943945E-4,
        -2.77777777730099687205E-3,
        8.33333333333331927722E-2};
        private static readonly double[] m_gamma_B = {
		-1.37825152569120859100E3,
		-3.88016315134637840924E4,
		-3.31612992738871184744E5,
		-1.16237097492762307383E6,
		-1.72173700820839662146E6,
		-8.53555664245765465627E5};
        private static readonly double[] m_gamma_C = {
		/* 1.00000000000000000000E0, */
		-3.51815701436523470549E2,
		-1.70642106651881159223E4,
		-2.20528590553854454839E5,
		-1.13933444367982507207E6,
		-2.53252307177582951285E6,
		-2.01889141433532773231E6};
        private static readonly double[] m_ndtr_P = {
		2.46196981473530512524E-10,
		5.64189564831068821977E-1,
		7.46321056442269912687E0,
		4.86371970985681366614E1,
		1.96520832956077098242E2,
		5.26445194995477358631E2,
		9.34528527171957607540E2,
		1.02755188689515710272E3,
		5.57535335369399327526E2};
        private static readonly double[] m_ndtr_Q = {
		//1.0
		1.32281951154744992508E1,
		8.67072140885989742329E1,
		3.54937778887819891062E2,
		9.75708501743205489753E2,
		1.82390916687909736289E3,
		2.24633760818710981792E3,
		1.65666309194161350182E3,
		5.57535340817727675546E2};
        private static readonly double[] m_ndtr_R = {
		5.64189583547755073984E-1,
		1.27536670759978104416E0,
		5.01905042251180477414E0,
		6.16021097993053585195E0,
		7.40974269950448939160E0,
		2.97886665372100240670E0};
        private static readonly double[] m_ndtr_S = {
	    //1.00000000000000000000E0, 
	    2.26052863220117276590E0,
	    9.39603524938001434673E0,
	    1.20489539808096656605E1,
	    1.70814450747565897222E1,
	    9.60896809063285878198E0,
	    3.36907645100081516050E0};
        private static readonly double[] m_ndtr_T = {
		9.60497373987051638749E0,
        9.00260197203842689217E1,
        2.23200534594684319226E3,
        7.00332514112805075473E3,
        5.55923013010394962768E4	};
        private static readonly double[] m_ndtr_U = {
		//1.00000000000000000000E0,
		3.35617141647503099647E1,
		5.21357949780152679795E2,
		4.59432382970980127987E3,
		2.26290000613890934246E4,
		4.92673942608635921086E4};
        #endregion

        #region coefs used in BVND
        private static readonly double[] BVND_WN20 = {
            0.01761400713915212, 0.04060142980038694,
            0.06267204833410906, 0.08327674157670475,
            0.1019301198172404,  0.1181945319615184,
            0.1316886384491766,  0.1420961093183821,
            0.1491729864726037,  0.1527533871307259  };

        private static readonly double[] BVND_XN20 = {
            -0.9931285991850949, -0.9639719272779138,
            -0.9122344282513259, -0.8391169718222188,
            -0.7463319064601508, -0.6360536807265150,
            -0.5108670019508271, -0.3737060887154196,
            -0.2277858511416451, -0.07652652113349733};

        private static readonly double[] BVND_WN12 =     {
            0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
            0.2031674267230659,  0.2334925365383547, 0.2491470458134029};

        private static readonly double[] BVND_XN12 ={
            -0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
            -0.5873179542866171, -0.3678314989981802, -0.1252334085114692 };

        private static readonly double[] BVND_WN6 =
        {
            0.1713244923791705,
            0.3607615730481384,
            0.4679139345726904
        };

        private static readonly double[] BVND_XN6 =
        {
            -0.9324695142031522,
            -0.6612093864662647,
            -0.2386191860831970
        };

        private static readonly double[] ndtri_P0 = {
            -5.99633501014107895267E1,
            9.80010754185999661536E1,
            -5.66762857469070293439E1,
            1.39312609387279679503E1,
            -1.23916583867381258016E0,
        };
        private static readonly double[] ndtri_Q0 = {
            /* 1.00000000000000000000E0,*/
            1.95448858338141759834E0,
            4.67627912898881538453E0,
            8.63602421390890590575E1,
            -2.25462687854119370527E2,
            2.00260212380060660359E2,
            -8.20372256168333339912E1,
            1.59056225126211695515E1,
            -1.18331621121330003142E0,
        };

        private static readonly double[] ndtri_P1 = {
            4.05544892305962419923E0,
            3.15251094599893866154E1,
            5.71628192246421288162E1,
            4.40805073893200834700E1,
            1.46849561928858024014E1,
            2.18663306850790267539E0,
            -1.40256079171354495875E-1,
            -3.50424626827848203418E-2,
            -8.57456785154685413611E-4,
        };
        private static readonly double[] ndtri_Q1 = {
            /*  1.00000000000000000000E0,*/
            1.57799883256466749731E1,
            4.53907635128879210584E1,
            4.13172038254672030440E1,
            1.50425385692907503408E1,
            2.50464946208309415979E0,
            -1.42182922854787788574E-1,
            -3.80806407691578277194E-2,
            -9.33259480895457427372E-4,
        };

        private static readonly double[] ndtri_P2 = {
            3.23774891776946035970E0,
            6.91522889068984211695E0,
            3.93881025292474443415E0,
            1.33303460815807542389E0,
            2.01485389549179081538E-1,
            1.23716634817820021358E-2,
            3.01581553508235416007E-4,
            2.65806974686737550832E-6,
            6.23974539184983293730E-9,
        };
        private static readonly double[] ndtri_Q2 = {
            /*  1.00000000000000000000E0,*/
            6.02427039364742014255E0,
            3.67983563856160859403E0,
            1.37702099489081330271E0,
            2.16236993594496635890E-1,
            1.34204006088543189037E-2,
            3.28014464682127739104E-4,
            2.89247864745380683936E-6,
            6.79019408009981274425E-9,
        };

        #endregion


        /// <summary>
        /// chi-square function (left hand tail).
        /// </summary>
        public static double chisq(double df, double x)
        {
            if (x < 0.0 || df < 1.0)
            {
                return 0.0;
            }
            return igam( df/2.0, x/2.0 );
        }
        /// <summary>
        /// chi-square function (right hand tail)
        /// </summary>
        public static double chisqc(double df, double x)
        {
            if (x < 0.0 || df < 1.0) return 0.0;

            return igamc(df / 2.0, x / 2.0);
        }

        /// <summary>
        /// Returns the sum of the first k terms of the Poisson distribution
        /// </summary>
        /// <param name="k"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double poisson(int k, double x)
        {
            if (k < 0 || x < 0) 
                return 0.0;

            return igamc((double)(k + 1), x);
        }

        public static double poissonc(int k, double x)
        {
            if (k < 0 || x < 0) 
                return 0.0;

            return igam((double)(k + 1), x);
        }

        public static double fac(double x)
        {
            if (x < 0.0)
            {
                return Double.NaN;
            }
            else
            {
                return gamma(x + 1.0);
            }
        }

        #region gamma
        /// <summary>
        /// Returns the gamma function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double gamma(double x)
        {
            int sgngam;
            return gamma(x, out sgngam);
        }
        /// <summary>
        /// Returns the gamma function
        /// </summary>
        private static double gamma(double x, out int sgngam)
        {
            sgngam = 1;
            if (Double.IsNaN(x))
            {
                return Double.NaN;
            }

            if (Double.IsPositiveInfinity(x))
            {
                return Double.PositiveInfinity;
            }

            if (Double.IsNegativeInfinity(x))
            {
                return Double.NaN;
            }
            int i;
            double p, z;

            double q = System.Math.Abs(x);

            if (q > 33.0)
            {
                if (x < 0.0)
                {
                    p = System.Math.Floor(q);
                    if (p == q)
                    {
                        return Double.NaN; //gamma: overflow
                    }
                    i = (int)p;
                    if ((i & 1) == 0)
                    {
                        sgngam = -1;
                    }
                    z = q - p;
                    if (z > 0.5)
                    {
                        p += 1.0;
                        z = q - p;
                    }
                    z = q * System.Math.Sin(System.Math.PI * z);
                    if (z == 0.0)
                    {
                        return sgngam > 0 ? Double.PositiveInfinity : Double.NegativeInfinity;//"gamma: overflow";
                    }
                    z = System.Math.Abs(z);
                    z = System.Math.PI / (z * stirf(q));
                }
                else
                {
                    z = stirf(x);
                }
                return sgngam * z;
            }

            z = 1.0;
            while (x >= 3.0)
            {
                x -= 1.0;
                z *= x;
            }

            while (x < 0.0)
            {
                if (x == 0.0)
                {
                    return Double.PositiveInfinity;//"gamma: overflow");
                }
                else if (x > -1.0E-9)
                {
                    return (z / ((1.0 + 0.5772156649015329 * x) * x));
                }
                z /= x;
                x += 1.0;
            }

            while (x < 2.0)
            {
                if (x == 0.0)
                {
                    return Double.PositiveInfinity;//"gamma: singular"
                }
                else if (x < 1.0E-9)
                {
                    return (z / ((1.0 + 0.5772156649015329 * x) * x));
                }
                z /= x;
                x += 1.0;
            }

            if (x == 2.0)
            {
                return z;
            }

            x -= 2.0;
            p = polevl(x, m_gamma_P, 6);
            q = polevl(x, m_gamma_Q, 7);
            return z * p / q;
        }

        /// <summary>
        /// Gamma function computed by Stirling's formula, The polynomial STIR is valid for 33 <= x <= 172.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double stirf(double x)
        {
            const double MAXSTIR = 143.01608;

            double w = 1.0 / x;
            double y = exp(x);

            w = 1.0 + w * polevl(w, m_gamma_STIR, 4);

            if (x > MAXSTIR)
            {
                // Avoid overflow in Math.Pow() 
                double v = pow(x, 0.5 * x - 0.25);
                y = v * (v / y);
            }
            else
            {
                y = pow(x, x - 0.5) / y;
            }
            y = SQTPI * y * w;
            return y;
        }

        /// <summary>
        /// complemented incomplete gamma function.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double igamc(double a, double x)
        {           
            double ans, ax, c, yc, r, t, y, z;
            double pk, pkm1, pkm2, qk, qkm1, qkm2;

            if (x <= 0 || a <= 0)
            {
                return 1.0;
            }

            if (x < 1.0 || x < a)
            {
                return 1.0 - igam(a, x);
            }

            ax = a * log(x) - x - lgam(a);

            if (ax < -MAXLOG)
            {
                return 0.0;//UNDERFLOW 
            }

            ax = exp(ax);

            // continued fraction 
            y = 1.0 - a;
            z = x + y + 1.0;
            c = 0.0;
            pkm2 = 1.0;
            qkm2 = x;
            pkm1 = x + 1.0;
            qkm1 = z * x;
            ans = pkm1 / qkm1;

            do
            {
                c += 1.0;
                y += 1.0;
                z += 2.0;
                yc = y * c;
                pk = pkm1 * z - pkm2 * yc;
                qk = qkm1 * z - qkm2 * yc;
                if (qk != 0)
                {
                    r = pk / qk;
                    t = System.Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                    t = 1.0;

                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (fabs(pk) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
            } while (t > MACHEP);

            return ans * ax;
        }


        /// <summary>
        /// Incomplete gamma integral.
        /// </summary>
        public static double igam(double a, double x)
        {
            double ans, ax, c, r;

            if (x <= 0 || a <= 0)
            {
                return 0.0;
            }

            if (x > 1.0 && x > a)
            {
                return 1.0 - igamc(a, x);
            }

            // Compute  x**a * exp(-x) / gamma(a)  
            ax = a * log(x) - x - lgam(a);
            if (ax < -MAXLOG)
            {
                return (0.0);
            }

            ax = exp(ax);

            // power series 
            r = a;
            c = 1.0;
            ans = 1.0;

            do
            {
                r += 1.0;
                c *= x / r;
                ans += c;
            } while (c / ans > MACHEP);

            return ans * ax / a;
        }

        public static double lgam(double x)
        {
            int sgngam;
            return lgam(x, out sgngam);
        }
        /// <summary>
        /// Logarithm of gamma function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double lgam(double x, out int sgngam)
        {
            double p, q, u, w, z;
            int i;

            sgngam = 1;

            if (Double.IsNaN(x))
            {
                return Double.NaN;
            }

            if (Double.IsInfinity(x))
            {
                return Double.PositiveInfinity; 
            }

            if (x < -34.0)
            {
                q = -x;
                w = lgam(q, out sgngam);
                p = System.Math.Floor(q);
                if (p == q)
                {
                    return Double.PositiveInfinity; //lgam: Overflow
                }

                i = (int)p;
                if ((i & 1) == 0)
                    sgngam = -1;
                else
                    sgngam = 1;

                z = q - p;
                if (z > 0.5)
                {
                    p += 1.0;
                    z = p - q;
                }
                z = q * System.Math.Sin(System.Math.PI * z);
                if (z == 0.0) 
                {
                    return Double.PositiveInfinity;
                }
                z = LOGPI - log(z) - w;
                return z;
            }

            if (x < 13.0)
            {
                z = 1.0;
                p = 0.0;
                u = x;
                while (u >= 3.0)
                {
                    p -= 1.0;
                    u = x + p;
                    z *= u;
                }
                while (u < 2.0)
                {
                    if (x == 0.0)
                        return Double.PositiveInfinity;
                    z /= u;
                    p += 1.0;
                    u = x + p;
                }
                if (z < 0.0)
                {
                    sgngam = -1;
                    z = -z;
                }
                else
                {
                    sgngam = 1;
                }

                if (u == 2.0)
                {
                    return log(z);
                }

                p -= 2.0;
                x = x + p;
                p = x * polevl(x, m_gamma_B, 5) / p1evl(x, m_gamma_C, 6);

                return (log(z) + p);
            }

            if (x > MAXLGM)
            {
                return sgngam == 1 ? Double.PositiveInfinity : Double.NegativeInfinity;
            }

            q = (x - 0.5) * log(x) - x + LS2PI;
            if (x > 1.0e8)
            {
                return (q);
            }

            p = 1.0 / (x * x);
            if (x >= 1000.0)
            {
                q += ((7.9365079365079365079365e-4 * p
                    - 2.7777777777777777777778e-3) * p
                    + 0.0833333333333333333333) / x;
            }
            else
            {
                q += polevl(p, m_gamma_A, 4) / x;
            }
            return q;
        }
 
        /// <summary>
        /// Gamma distribution function
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double gdtr( double a, double b, double x )
        {
            if( x < 0.0 )
            {
                return 0.0;
            }
            return igam( b, a * x );
        }
        /// <summary>
        /// Returns the integral from x to infinity of the gamma probability density function:
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double gdtrc( double a, double b, double x )
        {
            if( x < 0.0 )
            {
                return 0.0;
            }
            return igamc( b, a * x );
        }
        #endregion 

        #region normal
        /// <summary>
        /// Returns the area under the Gaussian probability density, integrated from minus infinity to x.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double NormalCdf(double x)
        {
            return ndtr(x); 
        }

        public static double InverseNormalCdf(double x)
        {
            return ndtri(x);
        }
        public static double NormalPdf(double x)
        {
            return Special.SQT2PI_INV * System.Math.Exp(-0.5* x * x);
        }

        /// <summary>
        /// The cumulative bivariate normal distribution function,
        /// his function is based on the method described by
        ///  Drezner, Z and G.O. Wesolowsky, (1990),
        ///  On the computation of the bivariate normal integral,
        ///  Journal of Statist. Comput. Simul. 35, pp. 101-107,
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="rho"></param>
        /// <returns></returns>
        public static double CBND(double x, double y, double rho)
        {
            return BVND(-x, -y, rho);
        }
        /// <summary>
        ///  A function for computing bivariate normal probabilities. 
        ///   BVND calculates the probability that X > DH and Y > DK.
        ///   
        /// This method is based on the work done by Alan Genz, Department of 
        ///   Mathematics, Washington State University. Pullman, WA 99164-3113
        ///   Email: alangenz@wsu.edu. This work was shared under a 3-clause BSD
        ///   license. Please see source file for more details and the actual
        ///   license text.
        ///   
        /// </summary>
        /// <param name="dh"></param>
        /// <param name="dk"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static double BVND(double dh, double dk, double r)
        {           
            double[] x;
            double[] w;

            if (System.Math.Abs(r) < 0.3)
            {
                // Gauss Legendre Points and Weights N =  6
                x = BVND_XN6;
                w = BVND_WN6;
            }
            else if (System.Math.Abs(r) < 0.75)
            {
                // Gauss Legendre Points and Weights N =  12
                x = BVND_XN12;
                w = BVND_WN12;
            }
            else
            {
                // Gauss Legendre Points and Weights N =  20
                x = BVND_XN20;
                w = BVND_WN20;
            }

            double h = dh;
            double k = dk;
            double hk = h * k;
            double bvn = 0;

            if (System.Math.Abs(r) < 0.925)
            {
                if (System.Math.Abs(r) > 0)
                {
                    double sh = (h * h + k * k) / 2;
                    double asr = System.Math.Asin(r);

                    for (int i = 0; i < x.Length; i++)
                    {
                        for (int j = -1; j <= 1; j += 2)
                        {
                            double sn = System.Math.Sin(asr * (j * x[i] + 1) / 2);
                            bvn = bvn + w[i] * System.Math.Exp((sn * hk - sh) / (1 - sn * sn));
                        }
                    }
                    bvn = bvn * asr / (2 * TWOPI);
                }

                return bvn + NormalCdf(-h) * NormalCdf(-k);
            }


            if (r < 0)
            {
                k = -k;
                hk = -hk;
            }

            if (System.Math.Abs(r) < 1)
            {
                double sa = (1 - r) * (1 + r);
                double A = System.Math.Sqrt(sa);
                double sb = (h - k);
                sb = sb * sb;
                double c = (4 - hk) / 8;
                double d = (12 - hk) / 16;
                double asr = -(sb / sa + hk) / 2;

                if (asr > -100)
                    bvn = A * System.Math.Exp(asr) * (1 - c * (sb - sa) * (1 - d * sb / 5) / 3 + c * d * sa * sa / 5);

                if (-hk < 100)
                {
                    double B = System.Math.Sqrt(sb);
                    bvn = bvn - System.Math.Exp(-hk / 2) * System.Math.Sqrt(TWOPI) * NormalCdf(-B / A) * B
                              * (1 - c * sb * (1 - d * sb / 5) / 3);
                }

                A = A / 2;

                for (int i = 0; i < x.Length; i++)
                {
                    for (int j = -1; j <= 1; j += 2)
                    {
                        double xs = (A * (j * x[i] + 1));
                        xs = xs * xs;
                        double rs = System.Math.Sqrt(1 - xs);
                        asr = -(sb / xs + hk) / 2;

                        if (asr > -100)
                        {
                            bvn = bvn + A * w[i] * System.Math.Exp(asr)
                                * (System.Math.Exp(-hk * xs / (2 * (1 + rs) * (1 + rs))) / rs
                                - (1 + c * xs * (1 + d * xs)));
                        }
                    }
                }

                bvn = -bvn / TWOPI;
            }

            if (r > 0)
                return bvn + NormalCdf(-System.Math.Max(h, k));

            bvn = -bvn;

            if (k <= h)
                return bvn;

            if (h < 0)
                return bvn + NormalCdf(k) - NormalCdf(h);

            return bvn + NormalCdf(-h) - NormalCdf(-k);
        }


        /// <summary>
        /// Returns the area under the Gaussian probability density, integrated from minus infinity to a.
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double ndtr(double a)
        {
            double x, y, z;

            x = a * SQRTH;
            z = fabs(x);

            if (z < 1.0) //if( z < SQRTH  ) ?
            {
                y = 0.5 + 0.5 * erf(x);
            }
            else
            {
                y = 0.5 * erfc(z);
                if (x > 0)
                {
                    y = 1.0 - y;
                }
            }

            return y;
        }

        /// <summary>
        /// Inverse of Normal distribution function
        /// </summary>
        /// <param name="y0"></param>
        /// <returns></returns>
        public static  double ndtri(double y0)
        {
            double x, y, z, y2, x0, x1;
            int code;
            
            if(y0 <= 0.0 )
            {                
                return Double.NegativeInfinity;
            }
            if(y0 >= 1.0 )
            {
                
                return Double.PositiveInfinity;
            }
            
            code = 1;
            y = y0;
            if (y > (1.0 - 0.13533528323661269189)) //0.135... = exp(-2) 
            {
                y = 1.0 - y;
                code = 0;
            }
            
            if (y > 0.13533528323661269189)
            {
                y = y - 0.5;
                y2 = y * y;
                x = y + y * (y2 * polevl(y2, ndtri_P0, 4) / p1evl(y2, ndtri_Q0, 8));
                x = x * s2pi;
                return (x);
            }
            
            x = sqrt(-2.0 * log(y));
            x0 = x - log(x) / x;
            z = 1.0 / x;
            if (x < 8.0) /* y > exp(-32) = 1.2664165549e-14 */
                x1 = z * polevl(z, ndtri_P1, 8) / p1evl(z, ndtri_Q1, 8);
            else
                x1 = z * polevl(z, ndtri_P2, 8) / p1evl(z, ndtri_Q2, 8);
            x = x0 - x1;
            if (code != 0)
                x = -x;
            return (x);
        }

        /// <summary>
        /// Complementary error function
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double erfc(double a)
        {
            double x, y, z, p, q;

            if (a < 0.0) 
                x = -a;
            else 
                x = a;

            if (x < 1.0)
            {
                return 1.0 - erf(a);
            }

            z = -a * a;

            if (z < -MAXLOG)
            {
                if (a < 0) 
                    return 2.0;
                else 
                    return 0.0;
            }

            z =  exp(z); // Special.expx2(a, -1);

            if (x < 8.0)
            {
                p = polevl(x, m_ndtr_P, 8);
                q = p1evl(x, m_ndtr_Q, 8);
            }
            else
            {
                p = polevl(x, m_ndtr_R, 5);
                q = p1evl(x, m_ndtr_S, 6);
            }

            y = (z * p) / q;

            if (a < 0)
            {
                y = 2.0 - y;
            }

            if (y == 0.0)
            {
                if (a < 0) 
                    return 2.0;
                else 
                    return 0.0;
            }


            return y;
        }

        /// <summary>
        ///  Error function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double erf(double x)
        {
            double y, z;

            if (fabs(x) > 1.0)
            {
                return (1.0 - erfc(x));
            }
            z = x * x;
            y = x * polevl(z, m_ndtr_T, 4) / p1evl(z, m_ndtr_U, 5);
            return y;
        }

        /// <summary>
        /// Computes y = exp(x*x)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="sign"></param>
        /// <returns></returns>
        private static double expx2 (double x, int sign)
        {
            const double M = 128.0;
            const double MINV = 0.0078125;

            double u, u1, m, f;
            
            x = fabs (x);

            if (sign < 0)
            {
                x = -x;
            }
            
            // Represent x as an exact multiple of M plus a residual.
            //  M is a power of 2 chosen so that exp(m * m) does not overflow
            //  or underflow and so that |x - m| is small.  
            m = MINV * System.Math.Floor(M * x + 0.5);
            f = x - m;

            // x^2 = m^2 + 2mf + f^2 
            u = m * m;
            u1 = 2 * m * f  +  f * f;

            if (sign < 0)
            {
                u = -u;
                u1 = -u1;
            }

            if ((u + u1) > MAXLOG)
            {
                return System.Double.PositiveInfinity;
            }

            // u is exact, u1 is small.  
            u = exp(u) * exp(u1);
            return(u);
        }

        private static readonly double[] a_coefs_BSM = new double[] { -25.44106049637, 41.39119773534, -18.61500062529, 2.50662823884 };
        private static readonly double[] b_coefs_BSM = new double[] { 3.13082909833, -21.06224101826, 23.08336743743, -8.47351093090, 1 };
        private static readonly double[] c_coefs_BSM = new double[] { 0.0000003960315187, 0.0000002888167364, 0.0000321767881768, 0.0003951896511919, 0.0038405729373609, 0.0276438810333863, 0.1607979714918209, 0.9761690190917186, 0.3374754822726147 };

        /// <summary>
        /// % Beasley-Springer-Moro algorithm for approximating the inverse normal.
        /// Reference:
        /// Paul Glasserman, Monte Carlo methods in financial engineering, vol. 53 of Applications of Mathematics(New York),  Springer-Verlag, new York, 2004, p.67-68
        /// </summary>
        /// <returns></returns>
        public static double InverseNormalCdfBSM(double u)
        {
            if (u < 0 || u > 1)
            {
                return Double.NaN;
            }
            if (u == 0.0)
            {
                return Double.NegativeInfinity;
            }
            if (u == 1.0)
            {
                return Double.PositiveInfinity;
            }

            double y = u - 0.5;
            double abs_y = System.Math.Abs(y);
            double res;

            if (abs_y <= 0.42)
            {
                double r = abs_y * abs_y;
                res = y * polevl(r, a_coefs_BSM) / polevl(r, b_coefs_BSM);
            } else
            {               
                double r = System.Math.Min(u, 1 - u);
                r = System.Math.Log(-System.Math.Log(r));
                res = polevl(r, c_coefs_BSM) * System.Math.Sign(y);                
            }
            return res;

        }
        #endregion 

        #region incbet
        /// <summary>
        /// Incomplete beta integral
        /// </summary>
        /// <param name="aa"></param>
        /// <param name="bb"></param>
        /// <param name="xx"></param>
        /// <returns></returns>
        public static double incbet(double aa, double bb, double xx )
        {
            double a, b, t, x, xc, w, y;
            int flag;

            if( aa <= 0.0 || bb <= 0.0 )
                return( 0.0 );

            if( (xx <= 0.0) || ( xx >= 1.0) )
            {
                if( xx == 0.0 )
                    return(0.0);
                if( xx == 1.0 )
                    return( 1.0 );
                return( 0.0 );
            }
            
            flag = 0;
            if( (bb * xx) <= 1.0 && xx <= 0.95)
            {
                t = pseries(aa, bb, xx);
                return t;
            }

            w = 1.0 - xx;

            // Reverse a and b if x is greater than the mean. 
            if( xx > (aa/(aa+bb)) )
            {
                flag = 1;
                a = bb;
                b = aa;
                xc = xx;
                x = w;
            }
            else
            {
                a = aa;
                b = bb;
                xc = w;
                x = xx;
            }

            if( flag == 1 && (b * x) <= 1.0 && x <= 0.95)
            {
                t = pseries(a, b, x);
                if( t <= MACHEP )
                    t = 1.0 - MACHEP;
                else
                    t = 1.0 - t;
                return( t );
            }

            // Choose expansion for better convergence. 
            y = x * (a+b-2.0) - (a-1.0);
            if( y < 0.0 )
                w = incbcf( a, b, x );
            else
                w = incbd( a, b, x ) / xc;
            
            // Multiply w by the factor
            //a      b   _             _     _
            //x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   
            y = a * log(x);
            t = b * log(xc);
            if( (a+b) < MAXGAM && fabs(y) < MAXLOG && fabs(t) < MAXLOG )
            {
                t = pow(xc,b);
                t *= pow(x,a);
                t /= a;
                t *= w;
                t *= gamma(a+b) / (gamma(a) * gamma(b));
                
                if (flag == 1)
                {
                    if (t <= MACHEP)
                        t = 1.0 - MACHEP;
                    else
                        t = 1.0 - t;
                }
                return (t);

            }
            // Resort to logarithms. 
            y += t + lgam(a+b) - lgam(a) - lgam(b);
            y += log(w/a);
            if( y < MINLOG )
                t = 0.0;
            else
                t = exp(y);
            
            if( flag == 1 )
            {
                if( t <= MACHEP )
                    t = 1.0 - MACHEP;
                else
                    t = 1.0 - t;
            }
            return( t );
        }

        /// <summary>
        /// Continued fraction expansion #1 for incomplete beta integral
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double incbcf(double a,double b, double x )
        {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, thresh;
            int n;

            k1 = a;
            k2 = a + b;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = b - 1.0;
            k7 = k4;
            k8 = a + 2.0;

            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            ans = 1.0;
            r = 1.0;
            n = 0;
            thresh = 3.0 * MACHEP;
            do
            {	
                xk = -( x * k1 * k2 )/( k3 * k4 );
                pk = pkm1 +  pkm2 * xk;
                qk = qkm1 +  qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                xk = ( x * k5 * k6 )/( k7 * k8 );
                pk = pkm1 +  pkm2 * xk;
                qk = qkm1 +  qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                if( qk != 0 )
                    r = pk/qk;
                if( r != 0 )
                {
                    t = fabs( (ans - r)/r );
                    ans = r;
                }
                else
                    t = 1.0;

                if( t < thresh )
                   return ans;
                
                k1 += 1.0;
                k2 += 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 -= 1.0;
                k7 += 2.0;
                k8 += 2.0;

                if( (fabs(qk) + fabs(pk)) > big )
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
                if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
                {
                    pkm2 *= big;
                    pkm1 *= big;
                    qkm2 *= big;
                    qkm1 *= big;
                }
            }
            while( ++n < 300 );            
            return(ans);
        }
        /// <summary>
        /// Continued fraction expansion #2 for incomplete beta integral
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double incbd(double a, double b, double x )
        {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, z, thresh;
            int n;

            k1 = a;
            k2 = b - 1.0;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = a + b;
            k7 = a + 1.0;
            k8 = a + 2.0;

            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            z = x / (1.0-x);
            ans = 1.0;
            r = 1.0;
            n = 0;
            thresh = 3.0 * MACHEP;
            do
            {
                xk = -( z * k1 * k2 )/( k3 * k4 );
                pk = pkm1 +  pkm2 * xk;
                qk = qkm1 +  qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                xk = ( z * k5 * k6 )/( k7 * k8 );
                pk = pkm1 +  pkm2 * xk;
                qk = qkm1 +  qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                if( qk != 0 )
                    r = pk/qk;
                if( r != 0 )
                {
                    t = fabs( (ans - r)/r );
                    ans = r;
                }
                else
                    t = 1.0;

                if( t < thresh )
                    return (ans);

                k1 += 1.0;
                k2 -= 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 += 1.0;
                k7 += 2.0;
                k8 += 2.0;

                if( (fabs(qk) + fabs(pk)) > big )
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
                if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
                {
                    pkm2 *= big;
                    pkm1 *= big;
                    qkm2 *= big;
                    qkm1 *= big;
                }
            }
            while( ++n < 300 );
            
            return(ans);
        }      
        #endregion

        #region private functions
        private static double fabs(double x)
        {
            return System.Math.Abs(x);
        }
        private static double sqrt(double x)
        {
            return System.Math.Sqrt(x);
        }
        private static double exp(double x)
        {
            return System.Math.Exp(x);
        }
        private static double log(double x)
        {
            return System.Math.Log(x);
        }
        private static double pow(double x, double a)
        {            
            return System.Math.Pow(x, a);
        }
        ///Power series for incomplete beta integral.  Use when b*x is small and x not too close to 1.
        private static double pseries(double a, double b, double x)
        {
            double s, t, u, v, n, t1, z, ai;
            ai = 1.0 / a;
            u = (1.0 - b) * x;
            v = u / (a + 1.0);
            t1 = v;
            t = u;
            n = 2.0;
            s = 0.0;
            z = MACHEP * ai;
            while (fabs(v) > z)
            {
                u = (n - b) * x / n;
                t *= u;
                v = t / (a + n);
                s += v;
                n += 1.0;
            }
            s += t1;
            s += ai;

            u = a * log(x);
            if ((a + b) < MAXGAM && fabs(u) < MAXLOG)
            {
                t = gamma(a + b) / (gamma(a) * gamma(b));
                s = s * t * pow(x, a);
            }
            else
            {
                t = lgam(a + b) - lgam(a) - lgam(b) + u + log(s);
                if (t < MINLOG)
                    s = 0.0;
                else
                    s = exp(t);
            }
            return (s);
        }
        /// <summary>
        /// Evaluates polynomial of degree n
        /// </summary>
        private static double polevl(double x, double[] coef, int n)
        {
            double ans;

            ans = coef[0];

            for (int i = 1; i <= n; i++)
            {
                ans = ans * x + coef[i];
            }

            return ans;
        }

        private static double polevl(double x, double[] coef)
        {
            double ans;

            ans = coef[0];

            for (int i = 1; i < coef.Length; i++)
            {
                ans = ans * x + coef[i];
            }

            return ans;
        }

        /// <summary>
        /// Evaluates polynomial of degree n ( assumtion that coef[n] = 1.0)
        /// </summary>	
        private static double p1evl(double x, double[] coef, int n)
        {
            double ans;

            ans = x + coef[0];

            for (int i = 1; i < n; i++)
            {
                ans = ans * x + coef[i];
            }

            return ans;
        }
        #endregion
    }
}
