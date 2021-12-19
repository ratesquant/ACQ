using System;
using System.Collections.Generic;
using System.Text;

using static System.Math;

namespace ACQ.Math.Roots
{
    //One dimensional Root-Finding   
    /// 1. |a - b| < AbsTol + RelTol min(|a|,|b|)
    /// 2. |x_1 - x_0| < AbsTol + RelTol |x_1|
    /// 3. |f| < AbsTol 
    public enum ConvergenceTest { Interval, Delta, Residual};
    public enum RootAlgorithm { Bisection, FalsePosition, Brent, Secant, Ridders, Steffenson };
    public delegate double obj_function(double x);

    public class IterationResults
    {
        public IterationResults(int nIterations, double x_lower, double x_upper, double root, double f_res)
        {
            m_nIterations = nIterations;
            m_x_lower     = x_lower;
            m_x_upper     = x_upper;
            m_root        = root;
            m_f_res       = f_res;
        }
        int m_nIterations;
        double m_x_lower;
        double m_x_upper;
        double m_root;
        double m_f_res;

        public int Iterations
        {
            get { return m_nIterations; }
        }

        public double Residual
        {
            get { return m_f_res; }
        }
        public double Root
        {
            get { return m_root; }
        }
        public double LeftBracket
        {
            get { return m_x_lower; }
        }
        public double RightBracket
        {
            get { return m_x_upper; }
        }
    }

    public class Root
    {
        public static double Find(obj_function funct, double x_lower, double x_upper)
        {
            return Find(funct, x_lower, x_upper, RootAlgorithm.Brent).Root;
        }
        public static IterationResults Find(obj_function funct, double x_lower, double x_upper, RootAlgorithm algorithm)
        {
            RootBase root;

            switch (algorithm)
            {
                case RootAlgorithm.Bisection: 
                    root = new Bisection();
                    break;
                case RootAlgorithm.Brent:
                    root = new Brent();
                    break;
                case RootAlgorithm.FalsePosition:
                    root = new FalsePosition();
                    break;
                case RootAlgorithm.Secant:
                    root = new Secant();
                    break;
                case RootAlgorithm.Ridders:
                    root = new Ridders();
                    break;
                case RootAlgorithm.Steffenson:
                    root = new Steffenson();
                    break;
                default: 
                    root = new Brent();
                    break;
            }
            return root.Solve(funct, x_lower, x_upper);
        }
    }
    public abstract class RootBase
    {
        protected int m_nMaxIt;
        protected double m_dAbsTol;
        protected double m_dRelTol;       
        protected ConvergenceTest m_ConvTest;

        protected RootBase(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt)
        {
            if (RelTol < 0.0)
                throw new ArgumentException("Root: relative tolerance is negative");

            if (AbsTol < 0.0)
                throw new ArgumentException("Root: absolute tolerance is negative");

            if (MaxIt < 0)
                throw new ArgumentException("Root: maximum number of iterations is negative");                

            m_nMaxIt  = MaxIt;
            m_dAbsTol = AbsTol;
            m_dRelTol = RelTol;
            m_ConvTest = ConvTest; 
        }

        protected RootBase()
        {
            m_nMaxIt   = 100;
            m_dAbsTol  = 1e-8;
            m_dRelTol  = 1e-6;
            m_ConvTest = ConvergenceTest.Delta;
        }

        protected bool CheckConvergence(double f, double xl, double xu, double x1, double x0)
        {
            //check these first to avoid round-off error problems
            if (f == 0 || x1 == x0 || xl == xu) return true;            

            switch(m_ConvTest)
            {
                case ConvergenceTest.Interval: return (Abs(xl - xu) < m_dAbsTol + m_dRelTol * Min(Abs(xl), Abs(xu)));
                case ConvergenceTest.Residual: return (Abs(f) < m_dAbsTol);
                case ConvergenceTest.Delta: return (Abs(x1 - x0) < m_dAbsTol + m_dRelTol * Abs(xl));
                default: return (Abs(f) < m_dAbsTol);
            }            
        }
        public abstract IterationResults Solve(obj_function funct, double x_lower, double x_upper);
        public double SolveSimple(obj_function funct, double x_lower, double x_upper)
        {
            return Solve(funct, x_lower, x_upper).Root;
        }
     }

    public class Bisection : RootBase 
    {
        public Bisection(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt) : base(AbsTol, RelTol, ConvTest, MaxIt) { }
        public Bisection() : base() { }

        public override IterationResults Solve(obj_function funct, double x_lower, double x_upper)
        {
            int nIteration = 0;
            double root;
            double f_lo, f_up, f_mid;            
            double x_lo, x_up, x_mid;

            x_lo = x_lower;
            x_up = x_upper;            

            f_lo = funct(x_lo);
            f_up = funct(x_up);            

            if ((f_lo < 0.0 && f_up < 0.0) || (f_lo > 0.0 && f_up > 0.0))
            {
                throw new Exception("Bisection: root is not bracketed");
            }

            if (f_lo == 0.0) return new IterationResults(nIteration, x_lo, x_lo, x_lo, f_lo);
            if (f_up == 0.0) return new IterationResults(nIteration, x_up, x_up, x_up, f_up);

            root = 0.5 * (x_lo + x_up);
            do
            {
                nIteration++;                

                x_mid = root;
                f_mid = funct(x_mid);

                if (f_mid == 0.0) return new IterationResults(nIteration, x_lo, x_up, x_mid, f_mid);
                
                if ((f_lo > 0.0 && f_mid < 0.0) || (f_lo < 0.0 && f_mid > 0.0))
                {
                    x_up = x_mid;
                    f_up = f_mid;
                }
                else
                {
                    x_lo = x_mid;
                    f_lo = f_mid;
                }
                root = 0.5 * (x_lo + x_up);

            } while (nIteration < m_nMaxIt && !CheckConvergence(f_mid, x_lo, x_up, root, x_mid));

            return new IterationResults(nIteration, x_lo, x_up, root, f_mid);
        }
    }

    public class FalsePosition : RootBase
    {
        public FalsePosition(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt) : base(AbsTol, RelTol, ConvTest, MaxIt) { }
        public FalsePosition() : base() { }

        public override IterationResults Solve(obj_function funct, double x_lower, double x_upper)
        {
            int nIteration = 0;
            double root;
            double f_lo, f_up, f_pos;
            double x_lo, x_up, x_pos;
            int side = 0;

            x_lo = x_lower;
            x_up = x_upper;

            f_lo = funct(x_lo);
            f_up = funct(x_up);

            if ((f_lo < 0.0 && f_up < 0.0) || (f_lo > 0.0 && f_up > 0.0))
            {
                throw new Exception("FalsePosition: root is not bracketed");
            }

            if (f_lo == 0.0) return new IterationResults(nIteration, x_lo, x_lo, x_lo, f_lo);
            if (f_up == 0.0) return new IterationResults(nIteration, x_up, x_up, x_up, f_up);

            root = x_up - (f_up * (x_up - x_lo) / (f_up - f_lo));
            do
            {
                nIteration++;
                
                x_pos = root;
                f_pos = funct(x_pos);

                if (f_pos == 0.0) return new IterationResults(nIteration, x_pos, x_pos, x_pos, f_pos);

                if ((f_lo > 0.0 && f_pos < 0.0) || (f_lo < 0.0 && f_pos > 0.0))
                {                    
                    x_up = x_pos;
                    f_up = f_pos;                    
                    if (side == -1) f_lo = 0.5 * f_lo;
                    side = -1;
                }
                else
                {                    
                    x_lo = x_pos;
                    f_lo = f_pos;                                        
                    if (side == 1) f_up = 0.5 * f_up;
                    side = 1;
                }
                root = x_up - (f_up * (x_up - x_lo) / (f_up - f_lo));
                
            } while (nIteration < m_nMaxIt && !CheckConvergence(f_pos, x_lo, x_up, root, x_pos));

            return new IterationResults(nIteration, x_lo, x_up, root, f_pos);
        }
    }
    
    public class Brent : RootBase
    {
        public Brent(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt) : base(AbsTol, RelTol, ConvTest, MaxIt) { }
        public Brent() : base() { }

        public override IterationResults Solve(obj_function funct, double x_lower, double x_upper)
        {
            int nIteration = 0;
            double root, prev_root;
            double f_lo, f_up;
            double x_lo, x_up;

            double  a,  b,  c;
            double fa, fb, fc;
            double d, e;
            double tol, m;

            x_lo = x_lower;
            x_up = x_upper;

            f_lo = funct(x_lo);
            f_up = funct(x_up);

            if ((f_lo < 0.0 && f_up < 0.0) || (f_lo > 0.0 && f_up > 0.0))
            {
                throw new Exception("Brent: root is not bracketed");
            }

            if (f_lo == 0.0) return new IterationResults(nIteration, x_lo, x_lo, x_lo, f_lo);
            if (f_up == 0.0) return new IterationResults(nIteration, x_up, x_up, x_up, f_up);

            a  = x_lo;
            fa = f_lo;

            b = x_up;
            fb = f_up;

            c = x_up;
            fc = f_up;

            d = x_up - x_lo;
            e = x_up - x_lo;

            root = 0.5 * (x_lo + x_up); 
            do
            {
                nIteration++;

                int ac_equal = 0;

                //Insure that b is the best result so far, a is the previous
                //value of b, and c is on the opposite of the zero from b.
                if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
                {
                    ac_equal = 1;
                    c = a;
                    fc = fa;
                    d = e = b - a;                    
                }

                if (Abs(fc) <Abs(fb))
                {
                    ac_equal = 1;
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }

                tol = 2.0 * Const.epsilon * Max(1.0, Abs(b));
                m = 0.5 * (c - b);

                if (fb == 0)
                {
                    return new IterationResults(nIteration, b, b, b, fb);
                }

                if (Abs(m) <= tol)
                {
                    root = b;

                    if (b < c)
                    {
                        x_lo = b;
                        x_up = c;
                    }
                    else
                    {
                        x_lo = c;
                        x_up = b;
                    }
                    return new IterationResults(nIteration, x_lo, x_up, root, fb);
                }

                // use bisection
                if (Abs(e) < tol || Abs(fa) <= Abs(fb))
                {
                    d = m;             
                    e = m;
                }
                else
                {
                    // use inverse quadratic interpolation 
                    double p, q, r;   
                    double s = fb / fa;

                    if (ac_equal==1)
                    {
                        p = 2 * m * s;
                        q = 1 - s;
                    }
                    else
                    {
                        q = fa / fc;
                        r = fb / fc;
                        p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                        q = (q - 1) * (r - 1) * (s - 1);
                    }

                    if (p > 0)
                    {
                        q = -q;
                    }
                    else
                    {
                        p = -p;
                    }

                    if (2 * p < Min(3 * m * q - Abs(tol * q), Abs(e * q)))
                    {
                        // interpolation
                        e = d;
                        d = p / q;
                    }
                    else
                    {
                        // interpolation failed, use bisection

                        d = m;
                        e = m;
                    }
                }

                a = b;
                fa = fb;

                if (Abs(d) > tol)
                {
                    b += d;
                }
                else
                {
                    b += (m > 0 ? +tol : -tol);
                }

                fb = funct(b);

                // Update the best estimate of the root and bounds on each iteration
                prev_root = root;
                root = b;

                if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
                {
                    c = a;
                }

                if (b < c)
                {
                    x_lo = b;
                    x_up = c;
                }
                else
                {
                    x_lo = c;
                    x_up = b;
                }

            } while (nIteration < m_nMaxIt && !CheckConvergence(fb, x_lo, x_up, root, prev_root));

            return new IterationResults(nIteration, x_lo, x_up, root, fb);
        }
    }
    
    public class Secant : RootBase
    {
        public Secant(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt) : base(AbsTol, RelTol, ConvTest, MaxIt) { }
        public Secant() : base() { }

        public override IterationResults Solve(obj_function funct, double x_lower, double x_upper)
        {
            int nIteration = 0;
            double root;
            double f_lo, f_up;
            double x_lo, x_up;
            double x, x_new;
            double f, df, df_new, f_new;

            x_lo = x_lower;
            x_up = x_upper;

            f_lo = funct(x_lo);
            f_up = funct(x_up);

            if (f_lo == 0.0) return new IterationResults(nIteration, x_lo, x_lo, x_lo, f_lo);
            if (f_up == 0.0) return new IterationResults(nIteration, x_up, x_up, x_up, f_up);

            if (Abs(f_lo) < Abs(f_up))
            {
                root = x_lo;
                f = f_lo;
            }
            else
            {
                root = x_up;
                f = f_up;
            }

            df = (f_lo - f_up) / (x_lo - x_up);
            do
            {
                nIteration++;

                x = root;
                x_new = x - (f / df);

                f_new = funct(x_new);
                df_new = (f_new - f) / (x_new - x);

                root = x_new;
                f = f_new;
                df = df_new;

                if (x_new < x)
                {
                    x_lo = x_new;
                    x_up = x;
                }
                else
                {
                    x_lo = x;
                    x_up = x_new;
                }
            } while (nIteration < m_nMaxIt && !CheckConvergence(f_new, x_lo, x_up, root, x));

            return new IterationResults(nIteration, x_lo, x_up, root, f_new);
        }
    }

    public class Ridders : RootBase
    {
        public Ridders(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt) : base(AbsTol, RelTol, ConvTest, MaxIt) { }
        public Ridders() : base() { }
        
        double SIGN(double a, double b)
        {
            return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
        }

        public override IterationResults Solve(obj_function funct, double x_lower, double x_upper)
        {
            int nIteration = 0;
            double root;
            double f_lo, f_up, f_mid;
            double x_lo, x_up, x_mid;
            double s, x_new, f_new;


            x_lo = x_lower;
            x_up = x_upper;

            f_lo = funct(x_lo);
            f_up = funct(x_up);

            if ((f_lo < 0.0 && f_up < 0.0) || (f_lo > 0.0 && f_up > 0.0))
            {
                throw new Exception("Ridders: root is not bracketed");
            }

            if (f_lo == 0.0) return new IterationResults(nIteration, x_lo, x_lo, x_lo, f_lo);
            if (f_up == 0.0) return new IterationResults(nIteration, x_up, x_up, x_up, f_up);
            
            do
            {
                nIteration++;

                x_mid = 0.5 * (x_lo + x_up);
                f_mid = funct(x_mid);

                s = Sqrt(f_mid * f_mid - f_lo * f_up);

                if (s == 0.0) return new IterationResults(nIteration, x_lo, x_up, x_mid, f_mid);

                x_new = x_mid + (x_mid - x_lo) * ((f_lo >= f_up ? 1.0 : -1.0) * f_mid / s);

                root = x_new;
                f_new = funct(root);

                if (f_new == 0.0) return new IterationResults(nIteration, x_lo, x_up, root, f_new);

                if (SIGN(f_mid, f_new) != f_mid)
                {
                    x_lo = x_mid;
                    f_lo = f_mid;
                    x_up = root;
                    f_up = f_new;
                }
                else if(SIGN(f_lo,f_new) != f_lo)
                {
                    x_up = root;
                    f_up = f_new;
                }
                else if (SIGN(f_up, f_new) != f_up)
                {
                    x_lo = root;
                    f_lo = f_new;
                }
                else
                {
                    throw new Exception("Ridders: root is not bracketed");
                }               

            } while (nIteration < m_nMaxIt && !CheckConvergence(f_new, x_lo, x_up, root, x_mid));

            return new IterationResults(nIteration, x_lo, x_up, root, f_new);
        }
    }

    public class Steffenson : RootBase
    {
        public Steffenson(double AbsTol, double RelTol, ConvergenceTest ConvTest, int MaxIt) : base(AbsTol, RelTol, ConvTest, MaxIt) { }
        public Steffenson() : base() { }

        public override IterationResults Solve(obj_function funct, double x_lower, double x_upper)
        {
            int nIteration = 0;
            double root;
            double f_lo, f_up;
            double x_lo, x_up;
            double x, x_new, x_prev;
            double f, df, df_new, f_new;

            x_lo = x_lower;
            x_up = x_upper;

            f_lo = funct(x_lo);
            f_up = funct(x_up);

            if (f_lo == 0.0) return new IterationResults(nIteration, x_lo, x_lo, x_lo, f_lo);
            if (f_up == 0.0) return new IterationResults(nIteration, x_up, x_up, x_up, f_up);

            if (Abs(f_lo) < Abs(f_up))
            {
                root = x_lo;                
                f = f_lo;
            }
            else
            {
                root = x_up;                
                f = f_up;
            }

            x = root;

            df = (f_lo - f_up) / (x_lo - x_up);
            do
            {
                nIteration++;

                x_prev = x;
                x = root;
                x_new = x - (f / df);

                f_new = funct(x_new);
                df_new = (f_new - f) / (x_new - x);

                if (nIteration < 3)
                    root = x_new;
                else
                {
                    double u = (x - x_prev);
                    double v = (x_new - 2 * x + x_prev);

                    if (v == 0)
                        root = x_new;  // avoid division by zero
                    else
                        root = x_prev - u * u / v;  // accelerated value 

                    f_new = funct(root);
                    df_new = (f_new - f) / (root - x_new);
                }

                f = f_new;
                df = df_new;

                if (x_new < x)
                {
                    x_lo = x_new;
                    x_up = x;
                }
                else
                {
                    x_lo = x;
                    x_up = x_new;
                }
            } while (nIteration < m_nMaxIt && !CheckConvergence(f_new, x_lo, x_up, root, x));

            return new IterationResults(nIteration, x_lo, x_up, root, f_new);
        }
    }
}
