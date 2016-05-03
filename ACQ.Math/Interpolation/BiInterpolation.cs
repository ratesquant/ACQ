using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    /// <summary>
    /// 2D interpolation on the rectangular grid, done using specified 1D interpolator in each dimension
    /// the method is general, but not efficient
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class BiInterpolation<T> : InterpolationBase2D where T : InterpolationInterface
    {
        public BiInterpolation(double[] x1, double[] x2, double[,] y)
            : base(x1, x2, y, true)
        {
        }
        public BiInterpolation(double[] x1, double[] x2, double[,] y, bool copyData)
            : base(x1, x2, y, copyData)
        {
        }

        /// <summary>
        /// returns number of support nodes,  
        /// </summary>
        /// <returns></returns>
        protected virtual int SupportSize
        {
            get
            {
                return 0;
            }
        }

        public double EvalRef(double x1, double x2)
        {
            double value = Double.NaN;

            int index1, index2;


            FindIndex(x1, x2, out index1, out index2);

            if (index1 > 0 && index2 > 0)
            {
                //slow reference method

                int n1 = m_x1.Length;
                int n2 = m_x2.Length;

                double[] yt = new double[n1];
                double[] y2 = new double[n2];

                InterpolationInterface interpolator;
                Type interpolator_type = typeof(T);

                for(int i=0; i<n2; i++)
                {
                    for (int j = 0; j < n1; j++)
                    {
                        yt[j] = m_y[i, j];
                    }
                    interpolator = Activator.CreateInstance(interpolator_type, m_x1, yt) as InterpolationInterface;
                    interpolator.Bounds = false;

                    y2[i] = interpolator.Eval(x1);
                }

                interpolator = Activator.CreateInstance(interpolator_type, m_x2, y2) as InterpolationInterface;

                return interpolator.Eval(x2);
    
            }

            return value;
        }

        public override double Eval(double x1, double x2)
        {
            double value = Double.NaN;

            int index1, index2;


            FindIndex(x1, x2, out index1, out index2);

            if (index1 > 0 && index2 > 0)
            {
                int sn = this.SupportSize;
                int n1 = m_x1.Length;
                int n2 = m_x2.Length;

                int j0, j1, i0, i1;
                //determine interpolation window
                //use all grind 
                if (sn <= 0)
                {
                    j0 = i0 = 0;
                    j1 = n1;
                    i1 = n2;
                }
                else
                {
                    j0 = System.Math.Max(0,  index1 - sn);
                    j1 = System.Math.Min(n1, index1 + sn);

                    i0 = System.Math.Max(0,  index2 - sn);
                    i1 = System.Math.Min(n2, index2 + sn);
                }



                double[] x1t = new double[j1 - j0];
                double[] yt = new double[j1 - j0];
                double[] x2t = new double[i1 - i0];
                double[] y2 = new double[i1 - i0];

                InterpolationInterface interpolator;
                Type interpolator_type = typeof(T);

                for (int j = j0; j < j1; j++)
                {
                    x1t[j - j0] = m_x1[j];
                }

                for (int i = i0; i < i1; i++)
                {
                    for (int j = j0; j < j1; j++)
                    {
                        yt[j - j0] = m_y[i, j];
                    }
                    interpolator = Activator.CreateInstance(interpolator_type, x1t, yt) as InterpolationInterface;
                    interpolator.Bounds = false;

                    y2[i - i0] = interpolator.Eval(x1);
                    x2t[i - i0] = m_x2[i];
                }

                interpolator = Activator.CreateInstance(interpolator_type, x2t, y2) as InterpolationInterface;

                return interpolator.Eval(x2);
            }

            return value;
        }
    }
}
