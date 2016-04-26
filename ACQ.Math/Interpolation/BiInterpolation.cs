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

        public override double Eval(double x1, double x2)
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

    }
}
