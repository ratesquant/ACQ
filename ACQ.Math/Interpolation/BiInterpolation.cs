using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class BiInterpolation<T> : InterpolationBase2D where T : InterpolationInterface
    {
        public BiInterpolation(double[] x1, double[] x2, double[,] y)
            : base(x1, x2, y, false)
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

                bool bounds = false; //this is only helpful for 1D

                InterpolationInterface interpolator;

                for(int i=0; i<n2; i++)
                {
                    for (int j = 0; j < n1; j++)
                    {
                        yt[j] = m_y[i, j];
                    }
                    interpolator = Activator.CreateInstance(typeof(T), m_x1, yt, bounds) as InterpolationInterface;

                    y2[i] = interpolator.Eval(x1);
                }

                interpolator = Activator.CreateInstance(typeof(T), m_x2, y2, bounds) as InterpolationInterface;

                return interpolator.Eval(x2);
    
            }

            return value;
        }

    }
}
