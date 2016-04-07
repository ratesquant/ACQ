using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Random
{
    public class BoxMuller
    { 
        RandomBase m_rng;
        double m_spareValue = Double.NaN;

        public BoxMuller(RandomBase rng)            
        {
            m_rng = rng;
        }

        public double Next()
        {
            if (!Double.IsNaN(m_spareValue))
            {
                double tmp = m_spareValue;
                m_spareValue = Double.NaN;
                return tmp;
            }

            // Generate two new gaussian values.
            double x, y, sqr;

            // We need a non-zero random point inside the unit circle.
            do
            {
                x = 2.0 * m_rng.NextDouble() - 1.0;
                y = 2.0 * m_rng.NextDouble() - 1.0;
                sqr = x * x + y * y;
            }
            while (sqr > 1.0 || sqr == 0);

            // Make the Box-Muller transformation.
            double fac = System.Math.Sqrt(-2.0 * System.Math.Log(sqr) / sqr);

            m_spareValue = x * fac;
            return y * fac;
        }

        public void Next(double[] sample)
        {
            int index = 0;

            while (true)
            {
                // Generate two new gaussian values.
                double x, y, sqr;

                // We need a non-zero random point inside the unit circle.
                do
                {
                    x = 2.0 * m_rng.NextDouble() - 1.0;
                    y = 2.0 * m_rng.NextDouble() - 1.0;
                    sqr = x * x + y * y;
                }
                while (sqr > 1.0 || sqr == 0);

                // Make the Box-Muller transformation.
                double fac = System.Math.Sqrt(-2.0 * System.Math.Log(sqr) / sqr);

                if (index < sample.Length)
                {
                    sample[index++] = x * fac;

                    if (index < sample.Length)
                        sample[index++] = y * fac;
                    else
                        break;
                }else
                    break;
            }
        }
    }

}
