using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class ExcelRandom
    {
        [ExcelFunction(Description = "Create random vector, each number [0,1) (mersenne twister)", Category = AddInInfo.Category)]
        public static object acq_random_vector(int seed, int size)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelVector.Tag, new object[] { seed, size, "acq_random_vector" },
                    (objectType, parameters) =>
                    {

                        ACQ.Math.Random.MersenneTwister rng = new Math.Random.MersenneTwister((uint)seed);
                        double[] v = new double[size];

                        for (int i = 0; i < size; i++)
                        {
                            v[i] = rng.NextDouble();
                        }
                        return v;

                    });

            }
        }

        [ExcelFunction(Description = "Create random vector with array seed (uses mersenne twister)", Category = AddInInfo.Category)]
        public static object acq_random_vector_ex(double[] seed, int size)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelVector.Tag, new object[] { seed, size, "acq_random_vector" },
                    (objectType, parameters) =>
                    {
                        uint[] useed = new uint[seed.Length];

                        for (int i = 0; i < seed.Length; i++)
                        {
                            useed[i] = (uint)seed[i];
                        }
                        ACQ.Math.Random.MersenneTwister rng = new Math.Random.MersenneTwister(useed);
                        double[] v = new double[size];

                        for (int i = 0; i < size; i++)
                        {
                            v[i] = rng.NextDouble();
                        }
                        return v;

                    });

            }
        }
    }

}
