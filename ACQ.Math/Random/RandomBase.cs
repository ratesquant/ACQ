using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Random
{
    public interface IRandomGenerator
    {
        /// <summary>
        /// Returns a random number in the interval [0, 1)
        /// </summary>
        /// <returns></returns>
        double NextDouble();
        int Next();
    }

    public abstract class RandomBase : IRandomGenerator
    {
        public abstract double NextDouble();
        public abstract int Next();
    }
}
