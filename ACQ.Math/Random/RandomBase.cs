using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Random
{
    public abstract class RandomBase
    {
        /// <summary>
        /// Returns a random number in the interval [0, 1)
        /// </summary>
        /// <returns></returns>
        public abstract double NextDouble();
        public abstract int Next();
    }
}
