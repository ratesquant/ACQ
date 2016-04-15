using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Random
{
    #region MTNInfo
    /*
       A C-program for MT19937, with initialization improved 2002/1/26.
       Coded by Takuji Nishimura and Makoto Matsumoto.

       Before using, initialize the state by using init_genrand(seed)  
       or init_by_array(init_key, key_length).

       Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
       All rights reserved.                          

       Redistribution and use in source and binary forms, with or without
       modification, are permitted provided that the following conditions
       are met:

         1. Redistributions of source code must retain the above copyright
            notice, this list of conditions and the following disclaimer.

         2. Redistributions in binary form must reproduce the above copyright
            notice, this list of conditions and the following disclaimer in the
            documentation and/or other materials provided with the distribution.

         3. The names of its contributors may not be used to endorse or promote 
            products derived from this software without specific prior written 
            permission.

       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
       "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
       LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
       A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
       CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


       Any feedback is very welcome.
       http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
       email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
    */
    #endregion
    /// <summary>
    /// Generates pseudo-random numbers using the Mersenne Twister algorithm (code adapted from mt19937ar).   
    /// </summary>
    /// <remarks>
    /// See <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
    /// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html</a> for details
    /// on the algorithm.
    /// </remarks>
    public class MersenneTwister : RandomBase
    {
        // Period parameters 
        private const int N = 624;
        private const int M = 397;
        private const uint MatrixA = 0x9908b0df; // constant vector a
        private const uint UpperMask = 0x80000000; // most significant w-r bits
        private const uint LowerMask = 0x7fffffff; // least significant r bits 

        // Tempering parameters
        private const uint TemperingMaskB = 0x9d2c5680U;
        private const uint TemperingMaskC = 0xefc60000U;

        private readonly uint[] m_mt = new uint[N]; // the array for the state vector  
        private int m_mti;

        private static readonly uint[] m_mag01 = { 0x0, MatrixA };

        private const double m_intNormInc = (1.0/4294967295.0); //1.0 / 2^32-1 
        private const double m_intNormExc = (1.0/4294967296.0); //1.0 / 2^32
        private const double m_scale = (1.0/9007199254740992.0);


        /// <summary>
        /// Creates a new pseudo-random number generator with a given seed.
        /// </summary>
        /// <param name="seed">A value to use as a seed.</param>
        public MersenneTwister(uint seed = 5489U)
        {
            init_genrand(seed);
        }

        public MersenneTwister(uint[] seed)
        {
            init_by_array(seed); //{0x123U, 0x234U, 0x345U, 0x456U}
        }

        public override double NextDouble()
        {
            return genrand_real2();
        }

        public override int Next()
        {
            int value = genrand_int31();

            if (value == 0x7fffffff) // to make it consistent with System.Random
                return Next();
            else
                return value;
        }

        /// <summary>
        /// initializes mt[N] with a seed 
        /// </summary>
        /// <param name="seed"></param>
        private void init_genrand(uint seed)
        {
            m_mt[0] = seed & 0xffffffffU;

            for (m_mti = 1; m_mti < N; m_mti++)
            {
                m_mt[m_mti] = (uint)(1812433253U * (m_mt[m_mti - 1] ^ (m_mt[m_mti - 1] >> 30)) + m_mti);
                // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
                // In the previous versions, MSBs of the seed affect   
                // only MSBs of the array _mt[].                        
                // 2002/01/09 modified by Makoto Matsumoto             
                m_mt[m_mti] &= 0xffffffffU;
                // for >32 bit machines
            }
        }

        private void init_by_array(uint[] key)
        {
            int i, j, k;
            init_genrand(19650218U);

            Int32 keyLength = key.Length;
            i = 1; j = 0;
            k = (N > keyLength ? N : keyLength);

            for (; k > 0; k--)
            {
                m_mt[i] = (uint)((m_mt[i] ^ ((m_mt[i - 1] ^ (m_mt[i - 1] >> 30)) * 1664525U)) + key[j] + j); /* non linear */
                m_mt[i] &= 0xffffffffU; // for WORDSIZE > 32 machines
                i++; j++;
                if (i >= N) { m_mt[0] = m_mt[N - 1]; i = 1; }
                if (j >= keyLength) j = 0;
            }

            for (k = N - 1; k > 0; k--)
            {
                m_mt[i] = (uint)((m_mt[i] ^ ((m_mt[i - 1] ^ (m_mt[i - 1] >> 30)) * 1566083941U)) - i); /* non linear */
                m_mt[i] &= 0xffffffffU; // for WORDSIZE > 32 machines
                i++;

                if (i >= N) { m_mt[0] = m_mt[N - 1]; i = 1; }
            }

            m_mt[0] = 0x80000000U; // MSB is 1; assuring non-zero initial array
        }

        /// <summary>
        /// generates a random number on [0,0xffffffff(4294967295)]-interval
        /// </summary>
        /// <returns></returns>
        private uint genrand_int32()
        {
            uint y;
            // mag01[x] = x * MATRIX_A  for x=0,1 

            if (m_mti >= N) 
            { 
                //generate N words at one time 
                int kk;

                for (kk=0; kk<N-M; kk++) 
                {
                    y = (m_mt[kk]&UpperMask)|(m_mt[kk+1]&LowerMask);
                    m_mt[kk] = m_mt[kk+M] ^ (y >> 1) ^ m_mag01[y & 0x1UL];
                }
                for (; kk<N-1; kk++) 
                {
                    y = (m_mt[kk]&UpperMask)|(m_mt[kk+1]&LowerMask);
                    m_mt[kk] = m_mt[kk+(M-N)] ^ (y >> 1) ^ m_mag01[y & 0x1UL];
                }
                y = (m_mt[N-1]&UpperMask)|(m_mt[0]&LowerMask);
                m_mt[N-1] = m_mt[M-1] ^ (y >> 1) ^ m_mag01[y & 0x1UL];

                m_mti = 0;
            }
  
            y = m_mt[m_mti++];

            // Tempering
            y ^= (y >> 11);
            y ^= (y << 7) & TemperingMaskB;
            y ^= (y << 15) & TemperingMaskC;
            y ^= (y >> 18);

            return y;
        }

        /// <summary>
        /// generates a random number on [0,0x7fffffff(2147483647)]-interval
        /// </summary>
        /// <returns></returns>
        private int genrand_int31()
        {
            return (int)(genrand_int32()>>1);
        }

        /// <summary>
        ///  generates a random number on [0,1]-real-interval
        /// </summary>
        /// <returns></returns>
        public double genrand_real1()
        {
            return genrand_int32()*m_intNormInc;
        }
        
        /// <summary>
        /// generates a random number on [0,1)-real-interval
        /// </summary>
        /// <returns></returns>
        public double genrand_real2()
        {
            return genrand_int32()*m_intNormExc;
        }

        /// <summary>
        /// generates a random number on (0,1)-real-interval
        /// </summary>
        /// <returns></returns>
        public double genrand_real3()
        {
            return (((double)genrand_int32()) + 0.5)*m_intNormExc;  //divided by 2^32
        }

        /// <summary>
        /// generates a random number on [0,1) with 53-bit resolution
        /// </summary>
        /// <returns></returns>
        public double genrand_res53() 
        { 
            uint a=genrand_int32()>>5, b=genrand_int32()>>6;
            return (a * 67108864.0 + b) * m_scale; 
        } 
    }
}
