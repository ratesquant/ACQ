using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math
{
    public class Primes
    {
        public static bool is_prime(int n)
        {
            bool result = true;
            if (n == 0 || n == 1) //not primes
            {
                result = false;
            }
            else
            {
                for (int i = 2; i <= n / 2; ++i)
                {
                    if (n % i == 0)
                    {
                        result = false;
                        break;
                    }
                }
            }
            return result;
        }
    }
}
