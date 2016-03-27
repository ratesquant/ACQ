using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;


namespace ACQ.Excel
{
    public class Examples
    {
        [ExcelFunction(Description = "acq_sum")]
        public static double acq_sum(double a, double b)
        {
            return a + b;
        }
    }
}
