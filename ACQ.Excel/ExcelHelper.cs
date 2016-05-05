using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    /*
     * ExcelError.ExcelErrorNum -  Problem with a number in the formula
     * ExcelError.ExcelErrorNull - Null range is specified as an input
     * ExcelError.ExcelErrorRef -  	Invalid cell reference
     * ExcelError.ExcelErrorValue - Wrong type of argument in a function or wrong type of operator
     * ExcelError.ExcelErrorNA - No value available (use this for Nan, excel plot ExcelErrorNum as zeros)
     * 
     */
    static class ExcelHelper
    {
        internal static object CheckNan(double value)
        {
            object result;

            if (Double.IsNaN(value) || Double.IsInfinity(value))
            {
                result = ExcelError.ExcelErrorNA;
            }
            else
            {
                result = value;
            }

            return result;
        }


        internal static bool IsMissingOrEmpty(object value)
        {
            return value is ExcelMissing || value is ExcelEmpty;
        }

        internal static T CheckValue<T>(object value, T defaultValue) where T : struct
        {
            T result = value is T ? (T)value : defaultValue;
           
            return result;
        }

        internal static T Check<T>(object value, T defaultValue) where T : class
        {
            T result = value is T ? value as T : defaultValue;

            return result;
        }

        internal static T CheckEnum<T>(object value, T defaultValue) where T : struct, IConvertible
        {
            T result = defaultValue;

            if (value is string)
            {
                string string_value = value as string;

                if (!String.IsNullOrWhiteSpace(string_value))
                {
                    T parsedValue;

                    if (Enum.TryParse<T>(value as string, true, out parsedValue))
                    {
                        result = parsedValue;
                    }
                }
            }
           
            return result;
        }

        internal static object[,] CreateArray(int size1, int size2, object value)
        {
            object[,] result = new object[size1, size2];

            for (int i = 0; i < size1; i++)
            {
                for (int j = 0; j < size2; j++)
                {
                    result[i, j] = value;
                }
            }

            return result;
        }

        internal static object[,] BoxArray<T>(T[,] a)
        {
            object[,] result = new object[a.GetLength(0), a.GetLength(1)];

            for (int i = 0; i < a.GetLength(0); i++)
            {
                for (int j = 0; j < a.GetLength(1); j++)
                {
                    result[i, j] = a[i, j];
                }
            }

            return result;
        }
    }
}
