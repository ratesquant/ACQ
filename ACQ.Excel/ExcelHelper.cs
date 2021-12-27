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
    [Flags]
    internal enum XlType : int
    {
        XlTypeNumber = 0x0001,
        XlTypeString = 0x0002,
        XlTypeBoolean = 0x0004,
        XlTypeReference = 0x0008,
        XlTypeError = 0x0010,
        XlTypeArray = 0x0040,
        XlTypeMissing = 0x0080,
        XlTypeEmpty = 0x0100,
        XlTypeInt = 0x0800,     // int16 in XlOper, int32 in XlOper12, never passed into UDF
    }    

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

        public static bool IsInteger(object item, out int result)
        {
            bool is_integer = false;
            result = 0; //only valid if is_integer == true

            if (item == null || item is ExcelError || item is ExcelMissing || item is ExcelEmpty)
                is_integer = false;
            else if (item is string)
            {
                double float_number;
                if (int.TryParse(item as string, out result))
                    is_integer = true;
                else if (double.TryParse(item as string, out float_number) && (float_number - System.Math.Floor(float_number)) == 0)
                {
                    result = (int)System.Math.Floor(float_number);
                    is_integer = true;
                }
            }
            else if (item is double)
            {
                double float_number = (double) item;
                if ((float_number - System.Math.Floor(float_number)) == 0)
                {
                    result = (int)System.Math.Floor(float_number);
                    is_integer = true;
                }
            }
            else
                is_integer = false;

            return is_integer;
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

        internal static T[] CheckArray<T>(object value) where T : struct
        {
            T[] result = null;

            object[,] array = value as object[,];

            if (array != null)
            {
                try
                {
                    int n = array.GetLength(0);
                    int m = array.GetLength(1);

                    if (n > 0)
                    {
                        var temp = new List<T>(n); //TODO: convert first column, think about how to do it in more general way

                        for (int i = 0; i < n; i++)
                        {
                            object item = array[i, 0];
                            if (!IsMissingOrEmpty(item)) //this is what excel does, it excludes empty cell from the selected range
                            {
                                temp.Add((T)item);
                            }
                        }
                        result = temp.ToArray();
                    }
                }
                catch(Exception e) //return null
                {
                    System.Diagnostics.Debug.Write(e.Message);
                }
            }

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
