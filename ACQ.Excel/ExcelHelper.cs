using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    static class ExcelHelper
    {
        internal static object CheckNan(double value)
        {
            object result;

            if (Double.IsNaN(value))
            {
                result = ExcelError.ExcelErrorNum;
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
    }
}
