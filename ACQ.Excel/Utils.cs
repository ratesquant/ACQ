using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    public static class Utils
    {
        private static string Excel_ToString(object item, bool ignore_errors = true)
        {
            string item_str = "";
      
            if (item == null)
                item_str = ignore_errors ? "" : "#NULL";
            else if (item is string)
                item_str = item as string;
            else if (item is ExcelError)
                item_str = ignore_errors ? "" : item.ToString().Replace("ExcelError", "#");
            else if (item is ExcelMissing)
                item_str = ignore_errors ? "" : "#Missing";
            else if (item is ExcelEmpty)
                item_str = ignore_errors ? "" : "#Empty";
            else
                item_str = item.ToString();      

            return item_str;
        }

        [ExcelFunction(Description = "Count number of unique elements in the array", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_count_unique(object[] x,
             [ExcelArgument(Description = "Optional flag if TRUE (default) ignores differences between error codes")] object ignore_error_codes)
        {
            object result = ExcelError.ExcelErrorNA;

            if (x == null )
            {
                result = ExcelError.ExcelErrorNA;
            }
            else
            {
                bool ignore_errors = ExcelHelper.CheckValue<bool>(ignore_error_codes, true);

                HashSet<string> unique_values = new HashSet<string>();
                for (int i = 0; i < x.Length; i++)
                {
                    unique_values.Add(Excel_ToString(x[i], ignore_errors));
                }
                result = unique_values.Count;
            }
            return result;
        }

        [ExcelFunction(Description = "Concatenates the elements of a specified array using the specified separator (optional) between each element", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_join(object[] x,
            [ExcelArgument(Description = "Optional separator: default is comma")] object separator)
        {
            object result = ExcelError.ExcelErrorNA;

            if (x == null)
            {
                result = ExcelError.ExcelErrorNA;
            }
            else
            {
                string sep = ExcelHelper.Check<string>(separator, ",");

                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < x.Length; i++)
                {
                    if (i != 0)
                    {
                        sb.Append(sep);
                    }
                    string item_str = Excel_ToString(x[i]);
                    
                    sb.Append(item_str);                    
                }
                result = sb.ToString();
            }
            return result;
        }

        [ExcelFunction]
        public static string acq_tostring(object value)
        {
            object result;
            var retVal = XlCall.TryExcel(XlCall.xlCoerce, out result, value, (int)XlType.XlTypeString);
            if (retVal == XlCall.XlReturn.XlReturnSuccess)
            {
                return (string)result;
            }
            else
                return "";

        }
    }
}