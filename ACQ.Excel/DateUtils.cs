using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    public static class DateUtils
    {
        [ExcelFunction(Description = "acq_convert_todate", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_convert_todate(object date)
        {
            if (date is ExcelMissing || date is ExcelEmpty)
            {
                return ExcelError.ExcelErrorRef;
            }
            else if (date is string)
            {
                string date_string = date as string;

                if (!String.IsNullOrWhiteSpace(date_string))
                {
                    DateTime value;
                    int iso_date;
                    double oad_date;
                    if (DateTime.TryParse(date as string, out value))
                    {
                        return value;
                    }
                    else if (Int32.TryParse(date_string, out iso_date))
                    {
                        return IntToDate(iso_date);
                    }
                    else if (Double.TryParse(date_string, out oad_date))
                    {
                        return DateTime.FromOADate((double)date);
                    }
                }
            }
            else if (date is int)
            {
                //assume it is ISO date YYYYMMDD
                return IntToDate((int)date);
            }
            else if (date is double)
            {
                if (System.Math.Abs((double)date - System.Math.Round((double)date, 0)) < 1e-6)
                    return IntToDate((int)System.Math.Round((double)date, 0));
                else
                    return DateTime.FromOADate((double)date);
            }

            return ExcelError.ExcelErrorValue;
        }

        public static DateTime IntToDate(int value)
        {
            int year = value / 10000;
            int day = value % 100;
            int month = (value - year * 10000)/100;

            return new DateTime(year, month, day); 
        }
    }
}
