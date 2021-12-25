using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel
{
    public static class DateUtils
    {
        /// <summary>
        /// Checks if a the year is a leap year (days in year is about 365.2425)
        /// Any year that is evenly divisible by 4 is a leap year, year that is evenly divisible by 100 (for example, 1900) is a leap year only if it is also evenly divisible by 400
        /// </summary>
        /// <param name="year"></param>
        /// <returns></returns>
        public static bool IsLeapYear(int year)
        {
            bool isLeap = (year % 400 == 0) || (year % 4 == 0 && year % 100 != 0);

            return isLeap;
        }

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

        [ExcelFunction(Description = "Returns next business day (doesn't use holiday calendar)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static DateTime acq_nextbusinessday(DateTime date)
        {
            DateTime businessday = date.AddDays(1);

            while (IsHoliday(businessday))
            {
                businessday = businessday.AddDays(1);
            }
            return businessday;
        }

        [ExcelFunction(Description = "Returns previous business day (doesn't use holiday calendar)", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static DateTime acq_prevbusinessday(DateTime date)
        {
            DateTime businessday = date.AddDays(-1);

            while (IsHoliday(businessday))
            {
                businessday = businessday.AddDays(-1);
            }
            return businessday;
        }

        [ExcelFunction(Description = "Adjust to business day", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static DateTime acq_adjustbusinessday(DateTime date,
            [ExcelArgument(Description = "Adjustment direction (-1, 1)")] object direction,
            [ExcelArgument(Description = "Adjusted day should be in the same month")] object same_month)
        {
            int step = (int)ExcelHelper.CheckValue<double>(direction, 1); //default is next
            bool mod = ExcelHelper.CheckValue<bool>(same_month, true); //default is true

            DateTime businessday = date;           

            while (IsHoliday(businessday))
            {
                businessday = businessday.AddDays(step);

                if (mod && businessday.Month != date.Month)
                {
                    step = -step;
                    businessday = date; //reset to start date, we know it is a holiday
                }
            }            

            return businessday;
        }

        [ExcelFunction(Description = "Check if Year is a leap year", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_isleap_year(object year)
        {
            object result = ExcelError.ExcelErrorValue;

            int n_year;
            if (ExcelHelper.IsInteger(year, out n_year))
            {
                result = IsLeapYear(n_year);
            }

            return result;
        }

        public static bool IsHoliday(DateTime date)
        {
            return date.DayOfWeek == DayOfWeek.Saturday || date.DayOfWeek == DayOfWeek.Sunday;
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
