using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

namespace ACQ.Excel
{
    public static class AddInInfo
    {
        public const string Category = "ACQ";
        public const string MenuName = "ACQ";
    }

    class AddInInit : IExcelAddIn
    {
        public void AutoOpen()
        {
            ExcelIntegration.RegisterUnhandledExceptionHandler(ErrorHandler);

            // Register Ctrl+Shift+H to show log window 
            XlCall.Excel(XlCall.xlcOnKey, "^H", "ShowLogWindow");

            // Register Ctrl+Shift+A to show introspection window 
            XlCall.Excel(XlCall.xlcOnKey, "^A", "ShowIntrospectionWindow"); 
        }

        public void AutoClose()
        { 

        }

        private object ErrorHandler(object exceptionObject)
        {
            if (exceptionObject != null)
            {
                LogDisplay.WriteLine("Exception: " + exceptionObject.ToString());
            }

            // return #VALUE into the cell.
            return ExcelError.ExcelErrorValue; 
        }
    }
}
