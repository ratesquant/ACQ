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
    }

    class AddInInit : IExcelAddIn
    {
        public void AutoOpen()
        {
            ExcelIntegration.RegisterUnhandledExceptionHandler(ErrorHandler);
        }

        public void AutoClose()
        { 

        }

        private object ErrorHandler(object exceptionObject)
        {
            ExcelReference caller = (ExcelReference)XlCall.Excel(XlCall.xlfCaller);

       
            LogDisplay.WriteLine("Error: " + exceptionObject.ToString());

            // return #VALUE into the cell.
            return ExcelError.ExcelErrorValue; 
        }
    }
}
