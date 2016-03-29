using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

namespace ACQ.Excel
{
    public static class Info
    {
        [ExcelFunction(Description = "acq_version", Category = AddInInfo.Category)]
        public static object acq_version()
        {
            return System.Reflection.Assembly.GetExecutingAssembly().GetName().Version.ToString();
        }
        [ExcelFunction(Description = "acq_excel_version", Category = AddInInfo.Category)]
        public static object acq_excel_version()
        {
            return ExcelDnaUtil.ExcelVersion;
        }

        [ExcelFunction(Description = "acq_addin_path", Category = AddInInfo.Category)]
        public static string acq_xllpath()
        {
            //ExcelDna.Integration.DnaLibrary.CurrentLibrary
            return ExcelDnaUtil.XllPath;// System.IO.Path.GetDirectoryName((string)XlCall.Excel(XlCall.xlGetName));
        }
        [ExcelCommand(MenuText = "Show Log Window", MenuName = AddInInfo.MenuName)]
        public static void ShowLogWindow()
        {
            LogDisplay.Show();
        }
    }
}
