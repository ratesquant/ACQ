using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class ExcelArray
    {
        private static readonly string m_tag = "#acqArray";

        public static string Tag
        {
            get
            {
                return m_tag;
            }
        }

        [ExcelFunction(Description = "Create Array Object", Category = AddInInfo.Category)]
        public static object acq_array_create(object[] x)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
            {
                return ExcelError.ExcelErrorRef;
            }
            else if (x == null)
            {
                return ExcelError.ExcelErrorRef;
            }
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelArray.Tag, new object[] { x, "acq_array_create" },
                    (objectType, parameters) =>
                    {
                        object[] a = (object[])x.Clone(); //shallow copy should be enought, becuase we are going to regenerate it every time argument changes 

                        return a;

                    });

            }
        }

        [ExcelFunction(Description = "Get array element", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_array_element(string handle, int index)
        {
            object[] array;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<object[]>(handle, out array))
            {
                if (array != null && index >= 0 && index < array.Length)
                {
                    return array[index];
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Get array size", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_array_size(string handle)
        {
            object[] array;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<object[]>(handle, out array))
            {
                if (array != null)
                {
                    return array.Length;
                }
            }
            return ExcelError.ExcelErrorRef;
        }
    }
}
