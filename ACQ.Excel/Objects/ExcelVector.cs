using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class ExcelVector
    {
        private static readonly string m_tag = "#Vector";

        public static string Tag
        {
            get
            {
                return m_tag;
            }
        }

        [ExcelFunction(Description = "Create Vector Object", Category = AddInInfo.Category)]
        public static object acq_vector_create(double[] x)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelVector.Tag, new object[] { x, "acq_vector_create" },
                    (objectType, paramaters) =>
                    {
                        return x.Clone();

                    });
 
            }
        }

        public static object acq_vector_get_element(string handle, int index)
        {
            double[] vector;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<double[]>(handle, out vector))
            {
                if (vector != null && index >= 0 && index < vector.Length)
                {
                    return vector[index];
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        public static object acq_vector_get_size(string handle)
        {
            double[] vector;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<double[]>(handle, out vector))
            {
                if (vector != null)
                {
                    return vector.Length;
                }
            }
            return ExcelError.ExcelErrorRef;
        }
    }
}
