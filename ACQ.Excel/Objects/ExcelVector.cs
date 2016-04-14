using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class ExcelVector
    {
        private static readonly string m_tag = "#acqVector";

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
                    (objectType, parameters) =>
                    {
                        return new ACQ.Math.Linalg.Vector(x);

                    });
 
            }
        }

        [ExcelFunction(Description = "Get vector element", Category = AddInInfo.Category, IsThreadSafe = false)]
        public static object acq_vector_element(string handle, int index)
        {
            ACQ.Math.Linalg.Vector vector;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Linalg.Vector>(handle, out vector))
            {
                if (vector != null && index >= 0 && index < vector.Size)
                {
                    return vector[index];
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Get vector size", Category = AddInInfo.Category, IsThreadSafe = false)]
        public static object acq_vector_size(string handle)
        {
            ACQ.Math.Linalg.Vector vector;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Linalg.Vector>(handle, out vector))
            {
                if (vector != null)
                {
                    return vector.Size;
                }
                else
                {
                    return ExcelError.ExcelErrorNA;
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Scale Vector", Category = AddInInfo.Category)]
        public static object acq_vector_scale(string handle, double scale)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                ACQ.Math.Linalg.Vector vector;

                if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Linalg.Vector>(handle, out vector))
                {
                    if (vector != null)
                    {
                        return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelVector.Tag, new object[] { handle, scale, "acq_vector_scale" },
                          (objectType, parameters) =>
                          {
                              return scale * vector;

                          });
                    }
                }
                return ExcelError.ExcelErrorRef;
            }
        }
    }
}
