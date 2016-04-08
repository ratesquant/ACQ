using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class ExcelMatrix
    {
        private static readonly string m_tag = "#acqMatrix";

        public static string Tag
        {
            get
            {
                return m_tag;
            }
        }

        [ExcelFunction(Description = "Create matrix object", Category = AddInInfo.Category)]
        public static object acq_matrix_create(double[,] x)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelMatrix.Tag, new object[] { x, "acq_matrix_create" },
                    (objectType, parameters) =>
                    {
                        return new ACQ.Math.Linalg.Matrix(x);

                    });

            }
        }

        [ExcelFunction(Description = "Get matrix element", Category = AddInInfo.Category)]
        public static object acq_matrix_element(string handle, int row, int column)
        {
            ACQ.Math.Linalg.Matrix matrix;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Linalg.Matrix>(handle, out matrix))
            {
                if (matrix != null && matrix.CheckIndexes(row, column))
                {
                    return matrix[row, column];
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Get matrix rows", Category = AddInInfo.Category)]
        public static object acq_matrix_rows(string handle)
        {
            ACQ.Math.Linalg.Matrix matrix;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Linalg.Matrix>(handle, out matrix))
            {
                if (matrix != null)
                {
                    return matrix.Rows;
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Get matrix columns", Category = AddInInfo.Category)]
        public static object acq_matrix_columns(string handle)
        {
            ACQ.Math.Linalg.Matrix matrix;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Linalg.Matrix>(handle, out matrix))
            {
                if (matrix != null)
                {
                    return matrix.Columns;
                }
            }
            return ExcelError.ExcelErrorRef;
        }
    }
}
