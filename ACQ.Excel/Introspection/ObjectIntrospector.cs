using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;
using System.Data;

namespace ACQ.Excel.Introspection
{
    public class ObjectIntrospector
    {
        object m_acq_object;
        MethodInfo m_toDataTable = null;

        private readonly static Dictionary<Type, MethodInfo> m_toTableMethods = new Dictionary<Type, MethodInfo>();

        static ObjectIntrospector()
        {
            //init ToDataTable methods
            foreach (MethodInfo method in typeof(ObjectIntrospector).GetMethods(BindingFlags.Static | BindingFlags.NonPublic))
            {
                ParameterInfo[] paramInfo =  method.GetParameters();
                if (paramInfo.Length == 1 && method.ReturnType.Equals(typeof(DataTable))) //method.Name == ToTable ?
                {   
                    Type arg = paramInfo[0].ParameterType;
                    m_toTableMethods[arg] = method;
                }
            }
        }

        public ObjectIntrospector(object acq_object)
        {
            m_acq_object = acq_object;

            if (m_acq_object != null)
            {
                Type arg_type = m_acq_object.GetType();

                //unfortunatly we have to iterate over all elements because types dont match exactly due to inheritance  
                foreach (Type base_type in m_toTableMethods.Keys)
                {
                    if (base_type.IsAssignableFrom(arg_type))
                    {
                        m_toDataTable = m_toTableMethods[base_type];
                        break;
                    }
                }
            }
        }

        public string Name
        {
            get
            {
                string name = "NULL";

                if (m_acq_object != null)
                {
                    name = m_acq_object.GetType().ToString(); 
                }

                return name;
            }
        }

        public override string ToString()
        {
            string name = "NULL";

            if (m_acq_object != null)
            {
                name = m_acq_object.ToString();
            }

            return name;
        }

        public DataTable Data
        {
            get
            {
                DataTable table = new DataTable();
 
                if (m_acq_object != null)
                {                  
                    if (m_toDataTable != null)
                    {
                        table = m_toDataTable.Invoke(null, new object[] {m_acq_object}) as DataTable;
                    }                 
                }

                return table;
 
            }
        }

        public object SelectedObject
        {
            get
            {
                return m_acq_object;
            }
        }

        public bool IsDataTableConvertable
        {
            get
            {
                return m_toDataTable != null;
            }
        }

        #region ToTable Methods
        private static DataTable ToDataTable(ACQ.Math.Linalg.Vector vector)
        {
            DataTable table = new DataTable(vector.GetType().ToString());
         
            DataColumn column;
            DataRow row;

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.Double");
            column.ColumnName = "value";
            column.ReadOnly = true;
            table.Columns.Add(column);
         
            // Create three new DataRow objects and add 
            // them to the DataTable
            for (int i = 0; i < vector.Size; i++)
            {
                row = table.NewRow();
                row["value"] = vector[i];
                table.Rows.Add(row);
            }
            return table;
        }

        private static DataTable ToDataTable(ACQ.Math.Interpolation.InterpolationBase interpolator)
        {
            DataTable table = new DataTable(interpolator.GetType().ToString());
         
            DataColumn column;
            DataRow row;

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.Double");
            column.ColumnName = "x";
            column.ReadOnly = true;
            table.Columns.Add(column);

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.Double");
            column.ColumnName = "y";
            column.ReadOnly = true;
            table.Columns.Add(column);
         
            // Create three new DataRow objects and add 
            // them to the DataTable
            for (int i = 0; i < interpolator.Size; i++)
            {
                Tuple<double, double> node = interpolator.GetNode(i);

                row = table.NewRow();
                row["x"] = node.Item1;
                row["y"] = node.Item2;
                table.Rows.Add(row);
            }
            return table;
        }

        private static DataTable ToDataTable(ACQ.Math.Linalg.Matrix matrix)
        {
            DataTable table = new DataTable(matrix.GetType().ToString());

            DataColumn column;
            DataRow row;

            string[] cnames = new string[matrix.Columns];
            // Create new DataColumn, set DataType, 
            // ColumnName and add to DataTable.    
            for (int i = 0; i < matrix.Columns; i++)
            {
                cnames[i] = i.ToString();
                column = new DataColumn();
                column.DataType = System.Type.GetType("System.Double");
                column.ColumnName = cnames[i];
                column.ReadOnly = true;
                table.Columns.Add(column);
            }

            // Create three new DataRow objects and add 
            // them to the DataTable
            for (int i = 0; i < matrix.Rows; i++)
            {
                row = table.NewRow();
                for (int j = 0; j < matrix.Columns; j++)
                {
                    row[cnames[j]] = matrix[i, j];
                }
                table.Rows.Add(row);
            }
            return table;
        }

        private static DataTable ToDataTable(Hashtable ht)
        {
            DataTable table = new DataTable(ht.GetType().ToString());

            DataColumn column;
            DataRow row;

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.String");
            column.ColumnName = "key";
            column.ReadOnly = true;
            table.Columns.Add(column);

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.String");
            column.ColumnName = "value";
            column.ReadOnly = true;
            table.Columns.Add(column);

            // Create three new DataRow objects and add 
            // them to the DataTable
            foreach (object key in ht.Keys)
            {
                row = table.NewRow();
                row["key"] = key.ToString();
                row["value"] = ht[key].ToString();
                table.Rows.Add(row);
            }
            return table;
        }

        private static DataTable ToDataTable(Dictionary<string, object> dict)
        {
            DataTable table = new DataTable(dict.GetType().ToString());

            DataColumn column;
            DataRow row;

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.String");
            column.ColumnName = "key";
            column.ReadOnly = true;
            table.Columns.Add(column);

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.String");
            column.ColumnName = "value";
            column.ReadOnly = true;
            table.Columns.Add(column);

            // Create three new DataRow objects and add 
            // them to the DataTable
            foreach (KeyValuePair<string, object> pair in dict)
            {
                row = table.NewRow();
                row["key"] = pair.Key;
                row["value"] = pair.Value!=null ? pair.Value.ToString() : "NULL";
                table.Rows.Add(row);
            }
            return table;
        }

        private static DataTable ToDataTable(ACQ.Math.Regression.IRegressionSummary regression)
        {
            return ToDataTable(regression.Summary);
        }

        private static DataTable ToDataTable(Dictionary<string, double> dict)
        {
            DataTable table = new DataTable(dict.GetType().ToString());

            DataColumn column;
            DataRow row;

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.String");
            column.ColumnName = "key";
            column.ReadOnly = true;
            table.Columns.Add(column);

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.Double");
            column.ColumnName = "value";
            column.ReadOnly = true;
            table.Columns.Add(column);

            // Create three new DataRow objects and add 
            // them to the DataTable
            foreach (KeyValuePair<string, double> pair in dict)
            {
                row = table.NewRow();
                row["key"] = pair.Key;
                row["value"] = pair.Value;
                table.Rows.Add(row);
            }
            return table;
        }

        private static DataTable ToDataTable(object[] array)
        {
            DataTable table = new DataTable(array.GetType().ToString());

            DataColumn column;
            DataRow row;

            column = new DataColumn();
            column.DataType = System.Type.GetType("System.String");
            column.ColumnName = "value";
            column.ReadOnly = true;
            table.Columns.Add(column);

            // Create three new DataRow objects and add 
            // them to the DataTable
            for (int i = 0; i < array.Length; i++)
            {
                row = table.NewRow();
                row["value"] = array[i].ToString();
                table.Rows.Add(row);
            }
            return table;
        }
        #endregion
    }
}
