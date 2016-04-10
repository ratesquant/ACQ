using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;

namespace ACQ.Excel.Introspection
{
    public class ObjectIntrospector
    {
        object m_acq_object;

        public ObjectIntrospector(object acq_object)
        {
            m_acq_object = acq_object;
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
                    //object specific code goes here, all new objects need to have ToTable function for introspection to work
                    if (m_acq_object is ACQ.Math.Linalg.Vector)
                    {
                        table = VectorToTable(m_acq_object as ACQ.Math.Linalg.Vector);
                    }
                    else if (m_acq_object is ACQ.Math.Linalg.Matrix)
                    {
                        table = MatrixToTable(m_acq_object as ACQ.Math.Linalg.Matrix);
                    }
                    else if (m_acq_object is Hashtable)
                    {
                        table = HashtableToTable(m_acq_object as Hashtable);
                    }
                    else if (m_acq_object is object[])
                    {
                        table = ArrayToTable(m_acq_object as object[]);
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

        public bool IsDataObject
        {
            get
            {
                return m_acq_object is ACQ.Math.Linalg.Vector ||
                       m_acq_object is ACQ.Math.Linalg.Matrix ||
                       m_acq_object is System.Collections.Hashtable ||
                       m_acq_object is object[];
            }
        }

        #region Private Convertion Methods
        private DataTable VectorToTable(ACQ.Math.Linalg.Vector vector)
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

        private DataTable MatrixToTable(ACQ.Math.Linalg.Matrix matrix)
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

        private DataTable HashtableToTable(Hashtable ht)
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

        private DataTable ArrayToTable(object[] array)
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
