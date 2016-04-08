using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class ExcelHashtable
    {
        private static readonly string m_tag = "#acqHashtable";

        public static string Tag
        {
            get
            {
                return m_tag;
            }
        }

        [ExcelFunction(Description = "Create Hashtable Object", Category = AddInInfo.Category)]
        public static object acq_hashtable_create(object[] keys, object[] values)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
            {
                return ExcelError.ExcelErrorRef;
            }
            else if (keys == null || values == null)
            {
                return ExcelError.ExcelErrorRef;
            }
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(ExcelHashtable.Tag, new object[] { keys, values, "acq_hashtable_create" },
                    (objectType, parameters) =>
                    {
                        Hashtable htable = new Hashtable();

                        for (int i = 0; i < keys.Length; i++)
                        {
                            htable[keys[i]] = values[i];
                        }

                        return htable;

                    });

            }
        }

        [ExcelFunction(Description = "Get hashtable element", Category = AddInInfo.Category)]
        public static object acq_hashtable_element(string handle, object key)
        {
            Hashtable htable;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<Hashtable>(handle, out htable))
            {
                if (htable != null && htable.ContainsKey(key))
                {
                    return htable[key];
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Get hashtable size", Category = AddInInfo.Category)]
        public static object acq_hashtable_size(string handle)
        {
            Hashtable htable;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<Hashtable>(handle, out htable))
            {
                if (htable != null)
                {
                    return htable.Count;
                }
            }
            return ExcelError.ExcelErrorRef;
        }
    }
}
