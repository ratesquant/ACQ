using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;

namespace ACQ.Excel.Objects
{
    public class NumDictionary : Dictionary<string, double> { };

    /// <summary>
    /// Excel Numeric dictionary: keys -  string, values - double
    /// </summary>
    public class ExcelNumDictionary
    {
        private static readonly string m_tag = "#acqNumDict";

        public static string Tag
        {
            get
            {
                return m_tag;
            }
        }

        [ExcelFunction(Description = "Create Numeric Dictionary Object", Category = AddInInfo.Category)]
        public static object acq_numdict_create(string[] keys, double[] values)
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
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(Tag, new object[] { keys, values, "acq_numdict_create" },
                    (objectType, parameters) =>
                    {
                        NumDictionary dict = new NumDictionary();

                        for (int i = 0; i < keys.Length; i++)
                        {
                            dict[keys[i]] = values[i];
                        }

                        return dict;

                    });

            }
        }

        [ExcelFunction(Description = "Get Dictionary element", Category = AddInInfo.Category)]
        public static object acq_dict_element(string handle, string key)
        {
            NumDictionary dict;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<NumDictionary>(handle, out dict))
            {
                if (dict != null && dict.ContainsKey(key))
                {
                    return dict[key];
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Get hashtable size", Category = AddInInfo.Category)]
        public static object acq_dict_size(string handle)
        {
            NumDictionary dict;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<NumDictionary>(handle, out dict))
            {
                if (dict != null)
                {
                    return dict.Count;
                }
            }
            return ExcelError.ExcelErrorRef;
        }
    }
}
