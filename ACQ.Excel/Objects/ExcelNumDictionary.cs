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
        public static object acq_numdict_create(object[] keys, double[] values)
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
                            dict[keys[i].ToString()] = values[i];
                        }

                        return dict;

                    });

            }
        }
        [ExcelFunction(Description = "Create Numeric Dictionary Object", Category = AddInInfo.Category)]
        public static object acq_numdict_create_fake_data(string name)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
            {
                return ExcelError.ExcelErrorRef;
            }
            else if (name == null)
            {
                return ExcelError.ExcelErrorRef;
            }
            else
            {
               return ACQ.Excel.Handles.GlobalCache.CreateHandle(Tag, new object[] { name, "acq_numdict_create_fake_data" },
                    (objectType, parameters) =>
                    {
                        NumDictionary dict = new NumDictionary();

                        System.Threading.Thread.Sleep(2000); //simulate loading from DB

                        for (int i = 0; i < 100; i++)
                        {

                            dict[String.Format("{0}-{1}", name, i)] = i;
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
