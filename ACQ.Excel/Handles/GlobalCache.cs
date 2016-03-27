using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using ExcelDna.Integration;

namespace ACQ.Excel.Handles
{
    class GlobalCache
    {
        private static HandleStorage m_storage = new HandleStorage();

        internal static object CreateHandle(string objectType, object[] parameters, Func<string, object[], object> maker)
        {
            return m_storage.CreateHandle(objectType, parameters, maker); 
        }

        internal static bool TryGetObject<T>(string name, out T value) where T : class
        {
            return m_storage.TryGetObject<T>(name, out value);
        }
    }
}
