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

        internal static object CreateHandleAsync(string objectType, object[] parameters, Func<string, object[], Task<object>> maker)
        {
            return m_storage.CreateHandleAsync(objectType, parameters, maker);
        }

        internal static bool TryGetObject<T>(string name, out T value) where T : class
        {
            return m_storage.TryGetObject<T>(name, out value);
        }

        internal static Tuple<bool, TResult> TryReadObject<T, TResult>(string name, Func<T, TResult> reader) where T : class
        {
            return m_storage.TryReadObject<T, TResult>(name, reader);
        }

        /// <summary>
        /// read object with argument 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <typeparam name="TResult"></typeparam>
        /// <typeparam name="TArg"></typeparam>
        /// <param name="name"></param>
        /// <param name="reader"></param>
        /// <param name="argument"></param>
        /// <returns></returns>
        internal static Tuple<bool, TResult> TryReadObject<T, TResult, TArg>(string name, Func<T, TArg, TResult> reader, TArg argument) where T : class
        {
            return m_storage.TryReadObject<T, TResult, TArg>(name, reader, argument); 
        }
    }
}
