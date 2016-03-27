using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using ExcelDna.Integration;

namespace ACQ.Excel.Handles
{
    class HandleStorage
    {
        private ReaderWriterLockSlim m_lock = new ReaderWriterLockSlim();
        private Dictionary<string, Handle> m_storage = new Dictionary<string, Handle>();

        internal object CreateHandle(string objectType, object[] parameters, Func<string, object[], object> maker)
        {
            return ExcelAsyncUtil.Observe(objectType, parameters, () =>
            {
                var value = maker(objectType, parameters);
                var handle = new Handle(this, objectType, value);

                m_lock.EnterWriteLock();

                try
                {
                    m_storage.Add(handle.Name, handle);
                }
                finally
                {
                    m_lock.ExitWriteLock();
                }
                return handle;

            });

        }

        internal bool TryGetObject(string name, out object value)
        {
            bool found = false;

            value = null;

            m_lock.EnterReadLock();

            try
            {
                Handle handle;

                if (m_storage.TryGetValue(name, out handle))
                {
                    value = handle.Value;
                    found = true;
                }
            }
            finally
            {
                m_lock.ExitReadLock();
            }
            return found;
        }

        internal bool TryGetObject<T>(string name, out T value) where T : class
        {
            bool found = false;

            value = default(T);

            m_lock.EnterReadLock();

            try
            {
                Handle handle;

                if (m_storage.TryGetValue(name, out handle))
                {
                    if (handle.Value is T)
                    {
                        value = handle.Value as T;
                        found = true;
                    }
                }
            }
            finally
            {
                m_lock.ExitReadLock();
            }
            return found;
        }

        internal void Remove(Handle handle)
        {
            object value;

            if (TryGetObject(handle.Name, out value))
            {
                m_lock.EnterWriteLock();

                try
                {
                    m_storage.Remove(handle.Name);
                }
                finally
                {
                    m_lock.ExitWriteLock();
                }

                var disp = value as IDisposable;

                if (disp != null)
                {
                    disp.Dispose();
                }
            }
        }
    }
}
