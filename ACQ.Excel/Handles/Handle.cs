using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using ExcelDna.Integration;

namespace ACQ.Excel.Handles
{
    interface IHandle
    {
        string Name { get; }
        object Value { get; }        
    }
    class Handle : IExcelObservable, IDisposable, IHandle
    {
        private static readonly object m_lock = new object();
        private static int m_index;

        private readonly HandleStorage m_storage;
        protected IExcelObserver m_observer;
        private readonly string m_name;
        protected object m_value;

        public Handle(HandleStorage storage, string tag, object value)
        {
            m_storage = storage;
            m_value = value;

            lock (m_lock)
            {
                m_name = String.Format("{0}:{1}", tag, m_index++);
            }
        }

        public virtual IDisposable Subscribe(IExcelObserver observer)
        {
            m_observer = observer;
            m_observer.OnNext(m_name);
            return this;
        }

        public void Dispose()
        {
            m_storage.Remove(this); 
        }

        public string Name
        {
            get
            {
                return m_name;
            }
        }

        public object Value
        {
            get 
            {
                return m_value;
            }
        }
    }

    internal class HandleAsync : Handle, IExcelObservable, IDisposable
    {
        private static readonly object m_sync = new object();
        private Task<object> m_valueAsync;

        public HandleAsync(HandleStorage storage, string tag, Task<object> valueAsync) : base(storage, tag, null)
        {
            m_valueAsync = valueAsync;
            valueAsync.ContinueWith(ProcessResult);
        }
        void ProcessResult(Task<object> valueAsync) 
        {
            lock (m_sync) 
            {
                try
                {
                    m_value = valueAsync.Result;
                    m_observer?.OnNext(Name);
                }
                catch (Exception ex)
                {
                    m_observer?.OnError(ex);
                }
            }
        }

        public override IDisposable Subscribe(IExcelObserver observer) 
        {
            lock (m_sync)
            { 
            }
            m_observer = observer;
            if (m_valueAsync.IsCompleted && m_valueAsync.Status == TaskStatus.RanToCompletion)
            {
                m_observer.OnNext(Name);
            }
            else if (m_valueAsync.IsFaulted)
            {
                m_observer.OnError(m_valueAsync.Exception);
            }
            return this;
        }
    }

}
