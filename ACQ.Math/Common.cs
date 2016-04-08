using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;

namespace ACQ.Math
{
    public class Common
    {
        public static Type[] GetClassTypes(Assembly assembly, string nameSpace)
        {
            return assembly.GetTypes().Where(t => t.IsClass && String.Equals(t.Namespace, nameSpace)).ToArray();
        }
    }
}
