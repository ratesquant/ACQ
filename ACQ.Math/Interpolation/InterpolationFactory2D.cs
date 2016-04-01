using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;

namespace ACQ.Math.Interpolation
{
    public class InterpolationFactory2D
    {
        private static Dictionary<string, Type> m_interpolation_types = new Dictionary<string, Type>();

        static InterpolationFactory2D()
        {
            Type[] types = GetClassTypes(System.Reflection.Assembly.GetExecutingAssembly(), typeof(InterpolationFactory2D).Namespace);

            Type base_type = typeof(InterpolationInterface2D);

            foreach(Type t in types)
            {
                if (!t.IsAbstract && base_type.IsAssignableFrom(t))
                {
                    m_interpolation_types[t.FullName.ToLower()] = t;
                }
            }
        }

        private static Type[] GetClassTypes(Assembly assembly, string nameSpace)
        {
            return assembly.GetTypes().Where(t => t.IsClass && String.Equals(t.Namespace, nameSpace)).ToArray();
        }

        public static Type GetInterpolationType(string method)
        {
            string name = String.Format("ACQ.Math.Interpolation.{0}Interpolation", method).ToLower();

            Type result = null;

            if (m_interpolation_types.ContainsKey(name))
            {
                result = m_interpolation_types[name];
            }
            return result;
        }

        public static InterpolationInterface2D GetInterpolator(string method, params object[] arguments)
        {
            InterpolationInterface2D interpolator = null;

            Type interpolator_type = GetInterpolationType(method);

            if (interpolator_type != null)
            {
                interpolator = Activator.CreateInstance(interpolator_type, arguments) as InterpolationInterface2D;
            }

            return interpolator;
        }
    }
}
