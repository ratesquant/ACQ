using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public class InterpolationFactory
    {
        private static Dictionary<string, Type> m_interpolation_types = new Dictionary<string, Type>();

        static InterpolationFactory()
        {
            Type base_type = typeof(InterpolationInterface);

            Type[] types = Common.GetClassTypes(System.Reflection.Assembly.GetExecutingAssembly(), base_type.Namespace);

            foreach(Type t in types)
            {
                if (!t.IsAbstract && base_type.IsAssignableFrom(t))
                {
                    m_interpolation_types[t.FullName.ToLower()] = t;
                }
            }
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
    //        return Type.GetType(String.Format("ACQ.Math.Interpolation.{0}Interpolation", method), false, true); //return null if missing, ignore case
        }

        public static InterpolationInterface GetInterpolator(Type type, double[] x, double[] y)
        {
            InterpolationInterface interpolator = Activator.CreateInstance(type, x, y) as InterpolationInterface;

            return interpolator; 
        }

        public static InterpolationInterface GetInterpolator(Type type, params object[] arguments)
        {
            InterpolationInterface interpolator = Activator.CreateInstance(type, arguments) as InterpolationInterface;

            return interpolator;
        }

        public static InterpolationInterface GetInterpolator(string method, params object[] arguments)
        {
            InterpolationInterface interpolator = null;

            Type interpolator_type = GetInterpolationType(method);

            if (interpolator_type != null)
            {
                interpolator = Activator.CreateInstance(interpolator_type, arguments) as InterpolationInterface;
            }

            return interpolator;
        }
    }
}
