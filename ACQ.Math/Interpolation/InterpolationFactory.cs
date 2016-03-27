using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Interpolation
{
    public enum enInterpolationMethod
    {
        Linear,
        Cubic,
        Hermit,
        Akima
    }

    public class InterpolationFactory
    {
        public static Type GetInterpolationType(enInterpolationMethod method)
        {
            Type interpolation_type = GetInterpolationType(method.ToString());

            return interpolation_type;
        }

        public static Type GetInterpolationType(string method)
        {
            return Type.GetType(String.Format("ACQ.Math.Interpolation.{0}Interpolation", method), false, true); //return null if missing, ignore case
        }

        public static InterpolationInterface GetInterpolator(Type type, double[] x, double[] y, bool bounds)
        {
            InterpolationInterface interpolator = Activator.CreateInstance(type, x, y, bounds) as InterpolationInterface;

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
