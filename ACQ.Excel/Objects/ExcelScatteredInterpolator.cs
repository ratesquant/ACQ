using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ExcelDna.Integration;
using ExcelDna.Logging;

namespace ACQ.Excel.Objects
{
    public class ExcelScatteredInterpolator
    {
        private static readonly string m_tag = "#acqScatteredInterpolator";
   
        [ExcelFunction(Description = "Create Interpolator object", Category = AddInInfo.Category)]
        public static object acq_interpolator_scattered_create(
            [ExcelArgument(Description = "nodes")] double[,] x,
            [ExcelArgument(Description = "function values ")]  double[] y,
            [ExcelArgument(Description = "scales")]  object scales,
            [ExcelArgument(Description = "multiquadrics")] object method)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                return ACQ.Excel.Handles.GlobalCache.CreateHandle(m_tag, new object[] { x, y, scales, method, "acq_interpolator_scattered_create" },
                    (objectType, parameters) =>
                    {
                        ACQ.Math.Interpolation.ScatteredInterpolationInterface interpolator = construct_interpolator(x, y, scales, method);

                        if (interpolator == null)
                            return ExcelError.ExcelErrorNull;
                        else
                            return interpolator;
                    });
            }
        }


        [ExcelFunction(Description = "Evaluate Scattered interpolation at specified point", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_interpolator_scattered_eval(
            [ExcelArgument(Description = "Interpolator object")] string handle,
            [ExcelArgument(Description = "Interpolation point")] double[] x)
        {
            ACQ.Math.Interpolation.ScatteredInterpolationInterface interpolator;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Interpolation.ScatteredInterpolationInterface>(handle, out interpolator))
            {
                if (interpolator != null)
                {
                    return ExcelHelper.CheckNan(interpolator.Eval(x));
                }
            }
            return ExcelError.ExcelErrorRef;
        }

        [ExcelFunction(Description = "Evaluate Scattered interpolation at specified point", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_interpolator_scattered_eval_x5(
            [ExcelArgument(Description = "Interpolator object")] string handle,
            [ExcelArgument(Description = "Interpolation point")] object x1,
            [ExcelArgument(Description = "Interpolation point")] object x2,
            [ExcelArgument(Description = "Interpolation point")] object x3,
            [ExcelArgument(Description = "Interpolation point")] object x4,
            [ExcelArgument(Description = "Interpolation point")] object x5)
        {
            ACQ.Math.Interpolation.ScatteredInterpolationInterface interpolator;

            if (ACQ.Excel.Handles.GlobalCache.TryGetObject<ACQ.Math.Interpolation.ScatteredInterpolationInterface>(handle, out interpolator))
            {
                if (interpolator != null)
                {
                    //i am not sure how to support arbitrary number of optional arguments, so limit it to 5
                    object[] xn = new object[] { x1, x2, x3, x4, x5 };
                    List<double> xp = new List<double>();

                    for (int i = 0; i < xn.Length; i++)
                    {
                        if (xn[i] is ExcelMissing || !(xn[i] is double))
                        {
                            break;
                        }
                        else
                        {
                            xp.Add((double)xn[i]);
                        }
                    }

                   return ExcelHelper.CheckNan(interpolator.Eval(xp.ToArray()));
                }
            }
            return ExcelError.ExcelErrorRef;
        }


        private static ACQ.Math.Interpolation.ScatteredInterpolationInterface construct_interpolator(double[,] x, double[] y, object input_scales, object method)
        {
            ACQ.Math.Interpolation.ScatteredInterpolationInterface interpolator = null;
            try
            {
                ACQ.Math.Interpolation.enRadialBasisFunction rbf_function = ExcelHelper.CheckEnum<ACQ.Math.Interpolation.enRadialBasisFunction>(method, ACQ.Math.Interpolation.enRadialBasisFunction.Linear);

                double[] scales = null; //null argument is ok

                if (!(input_scales is ExcelMissing) && input_scales != null)
                {
                    if (input_scales is double)
                    {
                        double const_scale = (double)input_scales;
                        scales = new double[x.GetLength(1)];

                        for (int i = 0; i < scales.Length; i++)
                        {
                            scales[i] = const_scale;
                        }
                    }
                    else if (input_scales is object[])
                    {
                        object[] vect_scale = input_scales as object[];
                        scales = new double[x.GetLength(1)];

                        for (int i = 0; i < scales.Length; i++)
                        {
                            if (i < vect_scale.Length && vect_scale[i] is double)
                                scales[i] = (double)vect_scale[i];
                            else
                                scales[i] = 1.0;

                        }
 
                    }
                }

                interpolator = new ACQ.Math.Interpolation.RbfInterpolation(x, y, rbf_function, scales, 0.0);
            }
            catch (Exception ex)
            {
                LogDisplay.WriteLine("Error: " + ex.ToString());
            }
            return interpolator;
        }
    }
}
