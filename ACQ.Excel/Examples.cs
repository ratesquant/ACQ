using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


using ExcelDna.Integration;

namespace ACQ.Excel
{
    public static class Examples
    {
        #region Sudoku
        [ExcelFunction(Description = "Solve Sudoku puzzle", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object[,] acq_sudoku_generate(object seed)
        {
            int[,] grid;

            if (seed != null && seed is double)
            {
                grid = ACQ.Math.Sudoku.Generate((int)(double)seed);
            }
            else
            {
                grid = ACQ.Math.Sudoku.Generate();
            }
            
            object[,] result = new object[grid.GetLength(0), grid.GetLength(1)];

            for (int i = 0; i < grid.GetLength(0); i++)
            {
                for (int j = 0; j < grid.GetLength(1); j++)
                {
                    if (grid[i, j] == 0)
                        result[i, j] = "";//ExcelEmpty.Value;
                    else
                        result[i, j] = grid[i, j];
                }
            }

            return result;
        }

        [ExcelFunction(Description = "Solve Sudoku puzzle", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object[,] acq_sudoku_solve(object[,] grid)
        {
            const int size = 9;

            int[,] sudoku = acq_sudoku_convertgrid(grid, size);

            if (sudoku != null)
            {
                List<int[,]> solutions = new List<int[,]>();

                ACQ.Math.Sudoku.Solve(sudoku, solutions, 1);

                if (solutions.Count > 0)
                {
                    return ExcelHelper.BoxArray(solutions[0]);
                }
            }
            
            return ExcelHelper.CreateArray(1, 1, ExcelError.ExcelErrorNull);
        }

        [ExcelFunction(Description = "Number of possible solutions", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static int acq_sudoku_solution_count(object[,] grid)
        {
            const int size = 9;
            const int max_count = 1024;

            int[,] sudoku = acq_sudoku_convertgrid(grid, size);

            int count = 0;

            if (sudoku != null)
            {
                List<int[,]> solutions = new List<int[,]>();

                ACQ.Math.Sudoku.Solve(sudoku, solutions, max_count);

                count = solutions.Count;                
            }

            return count;
        }

        private static int[,] acq_sudoku_convertgrid(object[,] grid, int size)
        {            
            int[,] sudoku = null;

            if (grid != null && grid.GetLength(0) == size && grid.GetLength(1) == size)
            {
                sudoku = new int[size, size];

                for (int i = 0; i < size; i++)
                {
                    for (int j = 0; j < size; j++)
                    {
                        object v = grid[i, j];
                        int sv = 0;

                        if (v is double)
                        {
                            sv = (int)((double)v);
                        }
                        else if (v is string)
                        {
                            double res;

                            if (Double.TryParse(v as string, out res))
                            {
                                sv = (int)res;
                            }
                        }
                        sudoku[i, j] = System.Math.Min(sv, size);
                    }
                }                
            }
            return sudoku;
        }
        #endregion

        [ExcelFunction(Description = "acq_async_load_example", Category = AddInInfo.Category)]
        public static object acq_async_load_example(string name)
        {
            object asyncResult = ExcelAsyncUtil.Run("AsyncLoadExample", name,
                delegate()
                {
                    System.Threading.Thread.Sleep(2000);
                    return String.Format("{0} was loaded",name);
                });

            if (asyncResult.Equals(ExcelError.ExcelErrorNA))
            {
                return ExcelError.ExcelErrorGettingData;
            }
            return asyncResult;
        }

        [ExcelFunction(Description = "acq_async_load_example2", Category = AddInInfo.Category)]
        public static void acq_async_load_example2(string name, ExcelAsyncHandle asyncHandle)
        {
            int managedThreadId = System.Threading.Thread.CurrentThread.ManagedThreadId;
            System.Threading.ThreadPool.QueueUserWorkItem(delegate(object state)
            {
                System.Threading.Thread.Sleep(2000);
                int completedThreadId = System.Threading.Thread.CurrentThread.ManagedThreadId;

                asyncHandle.SetResult(String.Format("{0} was loaded", name));
            });
        }

        /*
        [ExcelFunction(Description = "acq_async_load_example3", Category = AddInInfo.Category)]
        public static object acq_async_load_example3(string name)
        {
            var callingCell = (ExcelReference)XlCall.Excel(XlCall.xlfCaller);

            return ExcelAsyncUtil.Observe("acq_async_load_example3", new object[] { name },
                () =>
                {
                    //var longTask = Task.Run(delegate
                    var longTask = System.Threading.Tasks.Task.Factory.StartNew(()=>
                    {
                        System.Threading.Thread.Sleep(3000);
                        return String.Format("{0} was loaded", name);
                    });
                    return longTask.ToExcelObservable();
                });
        }*/

        /*
        [ExcelFunction(Description = "acq_sum", Category = AddInInfo.Category)]
        public static double acq_sum(double a, double b)
        {
            return a + b;
        }         

        [ExcelFunction(Description = "Compute interpolated value", Category = AddInInfo.Category, IsThreadSafe = true)]
        public static object acq_linear_interpolation(double xi, double[] x, double[] y, object bounds)
        {
            if (ExcelDnaUtil.IsInFunctionWizard())
                return ExcelError.ExcelErrorRef;
            else
            {
                ACQ.Math.Interpolation.InterpolationInterface interpolator = null;
                try
                {
                    interpolator = new ACQ.Math.Interpolation.LinearInterpolation(x, y); //create linear interpolator 
                    interpolator.Bounds = ExcelHelper.CheckValue(bounds, true);
                }
                catch (Exception ex)
                {
                    //LogDisplay.WriteLine("Error: " + ex.ToString());                    
                }

                if (interpolator != null)
                    return ExcelHelper.CheckNan(interpolator.Eval(xi)); //interpolate
                else
                    return ExcelError.ExcelErrorNA;
            }
        }
         */
    }
}
