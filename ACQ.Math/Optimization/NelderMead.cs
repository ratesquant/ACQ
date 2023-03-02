using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math.Optimization
{
	/// <summary>
	/// Gao, F. and Han, L.	Implementing the Nelder-Mead simplex algorithm with adaptive	parameters. 2012. Computational Optimization and Applications.  51:1, pp. 259-277
	/// https://github.com/yw5aj/nelder_mead_scipy/blob/master/nelder_mead_scipy.py
	/// https://www.mathworks.com/help/matlab/ref/fminsearch.html
	/// https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
	/// </summary>
	public class NelderMead
    {
		bool m_adaptive = false;
		double m_fatol = 1e-8;
		double m_xatol = 1e-4;

		int m_maxItPerDim = 200;
		double m_alpha = 1d; // alpha > 0  (reflection)
		double m_gamma = 2d; // gamma > 1 (expansion)
		double m_rho = 0.5;
		double m_sigma = 0.5;

		public delegate double ObjectiveFunction(double[] x);

		public NelderMead(bool adaptive = false)
		{
			m_adaptive = adaptive;			
		}

		private double MaxSimplexSize(double[][] simplex, int from_index, int n, int dim)
		{
			double max_size = 0;

			for (int i = 0; i < n; i++)
			{
				if (i != from_index)
				{
					double size = 0;
					for (int j = 0; j < dim; j++)
					{
						size += System.Math.Abs(simplex[from_index][j] - simplex[i][j]);
					}

					if (size > max_size)					
					{
						max_size = size;
					}					
				}
			}
			return max_size;
		}

		public double[] FindMinimum(double[] x0, ObjectiveFunction obj)
		{
			int N = x0.Length;

			if (m_adaptive)
			{
				m_alpha = 1d; //alpha, reflection
				m_gamma = 1 + 2d / N; //beta, Expansion
				m_rho = 0.75 - 1d / (2d * N); //gamma, contraction 
				m_sigma = 1 - 1d / N; //delta, Shrink
			}


			double[][] simplex = null;

			int maxiter = N * m_maxItPerDim;

			double nonzdelt = 0.05;
			double zdelt = 0.00025;

			if (simplex == null) // no initial_simplex
			{
				simplex = new double[N + 1][];

				for (int j = 0; j <= N; j++)
				{
					simplex[j] = new double[N];
				}

				for (int i = 0; i < N; i++)
				{
					simplex[0][i] = x0[i];
				}

				for (int k = 0; k < N; k++)
				{
					for (int i = 0; i < N; i++)
					{
						if (i != k)
						{
							simplex[k + 1][i] = x0[i];
						} else
						{
							simplex[k + 1][i] = (1 + nonzdelt) * x0[k] + ((x0[k] == 0) ? zdelt : 0d);							
						}
					}
				}
			}

			double[] fsim = new double[N + 1];
			int[] indices = new int[N + 1];

			for (int k = 0; k < N + 1; k++)
			{
				fsim[k] = obj(simplex[k]);
				indices[k] = k;
			}

			double[] centroid = new double[N];
			double[] x_reflection = new double[N];
			double[] x_expansion = new double[N];
			double[] x_contraction = new double[N];

			for(int iter = 0; iter < maxiter; iter++)			
			{
				Array.Sort(fsim, indices);

				//Console.WriteLine("{0}: {1}", iter, String.Join(",", fsim.Select(p => p.ToString()).ToArray()));

				double[] x_worst = simplex[indices[N]];

				// Check for convergence (check distance from first vertex and function delta)
				if ( System.Math.Abs(fsim[0] - fsim[N] ) < m_fatol || MaxSimplexSize(simplex, indices[0], N+1, N) < m_xatol)
						break;
				
				// Find centroid of the simplex excluding the vertex with highest function value (x_worst).			
				for (int i = 0; i < N; i++)
				{
					centroid[i] = 0;
					for (int j = 0; j < N; j++)
					{
						centroid[i] += simplex[indices[j]][i];						
					}
					centroid[i] /= N;
				}

				//Reflection				
				for (int i = 0; i < N; i++)
				{
					x_reflection[i] = centroid[i] + m_alpha * (centroid[i] - x_worst[i]);					
				}
				double f_reflection = obj(x_reflection);

				bool doshrink = false;

				// Expansion
				if (f_reflection < fsim[0])
				{					
					for (int i = 0; i < N; i++)
					{
						x_expansion[i] = centroid[i] + m_gamma * (x_reflection[i] - centroid[i]);
					}
					double f_expansion = obj(x_expansion);

					if (f_expansion < f_reflection)
					{
						fsim[N] = f_expansion;
						for (int i = 0; i < N; i++)
						{
							x_worst[i] = x_expansion[i];
						}
					}
					else
					{
						fsim[N] = f_reflection;
						for (int i = 0; i < N; i++)
						{
							x_worst[i] = x_reflection[i];
						}						
					}
				}
				else //f_reflection >= fsim[0]
				{
					if (f_reflection < fsim[N - 1])
					{						
						for (int i = 0; i < N; i++)
						{
							x_worst[i] = x_reflection[i];
						}

						fsim[N] = f_reflection;
					}
					else//f_reflection >= fsim[N - 1]
					{
						if (f_reflection < fsim[N])
						{
							// Contraction							
							for (int i = 0; i < N; i++)
							{
								x_contraction[i] = centroid[i] + m_rho * (x_reflection[i] - centroid[i]);
							}

							double f_contraction = obj(x_contraction);

							if (f_contraction <= f_reflection)
							{								
								fsim[N] = f_contraction;
								for (int i = 0; i < N; i++)
								{
									x_worst[i] = x_contraction[i];
								}
							}
							else
							{
								doshrink = true;
							}
						}else
						{
							// Inside Contraction							
							for (int i = 0; i < N; i++)
							{
								x_contraction[i] = centroid[i] - m_rho * (x_reflection[i] - centroid[i]);
							}

							double f_contraction = obj(x_contraction);

							if (f_contraction < fsim[indices[N]])
							{								
								fsim[N] = f_contraction;
								for (int i = 0; i < N; i++)
								{
									x_worst[i] = x_contraction[i];
								}
							}
							else
							{
								doshrink = true;
							}
						}
					}
				}


				//Shrink
				if (doshrink)
				{
					//Shrink
					double[] best_point = simplex[indices[0]];
					for (int j = 1; j <= N; j++)
					{
						double[] vertex = simplex[indices[j]];
						for (int i = 0; i < N; i++)
						{

							vertex[i] = best_point[i] + m_sigma * (vertex[i] - best_point[i]);
						}

						fsim[j] = obj(vertex);
					}
				}
			}

			return simplex[indices[0]];

		}

		public double[] FindMinimumEx(double[] x0, ObjectiveFunction obj)
		{
			int N = x0.Length;

			double alpha, gamma, rho, sigma; //α, β, γ , δ

			double[][] simplex = null;

			if (m_adaptive)
			{
				alpha = 1d; //alpha, reflection
				gamma = 1 + 2d / N; //beta, Expansion
				rho = 0.75 - 1d / (2d * N); //gamma, contraction 
				sigma = 1 - 1d / N; //delta, Shrink
			}
			else
			{ //alpha, beta gamma, delta paramaters
				alpha = 1d;
				gamma = 2d;
				rho = 0.5;
				sigma = 0.5;
			}

			int maxiter = N * 200;

			double nonzdelt = 0.05;
			double zdelt = 0.00025;

			if (simplex == null) // no initial_simplex
			{
				simplex = new double[N + 1][];

				for (int i = 0; i < N + 1; i++)
				{
					simplex[i] = new double[N];
				}

				for (int i = 0; i < N; i++)
				{
					simplex[0][i] = x0[i];
				}

				for (int k = 0; k < N; k++)
				{
					for (int i = 0; i < N; i++)
					{
						if (i != k)
						{
							simplex[k + 1][i] = x0[i];
						}
						else
						{
							if (x0[k] != 0)
								simplex[k + 1][i] = (1 + nonzdelt) * x0[k];
							else
							{
								simplex[k + 1][i] = zdelt;
							}
						}
					}
				}
			}

			// Evaluation
			double[] functionValues = new double[N + 1];
			int[] indices = new int[N + 1];
			for (int vertex_of_simplex = 0; vertex_of_simplex <= N; vertex_of_simplex++)
			{
				functionValues[vertex_of_simplex] = obj(simplex[vertex_of_simplex]);
				indices[vertex_of_simplex] = vertex_of_simplex;
			}

			// Infinite loop until convergence
			int iterations = 1;

			while (iterations < maxiter)
			{
				iterations++;

				// Order
				Array.Sort(functionValues, indices);
				
				Console.WriteLine("{0}", String.Join(",", functionValues.Select(p => p.ToString()).ToArray()));

				// Check for convergence
				if (System.Math.Abs(functionValues[0] - functionValues[N]) < m_fatol)					
				{
					break;
				}

				// Find centroid of the simplex excluding the vertex with highest functionvalue.
				double[] centroid = new double[N];

				for (int index = 0; index < N; index++)
				{
					centroid[index] = 0;
					for (int vertex_of_simplex = 0; vertex_of_simplex <= N; vertex_of_simplex++)
					{
						if (vertex_of_simplex != indices[N])
						{
							centroid[index] += simplex[vertex_of_simplex][index] / N;
						}
					}
				}

				//Reflection
				double[] reflection_point = new double[N];
				for (int index = 0; index < N; index++)
				{
					reflection_point[index] = centroid[index] + alpha * (centroid[index] - simplex[indices[N]][index]);
				}

				double reflection_value = obj(reflection_point);

				if (reflection_value >= functionValues[0] & reflection_value < functionValues[N - 1])
				{
					simplex[indices[N]] = reflection_point;
					functionValues[N] = reflection_value;
					continue;
				}


				// Expansion
				if (reflection_value < functionValues[0])
				{
					double[] expansion_point = new double[N];
					for (int index = 0; index < N; index++)
					{
						expansion_point[index] = centroid[index] + gamma * (reflection_point[index] - centroid[index]);
					}
					double expansion_value = obj(expansion_point);

					if (expansion_value < reflection_value)
					{
						simplex[indices[N]] = expansion_point;
						functionValues[N] = expansion_value;
					}
					else
					{
						simplex[indices[N]] = reflection_point;
						functionValues[N] = reflection_value;
					}
					continue;
				}

				// Contraction
				double[] contraction_point = new double[N];
				for (int index = 0; index < N; index++)
				{
					contraction_point[index] = centroid[index] + rho * (simplex[indices[N]][index] - centroid[index]);
				}

				double contraction_value = obj(contraction_point);

				if (contraction_value < functionValues[N])
				{
					simplex[indices[N]] = contraction_point;
					functionValues[N] = contraction_value;
					continue;
				}

				//Shrink
				double[] best_point = simplex[indices[0]];
				for (int vertex_of_simplex = 0; vertex_of_simplex <= N; vertex_of_simplex++)
				{
					for (int index = 0; index < N; index++)
					{

						simplex[indices[vertex_of_simplex]][index] = best_point[index] + sigma * (simplex[vertex_of_simplex][index] - best_point[index]);						
					}

					functionValues[vertex_of_simplex] = obj(simplex[indices[vertex_of_simplex]]);
				}

			}

			return simplex[0];

		}
	}
}
