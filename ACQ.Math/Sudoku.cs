using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math
{
    public class Sudoku : ICloneable
    {        
        private const int m_size = 9;
        private const int m_blocksize = 3;
        private int[,] m_grid; //from 1 - 9  (0-empty)        

        public Sudoku()
        {
            m_grid = Generate();
        }

        public object Clone() // ICloneable implementation
        {
            Sudoku mc = this.MemberwiseClone() as Sudoku;
            mc.m_grid = m_grid.Clone() as int[,];
            return mc;
        }

        public int Size
        {
            get
            {
                return m_size;
            }
        }

        public int BlockSize
        {
            get
            {
                return m_blocksize;
            }
        }

        public static int[,] Generate()
        {
            int seed = unchecked((int)DateTime.Now.Ticks);

            int[,] grid = Generate(seed);

            return grid;
        }

        public static int[,] Generate(int seed)
        {
            int i, k;   

            int[,] grid = new int[m_size, m_size];

            System.Random rnd = new System.Random(seed);

            int[] digits = new int[m_size];

            for (i = 0; i < m_size; i++)
            {
                digits[i] = i + 1;
            }

            //Initialize diagonal blocks with random numbers 
            for (i = 0; i < m_size; i += m_blocksize)
            {
                ShuffleArray(digits, rnd);

                for (k = 0; k < m_size; k++)
                {
                    int ik = k % m_blocksize;
                    int jk = k / m_blocksize;

                    grid[i + ik, i + jk] = digits[k];
                }
            }

            List<int[,]> solutions = new List<int[,]>();

            //Solve this puzzle (we only need one solution)
            Solve(grid, solutions, 1);

            grid = solutions[0];

            GeneratePuzzle(grid, rnd);

            return grid;
        }

        /// <summary>
        /// Remove random elements from solved puzzle while solution is unique 
        /// (if we remove to many elements solution will not be unique)
        /// </summary>
        private static void GeneratePuzzle(int[,] grid, System.Random rnd)
        {            
            int[] vIndex = new int[m_size * m_size];

            for (int i = 0; i < vIndex.Length; i++)
                vIndex[i] = i;

            ShuffleArray(vIndex, rnd);

            int x, y;

            List<int[,]> solutions = new List<int[,]>();

            for (int i = 0; i < vIndex.Length; i++)
            {
                x = vIndex[i] % m_size;
                y = vIndex[i] / m_size;

                int nDigit = grid[x, y];

                grid[x, y] = 0;

                Solve(grid, solutions, 2);

                //if puzzle has more than one solution, restore removed digit
                if (solutions.Count > 1)
                    grid[x, y] = nDigit;
            }

            Solve(grid, solutions, 2);
        }

        /// <summary>
        /// Fisher-Yates shuffle, also known as the Knuth shuffle
        /// </summary>
        /// <param name="array"></param>
        public static void ShuffleArray(int[] array, System.Random rnd)
        {
            int n = array.Length;
            while (--n > 0)
            {
                int k = rnd.Next(n + 1);  // 0 <= k <= n (!)
                int temp = array[n];
                array[n] = array[k];
                array[k] = temp;
            }
        }

        public int this[int i, int j]
        {
            get
            {
                return m_grid[i, j];
            }
            set
            {
                m_grid[i, j] = value;
            }
        }

        public bool Valid
        {
            get
            {
                return IsValid(m_grid);
            }
        }

        private static bool IsValid(int[,] grid)
        {
            int i, j, k, ii, jj;
            int[] digits = new int[m_size + 1];


            //check all rows
            for (i = 0; i < m_size; i++)
            {                
                for (j = 0; j < m_size; j++)
                {
                    digits[grid[i, j]]++;
                }

                for (k = 1; k <= m_size; k++)
                {
                    if (digits[k] != 1) //there should be exactly one of each digit
                    {
                        return false;
                    }
                    digits[k] = 0; //init for the next row
                }
            }

            //check all cols
            for (i = 0; i < m_size; i++)
            {                
                for (j = 0; j < m_size; j++)
                {
                    digits[grid[j, i]]++;
                }

                for (k = 1; k <= m_size; k++)
                {
                    if (digits[k] != 1)
                    {
                        return false;
                    }
                    digits[k] = 0; //init for the next row
                }
            }

            //check all blocks
            int blocks = m_size / m_blocksize;

            for (i = 0; i < blocks; i++)
            {
                for (j = 0; j < blocks; j++)
                {
                    for (ii = 0; ii < m_blocksize; ii++)
                        for (jj = 0; jj < m_blocksize; jj++)
                            digits[grid[ii + i * m_blocksize, jj + j * m_blocksize]]++;

                    for (k = 1; k <= m_size; k++)
                    {
                        if (digits[k] != 1)
                        {
                            return false;
                        }
                        digits[k] = 0;
                    }
                }
            }
            return true;
        }

        public bool IsSolved
        {
            get
            {
                return CheckSolution(m_grid);
            }
        }

        private static bool CheckSolution(int[,] grid)
        {
            for (int i = 0; i < m_size; i++)
            {
                for (int j = 0; j < m_size; j++)
                {
                    if (grid[i, j] == 0)
                        return false;
                }
            }
            return IsValid(grid);
        }
        /// <summary>
        /// Checks conflict between original and this sudoki game
        /// </summary>
        /// <param name="original"></param>
        /// <returns></returns>
        public bool IsValid(int index_i, int index_j)
        {
            bool bValid = true;

            if (m_grid[index_i, index_j] != 0)
            {
                int i, j;

                int k = m_grid[index_i, index_j];

                for (i = 0; i < m_size; i++)
                {
                    if ((m_grid[i, index_j] == k && i != index_i) ||
                        (m_grid[index_i, i] == k && i != index_j))
                    {
                        bValid = false;
                        break;
                    }
                }

                int ii = (index_i / m_blocksize) * m_blocksize;
                int jj = (index_j / m_blocksize) * m_blocksize;

                for (i = 0; i < m_blocksize; i++)
                {
                    for (j = 0; j < m_blocksize; j++)
                    {
                        if (m_grid[ii + i, jj + j] == k && (ii + i) != index_i && (jj + j) != index_j)
                        {
                            bValid = false;
                            break;
                        }
                    }
                }
            }
            return bValid;
        }

        public int ValidDigits(int index_i, int index_j, List<int> vValid)
        {
            return ValidDigits(m_grid, index_i, index_j, vValid);
        }

        private static int ValidDigits(int[,] grid, int index_i, int index_j, List<int> valid_digits)
        {
            int i, j, ii, jj;

            valid_digits.Clear();

            if (grid[index_i, index_j] == 0)
            {
                int[] digits = new int[m_size + 1]; //digits that are already there

                //check horizontal and vertical lines 
                for (i = 0; i < m_size; i++)
                {
                    digits[grid[i, index_j]] = 1;
                    digits[grid[index_i, i]] = 1;
                }

                //check block
                ii = (index_i / m_blocksize) * m_blocksize;
                jj = (index_j / m_blocksize) * m_blocksize;

                for (i = 0; i < m_blocksize; i++)
                {
                    for (j = 0; j < m_blocksize; j++)
                    {
                        digits[grid[ii + i, jj + j]] = 1;
                    }
                }

                //add digits that are not on the invalid digits list 
                for (i = 0; i < m_size; i++)
                {
                    if (digits[i + 1] == 0)
                    {
                        valid_digits.Add(i + 1);
                    }
                }

            }
            return valid_digits.Count;
        }

        public List<int[,]> SolvePuzzle(int nMaxSolutions)
        {
            List<int[,]> solutions = new List<int[,]>();

            Solve(m_grid, solutions, nMaxSolutions);

            return solutions;
        }

        public static void Solve(int[,] grid, List<int[,]> vSolutions, int nMaxSol)
        {            
            vSolutions.Clear();

            bool allSoultionsFound = SolveRecursive(grid, vSolutions, nMaxSol);
        }
        /// <summary>
        /// This function is called recursively
        /// </summary>
        /// <param name="grid"></param>
        /// <param name="solutions"></param>
        /// <param name="nMaxSol"></param>
        /// <returns>False is not all solutions have been found</returns>
        private static bool SolveRecursive(int[,] grid, List<int[,]> solutions, int nMaxSol)
        {
            if (solutions.Count >= nMaxSol)
                return false;

            int i, j, k;

            int not_set = -1;
            int index_i = not_set;
            int index_j = not_set;

            int min_digits = Int32.MaxValue;
            List<int> valid_digits = new List<int>(m_size);
            List<int> temp = new List<int>(m_size);

            //find the empty spot with fewer possible digits
            for (i = 0; i < m_size; i++)
            {
                for (j = 0; j < m_size; j++)
                {
                    if (grid[i, j] == 0)
                    {
                        ValidDigits(grid, i, j, temp);

                        int nCount = temp.Count;
                        if (nCount < min_digits)
                        {
                            index_i = i;
                            index_j = j;
                            min_digits = nCount;

                            valid_digits.Clear();
                            foreach (int item in temp)
                            {
                                valid_digits.Add(item);
                            }

                            if (nCount == 1) 
                                break;
                        }
                    }
                }
            }

            //grid is full 
            if (index_i == not_set)
            {
                if (CheckSolution(grid))
                {
                    solutions.Add((int[,])grid.Clone());
                }                
                return true;
            }
            else
            {
                //loop through all possibilities
                for (i = 0; i < valid_digits.Count; i++)
                {
                    k = valid_digits[i];

                    int digit = grid[index_i, index_j]; //save previous digit
                    grid[index_i, index_j] = k;

                    //call this function recursively!
                    SolveRecursive(grid, solutions, nMaxSol);

                    grid[index_i, index_j] = digit;//restore previous digit

                }
                return true;
            }
        }
    }
}
