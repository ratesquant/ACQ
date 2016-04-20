using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ACQ.Math
{
    /// <summary>
    /// Collection of sorting algoritms. 
    /// </summary>
    public static class Sort
    {
        #region Swap function
        private static void Swap<T>(T[] a, int i, int j)
        {
            T temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }
        #endregion

        #region O(n^2) Sorting Algorithms
        public static void Bubble<T>(T[] a) where T : IComparable<T>
        {
            if (a == null)
                throw new ArgumentNullException();

            int n = a.Length;
            bool swapped = true;

            do
            {
                swapped = false;
                n = n - 1;
                for (int i = 0; i < n; i++)
                {
                    if (a[i].CompareTo(a[i + 1]) > 0)
                    {
                        Swap(a, i, i + 1);
                        swapped = true;
                    }
                }
            } while (swapped);
        }

        /// <summary>
        /// A bi-directional bubble sort demonstration algorithm 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="a"></param>
        public static void Cocktail<T>(T[] a) where T : IComparable<T>
        {
            if (a == null)
                throw new ArgumentNullException();

            int j;
            int limit = a.Length;
            int st = -1;
            T temp;

            while (st < limit)
            {
                bool swapped = false;
                st++;
                limit--;
                for (j = st; j < limit; j++)
                {
                    if (a[j].CompareTo(a[j + 1]) > 0)
                    {
                        temp = a[j];
                        a[j] = a[j + 1];
                        a[j + 1] = temp;
                        swapped = true;
                    }
                }
                if (!swapped)
                {
                    return;
                }
                for (j = limit; --j >= st; )
                {
                    if (a[j].CompareTo(a[j + 1]) > 0)
                    {
                        temp = a[j];
                        a[j] = a[j + 1];
                        a[j + 1] = temp;
                        swapped = true;
                    }
                }

                if (!swapped)
                {
                    return;
                }
            }
        }
        /// <summary>
        /// n^2 sort that is more efficient than other n^2 sorts such as selection sort or bubble sort
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="a"></param>
        public static void Insertion<T>(T[] a) where T : IComparable<T>
        {
            int n = a.Length;
            int j;
            T temp;

            for (int i = 1; i < n; i++)
            {
                temp = a[i];

                j = i - 1;

                while (j >= 0 && a[j].CompareTo(temp) > 0)
                {
                    a[j + 1] = a[j];
                    j--;
                }
                a[j + 1] = temp;
            }
        }

        public static void Selection<T>(T[] a) where T : IComparable<T>
        {
            int min_index;
            int n = a.Length;

            for (int i = 0; i < n - 1; i++)
            {
                min_index = i;

                for (int j = i + 1; j < n; j++)
                {
                    if (a[j].CompareTo(a[min_index]) < 0)
                    {
                        min_index = j;
                    }
                }
                Swap(a, i, min_index);
            }
        }
        #endregion

        #region Shell Sort
        public static void Shell<T>(T[] a) where T : IComparable<T>
        {
            if (a == null)
                throw new ArgumentNullException();

            int n = a.Length;
            int i, j, inc;
            T temp;

            inc = n / 2;

            while (inc > 0)
            {
                for (i = inc; i < n; i++)
                {
                    j = i;
                    temp = a[i];
                    while ((j >= inc) && (a[j - inc].CompareTo(temp) > 0))
                    {
                        a[j] = a[j - inc];
                        j = j - inc;
                    }
                    a[j] = temp;
                }


                if (inc == 2)
                    inc = 1;
                else
                    inc = inc * 5 / 11;
            }
        }
        #endregion

        #region MergeSort
        public static void Merge<T>(T[] a) where T : IComparable<T>
        {
            MergeSort(a, 0, a.Length - 1);
        }
        private static void MergeSort<T>(T[] a, int left, int right) where T : IComparable<T>
        {
            int lo = left; int hi = right;

            if (lo >= hi)
            {
                return;
            }

            int mid = (lo + hi) / 2;
            // Partition the list into two lists and sort them recursively
            MergeSort(a, lo, mid);
            MergeSort(a, mid + 1, hi);

            //Merge the two sorted lists
            int end_lo = mid;
            int start_hi = mid + 1;
            while ((lo <= end_lo) && (start_hi <= hi))
            {
                if (a[lo].CompareTo(a[start_hi]) < 0)
                {
                    lo++;
                }
                else
                {
                    //a[lo] >= a[start_hi] The next element comes from the second list,
                    //move the a[start_hi] element into the next
                    //position and shuffle all the other elements up. 
                    T temp = a[start_hi];
                    for (int k = start_hi - 1; k >= lo; k--)
                    {
                        a[k + 1] = a[k];
                    }
                    a[lo] = temp;
                    lo++;
                    end_lo++;
                    start_hi++;
                }
            }
        }
        #endregion

        #region HeapSort
        public static void Heap<T>(T[] a) where T : IComparable<T>
        {
            int n = a.Length;
            T temp;

            for (int k = n / 2; k > 0; k--)
            {
                HeapDown(a, k, n);
            }
            do
            {
                temp = a[0];
                a[0] = a[n - 1];
                a[n - 1] = temp;
                n = n - 1;
                HeapDown(a, 1, n);
            }
            while (n > 1);
        }
        private static void HeapDown<T>(T[] a, int k, int n) where T : IComparable<T>
        {
            T temp = a[k - 1];
            while (k <= n / 2)
            {
                int j = k + k;
                if ((j < n) && (a[j - 1].CompareTo(a[j]) < 0))
                {
                    j++;
                }
                if (temp.CompareTo(a[j - 1]) >= 0)
                {
                    break;
                }
                else
                {
                    a[k - 1] = a[j - 1];
                    k = j;
                }
            }
            a[k - 1] = temp;
        }
        #endregion

        #region QuickSort
        public static void Quick<T>(T[] a) where T : IComparable<T>
        {
            if (a == null)
                throw new ArgumentNullException();

            QuickSort(a, 0, a.Length - 1);
        }

        private static void QuickSort<T>(T[] a, int left, int right) where T : IComparable<T>
        {
            if (left >= right)
                return;

            int index = QuickPartitionSimple(a, left, right);
            QuickSort(a, left, index - 1);
            QuickSort(a, index + 1, right);
        }



        //Partition the array into two halves and return the
        //index about which the array is partitioned
        private static int QuickPartition<T>(T[] a, int left, int right) where T : IComparable<T>
        {
            int l = left;
            int r = right;
            T temp;

            int pivotIndex = (right + left) / 2;
            T pivot = a[pivotIndex];

            temp = a[left];
            a[left] = pivot;
            a[pivotIndex] = temp;

            for (; ; )
            {
                //find element larger than pivot
                while (a[l].CompareTo(pivot) <= 0 && l < r)
                    l++;
                //find element smaller than pivot
                while (a[r].CompareTo(pivot) > 0 && l <= r)
                    r--;

                if (l >= r) break;

                //swap
                temp = a[l];
                a[l] = a[r];
                a[r] = temp;
            }

            a[left] = a[r];
            a[r] = pivot;

            // elements a[left] to a[l] are less than or equal to pivot
            // elements a[r] to a[right] are greater than pivot

            return r;
        }
        /// <summary>
        /// Simple version of partitioning algorithm that does a lot of unnessesary swaps 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="a"></param>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        private static int QuickPartitionSimple<T>(T[] a, int left, int right) where T : IComparable<T>
        {
            int pivotIndex = (right + left) / 2;
            int i, index;

            // Tri-Median Methode! 
            //if (a[left].CompareTo(a[pivotIndex])>0) Swap(a, left, pivotIndex);
            //if (a[left].CompareTo(a[right]) > 0) Swap(a, left, right);
            //if (a[pivotIndex].CompareTo(a[right])>0) Swap(a, pivotIndex, right); 

            T pivot = a[pivotIndex];
            T temp;

            temp = a[right];
            a[right] = pivot;
            a[pivotIndex] = temp;

            for (i = left, index = left; i < right; i++)
            {
                //find elements that are less than pivot
                if (a[i].CompareTo(pivot) < 0)
                {
                    //swap
                    temp = a[i];
                    a[i] = a[index];
                    a[index] = temp;

                    index++;
                }
            }

            temp = a[index];
            a[index] = pivot;
            a[right] = temp;

            return index;
        }
        #endregion

        #region QuickerSort (uses left element as pivot)
        public static void Quicker<T>(T[] a) where T : IComparable<T>
        {
            QuickSort(a, 0, a.Length - 1);
        }

        private static void QuickerSort<T>(T[] a, int i, int j) where T : IComparable<T>
        {
            if (i < j)
            {
                int q = Partition(a, i, j);
                QuickSort(a, i, q);
                QuickSort(a, q + 1, j);
            }
        }

        private static int Partition<T>(T[] a, int p, int r) where T : IComparable<T>
        {
            T x = a[p];
            int i = p - 1;
            int j = r + 1;
            T tmp;
            while (true)
            {
                do
                {
                    j--;
                } while (a[j].CompareTo(x) > 0);
                do
                {
                    i++;
                } while (a[i].CompareTo(x) < 0);
                if (i < j)
                {
                    tmp = a[i];
                    a[i] = a[j];
                    a[j] = tmp;
                }
                else return j;
            }
        }
        #endregion
    }

    /// <summary>
    /// Partial Sorting (for doubles only), adapted from R
    /// </summary>
    public static class PSort
    {
        /// <summary>
        /// Partial sort so that x[k] is in the correct place, smaller to left, larger to right
        /// </summary>
        /// <param name="x"></param>
        /// <param name="k"></param>
       public static void Psort(double[] x, int k)
       {
            rPsort2(x, 0, x.Length - 1, k); //dont check range for k here, let c# take care of this 
       }

        private static int rcmp(double x, double y, bool nalast)
        {
            bool nax = Double.IsNaN(x);
            bool nay = Double.IsNaN(y);
            if (nax && nay) return 0;
            if (nax) return nalast ? 1 : -1;
            if (nay) return nalast ? -1 : 1;
            if (x < y) return -1;
            if (x > y) return 1;
            return 0;
        }

        private static void rPsort2(double[] x, int lo, int hi, int k)
        {
            double v, w;
            
            bool nalast=true; //put nan to the end
            int L, R, i, j;
            
            for (L = lo, R = hi; L < R; )
            {
                v = x[k];	
                for(i = L, j = R; i <= j;) 
                {
                    while (rcmp(x[i], v, nalast) < 0) i++;
                    while (rcmp(v, x[j], nalast) < 0) j--;
                    if (i <= j) 
                    {
                        w = x[i]; 
                        x[i++] = x[j]; 
                        x[j--] = w;
                    }
                }
                if (j < k) L = i;	
                if (k < i) R = j;	
            }
        } 
    }
}
