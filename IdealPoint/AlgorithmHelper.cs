#define NOCHECKER

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace Algorithms
{
    /// <summary>
    /// Вспомогательные функции, например выбор ближайшего по модулю или рандом рулеткой
    /// Random потокобезопасен
    /// </summary>
    public class AlgorithmHelper
    {
        //private static Random _rnd = new Random();

        private static int seed = Environment.TickCount;
        private static ThreadLocal<Random> _threadLocal = new ThreadLocal<Random>(() => new Random(Interlocked.Increment(ref seed)));
        public static Random RandomInstance
        {
            get
            {
                return _threadLocal.Value;
            }
        }

        public static double NearestByModul(List<double> arr, double val)
        {
            #if (CHECKER)
            if (arr == null || arr.Count() == 0) throw new ArgumentException();
            #endif

            int idx = 0;
            double min = Math.Abs(val - arr[0]);
            double currentMin;
            for (int i = 1; i < arr.Count(); i++)
            {
                currentMin = Math.Abs(val - arr[i]);
                if (currentMin < min)
                {
                    min = currentMin;
                    idx = i;
                }
            }

            return arr[idx];
        }

        public static double[] NearestByDistance(List<double[]> arr, double[] val)
        {
#if (CHECKER)
            if (arr == null || arr.Count() == 0 || val == null) throw new ArgumentNullException();
            if (arr.Any(x => x == null || x.Count() == 0)) throw new ArgumentNullException();
            if (arr.Any(x => x.Count() != val.Count())) throw new ArgumentException();
#endif

            double min = GetDistance(arr[0], val);
            int idx = 0;
            double currentMin;
            for (int i = 1; i < arr.Count(); i++)
            {
                currentMin = GetDistance(arr[i], val);
                if (currentMin < min)
                {
                    min = currentMin;
                    idx = i;
                }
            }

            return arr[idx];
        }
        
        public static double? NearestWithUpperBound(List<double> arr, double val, double upperBound)
        {
#if (CHECKER)
            if (arr == null || arr.Count() == 0) throw new ArgumentException();
#endif

            arr = arr.Where(x => x <= upperBound).ToList();

            if (arr.Count() == 0)
                return null;

            return NearestByModul(arr, val);
        }

        public static double? NearestWithLowerBound(List<double> arr, double val, double lowerBound)
        {
#if (CHECKER)
            if (arr == null || arr.Count() == 0) throw new ArgumentException();
#endif

            arr = arr.Where(x => x >= lowerBound).ToList();

            if (arr.Count() == 0)
                return null;

            return NearestByModul(arr, val);
        }

        public static IEnumerable<IEnumerable<T>> CartesianProduct<T>
            (IEnumerable<IEnumerable<T>> sequences)
        {
            IEnumerable<IEnumerable<T>> emptyProduct = new[] { Enumerable.Empty<T>() };
            return sequences.Aggregate(emptyProduct, (accumulator, sequence) =>
               from accseq in accumulator
               from item in sequence
               select accseq.Concat(new[] { item }));
        }

        public static int GetRouletIndex(List<double> probabilities)
        {
#if (CHECKER)
            if (probabilities == null || probabilities.Count() == 0) throw new ArgumentException();
#endif

            double p = RandomInstance.NextDouble();
            double lower = 0, upper = 0;
            for (int i = 0; i < probabilities.Count(); i++)
            {
                upper += probabilities[i];

                if (p >= lower && p < upper)
                    return i;

                lower = upper;
            }

            return -1;
        }

        public static int GetRouletIndexCloserBetter(List<double> arr, double value, double irrelevant = 0.1)
        {
#if (CHECKER)
            if (arr == null || arr.Count() == 0) throw new ArgumentException();
#endif

            List<double> probability = arr.Select(x => 1 / (irrelevant + Math.Abs(value - x))).ToList();
            double sum = probability.Sum();
            probability = probability.Select(x => x / sum).ToList();

            return GetRouletIndex(probability);
        }

        public static double GetDistance(double[] v1, double[] v2)
        {
#if (CHECKER)
            if (v1 == null || v2 == null || v1.Length == 0 || v2.Length == 0 || v1.Length != v2.Length) throw new ArgumentException();
#endif

            return Math.Sqrt(v1.Zip(v2, (x, y) => { double k = x - y; return k * k; }).Sum());
        }

        public static double GetLength(double[] v1)
        {
#if (CHECKER)
            if (v1 == null || v1.Length == 0) throw new ArgumentException();
#endif

            return Math.Sqrt(v1.Sum(x => x * x));
        }

        public static bool ListsEquality(List<double[]> list1, List<double[]> list2)
        {
            if (list1 == null || list2 == null) throw new ArgumentNullException();

            return list1.Count() == list2.Count() &&
                !list1.Zip(list2, (x, y) => x.SequenceEqual(y)).Any(x => x == false);
        }

        public static bool ListsEqualityNoOrder<T>(List<T> list1, List<T> list2)
        {
#if (CHECKER)
            if (list1 == null || list2 == null) throw new ArgumentNullException();
#endif

            return list1.Count() == list2.Count() && list1.All(x => list2.Contains(x));
        }

        public static bool ListHasValue(List<double[]> list, double[] value)
        {
            return list.FindIndex(x => x.SequenceEqual(value)) != -1;
        }

        public static bool IsArrayLower(double[] arrLower, double[] arrUpper)
        {
            for (int i = 0; i < arrLower.Count(); i++)
            {
                if (arrLower[i] > arrUpper[i])
                    return false;
            }

            return true;
        }

        public static double[] ListToVector(List<double[]> list)
        {
            return list.SelectMany(x => x).ToArray();
        }

        public static List<T[]> RemoveDuplicates<T>(List<T[]> list)
        {
#if (CHECKER)
            if (list == null) throw new ArgumentNullException();
            if (list.Count() == 0) return new List<double[]>();
#endif

            List<T[]> result = new List<T[]>();
            foreach(var el in list)
            {
                if (result.FindIndex(x => x.SequenceEqual(el)) == -1)
                    result.Add(el);
            }

            return result;
        }

        public static double[] GetSumOnInterval(List<double[]> list, int start, int end)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0 || list.Any(x => x == null)) throw new ArgumentNullException();
            if (list.Any(x => x.Length != list[0].Length)) throw new ArgumentException();
            if (start >= end) throw new ArgumentException();
#endif

            double[] result = list[0].Select(x => 0.0).ToArray();
            for (int i = start; i < end; i++)
            {
                result = result.Zip(list[i], (x, y) => x + y).ToArray();
            }

            return result;
        }

        public static double GetSumOnInterval(List<double> list, int start, int end)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0) throw new ArgumentNullException();
            if (start >= end) throw new ArgumentException();
#endif

            double result = 0;
            for (int i = start; i < end; i++)
            {
                result += list[i];
            }

            return result;
        }

        public static List<T[]> CreateListOfArrays<T>(int listLength, int arrayLength, T defaultValue)
        {
#if (CHECKER)
            if (listLength < 1 || arrayLength < 1) throw new ArgumentException();
#endif
            T[][] result = new T[listLength][];
            for (int i = 0; i < listLength; i++)
            {
                result[i] = new T[arrayLength];
                for (int j = 0; j < arrayLength; j++)
                    result[i][j] = defaultValue;
            }
            return result.ToList();
        }

        public static List<T> CreateListOfElements<T>(int listLength, T defaultValue)
        {
#if (CHECKER)
            if (listLength < 1) throw new ArgumentException();
#endif

            List<T> result = new List<T>();
            for (int i = 0; i < listLength; i++)
            {
                result.Add(defaultValue);
            }
            return result;
        }

        public static List<double[]> GetDeviations(List<double[]> list1, List<double[]> list2)
        {
#if (CHECKER)
            if (list1 == null || list1.Count() == 0 || list1.Any(x => x == null) || list2 == null || list2.Count() == 0 || list2.Any(x => x == null)) throw new ArgumentNullException();
            if (list1.Any(x => x.Length != list1[0].Length) || list2.Any(x => x.Length != list2[0].Length)) throw new ArgumentException();
            if (list1.Count() != list2.Count() || list1[0].Count() != list2[0].Count()) throw new ArgumentException();
#endif

            return list1.Zip(list2, (el1, el2) => el1.Zip(el2, (x,y) => x - y).ToArray()).ToList();
        }

        public static double[] GetSumOfVecrots(double[] v1, double[] v2)
        {
#if (CHECKER)
            if (v1 == null || v2 == null || v1.Length == 0 || v2.Length == 0 || v1.Length != v2.Length) throw new ArgumentException();
#endif

            return v1.Zip(v2, (x, y) => x + y).ToArray();
        }

        public static double[] GetDifOfVectors(double[] v1, double[] v2)
        {
#if (CHECKER)
            if (v1 == null || v2 == null || v1.Length == 0 || v2.Length == 0 || v1.Length != v2.Length) throw new ArgumentException();
#endif

            return v1.Zip(v2, (x, y) => x - y).ToArray();
        }

        public static List<double>[] DivideIntoLists(List<double[]> list)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0 || list.Any(x => x == null)) throw new ArgumentNullException();
            if (list.Any(x => x.Length != list[0].Length)) throw new ArgumentException();
#endif

            int dim = list[0].Count(),
                period = list.Count();
            List<double>[] result = list[0].Select(x => new List<double>()).ToArray();

            for (int i = 0; i < period; i++)
            {
                for (int j = 0; j < dim; j ++)
                {
                    result[j].Add(list[i][j]);
                }
            }

            return result;
        }

        public static List<double[]> CombineIntoList(List<double>[] lists)
        {
#if (CHECKER)
            if (lists == null || lists.Any(x => x == null || x.Count() == 0)) throw new ArgumentNullException();
            if (lists.Any(x => x.Count() != lists[0].Count())) throw new ArgumentException();
#endif

            int dim = lists.Count();
            int period = lists[0].Count();
            List<double[]> result = new List<double[]>();

            for (int i = 0; i < period; i++)
            {
                result.Add(lists.Select(x => x[i]).ToArray());
            }

            return result;
        }

        public static List<double> NormalizeList(List<double> list)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0) throw new ArgumentNullException();
#endif

            double min = list.Min(), max = list.Max();

            if (min != max)
            {
                return list.Select(x => (x - min) / (max - min)).ToList();
            }
            else
            {
                return list.Select(x => 1.0).ToList();
            }
        }

        public static List<double[]> NormalizeByComponents(List<double[]> list)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0 || list.Any(x => x == null)) throw new ArgumentNullException();
            if (list.Any(x => x.Length != list[0].Length)) throw new ArgumentException();
#endif
            List<double>[] temp = DivideIntoLists(list);

            for (int i = 0; i < temp.Count(); i++)
            {
                temp[i] = NormalizeList(temp[i]);
            }

            return CombineIntoList(temp);
        }

        public static double[] GetMaxInListByComponents(List<double[]> list)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0 || list.Any(x => x == null)) throw new ArgumentNullException();
            if (list.Any(x => x.Length != list[0].Length)) throw new ArgumentException();
#endif

            double[] max = list[0].Select(x => double.NegativeInfinity).ToArray();

            list.ForEach((el) =>
            {
                max = max.Zip(el, (x, y) => x > y ? x : y).ToArray();
            });

            return max;
        }

        public static double[] GetMinInListByComponents(List<double[]> list)
        {
#if (CHECKER)
            if (list == null || list.Count() == 0 || list.Any(x => x == null)) throw new ArgumentNullException();
            if (list.Any(x => x.Length != list[0].Length)) throw new ArgumentException();
#endif

            double[] min = list[0].Select(x => double.PositiveInfinity).ToArray();

            list.ForEach((el) =>
            {
                min = min.Zip(el, (x, y) => x > y ? y : x).ToArray();
            });

            return min;
        }

        public static List<double> GetIntegralDifference(double startValue, List<double> positiveList, List<double> negativeList, List<double> neutralList)
        {
            List<double> result = new List<double>() { startValue };
            int count = positiveList.Count();

            for (int i = 0;  i < count; i++)
            {
                double prevVal = result.Last();
                result.Add(prevVal + positiveList[i] - negativeList[i] + neutralList[i]);
            }

            return result;
        }

        public static List<int> GetGridOnInterval(int intervalsStart, int intervalEnd, int gridSize)
        {
            List<int> result = new List<int>();
            int intervalLength = intervalEnd - intervalsStart + 1;
            if (intervalLength <= gridSize)
            {
                for(int i = intervalsStart; i <= intervalEnd; i++)
                {
                    result.Add(i);
                }
            }
            else
            {
                double delta = (double)intervalLength / (gridSize - 1);
                for (int i = 0; i < gridSize - 1; i++ )
                {
                    result.Add((int)Math.Round(i * delta) + intervalsStart);
                }
                result.Add(intervalEnd);
            }

            return result;
        }

        public static T[] GetArrayByMask<T>(T[] source, bool[] mask)
        {
            if (source.Count() != mask.Count()) throw new ArgumentException();

            return source.Where((x, i) => mask[i]).ToArray();
        }

        public static List<T[]> MaskListOfArrays<T>(List<T[]> list, bool[] mask)
        {
            return list.Select(x => GetArrayByMask(x, mask)).ToList();
        }

        public static string GetStringToClipboard<T>(IEnumerable<IEnumerable<T>> list)
        {
            var toStrings = list.Select(arr => String.Join(" ", arr)).ToList();            
            return String.Join(System.Environment.NewLine, toStrings);
        }


        public static T[][] ToJuggedArray<T>(T[,] matrix)
        {
            int rows = matrix.GetLength(0), cols = matrix.GetLength(1);
            T[][] result = new T[rows][];
            for (int i = 0; i < rows; i++)
            {
                result[i] = new T[cols];
                for (int j = 0; j < cols; j++)
                    result[i][j] = matrix[i, j];
            }
            return result;
        }
    }
}
