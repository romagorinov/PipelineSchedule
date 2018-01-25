using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using Accord.Math;

namespace Algorithms
{
    public class SectionMathModel
    {
        #region Поля

        List<Tuple<double, double[][]>> _regimes;
        double[] _pumpSigns;
        List<Tuple<double, double[], double>> _maxFlows;
        int _period;
        int _pumpCount;

        List<List<Tuple<double, double[][]>>> _avalRegimesOnIntervals;

        #endregion

        #region Свойства

        #endregion

        #region Конструкторы

        public SectionMathModel(List<Tuple<double, double[][]>> regimes, double[] pumpSigns, List<Tuple<double, double[], double>> maxFlows)
        {
            _regimes = regimes;
            _pumpSigns = pumpSigns;
            _maxFlows = maxFlows;
            _period = maxFlows.Count();
            _pumpCount = pumpSigns.Count();

            _avalRegimesOnIntervals = new List<List<Tuple<double, double[][]>>>();
            for (int i = 0; i < _period; i++)
            {
                var avalRegimesOnInterval = new List<Tuple<double, double[][]>>();
                var maxFlow = _maxFlows[i];

                foreach(var regime in _regimes)
                {
                    // Проверяем вход ТУ
                    if (regime.Item1 > maxFlow.Item1)
                        continue;

                    // Проверяем все подкачки
                    if (regime.Item2.Zip(maxFlow.Item2, (r, m) => r[0] > m).Any(x => x))
                        continue;

                    // Првоеряем выход ТУ
                    if (GetOut(regime.Item1, regime.Item2.Select(x => x[0]).ToArray()) > maxFlow.Item3)
                        continue;

                    /*var pumpVolumes = regime.Item2.Zip(pumpSigns, (r, s) => r[0] * s).Aggregate(new List<double> { regime.Item1 }, (total,current) =>
                    {
                        total.Add(total.Last() + current);
                        return total;
                    });
                    if (pumpVolumes.Any(vol => vol > maxFlow.Item1))
                        continue;*/

                    avalRegimesOnInterval.Add(regime);    
                }

                _avalRegimesOnIntervals.Add(avalRegimesOnInterval);
            }
        }

        #endregion

        #region Методы

        private double GetOut(double input, double[] pumps)
        {
            return pumps.Count() == 0 ? input : pumps.Zip(_pumpSigns, (r, s) => r * s).Sum() + input;
        }

        public List<Tuple<double, double[]>> DecomposeVolume(Tuple<double, double[]> volume, int start, int end)
        {
            var maxFlows = _maxFlows.GetRange(start, end - start).ToList();
            var avalRegimesOnIntervals = _avalRegimesOnIntervals.GetRange(start, end - start).ToList();
            List<Tuple<double, double[]>> volumes

        }

        #endregion
    }

    //public class SectionMathModel
    //{
    //    #region Поля

    //    public const double EPS = 1E-10;

    //    List<double[]> _regimes;

    //    double[] _normCoefficients;

    //    List<double[]> _normRegimes;

    //    List<double[]> _maxFlows;

    //    List<double[]> _normMaxFlows;

    //    #endregion

    //    #region Свойства

    //    public List<double[]> Regimes => _regimes;

    //    public List<double[]> NormRegimes => _normRegimes;

    //    public double[] NormCoefficients => _normCoefficients;

    //    public int Dimension => _regimes[0].Count();

    //    #endregion

    //    #region Конструкторы

    //    public SectionMathModel(List<double[]> regimes, List<double[]> maxFlows)
    //    {
    //        _regimes = regimes;
    //        _normCoefficients = AlgorithmHelper.GetMaxInListByComponents(regimes);
    //        _normRegimes = AlgorithmHelper.NormalizeByComponents(regimes);
    //        _maxFlows = maxFlows;
    //        _normMaxFlows = _maxFlows.Select(maxFlow => Normalize(maxFlow)).ToList();
    //    }

    //    #endregion

    //    #region Методы

    //    private double[] RegimeToVolume(double[] regime, int period)
    //    {
    //        return regime.Select(x => x * period).ToArray();
    //    }

    //    private double[] VolumeToRegime(double[] volume, int period)
    //    {
    //        return volume.Select(x => x / period).ToArray();
    //    }

    //    private double[] Normalize(double[] val)
    //    {
    //        return val.Zip(_normCoefficients, (x, y) => x / y).ToArray();
    //    }

    //    private double[] Denormalize(double[] val)
    //    {
    //        return val.Zip(_normCoefficients, (x, y) => x * y).ToArray();
    //    }

    //    private double GetRegimesSimilarity(double[] regime1, double[] regime2)
    //    {
    //        return AlgorithmHelper.GetDistance(regime1, regime2);
    //    }

    //    private double GetCombinationPreference(double[][] regimes, double[] targetRegime)
    //    {
    //        return regimes.Select(regime => GetRegimesSimilarity(regime, targetRegime)).Sum();
    //    }

    //    private double[] GetBestDecomposition(double[] volume, int period, double[][] regimes, double[] weights)
    //    {
    //        int regimesCount = regimes.Count();

    //        // Ограничения
    //        double[,] A = new double[regimesCount + Dimension, regimesCount + 1];
    //        //double[] b = new double[regimesCount + Dimension];

    //        //b[0] = period;
    //        for (int i = 0; i < regimesCount; i++)
    //            A[0, i] = 1.0;
    //        A[0, regimesCount] = period;

    //        for (int i = 0; i < regimesCount; i++)
    //        {
    //            //b[i + 1] = 0.0;
    //            A[i + 1, i] = 1.0;
    //            A[i + 1, regimesCount] = 0.0;
    //        }

    //        for (int i = 1; i < Dimension; i++)
    //        {
    //            //b[i + regimesCount] = - volume[i];
    //            for (int j = 0; j < regimesCount; j++)
    //                A[i + regimesCount, j] = - regimes[j][i];
    //            A[i + regimesCount, regimesCount] = -volume[i];
    //        }

    //        // Функция оптимизации
    //        double[,] Q = new double[regimesCount, regimesCount];
    //        double[] d = new double[regimesCount];

    //        for (int i = 0; i < Dimension; i++)
    //        {
    //            for (int j = 0; j < regimesCount; j++)
    //            {
    //                d[j] += -2 * weights[i] * volume[i] * regimes[j][i];
    //                for (int k = 0; k < regimesCount; k++)
    //                {
    //                    Q[j, k] += 2 * weights[i] * regimes[j][i] * regimes[k][i];
    //                }
    //            }
    //        }

    //        AlgorithmHelper.GetStringToClipboard(AlgorithmHelper.ToJuggedArray(Q));

    //        alglib.minqpstate state;
    //        alglib.minqpreport rep;
    //        alglib.minqpcreate(regimesCount, out state);
    //        alglib.minqpsetquadraticterm(state, Q);
    //        alglib.minqpsetlinearterm(state, d);
    //        int[] ct = AlgorithmHelper.CreateListOfElements(regimesCount + Dimension, 1).ToArray();
    //        ct[0] = 0;
    //        alglib.minqpsetlc(state, A, ct);
    //        double[] x = new double[regimesCount];
    //        alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
    //        alglib.minqpoptimize(state);
    //        alglib.minqpresults(state, out x, out rep);

    //        if (rep.terminationtype <= 4 && rep.terminationtype >= 1)
    //            return x;
    //        else
    //            return null;
    //    }

    //    public Tuple<double[], int>[] GetOptimalDecomposition(double[] volume, double[] minVolume, double[] maxVolume, int period, double[] weights)
    //    {
    //        var vols = Regimes.Select(regime => regime.Select(x => x * period).ToArray()).ToList();

    //        var normVolume = Normalize(volume);
    //        var normRegime = VolumeToRegime(normVolume, period);
    //        var normMinVolume = Normalize(minVolume);
    //        var normMaxVolume = Normalize(maxVolume);

    //        var lowerBound = GetQuadLowerBound(period, true);
    //        bool[] needZeroRegime = normVolume.Zip(lowerBound, (x, y) => x < y).ToArray();

    //        var regimesNumbers = _normRegimes.Select((x, i) => i).ToArray();
    //        Tuple<double[][], double[], double[], double> bestSolution = null;
    //        var combinations = Combinatorics.Combinations(regimesNumbers, Dimension + 1).ToList();
    //        ConcurrentBag<Tuple<double[][], double[], double[], double>> combinationParamsBag = new ConcurrentBag<Tuple<double[][], double[], double[], double>>();
    //        foreach (var combination in combinations)
    //        {
    //            var regimesCombination = combination.Select(x => _normRegimes[x]).ToArray();
    //            if (regimesCombination.Any(regime => regime.Zip(needZeroRegime, (x,y) => x == 0 && !y).Any(x => x))) continue;

    //            var optimalTimeCombination = GetBestDecomposition(normVolume, period, regimesCombination, weights);

    //            if (optimalTimeCombination == null) continue;

    //            // Перекачанные объемы
    //            double[] vol = new double[Dimension];
    //            for (int i = 0; i < Dimension + 1; i++)
    //            {
    //                for (int j = 0; j < Dimension; j++)
    //                {
    //                    vol[j] += regimesCombination[i][j] * optimalTimeCombination[i];
    //                }
    //            }

    //            combinationParamsBag.Add(new Tuple<double[][], double[], double[], double>(
    //                regimesCombination, optimalTimeCombination, vol, GetCombinationPreference(regimesCombination, normRegime)));
    //        }
    //        if (combinationParamsBag.Count() == 0)
    //            return null;

    //        // Сначала найдем удовлетворяющие нас решения
    //        var feasibles = combinationParamsBag.Where(x =>
    //        {
    //            for (int i = 0; i < Dimension; i++)
    //            {
    //                if (x.Item3[i] < normMinVolume[i] || x.Item3[i] > normMaxVolume[i])
    //                    return false;
    //            }

    //            return true;
    //        }).ToList();

    //        if (feasibles.Count() > 0)
    //        {
    //            feasibles.Sort((x1, x2) =>
    //            {
    //                if (x1.Item4 > x2.Item4) return 1;
    //                else if (x1.Item4 < x2.Item4) return -1;
    //                else return 0;
    //            });
    //            bestSolution = feasibles[0];
    //        }
    //        else
    //        {
    //            // Найдем лучшее решение
    //            var list = combinationParamsBag.ToList();
    //            list.Sort((x1, x2) =>
    //            {
    //                double val1 = x1.Item3.Zip(normVolume, (y1, y2) => y1 - y2).Select(y => y * y).Zip(weights, (y1, y2) => y1 * y2).Sum(),
    //                    val2 = x2.Item3.Zip(normVolume, (y1, y2) => y1 - y2).Select(y => y * y).Zip(weights, (y1, y2) => y1 * y2).Sum();
    //                if (val1 > val2) return 1;
    //                else if (val1 < val2) return -1;
    //                else return 0;
    //            });
    //            bestSolution = list[0];
    //        }

    //        if (bestSolution == null)
    //            return null;

    //        // Теперь нужо округлить
    //        var solution = bestSolution.Item2.Select(x => Math.Round(x)).ToArray();
    //        double periodModulo = period - solution.Sum();
    //        List<Tuple<double[], int>> sol = new List<Tuple<double[], int>>();
    //        for (int i = 0; i < bestSolution.Item1.Count(); i++)
    //        {
    //            if (solution[i] > Math.Abs(periodModulo))
    //            {
    //                solution[i] += periodModulo;
    //            }

    //            if (solution[i] != 0)
    //                sol.Add(new Tuple<double[], int>(bestSolution.Item1[i], (int)solution[i]));
    //        }

    //        return sol.Select(x => new Tuple<double[], int>(Denormalize(x.Item1), x.Item2)).ToArray();
    //    }

    //    private double[] GetQuadLowerBound(int period, bool norm = false)
    //    {
    //        double[] minWithoutZeros = AlgorithmHelper.CreateListOfElements(Dimension, double.PositiveInfinity).ToArray();
    //        for (int i = 0; i < _regimes.Count(); i++)
    //        {
    //            minWithoutZeros = minWithoutZeros.Zip(_regimes[i], (x, y) => y != 0 && y < x ? y : x).ToArray();
    //        }
    //        if (norm)
    //            return Normalize(minWithoutZeros.Select(x => x * period).ToArray());
    //        else
    //            return minWithoutZeros.Select(x => x * period).ToArray();
    //    }

    //    private double[] GetQuadUpperBound(int period, bool norm = false)
    //    {
    //        return AlgorithmHelper.GetMaxInListByComponents(norm ? _normRegimes : _regimes);
    //    }

    //    public bool CanPump(double[] volume)
    //    {
    //        return false;
    //    }

    //    #endregion

    //    #region Не нужно

    //    /// <summary>
    //    /// Получает все разбиения числа n на partsCount частей a[1] + a[2] + ... + a[m] = n;
    //    /// </summary>
    //    /// <param name="n"> Число, которое нужно разбиеть </param>
    //    /// <param name="m"> Количество чисел, на которые нужно разбиеть </param>
    //    /// <returns></returns>
    //    //private List<int[]> GetPartitions(int n, int m)
    //    //{
    //    //    List<int[]> result = new List<int[]>();

    //    //    // Инициализация
    //    //    int[] a = new int[m];
    //    //    a[0] = n - m + 1;
    //    //    for (int  i = 1; i < m; i++)
    //    //    {
    //    //        a[i] = 1;
    //    //    }
    //    //    int j, s, x;

    //    //    while(true)
    //    //    {
    //    //        while (true)
    //    //        {
    //    //            // Посещаем разбиение, получаем все перестановки
    //    //            result.AddRange(Combinatorics.Permutations(a));

    //    //            if (a[1] < a[0] - 1)
    //    //            {
    //    //                a[0] = a[0] - 1;
    //    //                a[1] = a[1] + 1;
    //    //            }
    //    //            else break;
    //    //        }

    //    //        if (m == 2)
    //    //            return result;

    //    //        j = 2;
    //    //        s = a[0] + a[1] - 1;

    //    //        while (a[j] >= a[0] - 1)
    //    //        {
    //    //            s += a[j];
    //    //            j++;
    //    //            if (j == m)
    //    //                return result;
    //    //        }

    //    //        x = a[j] + 1;
    //    //        a[j] = x;
    //    //        j--;

    //    //        for (int i = j; i > 0; i--)
    //    //        {
    //    //            a[i] = x;
    //    //            s -= x;
    //    //        }

    //    //        a[0] = s;
    //    //    }
    //    //}

    //    #endregion
    //}
}
