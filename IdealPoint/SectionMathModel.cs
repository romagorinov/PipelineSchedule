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

        List<RegimeMathModel> _regimes;
        List<RepairMathModel> _repairs;
        int _period;
        int _pumpsCount;
        double[] _pumpsSigns;

        List<ConvexHull> _convexHulls;
        List<List<RegimeMathModel>> _avalRegimes;
        Dictionary<RepairMathModel, List<int>> _sameRepairIntervals;

        #endregion

        #region Свойства

        public Tuple<double, double[]> NormCoefficients
        {
            get;
            set;
        }

        #endregion

        #region Конструкторы

        public SectionMathModel(List<Tuple<double, double[][]>> regimes, double[] pumpSigns, List<Tuple<double, double[], double>> maxFlows)
        {
            // Отыщем нормирующие коэффициенты
            Helpers.Pair<double, double[]> maxPair = new Helpers.Pair<double, double[]>(double.NegativeInfinity, pumpSigns.Select(x => double.NegativeInfinity).ToArray());
            for (int i = 0; i < regimes.Count(); i++)
            {
                if (regimes[i].Item1 > maxPair.Item1)
                    maxPair.Item1 = regimes[i].Item1;

                maxPair.Item2 = maxPair.Item2.Zip(regimes[i].Item2, (x, y) => x > y[1] ? x : y[1]).ToArray();
            }
            NormCoefficients = new Tuple<double, double[]>(maxPair.Item1, maxPair.Item2);

            _period = maxFlows.Count();
            _pumpsCount = pumpSigns.Count();
            _pumpsSigns = pumpSigns.ToList().ToArray();
            _regimes = regimes
                .Select(regime => new RegimeMathModel(
                    regime.Item1,
                    regime.Item2.Select(x => x[1]).ToArray(),
                    regime.Item2.Select(x => x[0]).ToArray(),
                    pumpSigns))
                .ToList();

            _convexHulls = new List<ConvexHull>();
            _avalRegimes = new List<List<RegimeMathModel>>();
            Dictionary<RepairMathModel, Tuple<ConvexHull, List<RegimeMathModel>>> trace = new Dictionary<RepairMathModel, Tuple<ConvexHull, List<RegimeMathModel>>>();
            _sameRepairIntervals = new Dictionary<RepairMathModel, List<int>>();
            _repairs = new List<RepairMathModel>();
            for (int i = 0; i < _period; i++)
            {
                var rep = new RepairMathModel(maxFlows[i].Item1, maxFlows[i].Item2, maxFlows[i].Item3);
                int idx = _repairs.FindIndex(r => r == rep);
                if (idx != -1)
                    _repairs.Add(_repairs[idx]);
                else
                    _repairs.Add(rep);

                if (!trace.ContainsKey(_repairs[i]))
                {
                    var avalRegimes = _regimes.Where(regime => regime.CanUse(_repairs[i])).ToList();
                    var convexHull = GetConvex(avalRegimes, _repairs[i]);
                    trace.Add(_repairs[i], new Tuple<ConvexHull, List<RegimeMathModel>>(convexHull, avalRegimes));
                    _convexHulls.Add(convexHull);
                    _avalRegimes.Add(avalRegimes);
                    _sameRepairIntervals.Add(_repairs[i], new List<int>() { i });
                }
                else
                {
                    _convexHulls.Add(trace[_repairs[i]].Item1);
                    _avalRegimes.Add(trace[_repairs[i]].Item2);
                    _sameRepairIntervals[_repairs[i]].Add(i);
                }
            }
        }

        #endregion

        #region Методы

        /// <summary>
        /// Получает выпуклую оболочку из точек режимов
        /// </summary>
        /// <param name="regimes">Режимы</param>
        /// <param name="repair">Ремонтная работа на интервале</param>
        /// <returns></returns>
        private static ConvexHull GetConvex(List<RegimeMathModel> regimes, RepairMathModel repair)
        {
            var regimeConvex = regimes
                .Select(regime => regime.GetSystemOfInequalities(repair))
                .Select(tuple => ConvexHull.CreateFromHPolytope(tuple.Item1, tuple.Item2))
                .Select((convex,i) => convex.ConvexHullPoints.Select(point => Convert(regimes[i].Gin, point)).ToList())
                .SelectMany(x => x).ToList();
            
            return new ConvexHull(regimeConvex);
        }
        
        public List<double[]> GetContinuousSchedule(double inVolume, double[] pumpsVolume, int start, int end, double[] weights = null)
        {
            int period = end - start + 1;
            var repairs = _repairs.GetRange(start, period).ToList();
            var convexHulls = _convexHulls.GetRange(start, period).ToList();
            var sameRepairsIntervals = _sameRepairIntervals.ToDictionary(kv => kv.Key, kv => kv.Value.Where(x => x >= start && x <= end).Select(x => x - start).ToList()).Where(x => x.Value.Count() > 0).ToDictionary(kv => kv.Key, kv => kv.Value);

            if (weights == null)
                weights = AlgorithmHelper.CreateListOfElements(_pumpsCount + 1, 1.0).ToArray();

            // Получим средний расход на весь период
            var schedule = new List<double[]>(period);
            double avgInVol = inVolume / period;
            double[] avgPumpsVol = pumpsVolume.Select(x => x / period).ToArray();
            double[] avgPoint = Convert(avgInVol, avgPumpsVol);
            for (int i = 0; i < period; i++)
            {
                schedule.Add(avgPoint.ToList().ToArray());
            }

            Dictionary<RepairMathModel, bool> isIntervalCorrected = sameRepairsIntervals.ToDictionary(kv => kv.Key, kv => false);
            int uncorrectedPeriod = period;
            while (true)
            {
                double[] sumCorrectedVolume = new double[_pumpsCount + 1];
                foreach(var kv in sameRepairsIntervals)
                {
                    int intervalLength = kv.Value.Count();
                    int intervalIdx = kv.Value[0];
                    var currentPoint = schedule[intervalIdx];
                    if (!isIntervalCorrected[kv.Key] && !convexHulls[intervalIdx].IsPointInConvexHull(currentPoint))
                    {
                        var nearestPoint = convexHulls[intervalIdx].FindNearestInteriorPoint(currentPoint, true, weights);
                        var correctedVolume = currentPoint.Subtract(nearestPoint);
                        kv.Value.ForEach(x => schedule[x] = nearestPoint);
                        sumCorrectedVolume = sumCorrectedVolume.Add(correctedVolume.Multiply(intervalLength));
                        uncorrectedPeriod -= intervalLength;
                        isIntervalCorrected[kv.Key] = true;
                    }
                }

                if (sumCorrectedVolume.All(x => x == 0.0))
                    break;

                double[] avgCorrectedVolume = sumCorrectedVolume.Select(x => x / uncorrectedPeriod).ToArray();

                foreach (var kv in sameRepairsIntervals)
                {
                    if (!isIntervalCorrected[kv.Key])
                        kv.Value.ForEach(x => schedule[x] = schedule[x].Add(avgCorrectedVolume));
                }
            }

            return schedule;
        }

        public List<Tuple<double, double[]>> GetRegimesDecomposition(double inVolume, double[] pumpsVolume, int start, int end, double[] weights = null)
        {
            if (weights == null)
                weights = AlgorithmHelper.CreateListOfElements(_pumpsCount + 1, 1.0).ToArray();

            var normWeights = GetNormWeights(weights);

            var normQuadraticWeights = weights;//GetNormQuadraticWeights(weights);

            double[] pumpsQuadraticWeights = weights;// normQuadraticWeights.ToList().GetRange(1, _pumpsCount).ToArray();

            int period = end - start + 1;
            var repairs = _repairs.GetRange(start, period).ToList();
            var avalRegimesOnIntervals = _avalRegimes.GetRange(start, period).ToList();
            var convexHulls = _convexHulls.GetRange(start, period).ToList();

            var avaliableContinuousSchedule = GetContinuousSchedule(inVolume, pumpsVolume, start, end, weights);
            var avaliableSum = Convert(AlgorithmHelper.GetSumOnInterval(avaliableContinuousSchedule, 0, avaliableContinuousSchedule.Count()));
            inVolume = avaliableSum.Item1;
            pumpsVolume = avaliableSum.Item2;            

            List<Tuple<double, double[]>> schedule = new List<Tuple<double, double[]>>();

            double leftoverInVolume = inVolume;
            double[] leftoverPumpsVolume = pumpsVolume;
            for (int i = 0; i < period; i++)
            {
                var repair = repairs[i];

                var avalRegimes = avalRegimesOnIntervals[i].Where(regime => !leftoverPumpsVolume.Zip(regime.GpumpMin, (x, y) => x == 0 && y > 0).Any(x => x) && !(leftoverInVolume == 0 && regime.Gin > 0)).ToList();
                if (avalRegimes.Count() == 0)
                    throw new Exception();

                var avalNonZeroRegimes = avalRegimes.Where(regime => !leftoverPumpsVolume.Zip(regime.GpumpMax, (x, y) => x > 0 && y == 0).Any(x => x) && !(leftoverInVolume > 0 && regime.Gin == 0)).ToList();
                if (avalNonZeroRegimes.Count() > 0)
                    avalRegimes = avalNonZeroRegimes;

                List<double[]> continuousSchedule = GetContinuousSchedule(leftoverInVolume, leftoverPumpsVolume, start + i, end, normQuadraticWeights);
                double currentInVolume = continuousSchedule[0][0];
                double[] currentPumpsVolume = continuousSchedule[0].ToList().GetRange(1, _pumpsCount).ToArray();

                Tuple<double, double[]> selectedRegime;

                // Выбираем режим, ближайший по дистанции
                var regimeBestStates = avalRegimes.Select(regime =>
                {
                    double[] pumpVol;
                    if (regime.CanPump(currentPumpsVolume, repair))
                        pumpVol = currentPumpsVolume;
                    else
                        pumpVol = regime.GetNearestPumpsPoit(currentPumpsVolume, repair, pumpsQuadraticWeights);
                    
                    return Convert(regime.Gin, pumpVol);
                }).ToList();
                var distances = regimeBestStates.Select(x => x.Multiply(normWeights)).ToList();
                int selectedIdx = AlgorithmHelper.NearestByDistance(
                    regimeBestStates.Select(x => x.Multiply(normWeights)).ToList(),
                    Convert(currentInVolume, currentPumpsVolume).Multiply(normWeights));
                selectedRegime = Convert(regimeBestStates[selectedIdx]);


                schedule.Add(selectedRegime);

                leftoverInVolume -= selectedRegime.Item1;
                if (leftoverInVolume < 0.0)
                    leftoverInVolume = 0.0;

                leftoverPumpsVolume = leftoverPumpsVolume.Subtract(selectedRegime.Item2);
                leftoverPumpsVolume = leftoverPumpsVolume.Select(x => x < 0.0 ? 0.0 : x).ToArray();
            }

            return schedule;
        }

        public double[] GetNormWeights(double[] weights)
        {
            double[] result = weights.Select(x => x).ToArray();
            result[0] /= NormCoefficients.Item1;
            for (int i = 1; i < _pumpsCount + 1; i++)
            {
                result[i] /= NormCoefficients.Item2[i - 1];
            }
            return result;
        }

        public double[] GetNormQuadraticWeights(double[] weights)
        {
            double[] result = weights.Select(x => Math.Sqrt(x)).ToArray();
            result[0] /= NormCoefficients.Item1;
            for (int i = 1; i < _pumpsCount + 1; i++)
            {
                result[i] /= NormCoefficients.Item2[i - 1];
            }
            return result;
        }

        public List<double[]> AddOutputElement(List<double[]> schedule)
        {
            var result = new List<double[]>();
            for (int i = 0; i < schedule.Count(); i++)
            {
                var list = schedule[i].ToList();
                list.Add(list.GetRange(1, list.Count() - 1).Zip(_pumpsSigns, (x, y) => x * y).Sum() + list[0]);
                result.Add(list.ToArray());
            }
            return result;
        }

        public static double[] Convert(double val, double[] vals)
        {
            double[] result = new double[vals.Length + 1];
            result[0] = val;
            Array.Copy(vals, 0, result, 1, vals.Length);
            return result;
        }

        public static double[] Convert(Tuple<double, double[]> val)
        {
            return Convert(val.Item1, val.Item2);
        }

        public static Tuple<double, double[]> Convert(double[] val)
        {
            return new Tuple<double, double[]>(val[0], val.ToList().GetRange(1, val.Count() - 1).ToArray());
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
