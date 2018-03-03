using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.Math.Optimization;

namespace Algorithms
{
    public class SectionWithPumpsMathModel : ISection
    {
        #region Поля

        #region Константы

        public const double EPS = 1E-8;
        public const int DIGITS = 8;

        #endregion

        #region Не расчетные

        List<RegimeMathModel> _regimes;
        List<RepairMathModel> _repairs;
        int _period;
        int _pumpsCount;
        double[] _pumpsSigns;
        double[] _weights;

        #endregion

        #region Расчетные

        public IntervalsParameters _currentIntervalsParameters;

        public Error err = new Error() { errorNumber = 1, msg = "Расчетов не производилось" };

        private RepairMathModel _notRepair;

        List<int> _reapirDayIntervals;

        List<Tuple<HashSet<RegimeMathModel>, ConvexHull>> _combinationsConvex;

        #endregion

        #endregion

        #region Свойства

        public double[] Weights
        {
            get => _weights;
            set
            {
                _weights = value ?? AlgorithmHelper.CreateListOfElements(_pumpsCount + 1, 1.0).ToArray();
            }
        }
        
        public int Dimension => _pumpsCount + 1;

        public List<List<int>> ControlAvaliableIntervals
        {
            get
            {
                return _currentIntervalsParameters.uniformityIntervals;
            }
        }

        #endregion

        #region Конструкторы

        public SectionWithPumpsMathModel(List<Tuple<double, double[][]>> regimes, double[] pumpSigns, List<Tuple<double, double[], double>> maxFlows, double[] weights = null)
        {
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
            _repairs = new List<RepairMathModel>();
            for (int i = 0; i < _period; i++)
            {
                var rep = new RepairMathModel(maxFlows[i].Item1, maxFlows[i].Item2, maxFlows[i].Item3);
                int idx = _repairs.FindIndex(r => r == rep);

                if (idx != -1)
                {
                    _repairs.Add(_repairs[idx]);
                }
                else
                {
                    _repairs.Add(rep);
                }
            }

            List<RepairMathModel> distinctRepairs = _repairs.Distinct().ToList();
            List<int> repCount = distinctRepairs.Select(repair => _repairs.Count(x => x == repair)).ToList();
            _notRepair = distinctRepairs[repCount.IndexOf(repCount.Max())];

            // Отыщем нормирующие коэффициенты
            Helpers.Pair<double, double[]> maxPair = new Helpers.Pair<double, double[]>(double.NegativeInfinity, pumpSigns.Select(x => double.NegativeInfinity).ToArray());
            for (int i = 0; i < regimes.Count(); i++)
            {
                if (regimes[i].Item1 > maxPair.Item1)
                    maxPair.Item1 = regimes[i].Item1;

                maxPair.Item2 = maxPair.Item2.Zip(regimes[i].Item2, (x, y) => x > y[1] ? x : y[1]).ToArray();
            }
            
            Weights = weights;

            HashSet<int> repairDays = new HashSet<int>();
            for (int i = 0; i < _period; i++)
            {
                int day = i / 24;
                if (_repairs[i] != _notRepair)
                    repairDays.Add(day);
            }

            _reapirDayIntervals = new List<int>();
            foreach (var day in repairDays.ToList())
                for (int hour = 0; hour < 24; hour++)
                    _reapirDayIntervals.Add(day * 24 + hour);

            _combinationsConvex = new List<Tuple<HashSet<RegimeMathModel>, ConvexHull>>();
            foreach (var regimesCombination in Combinatorics.Combinations(_regimes.ToArray(), Dimension + 1))
            {
                var list = regimesCombination.ToList();
                _combinationsConvex.Add(new Tuple<HashSet<RegimeMathModel>, ConvexHull>(new HashSet<RegimeMathModel>(list), GetConvex(list, _notRepair)));
            }
        }

        #endregion

        #region Вспомогтельные структуры и классы

        public struct Error
        {
            /**
             * 0 - Расчет успешно произведен
             * 1 - Расчетов не производилось
             * 2 - Невозможно составить непрерывное расписание
             * 3 - Невозможно составить дискретное расписание
             * 4 - Ошибка в параметрах ТУ
             **/
            public int errorNumber;
            public string msg;
        }
        
        public class IntervalsParameters
        {
            public Dictionary<int, HashSet<RegimeMathModel>> avaliableRegimesOnIntervals;
            public Dictionary<Tuple<HashSet<RegimeMathModel>, RepairMathModel>, List<int>> sameParamsIntervals;
            public Dictionary<Tuple<HashSet<RegimeMathModel>, RepairMathModel>, ConvexHull> convexHulls;
            public Dictionary<int, ConvexHull> convexHullOnIntervals;
            public List<List<int>> uniformityIntervals;

            public IntervalsParameters() { }
        }

        #endregion

        #region Методы

        #region Служебные функции

        private Dictionary<Tuple<HashSet<RegimeMathModel>, RepairMathModel>, List<int>> GetSameParamsIntervals(int[] indexes)
        {
            return _currentIntervalsParameters.sameParamsIntervals
                .ToDictionary(x => x.Key, x => x.Value.Where(y => indexes.Contains(y)).ToList())
                .Where(kv => kv.Value.Count() > 0)
                .ToDictionary(kv => kv.Key, kv => kv.Value);
        }

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
        
        /// <summary>
        /// Получить схожесть режимов
        /// </summary>
        /// <param name="regime1"></param>
        /// <param name="regime2"></param>
        /// <returns></returns>
        private double GetRegimesSimilarity(RegimeMathModel regime1, RegimeMathModel regime2)
        {
            return Math.Abs(regime1.Gin - regime2.Gin);
        }

        /// <summary>
        /// Получить предпочтительность комбинации
        /// </summary>
        /// <param name="regimes"></param>
        /// <param name="targetRegime"></param>
        /// <returns></returns>
        private double GetCombinationPreference(RegimeMathModel[] regimes, RegimeMathModel targetRegime)
        {
            return regimes.Select(regime => GetRegimesSimilarity(regime, targetRegime)).Sum();
        }

        #endregion

        #region Функции для расчетов 
        
        public List<Tuple<List<double[]>, List<int>>> DecomposeVolumes(List<Tuple<double[], int[]>> volumes)
        {
            var temp = volumes.SelectMany(x => x.Item2);
            if (temp.Count() != temp.Distinct().Count())
                throw new Exception();

            List<Tuple<List<double[]>, List<int>>> result = new List<Tuple<List<double[]>, List<int>>>();
            foreach(var tuple in volumes)
            {
                var volume = Convert(tuple.Item1);
                var indexes = tuple.Item2;

                var continuousSchedule = GetContinuousSchedule(volume.Item1, volume.Item2, indexes);
                if (continuousSchedule == null)
                {
                    err = new Error() { errorNumber = 2, msg = "На заданных режимах невозможно составить ПГДН для заданных объемов на стадии равномерного распределения объемов по часам." };
                    return null;
                }

                List<double[]> discreteSchedule = GetDiscreteSchedule(volume.Item1, volume.Item2, indexes);
                if (discreteSchedule == null)
                {
                    err = new Error() { errorNumber = 3, msg = "Не удается разложить объемы перекачки на дискретные режимы." };
                    return null;
                }

                var res = new List<Tuple<List<double[]>, List<int>>>() { new Tuple<List<double[]>, List<int>>(new List<double[]>(), new List<int>()) };

                for (int i = 0; i < indexes.Count(); i++)
                {
                    int idx = indexes[i];
                    if (_repairs[idx] != _notRepair)
                    {
                        res.Add(new Tuple<List<double[]>, List<int>>(new List<double[]>() { discreteSchedule[i] }, new List<int>() { idx }));
                    }
                    else
                    {
                        res[0].Item1.Add(discreteSchedule[i]);
                        res[0].Item2.Add(idx);
                    }
                }

                result.AddRange(res);
            }

            err = new Error() { errorNumber = 0, msg = "Расчет произведен без ошибок." };
            return result;
        }
        
        private List<double[]> GetContinuousSchedule(double inVolume, double[] pumpsVolume, int[] indexes)
        {
            int period = indexes.Count();
            var sameParamsIntervals = GetSameParamsIntervals(indexes);
            double[] volumePoint = Convert(inVolume, pumpsVolume);
            var normQuadraticWeights = volumePoint.Zip(_weights, (x, y) => Math.Sqrt(y) / x).Select(x => double.IsInfinity(x) ? 1 : x).ToArray();

            var str = AlgorithmHelper.GetStringToClipboard(sameParamsIntervals.ElementAt(0).Key.Item1.SelectMany(regime =>
            {
                return new List<double[]>() {
                    new double[] { regime.Gin, regime.GpumpMax[0] },
                    new double[] { regime.Gin, regime.GpumpMin[0] }
                };
            }));
            
            int dim = _pumpsCount + 1;
            int eqNumber = 0;
            int varsCount = dim * sameParamsIntervals.Count();
            List<double[]> constraintsA = new List<double[]>();
            List<double> constraintsB = new List<double>();
            
            for (int i = 0; i < dim; i++)
            {
                double[] left = new double[varsCount];
                for (int j = 0; j < sameParamsIntervals.Count(); j++)
                    left[j * dim + i] = sameParamsIntervals.ElementAt(j).Value.Count();
                double right = volumePoint[i];

                constraintsA.Add(left);
                constraintsB.Add(right);

                eqNumber++;
            }

            for (int i = 0; i < sameParamsIntervals.Count(); i++)
            {
                var kv = sameParamsIntervals.ElementAt(i);
                double[][] convA = _currentIntervalsParameters.convexHulls[kv.Key].Amatrix.ToJagged();
                double[] convB = _currentIntervalsParameters.convexHulls[kv.Key].Bvector;
                int constrNum = convB.Count();
                
                for (int j = 0; j < constrNum; j ++)
                {
                    double[] left = new double[varsCount];
                    Array.Copy(convA[j], 0, left, dim * i, dim);
                    left = left.Multiply(-1.0);
                    double right = - convB[j];

                    if (j >= _currentIntervalsParameters.convexHulls[kv.Key].EqNumber)
                    {
                        constraintsA.Add(left);
                        constraintsB.Add(right);
                    }
                    else
                    {
                        constraintsA.Insert(0, left);
                        constraintsB.Insert(0, right);
                        eqNumber++;
                    }
                }
            }

            double[,] A = constraintsA.ToArray().ToMatrix();
            double[] b = constraintsB.ToArray();

            double[,] Q = new double[varsCount, varsCount];
            double[] d = new double[varsCount];
            double[] avgPoint = volumePoint.Select(x => x / period).Select(x => x == 0 ? EPS : x).ToArray();
            for (int i = 0; i < varsCount; i++)
            {
                int curDim = i % dim;
                Q[i, i] =  2 * 1.0 * normQuadraticWeights[curDim];
                d[i] = -2 * avgPoint[curDim] * normQuadraticWeights[curDim];
            }

            str = AlgorithmHelper.GetStringToClipboard(A.ToJagged());

            var solution = AlgorithmHelper.SolveQP(Q, d, A, b, eqNumber);
            if (solution == null)
                return null;

            var solutionList = solution.ToList();
            var result = AlgorithmHelper.CreateListOfElements<double[]>(period, null);
            for (int i = 0; i < sameParamsIntervals.Count(); i++)
            {
                var val = solutionList.GetRange(i * dim, dim).ToArray();
                foreach (var idx in sameParamsIntervals.ElementAt(i).Value)
                    result[indexes.IndexOf(idx)] = val;
            }

            return result;
        }
        
        public List<double[]> GetDiscreteSchedule(double inVolume, double[] pumpsVolume, int[] indexes)
        {
            int period = indexes.Count();
            var sameParamsIntervals = GetSameParamsIntervals(indexes).ToList();
            sameParamsIntervals.Sort((el1, el2) =>
            {
                int regimeCount1 = el1.Key.Item1.Count(), regimeCount2 = el2.Key.Item1.Count();
                if (regimeCount1 > regimeCount2)
                    return 1;
                else if (regimeCount1 < regimeCount2)
                    return -1;
                else
                    return 0;
            });
            List<int> indexesList = indexes.ToList();
            
            List<double[]> schedule = new double[period][].ToList();
            pumpsVolume = pumpsVolume.Select(x => x).ToArray();
            // На ремонты ставим один и тот же режим
            for (int i = 0; i < sameParamsIntervals.Count(); i++)
            {
                var repair = sameParamsIntervals[i].Key.Item2;
                var continuousSchedule = GetContinuousSchedule(inVolume, pumpsVolume, indexesList.ToArray());
                if (continuousSchedule == null)
                    return null;
                var currentIndexes = sameParamsIntervals[i].Value;
                int currentPeriod = currentIndexes.Count();
                double[] sumOnCurrentInterval = new double[_pumpsCount + 1];
                for (int j = 0; j < currentIndexes.Count(); j++)
                {
                    sumOnCurrentInterval = sumOnCurrentInterval.Add(continuousSchedule[indexesList.IndexOf(currentIndexes[j])]);
                }
                List<Tuple<double[], int>> res = null;

                if (i < sameParamsIntervals.Count() - 1)
                {
                    double inAvgVolume = sumOnCurrentInterval[0] / currentPeriod;
                    double[] pumpsAvgVolume = Convert(sumOnCurrentInterval).Item2.Select(x => x / currentPeriod).ToArray();
                    var avalRegimes = sameParamsIntervals[i].Key.Item1.ToList();
                    // Смотрим, какими можно соблюдать подкачки
                    var avalPumpRegimes = avalRegimes.Where(x => x.CanPump(pumpsAvgVolume, repair)).ToList();

                    RegimeMathModel bestRegime = null;
                    if (avalPumpRegimes.Count() != 0)
                    {
                        // Если можно соблюсти перекачки, то просто отбираем лучший по главной линии режим
                        bestRegime = avalPumpRegimes[AlgorithmHelper.NearestByModul(avalPumpRegimes.Select(x => x.Gin).ToList(), inAvgVolume)];
                        res = new List<Tuple<double[], int>>() { new Tuple<double[], int>(Convert(bestRegime.Gin, pumpsAvgVolume), currentPeriod) };
                    }
                    else
                    {
                        // Отбираем лучший по главной линии и вычисоляем точку подкачек
                        bestRegime = avalRegimes[AlgorithmHelper.NearestByModul(avalRegimes.Select(x => x.Gin).ToList(), inAvgVolume)];
                        res = new List<Tuple<double[], int>>() { new Tuple<double[], int>(Convert(bestRegime.Gin, bestRegime.GetNearestPumpsPoit(pumpsAvgVolume, repair)), currentPeriod) };
                    }
                }
                else
                {
                    var allAvalRegimes = sameParamsIntervals[i].Key.Item1;
                    double[] targetPoint = Convert(inVolume, pumpsVolume).Select(x => x / currentPeriod).ToArray();
                    var regimesCombinations = _combinationsConvex.Where(combination =>
                    {
                        var combinationHash = combination.Item1;
                        if (allAvalRegimes.IsSupersetOf(combinationHash))
                        {
                            if (combination.Item2.IsPointInConvexHull(targetPoint))
                                return true;
                            else
                                return false;
                        }
                        else
                            return false;
                    }).ToList();
                    RegimeMathModel targetRegime = new RegimeMathModel(sumOnCurrentInterval[0] / currentPeriod, new double[_pumpsCount], new double[_pumpsCount], new double[_pumpsCount]);
                    regimesCombinations.Sort((el1, el2) =>
                    {
                        double val1 = GetCombinationPreference(el1.Item1.ToArray(), targetRegime), val2 = GetCombinationPreference(el2.Item1.ToArray(), targetRegime);
                        if (val1 > val2)
                            return 1;
                        else if (val1 < val2)
                            return -1;
                        else
                            return 0;
                    });
                    foreach(var combination in regimesCombinations)
                    {
                        res = Decompose(sumOnCurrentInterval.Select(x => Round(x)).ToArray(), combination.Item1.ToArray(), currentIndexes.Count(), sameParamsIntervals[i].Key.Item2);
                        if (res != null)
                            break;
                    }
                    if (res == null)
                        return null;
                }
                double[] temp = res.Aggregate(new double[_pumpsCount + 1], (total, current) => total.Add(current.Item1.Multiply(current.Item2)));
                inVolume -= temp[0];
                if (inVolume < EPS)
                    inVolume = 0.0;
                pumpsVolume = pumpsVolume.Subtract(Convert(temp).Item2).Select(x => x < EPS ? 0.0 : x).ToArray();
                indexesList = indexesList.Except(currentIndexes).ToList();

                int curIdx = 0;
                foreach(var t in res)
                {
                    for (int j = 0; j < t.Item2; j++)
                    {
                        schedule[indexes.IndexOf(currentIndexes[curIdx])] = t.Item1;
                        curIdx++;
                    }
                }
            }

            return schedule;
        }

        private List<Tuple<double[], int>> Decompose(double[] volumes, RegimeMathModel[] regimes, int period, RepairMathModel repair)
        {
            Tuple<double, double[]> volumesPair = Convert(volumes);
            double[] avgPumps = volumesPair.Item2.Select(x => x / period).ToArray();
            int rCount = regimes.Count();
            int varCount = rCount * (_pumpsCount + 1);

            /*
            // Нижние и верхние границы
            double[] uBound = new double[varCount], lBound = new double[varCount];
            for (int i = 0; i < rCount; i++)
            {
                uBound[i] = period;
                lBound[i] = 0;
            }
            for (int i = rCount; i < varCount; i++)
            {
                int pumpNumber = i / rCount - 1;
                int regimeNumber = i % rCount;
                uBound[i] = regimes[regimeNumber].GpumpMax[pumpNumber];
                lBound[i] = regimes[regimeNumber].GpumpMin[pumpNumber];
            }

            List<double[]> Cmatrix = new List<double[]>();
            List<int> Ct = new List<int>();
            // Линейные равенства
            Cmatrix.Add(new double[varCount + 1]);
            Cmatrix.Add(new double[varCount + 1]);
            for (int i = 0; i < rCount; i++)
            {
                Cmatrix[0][i] = regimes[i].Gin;
                Cmatrix[1][i] = 1.0;
            }
            Cmatrix[0][varCount] = volumes[0];
            Cmatrix[1][varCount] = period;
            Ct.Add(0);
            Ct.Add(0);

            // Линейные неравенства (ограничения подкачек)
            for (int i = 0; i < rCount; i++)
            {
                var pumpsSystem = regimes[i].GetSystemOfInequalities(repair);
                double[] tempRow = new double[varCount + 1];
                for (int j = 0; j < _pumpsCount; j++)
                {
                    int varNumber = rCount + j * rCount + i;
                    tempRow[varNumber] = -pumpsSystem.Item1[pumpsSystem.Item2.Count() - 1][j];
                }
                tempRow[varCount] = -pumpsSystem.Item2[pumpsSystem.Item2.Count() - 1];
                Cmatrix.Add(tempRow);
                Ct.Add(1);
            }

            // Нелинейные равенства (Якобиан)
            Action<double[], double[], double[,], object> Jac = (double[] x, double[] fi, double[,] jac, object obj) =>
            {
                // Целевая функция
                fi[0] = 1.0;
                for (int i = 0; i < varCount; i++)
                {
                    jac[0, i] = 0.0;
                }

                // Равенство подкачек
                for (int i = 0; i < _pumpsCount; i++)
                {
                    fi[i + 1] = 0.0;
                    for (int j = 0; j < rCount; j++)
                    {
                        fi[i + 1] += x[j] * x[rCount + rCount * i + j];
                        jac[i + 1, j] += x[rCount + rCount * i + j];
                        jac[i + 1, rCount + rCount * i + j] += x[j];
                    }
                }
            };

            double[] scaleCoef = new double[varCount];
            for (int i = 0; i < rCount; i++)
                scaleCoef[i] = period;
            for (int i = 0; i < _pumpsCount; i++)
                for (int j = 0; j < rCount; j++)
                    scaleCoef[rCount + i * rCount + j] = volumes[i + 1];
            double[] x0 = new double[varCount];
            int outerits = 5;
            double rho = 1000;
            alglib.minnlcstate state;
            alglib.minnlcreport rep;
            alglib.minnlccreate(varCount, x0, out state);
            alglib.minnlcsetalgoaul(state, rho, outerits);
            alglib.minnlcsetcond(state, 0, 0, 0, 0);
            alglib.minnlcsetscale(state, scaleCoef);
            alglib.minnlcsetbc(state, lBound, uBound);
            alglib.minnlcsetlc(state, Cmatrix.ToArray().ToMatrix(), Ct.ToArray());
            alglib.minnlcsetnlc(state, _pumpsCount, 0);
            alglib.minnlcoptimize(state, new alglib.ndimensional_jac(Jac), null, null);
            double[] x1 = new double[varCount];
            alglib.minnlcresults(state, out x1, out rep);
            */

            double avgTime = period / regimes.Count();
            var f = new NonlinearObjectiveFunction(
                numberOfVariables: varCount,
                function: (x) =>
                {
                    double val = 0.0;
                    for(int i = 0; i < rCount; i++)
                    {
                        double temp = x[i] - avgTime;
                        val += temp * temp;
                    }
                    return val;
                },
                gradient: (x) =>
                {
                    double[] grad = new double[varCount];
                    for (int i = 0; i < rCount; i++)
                    {
                        grad[i] = 2 * x[i] - 2 * avgTime;
                    }
                    return grad;
                });

            var constraints = new List<IConstraint>();

            // Равенство подкачек
            for (int j = 0; j < _pumpsCount; j++)
            {
                int currentPump = j;
                constraints.Add(new NonlinearConstraint(
                    numberOfVariables: varCount,
                    function: (x) =>
                    {
                        double val = 0;
                        for (int i = 0; i < rCount; i++)
                        {
                            val += x[i] * x[rCount + rCount * currentPump + i];
                        }
                        return val;
                    },
                    shouldBe: ConstraintType.EqualTo,
                    value: volumes[j + 1],
                    gradient: (x) =>
                    {
                        double[] grad = new double[varCount];
                        for (int i = 0; i < rCount; i++)
                        {
                            grad[i] = x[rCount + rCount * currentPump + i];
                            grad[rCount + rCount * currentPump + i] = x[i];
                        }
                        return grad;
                    },
                    withinTolerance: 1E-5
                    ));
            }

            List<double[]> Amatrix = new List<double[]>();
            List<double> Bmatrix = new List<double>();

            // Равенство объему на входе
            double[] tempRow = new double[varCount];
            for (int i = 0; i < rCount; i++)
            {
                tempRow[i] = regimes[i].Gin;
            }
            Amatrix.Add(tempRow);
            Bmatrix.Add(volumes[0]);

            // Равенство времени периоду
            tempRow = new double[varCount];
            for (int i = 0; i < rCount; i++)
            {
                tempRow[i] = 1.0;
            }
            Amatrix.Add(tempRow);
            Bmatrix.Add(period);

            // Время больше 0
            for (int i = 0; i < rCount; i++)
            {
                tempRow = new double[varCount];
                tempRow[i] = 1.0;
                Amatrix.Add(tempRow);
                Bmatrix.Add(0.0);
            }

            // Ограничения режима
            for (int i = 0; i < rCount; i++)
            {
                var pumpsSystem = regimes[i].GetSystemOfInequalities(repair);
                for (int k = 0; k < pumpsSystem.Item2.Count(); k++)
                {
                    tempRow = new double[varCount];
                    for (int j = 0; j < _pumpsCount; j++)
                    {
                        int varNumber = rCount + j * rCount + i;
                        tempRow[varNumber] = -pumpsSystem.Item1[k][j];
                    }
                    Amatrix.Add(tempRow);
                    Bmatrix.Add(-pumpsSystem.Item2[k]);
                }
            }

            List<LinearConstraint> linearConstraints = new List<LinearConstraint>();
            for (int i = 0; i < Amatrix.Count(); i++)
                linearConstraints.Add(new LinearConstraint(varCount)
                {
                    CombinedAs = Amatrix[i],
                    ShouldBe = i > 1 ? ConstraintType.GreaterThanOrEqualTo : ConstraintType.EqualTo,
                    Value = Bmatrix[i],
                    Tolerance = 1E-5
                });

            constraints.AddRange(linearConstraints);

            var solver = new AugmentedLagrangian(f, constraints);
            
            // Выберем точку где-то в середине
            solver.Minimize(new double[varCount]);

            var solution = solver.Solution.ToList();

            double solutionPeriod = solution.GetRange(0, rCount).Sum();
            double[] pumps = new double[_pumpsCount];
            if (solver.Status != AugmentedLagrangianStatus.Converged)
            {
                if (!((int)Math.Round(solutionPeriod) == period))
                    return null;

                double input = solution.GetRange(0, rCount).Zip(regimes, (x, y) => x * y.Gin).Sum();
                for (int i = 0; i < _pumpsCount; i++)
                {
                    pumps[i] = solution.GetRange(rCount + i * rCount, rCount).Zip(solution.GetRange(0, rCount), (x, y) => x * y).Sum();
                }
                var solutionPoint = Convert(input, pumps);
                if (solutionPoint.Zip(volumes, (x, y) => Math.Abs(x - y)).Any(x => x > 1))
                    return null;
            }

            // Округлим
            int sumPeriod = 0;
            int max = int.MinValue;
            int maxIdx = -1;
            for(int i = 0; i < rCount; i++)
            {
                solution[i] = Math.Round(solution[i]);
                int temp = (int)solution[i];
                if (max < temp)
                {
                    max = temp;
                    maxIdx = i;
                }
                sumPeriod += temp;
            }
            int difPeriod = period - sumPeriod;
            solution[maxIdx] += difPeriod;

            List<Tuple<double[], int>> result = new List<Tuple<double[], int>>();
            for (int i = 0; i < rCount; i++)
            {
                if (solution[i] != 0.0)
                {
                    double[] temp = new double[_pumpsCount + 1];
                    temp[0] = regimes[i].Gin;
                    for (int j = 0; j < _pumpsCount; j++)
                    {
                        temp[j + 1] = solution[rCount + j * rCount + i];
                    }
                    if (!regimes[i].CanPump(Convert(temp).Item2, repair))
                        return null;
                    result.Add(new Tuple<double[], int>(temp, (int)solution[i]));
                }
            }
            result.Sort((el1, el2) =>
            {
                if (el1.Item1[0] > el2.Item1[0])
                    return -1;
                else if (el1.Item1[0] < el2.Item1[0])
                    return 1;
                else
                    return 0;
            });
            return result;
        }

        #endregion

        #region Преобразования

        public double[] AddOutputElement(double[] val)
        {
            var list = val.ToList();
            list.Add(list.GetRange(1, list.Count() - 1).Zip(_pumpsSigns, (x, y) => x * y).Sum() + list[0]);
            return list.ToArray();
        }

        public static double[] RemoveOutputElement(double[] val)
        {
            var result = val.ToList();
            result.RemoveAt(result.Count() - 1);
            return result.ToArray();
        }

        public static T[] Convert<T>(T val, T[] vals)
        {
            T[] result = new T[vals.Length + 1];
            result[0] = val;
            Array.Copy(vals, 0, result, 1, vals.Length);
            return result;
        }

        public static T[] Convert<T>(Tuple<T, T[]> val)
        {
            return Convert(val.Item1, val.Item2);
        }

        public static Tuple<T, T[]> Convert<T>(T[] val)
        {
            return new Tuple<T, T[]>(val[0], val.ToList().GetRange(1, val.Count() - 1).ToArray());
        }
        
        public double Round(double val)
        {
            return Math.Round(val, DIGITS);
        }

        public List<double[]> AddOutputComponent(List<double[]> schedule)
        {
            return schedule.Select(x => AddOutputElement(x)).ToList();
        }

        #endregion

        #region Реализация интерфейсов

        public void CalcDefaultIntervalsParameters(List<Tuple<double[], int[]>> volumes)
        {
            var temp = volumes.SelectMany(x => x.Item2);
            if (temp.Count() != temp.Distinct().Count())
                throw new Exception();

            // Предрасчеты
            Dictionary<int, HashSet<RegimeMathModel>> avaliableRegimesOnIntervals = new Dictionary<int, HashSet<RegimeMathModel>>();
            foreach (var tuple in volumes)
            {
                var volumeIn = tuple.Item1[0];
                var volumePump = Convert(tuple.Item1).Item2;
                var indexes = tuple.Item2.ToList();

                // Нулевой режим
                if (volumeIn == 0.0)
                {
                    indexes.ForEach(idx => avaliableRegimesOnIntervals.Add(idx, new HashSet<RegimeMathModel>() { _regimes.First(regime => regime.Gin == 0.0) }));
                }
                else
                {
                    foreach (var idx in indexes)
                    {
                        var repair = _repairs[idx];

                        // Допустимы по ремонтам, обязательно
                        var canUseRegimes = _regimes.Where(regime => regime.CanUse(repair));

                        if (canUseRegimes.Count() == 0)
                            throw new Exception();

                        // Не качаем, когда не нужно, на подкачки, обязательно
                        canUseRegimes = canUseRegimes.Where(regime => !volumePump.Zip(regime.GpumpMin, (x, y) => x == 0 && y != 0).Any(x => x));

                        // Убираем нулевой режим
                        var technologyRegimes = canUseRegimes.Where(regime => regime.Gin != 0);

                        if (technologyRegimes.Count() == 0)
                            technologyRegimes = canUseRegimes;

                        // Качаем на подкачки, когда нужно
                        var needUseRegimes = technologyRegimes.Where(regime => !regime.GpumpMax.Zip(volumePump, (x, y) => x == 0 && y != 0).Any(x => x));

                        if (needUseRegimes.Count() == 0)
                            needUseRegimes = technologyRegimes;

                        avaliableRegimesOnIntervals.Add(idx, new HashSet<RegimeMathModel>(needUseRegimes));
                    }
                }
            }

            Dictionary<Tuple<HashSet<RegimeMathModel>, RepairMathModel>, List<int>> sameParamsIntervals = new Dictionary<Tuple<HashSet<RegimeMathModel>, RepairMathModel>, List<int>>();
            foreach (var kv in avaliableRegimesOnIntervals)
            {
                var idx = kv.Key;
                var hashSet = kv.Value;
                var repair = _repairs[idx];
                var existentKey = sameParamsIntervals.Keys.FirstOrDefault(x => x.Item1.SetEquals(hashSet) && x.Item2 == repair);

                if (existentKey == null)
                    sameParamsIntervals.Add(new Tuple<HashSet<RegimeMathModel>, RepairMathModel>(hashSet, repair), new List<int>() { idx });
                else
                    sameParamsIntervals[existentKey].Add(idx);
            }

            Dictionary<Tuple<HashSet<RegimeMathModel>, RepairMathModel>, ConvexHull> convexHulls = sameParamsIntervals
                .Select(x => new { KEY = x.Key, VALUE = GetConvex(x.Key.Item1.ToList(), x.Key.Item2) })
                .ToDictionary(x => x.KEY, x => x.VALUE);

            Dictionary<int, ConvexHull> convexHullOnIntervals = new Dictionary<int, ConvexHull>();
            foreach (var kv in sameParamsIntervals)
                foreach (var idx in kv.Value)
                    convexHullOnIntervals.Add(idx, convexHulls[kv.Key]);

            _currentIntervalsParameters = new IntervalsParameters()
            {
                avaliableRegimesOnIntervals = avaliableRegimesOnIntervals,
                convexHullOnIntervals = convexHullOnIntervals,
                convexHulls = convexHulls,
                sameParamsIntervals = sameParamsIntervals,
                uniformityIntervals = volumes.Select(kv => kv.Item2.Where(x => _repairs[x] == _notRepair).ToList()).ToList()
            };
            _currentIntervalsParameters.uniformityIntervals.ForEach(x => x.Sort());
        }

        public List<Tuple<List<double[]>, List<int>>> GetSchedule(List<Tuple<double[], int[]>> volumes)
        {
            return DecomposeVolumes(volumes);
        }

        public List<double[]> GetFullSchedule(List<double[]> schedule)
        {
            return AddOutputComponent(schedule);
        }

        #endregion

        #endregion
    }
}
