using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using Accord.Math;
using GeneticSharp;
using GeneticSharp.Domain.Fitnesses;
using GeneticSharp.Domain.Chromosomes;
using GeneticSharp.Domain.Mutations;
using Accord.Math.Optimization;
using System.Threading;

namespace Algorithms
{
    public class SectionMathModel
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
        bool _inZeroAvaliable;
        bool[] _pumpsZeroAvaliable;
        double[] _weights;

        #endregion

        #region Расчетные

        Dictionary<RepairMathModel, ConvexHull> _convexHulls;
        Dictionary<RepairMathModel, List<RegimeMathModel>> _avalRegimes;
        Dictionary<RepairMathModel, List<int>> _sameRepairIntervals;
        Tuple<double, double[]> _normCoefficients;

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

        public double[] UpperBound
        {
            get;
            set;
        }

        public double[] LowerBound
        {
            get;
            set;
        }

        #endregion

        #region Конструкторы

        public SectionMathModel(List<Tuple<double, double[][]>> regimes, double[] pumpSigns, List<Tuple<double, double[], double>> maxFlows, double[] weights = null, bool inZeroAvaliable = true, bool[] pumpsZeroAvaliable = null)
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
            _sameRepairIntervals = new Dictionary<RepairMathModel, List<int>>();
            _repairs = new List<RepairMathModel>();
            for (int i = 0; i < _period; i++)
            {
                var rep = new RepairMathModel(maxFlows[i].Item1, maxFlows[i].Item2, maxFlows[i].Item3);
                int idx = _repairs.FindIndex(r => r == rep);

                if (idx != -1)
                {
                    _repairs.Add(_repairs[idx]);
                    _sameRepairIntervals[_repairs[idx]].Add(i);
                }
                else
                {
                    _repairs.Add(rep);
                    _sameRepairIntervals.Add(rep, new List<int>() { i });
                }
            }

            // Отыщем нормирующие коэффициенты
            Helpers.Pair<double, double[]> maxPair = new Helpers.Pair<double, double[]>(double.NegativeInfinity, pumpSigns.Select(x => double.NegativeInfinity).ToArray());
            for (int i = 0; i < regimes.Count(); i++)
            {
                if (regimes[i].Item1 > maxPair.Item1)
                    maxPair.Item1 = regimes[i].Item1;

                maxPair.Item2 = maxPair.Item2.Zip(regimes[i].Item2, (x, y) => x > y[1] ? x : y[1]).ToArray();
            }
            _normCoefficients = new Tuple<double, double[]>(maxPair.Item1, maxPair.Item2);

            SetZeroAvaliable(inZeroAvaliable, pumpsZeroAvaliable);
            Weights = weights;
            UpperBound = Convert(1.0, _pumpsSigns.Select(x => 1.01).ToArray());
            LowerBound = Convert(1.0, _pumpsSigns.Select(x => 0.99).ToArray());
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
             **/
            public int errorNumber;
            public string msg;
        }

        public Error err = new Error() { errorNumber = 1, msg = "Расчетов не производилось" };
        
        #endregion

        #region Методы

        #region Служебные функции

        /// <summary>
        ///  Пересчет доступных режимов и выпуклых оболочек для периодов при заданных допустимых нулевых режимах 
        /// </summary>
        /// <param name="inZeroAvaliable"></param>
        /// <param name="pumpsZeroAvaliable"></param>
        public void SetZeroAvaliable(bool inZeroAvaliable, bool[] pumpsZeroAvaliable = null)
        {
            _inZeroAvaliable = inZeroAvaliable;
            _pumpsZeroAvaliable = pumpsZeroAvaliable ?? _pumpsSigns.Select(x => true).ToArray();

            _convexHulls = new Dictionary<RepairMathModel, ConvexHull>();
            _avalRegimes = new Dictionary<RepairMathModel, List<RegimeMathModel>>();
            foreach(var kv in _sameRepairIntervals)
            {
                var repair = kv.Key;

                var canUseRegimes = _regimes.Where(regime => regime.CanUse(repair)).ToList();
                if (canUseRegimes.Count() == 0) throw new Exception();
                var conditionRegimes = canUseRegimes.Where(regime => regime.MeetsZeroConditions(_inZeroAvaliable, _pumpsZeroAvaliable)).ToList();

                var avaliableRegimes = conditionRegimes.Count() == 0 ? canUseRegimes : conditionRegimes;
                _avalRegimes.Add(repair, avaliableRegimes);

                _convexHulls.Add(repair, GetConvex(avaliableRegimes, repair));
            }
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
        /// Получить интервалы с одинаковыми ремонтными работами
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        private Dictionary<RepairMathModel, List<int>> GetRepairIntervals(int[] indexes)
        {
            int period = indexes.Count();
            return _sameRepairIntervals
                .ToDictionary(kv => kv.Key, kv => kv.Value.Where(x => indexes.Contains(x)).ToList())
                .Where(x => x.Value.Count() > 0)
                .ToDictionary(kv => kv.Key, kv => kv.Value);
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

        public List<double[]> DecomposeVolumes(double inVolume, double[] pumpsVolume, int[] indexes)
        {
            var continuousSchedule = GetContinuousSchedule(inVolume, pumpsVolume, indexes);
            if (continuousSchedule == null)
            {
                err = new Error() { errorNumber = 2, msg = "На заданных режимах невозможно составить ПГДН для заданных объемов на стадии равномерного распределения объемов по месяцу." };
                return null;
            }

            List<double[]> schedule = GetDiscreteSchedule(inVolume, pumpsVolume, indexes);
            if (schedule == null)
            {
                err = new Error() { errorNumber = 3, msg = "Не удается разложить объемы перекачки на дискретные режимы." };
                return null;
            }

            err = new Error() { errorNumber = 0, msg = "Расчет произведен без ошибок." };
            return schedule;
        }

        public List<double[]> GetContinuousSchedule(double inVolume, double[] pumpsVolume, int[] indexes)
        {
            int period = indexes.Count();
            var sameRepairIntervals = GetRepairIntervals(indexes);
            double[] volumePoint = Convert(inVolume, pumpsVolume);
            var normQuadraticWeights = volumePoint.Zip(_weights, (x, y) => Math.Sqrt(y) / x).ToArray();

            int dim = _pumpsCount + 1;
            int eqNumber = 0;
            int varsCount = dim * sameRepairIntervals.Count();
            List<double[]> constraintsA = new List<double[]>();
            List<double> constraintsB = new List<double>();
            
            for (int i = 0; i < dim; i++)
            {
                double[] left = new double[varsCount];
                for (int j = 0; j < sameRepairIntervals.Count(); j++)
                    left[j * dim + i] = sameRepairIntervals.ElementAt(j).Value.Count();
                double right = volumePoint[i];

                constraintsA.Add(left);
                constraintsB.Add(right);

                eqNumber++;
            }

            for (int i = 0; i < sameRepairIntervals.Count(); i++)
            {
                var kv = sameRepairIntervals.ElementAt(i);
                double[][] convA = _convexHulls[kv.Key].Amatrix.ToJagged();
                double[] convB = _convexHulls[kv.Key].Bvector;
                int constrNum = convB.Count();
                
                for (int j = 0; j < constrNum; j ++)
                {
                    double[] left = new double[varsCount];
                    Array.Copy(convA[j], 0, left, dim * i, dim);
                    left = left.Multiply(-1.0);
                    double right = - convB[j];

                    if (j >= _convexHulls[kv.Key].EqNumber)
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
            double[] avgPoint = volumePoint.Select(x => x / period).ToArray();
            for (int i = 0; i < varsCount; i++)
            {
                int curDim = i % dim;
                Q[i, i] =  2 * 1.0 * normQuadraticWeights[curDim];
                d[i] = -2 * avgPoint[curDim] * normQuadraticWeights[curDim];
            }
            

            var solution = AlgorithmHelper.SolveQP(Q, d, A, b, eqNumber);
            if (solution == null)
                return null;

            var solutionList = solution.ToList();
            var result = AlgorithmHelper.CreateListOfElements<double[]>(period, null);
            for (int i = 0; i < sameRepairIntervals.Count(); i++)
            {
                var val = solutionList.GetRange(i * dim, dim).ToArray();
                foreach (var idx in sameRepairIntervals.ElementAt(i).Value)
                    result[indexes.IndexOf(idx)] = val;
            }

            return result;
        }
        
        public List<double[]> GetDiscreteSchedule(double inVolume, double[] pumpsVolume, int[] indexes)
        {
            int period = indexes.Count();
            var repairIntervals = GetRepairIntervals(indexes).ToList();
            repairIntervals.Sort((el1, el2) =>
            {
                int regimeCount1 = _avalRegimes[el1.Key].Count(), regimeCount2 = _avalRegimes[el2.Key].Count();
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
            for (int i = 0; i < repairIntervals.Count(); i++)
            {
                var repair = repairIntervals[i].Key;
                var continuousSchedule = GetContinuousSchedule(inVolume, pumpsVolume, indexesList.ToArray()).ToList();
                if (continuousSchedule == null)
                    return null;
                var currentIndexes = repairIntervals[i].Value;
                int currentPeriod = currentIndexes.Count();
                double[] sumOnCurrentInterval = new double[_pumpsCount + 1];
                for (int j = 0; j < currentIndexes.Count(); j++)
                {
                    sumOnCurrentInterval = sumOnCurrentInterval.Add(continuousSchedule[indexesList.IndexOf(currentIndexes[j])]);
                }
                List<Tuple<double[], int>> res = null;

                if (i < repairIntervals.Count() - 1)
                {
                    double inAvgVolume = sumOnCurrentInterval[0] / currentPeriod;
                    double[] pumpsAvgVolume = Convert(sumOnCurrentInterval).Item2.Select(x => x / currentPeriod).ToArray();
                    var avalRegimes = _avalRegimes[repair].ToList();
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
                    var regimesCombinations = Combinatorics.Combinations(_avalRegimes[repairIntervals[i].Key].ToArray(), _pumpsCount + 2).ToList();
                    RegimeMathModel targetRegime = new RegimeMathModel(sumOnCurrentInterval[0] / currentPeriod, new double[_pumpsCount], new double[_pumpsCount], new double[_pumpsCount]);
                    regimesCombinations = regimesCombinations.Where(combination =>
                    {
                        bool upper = false, lower = false;
                        double t = Round(targetRegime.Gin);
                        foreach (var regime in combination)
                        {
                            double r = Round(regime.Gin);
                            if (r < t)
                                lower = true;
                            else if (r > t)
                                upper = true;
                            else
                                lower = upper = true;
                        }
                        return lower && upper;
                    }).ToList();
                    regimesCombinations.Sort((el1, el2) =>
                    {
                        double val1 = GetCombinationPreference(el1, targetRegime), val2 = GetCombinationPreference(el2, targetRegime);
                        if (val1 > val2)
                            return 1;
                        else if (val1 < val2)
                            return -1;
                        else
                            return 0;
                    });
                    int counter = 0;
                    CancellationTokenSource tokenSource = new CancellationTokenSource();
                    CancellationToken token = tokenSource.Token;
                    Parallel.ForEach(regimesCombinations, (combination, state) =>
                    {
                        var curRes = Decompose(sumOnCurrentInterval.Select(x => Round(x)).ToArray(), combination, currentIndexes.Count(), repairIntervals[i].Key, token);
                        Interlocked.Increment(ref counter);
                        if (curRes != null)
                        {                            
                            res = curRes;
                            state.Stop();
                            tokenSource.Cancel();
                        }
                    });
                    if (res == null)
                        return null;
                }
                double[] temp = res.Aggregate(new double[_pumpsCount + 1], (total, current) => total.Add(current.Item1.Multiply(current.Item2)));
                inVolume -= temp[0];
                pumpsVolume = pumpsVolume.Subtract(Convert(temp).Item2);
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

        private List<Tuple<double[], int>> Decompose(double[] volumes, RegimeMathModel[] regimes, int period, RepairMathModel repair, CancellationToken token)
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

            solver.Token = token;

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
        
        #endregion

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
