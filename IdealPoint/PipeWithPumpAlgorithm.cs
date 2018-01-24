// 5.705
#define PARALLEL

//16.5

using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;

namespace Algorithms
{
    public class PipeWithPumpAlgorithm: ReservoirBalanceAlgorithm.TechnologicalSectionAlgorithm
    {
        #region Базовые показатели технологических участков

        // Список режимов
        private List<double[]> _regimes;

        // Список максимальных расходов по всем трубам (с учетом максимальных режимов работы на данном интервале времени)
        private List<double[]> _upperBounds;
        
        // Период
        private int _period;

        // Размерность
        private int _dimension;

        // Алгоритмы для каждой из труб
        private List<SinglePipeAlgorithm> _algorithms;

        // Показывает, по каким трубам нельзя переливать
        private bool[] _constraints;

        #endregion

        #region Дополнительные расчеты

        // Доступные режимы для каждого интервала
        private List<List<double[]>> _avaliableRegimesOnIntervals;

        // Максимальное и минимальное значение по трубам (без нулевого режима) - для нормирования
        private double[] _maxFlows;
        private double[] _minFlows;

        // Максимальный и минимальный объемы, которые можно перекачать (не включая нулевой режим)
        private double[] _maxVolumes;
        private double[] _minVolumes;

        private List<int[]> _combinations;
        private List<List<double[]>> _regimesCombinations;

        public List<List<double[]>> GetAvaliableRegimesOnIntervals
        {
            get
            {
                return _avaliableRegimesOnIntervals.Select(x => x.Select(y => y.Select(z => z).ToArray()).ToList()).ToList();
            }
        }

        #endregion

        #region Нормированные показатели

        // Нормированые режимы
        private List<double[]> _normRegimes;

        // Доступные нормированные режимы для каждого интервала
        private List<List<double[]>> _avaliableNormRegimesOnIntervals;

        private List<List<double[]>> _normRegimesCombinations;

        #endregion

        #region Конструкторы

        public PipeWithPumpAlgorithm(List<double[]> regimes, List<double[]> maxFlows, bool[] constraints)
        {
            _regimes = regimes.Select(x => x.Select(y => y).ToArray()).ToList() ?? throw new ArgumentNullException();
            _upperBounds = maxFlows.Select(x => x.Select(y => y).ToArray()).ToList() ?? throw new ArgumentNullException();
            _avaliableRegimesOnIntervals = _upperBounds.Select(maxFlow => _regimes.Where(regime => maxFlow.Zip(regime, (x, y) => x - y).Min() >= 0).ToList()).ToList();
            _constraints = constraints;

            PreCalculations();
        }

        public PipeWithPumpAlgorithm(List<double[]> regimes, List<List<double[]>> avaliableRegimesOnIntervals, bool[] constraints)
        {
            _regimes = regimes.Select(x => x.Select(y => y).ToArray()).ToList() ?? throw new ArgumentNullException();
            _avaliableRegimesOnIntervals = avaliableRegimesOnIntervals.Select(r => r.Select(regime => _regimes.First(x => x.SequenceEqual(regime))).ToList()).ToList();
            _upperBounds = _avaliableRegimesOnIntervals.Select(x => x[0].Select(y => 0.0).ToArray()).ToList();
            _constraints = constraints;

            PreCalculations();
        }

        private void PreCalculations()
        {
            _period = _upperBounds.Count();
            _dimension = _regimes[0].Count();

            _maxFlows = _regimes.Aggregate(_regimes[0].Select(x => 0.0).ToArray(), (total, current) => current.Zip(total, (x, y) => x > y ? x : y).ToArray());

            _minFlows = _regimes.Aggregate(_regimes[0].Select(x => Double.MaxValue).ToArray(), (total, current) => current.Zip(total, (x, y) => (x < y && x != 0) ? x : y).ToArray());

            _normRegimes = _regimes.Select(x => x.Zip(_maxFlows, (y, z) => y / z).ToArray()).ToList();

            _avaliableNormRegimesOnIntervals = _avaliableRegimesOnIntervals.Select(x => x.Select(y => _normRegimes[_regimes.IndexOf(y)]).ToList()).ToList();

            _maxVolumes = _regimes[0].Select(x => 0.0).ToArray();
            _minVolumes = _regimes[0].Select(x => 0.0).ToArray();
            for (int i = 0; i < _period; i++)
            {
                var r = _avaliableRegimesOnIntervals[i];
                double[] max = _maxVolumes.Select(x => 0.0).ToArray(),
                    min = _maxVolumes.Select(x => Double.MaxValue).ToArray();

                foreach (var regime in r)
                {
                    max = max.Zip(regime, (x, y) => x > y ? x : y).ToArray();
                    min = min.Zip(regime, (x, y) => (y == 0 || x < y) ? x : y).ToArray();
                }
                min = min.Select(x => x == Double.MaxValue ? 0 : x).ToArray();

                _maxVolumes = _maxVolumes.Zip(max, (x, y) => x + y).ToArray();
                _minVolumes = _minVolumes.Zip(min, (x, y) => x + y).ToArray();
                _upperBounds[i] = max;
            }

            _algorithms = new List<SinglePipeAlgorithm>();
            for (int i = 0; i < _dimension; i++)
            {
                List<List<double>> a = _avaliableRegimesOnIntervals.Select(x => x.Select(y => y[i]).Distinct().ToList()).ToList();
                List<double> r = _regimes.Select(x => x[i]).Distinct().ToList();
                _algorithms.Add(new SinglePipeAlgorithm(r, a));
            }

            //_combinations = Combinatorics.Combinations(_regimes.Select((x, i) => i).ToArray(), _dimension + 1).ToList();
            //_regimesCombinations = _combinations.Select(x => x.Select(y => _regimes[y]).ToList()).ToList();
            //_normRegimesCombinations = _combinations.Select(x => x.Select(y => _normRegimes[y]).ToList()).ToList();

        }
            
        #endregion

        #region Методы

        public List<double[]> FormSchedule(double[] targets, double[] weights, List<double[]> fixValues = null)
        {
            if (targets == null || weights == null)
            {
                throw new ArgumentNullException();
            }

            // Понижаем объемы по трубе до максимально возможных
            targets = targets.Zip(_maxVolumes, (x, y) => x > y ? y : x).ToArray();

            if (fixValues == null)
            {
                fixValues = AlgorithmHelper.CreateListOfArrays(_period, _dimension, -1.0);
            }

            List<double[]> notSortedSchedule = fixValues.Select(x => x.Select(y => y).ToArray()).ToList();
            for (int i = 0; i < _dimension; i++)
            {
                targets[i] -= fixValues.Sum(x => x[i] < 0 ? 0 : x[i]);
            }
            targets = targets.Select(x => x < 0 ? 0 : x).ToArray();

            //////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Получаем режимы, которыми доступно перекачать данный объем
            //////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Не запланированные объемы перекачки
            double[] currentTarget = targets.Select(x => x).ToArray();
            // Не запланированное время
            int currentPeriod = _period - notSortedSchedule.Count(x => x[0] > 0);
            // Начинаем заполнять интервалы
            for (int idx = 0; idx < _period; idx++)
            {
                if (notSortedSchedule[idx][0] >= 0)
                {
                    continue;
                }

                // Выберем доступные режимы
                List<double[]> strongConstrains = new List<double[]>(),
                    softConstrains = new List<double[]>();
                bool strong,
                    soft;
                double[] rgm;
                for (int i = 0; i < _avaliableNormRegimesOnIntervals[idx].Count(); i++)
                {
                    soft = true;
                    strong = true;
                    rgm = _avaliableNormRegimesOnIntervals[idx][i];
                    for (int j = 0; j < _dimension; j++)
                    {
                        // Если нельзя качать по трубе, а качаем, то режим не подходит
                        if (currentTarget[j] < _minFlows[j] && rgm[j] > 0 && _constraints[j])
                        {
                            soft = false;
                            break;
                        }

                        // Если нужно качать по трубе, но не качаем, то режим не подходит
                        if (strong && currentTarget[j] > _minFlows[j] && rgm[j] == 0)
                            strong = false;
                    }

                    if (soft)
                    {
                        softConstrains.Add(rgm);

                        if (strong)
                            strongConstrains.Add(rgm);
                    }
                }

                // Выбираем доступные для интервала режимы
                List<double[]> currentAvaliableRegimes = strongConstrains.Count() == 0 ? softConstrains : strongConstrains;

                // Нормированный режим для текущего объема перекачки
                //double[] normTargetRegime = currentTarget.Zip(_maxFlows, (x,y) => x / (currentPeriod * y)).ToArray();
                double[] targetRegime =  FormScheduleContinuousRegimeField(targets, notSortedSchedule, idx)[0].Zip(_maxFlows, (x,y) => x / y).ToArray();

                // Ищем ближайший к нему режим среди доступных
                double[] bestRegime = AlgorithmHelper.NearestByDistance(currentAvaliableRegimes, targetRegime);

                notSortedSchedule[idx] = _regimes[_normRegimes.IndexOf(bestRegime)];

                // Уменьшаем цель
                currentTarget = currentTarget.Zip(notSortedSchedule[idx], (x, y) => x - y).ToArray();

                // Уменьшаем время
                currentPeriod--;
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Расставим полученные режимы (только те, которые не фиксированы)
            //////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Сгруппируем режимы, получим количество каждого режима
            var grouped = notSortedSchedule
                .Where((x,i) => fixValues[i][0] < 0)
                .GroupBy(x => x).ToDictionary(x => x.Key, x => x.Count());
            // Получим доступные на интервалах режимы
            List<List<double[]>> avaliableRegimes = _avaliableRegimesOnIntervals.Select(x => x.Where(y => grouped.ContainsKey(y)).ToList()).ToList();

            // Пока просто отсортируем режимы в пределах интервалов с одинаковыми возможными режимами

            // Сгруппируем индексы по списку доступных режимов
            Dictionary<List<double[]>, List<int>> indexesOfSameRegimes = new Dictionary<List<double[]>, List<int>>();
            for (int i = 0; i < _period; i++)
            {
                if (fixValues[i][0] >= 0)
                {
                    continue;
                }

                List<double[]> cur = avaliableRegimes[i];
                List<double[]> key = null;
                foreach (var el in indexesOfSameRegimes)
                {
                    if (AlgorithmHelper.ListsEqualityNoOrder(cur, el.Key))
                    {
                        key = el.Key;
                        break;
                    }
                }
                if (key == null)
                {
                    indexesOfSameRegimes.Add(cur, new List<int>());
                    key = cur;
                }

                indexesOfSameRegimes[key].Add(i);
            }

            // Сортируем в пределах возможных интервалов (по дистанции между режимами и координатой (0,...,0) )

            List<double[]> sortedSchedule = fixValues.Select(x => x.Select(y => y).ToArray()).ToList();
            foreach (var el in indexesOfSameRegimes)
            {
                List<double[]> forSort = el.Value.Select((x) => _normRegimes[_regimes.IndexOf(notSortedSchedule[x])]).ToList();
                forSort.Sort((x, y) =>
                {
                    double dx = AlgorithmHelper.GetLength(x), dy = AlgorithmHelper.GetLength(y);

                    if (dx == dy) return 0;
                    else if (dx > dy) return -1;
                    else return 1;
                });
                for (int i = 0; i < el.Value.Count(); i++)
                {
                    sortedSchedule[el.Value[i]] = _regimes[_normRegimes.IndexOf(forSort[i])].Select(x => x).ToArray();
                }
            }

            return sortedSchedule;
        }

        public List<double[]> FormScheduleContinuousRegimeField(double[] targets, List<double[]> fixValues = null, int neededPoint = -1)
        {
            if (fixValues == null)
            {
                fixValues = AlgorithmHelper.CreateListOfArrays(_period, _dimension, -1.0);
            }

            List<double[]> result = AlgorithmHelper.CombineIntoList(_algorithms.Select((x, i) => x.FormScheduleContinuousRegimeField(targets[i], fixValues.Select(y => y[i]).ToList(), neededPoint)).ToArray());
            
            return result;

        }

        public List<double[]> FormScheduleByVectors(double[] targets, List<double[]> fixValues = null)
        {
            if (fixValues == null)
            {
                fixValues = AlgorithmHelper.CreateListOfArrays(_period, _dimension, -1.0);
            }

            if (targets == null || fixValues == null || fixValues.Count() != _period || fixValues.Any(x => x == null)) throw new ArgumentNullException();
            if (fixValues.Any(x => x.Count() != targets.Count())) throw new ArgumentException();

            double[] normTargets = targets.Zip(_maxFlows, (x, y) => x / y).ToArray();

            List<double[]> decompositions = _normRegimesCombinations.Select(x => GetDimensionDecomposition(normTargets, x, _period)).ToList();
            List<double[]> test = decompositions.Select((x, i) =>
            {
                if (x == null)
                    return null;
                var vol = _regimesCombinations[i].Zip(x, (comb, decomp) => comb.Select(c => c * decomp).ToArray()).ToList();
                return AlgorithmHelper.GetSumOnInterval(vol, 0, vol.Count());
            }).ToList();

            return null;
        }

        private double[] GetDimensionDecomposition(double[] targets, List<double[]> regimes, double period)
        {
            if (targets == null || regimes == null || regimes.Count() == 0 || regimes.Any(x => x == null )) throw new ArgumentNullException();
            if (regimes.Any(x => x.Count() != targets.Count())) throw new ArgumentNullException();
            if (regimes.Count() != targets.Count() + 1) throw new ArgumentException();

            int size = regimes.Count();

            double[][] A = new double[size][];
            for (int i = 0; i < size - 1; i++)
            {
                A[i] = new double[size];
                for (int  j = 0; j < size; j++)
                {
                    A[i][j] = regimes[j][i];
                }
            }
            A[size - 1] = regimes.Select(x => 1.0).ToArray();

            double[] B = new double[size];
            for (int i = 0; i < size - 1; i++)
            {
                B[i] = targets[i];
            }
            B[size - 1] = period;

            double[] X = A.Solve(B, true);

            if (X.Any(x => x < 0)) return null;
            else return X;            
        }

        #endregion

        #region Реализация интерфейса TechnologicalSectionAlgorithm

        public List<double[]> GetSchedule(double[] targets, List<double[]> fixValues = null)
        {
            if (fixValues == null)
            {
                fixValues = AlgorithmHelper.CreateListOfArrays(Period, Dimension, -1.0);
            }

            return FormSchedule(targets, targets.Select(x => 1.0).ToArray(), fixValues);
        }

        public List<double[]> GetContinuousSchedule(double[] targets, List<double[]> fixValues = null)
        {
            return FormScheduleContinuousRegimeField(targets, fixValues);
        }

        public bool IsRegimeAvaliableOnInterval(double[] regime, int interval)
        {
            for (int i = 0; i < _avaliableRegimesOnIntervals[interval].Count(); i++)
            {
                if (_avaliableRegimesOnIntervals[interval][i].SequenceEqual(regime))
                    return true;
            }

            return false;
        }

        public List<double[]> Regimes
        {
            get
            {
                return _regimes.Select(x => x.Select(y => y).ToArray()).ToList();
            }
        }

        public List<double[]> NormRegimes
        {
            get
            {
                return _normRegimes.Select(x => x.Select(y => y).ToArray()).ToList();
            }
        }

        public int Period
        {
            get { return _period; }
        }

        public int Dimension
        {
            get { return _dimension; }
        }

        public List<double[]> UpperBounds
        {
            get
            {
                return _upperBounds.Select(x => x.Select(y => y).ToArray()).ToList();
            }
        }

        public double[] MaxFlows
        {
            get
            {
                return _maxFlows.Select(x => x).ToArray();
            }
        }

        #endregion

        #region Ген. алг. Не нужен пока.
        //        /// <summary>
        //        /// Формирование расписания генетическим алгоритмом
        //        /// </summary>
        //        /// <param name="targets">Список целей по всем трубам</param>
        //        /// <param name="weights">Список весов важности труб</param>
        //        /// <returns></returns>
        //        public List<DiscreteSchedule> GetGeneticSchedule(double[] targets, double[] weights)
        //        {
        //            if (targets == null || weights == null)
        //            {
        //                throw new ArgumentNullException();
        //            }

        //            GeneticMachine gm = new GeneticMachine(targets, _avaliableRegimesOnIntervals, new GeneticMachine.GeneticParameters()
        //            {
        //                populationNumber = 50,
        //                childrenNumber = 50,
        //                mutateProbability = 0.4,
        //                childrenNextGenerationNumber = 20,
        //                eliteParentsNextGenerationNumber = 5,
        //                mutatsNextGenerationNumber = 15
        //            });

        //            return gm.Run(1000).GetSchedules();
        //        }

        //        #region Собственная реализация генетического алгоритма        

        //        private class TechnologicalSectionIndividual
        //        {
        //            private double[] _targetRegime;
        //            private List<List<double[]>> _avaliableRegimesOnIntervals;
        //            private List<List<double>> _probabilities;

        //            private int _maxChangesNumber;
        //            private double[] _maxDifference;
        //            private double[] _target;

        //            private Random _rnd = new Random();

        //            private List<double[]> _schedule;

        //            public double Fitness
        //            {
        //                get;
        //                private set;
        //            }

        //            public TechnologicalSectionIndividual(TechnologicalSectionIndividual i)
        //            {
        //                _targetRegime = i._targetRegime;
        //                _avaliableRegimesOnIntervals = i._avaliableRegimesOnIntervals;
        //                _probabilities = i._probabilities;
        //                _maxChangesNumber = i._maxChangesNumber;
        //                _maxDifference = i._maxDifference;
        //                _target = i._target;
        //            }

        //            // Базовый конструктор, осуществляющий все предрасчеты
        //            public TechnologicalSectionIndividual(double[] target, List<List<double[]>> avaliableRegimesOnIntervals)
        //            {
        //                _avaliableRegimesOnIntervals = avaliableRegimesOnIntervals;
        //                _target = target;
        //                double[] maxTarget = _avaliableRegimesOnIntervals.Aggregate(target.Select(x => 0.0), (totalSum, regimes) =>
        //                {
        //                    var max = regimes.Aggregate(target.Select(x => 0.0), (totalMax, regime) =>
        //                    {
        //                        return totalMax.Zip(regime, (x, y) => x > y ? x : y);
        //                    });
        //                    return totalSum.Zip(max, (x, y) => x + y);
        //                }).ToArray();
        //                _maxDifference = _target.Zip(maxTarget, (x, y) => x > y / 2 ? x : y - x).ToArray();
        //                _targetRegime = maxTarget.Zip(target, (x, y) => (x < y ? x : y) / avaliableRegimesOnIntervals.Count()).ToArray();

        //                double maxLength = Math.Sqrt(_targetRegime.Count());
        //                _probabilities = _avaliableRegimesOnIntervals.Select((regimesOnInterval) =>
        //                {
        //                    if (regimesOnInterval.Count() == 1)
        //                    {
        //                        return new List<double>() { 1 };
        //                    }

        //                    // Ищем максимум и минимум по каждой оси
        //                    double[] max = _targetRegime.Select(x => x).ToArray(), min = _targetRegime.Select(x => x).ToArray();
        //                    foreach (var regime in regimesOnInterval)
        //                    {
        //                        max = max.Zip(regime, (r, m) => r > m ? r : m).ToArray();
        //                        min = min.Zip(regime, (r, m) => r < m ? r : m).ToArray();
        //                    }

        //                    // Нормируем каждую ось и ищем расстояния от каждого режима до целевого режима
        //                    double[] normTarget = _targetRegime.Select((x, i) => (x - min[i]) / (max[i] - min[i])).ToArray();
        //                    List<double> result = regimesOnInterval
        //                        .Select(regime => regime.Select((x, i) => (x - min[i]) / (max[i] - min[i])))
        //                        .Select(normRegime => Math.Sqrt(normRegime.Zip(normTarget, (x, y) => x - y).Sum(x => x * x)) / maxLength)
        //                        .Select(length => maxLength * 1.001 - length)
        //                        .ToList();

        //                    return result.Select(x => x / result.Sum()).ToList();
        //                }).ToList();
        //                _maxChangesNumber = _avaliableRegimesOnIntervals.Count() * 3 / 24;
        //            }

        //            // Случайная генерация расписания
        //            private void GenerateSchedule()
        //            {
        //                int period = _avaliableRegimesOnIntervals.Count();
        //                _schedule = new List<double[]>();
        //                for (int i = 0; i < period; i++)
        //                {
        //                    _schedule.Add(_avaliableRegimesOnIntervals[i][AlgorithmHelper.GetRouletIndex(_probabilities[i])]);
        //                }                
        //            }

        //            // Расчет фитнесфункции
        //            private void CalcFitness()
        //            {
        //                double volumeDifferenceCriteria = _schedule
        //                    .Aggregate(_target.Select(x => 0.0), (total, current) => total.Zip(current, (x, y) => x + y))
        //                    .Zip(_target, (x, y) => Math.Abs(x - y))
        //                    .Zip(_maxDifference, (x, y) => x / y)
        //                    .Sum(x => x * x);
        //                volumeDifferenceCriteria = Math.Sqrt(volumeDifferenceCriteria) / Math.Sqrt(_target.Count());                    

        //                double changesCriteria = 0;
        //                for (int i = 1; i < _avaliableRegimesOnIntervals.Count(); i++)
        //                {
        //                    if (_schedule[i] != _schedule[i - 1])
        //                    changesCriteria++;
        //                }
        //                changesCriteria /= _maxChangesNumber;

        //                Fitness = Math.Max(volumeDifferenceCriteria, changesCriteria);
        //            }

        //            // Создает пустой клон, без расписания
        //            private TechnologicalSectionIndividual CreateEmptyClone()
        //            {
        //                return new TechnologicalSectionIndividual(this);
        //            }

        //            // Создает полный клон
        //            private TechnologicalSectionIndividual CreateFullClone()
        //            {
        //                TechnologicalSectionIndividual individ = CreateEmptyClone();
        //                individ._schedule = _schedule.Select(x => x).ToList();
        //                individ.Fitness = Fitness;
        //                return individ;
        //            }

        //            // Возврящает мутанта из данного индивида (мутирует отрезок расписания)
        //            public TechnologicalSectionIndividual Mutate()
        //            {
        //                TechnologicalSectionIndividual mutant = CreateFullClone();
        //                int period = _avaliableRegimesOnIntervals.Count();
        //                int start = _rnd.Next(period - 1), end = _rnd.Next(start + 1, period + 1);
        //                double[] value = _avaliableRegimesOnIntervals[start][0].Select(x => -1.0).ToArray();
        //                for (int i = start; i < end; i++)
        //                {
        //                    if (_avaliableRegimesOnIntervals[i].IndexOf(value) == -1)
        //                    {
        //                        value = _avaliableRegimesOnIntervals[i][_rnd.Next(_avaliableRegimesOnIntervals[i].Count())];
        //                    }
        //                    mutant._schedule[i] = value;
        //                }
        //                mutant.CalcFitness();
        //                return mutant;
        //            }

        //            // Возвращает случайную особь
        //            public TechnologicalSectionIndividual CreateRandomIndividual()
        //            {
        //                TechnologicalSectionIndividual individ = CreateEmptyClone();
        //                individ.GenerateSchedule();
        //                individ.CalcFitness();
        //                return individ;
        //            }

        //            // Скрещивает две особи (от 1 до 3 точек случайным образом)
        //            public TechnologicalSectionIndividual Crossingover(TechnologicalSectionIndividual pair)
        //            {
        //                TechnologicalSectionIndividual individ = CreateEmptyClone();
        //                int numberOfPoints = _rnd.Next(1, (int)(0.1 * _avaliableRegimesOnIntervals.Count())) + 2;
        //                int period = _avaliableRegimesOnIntervals.Count();
        //                List <int> points = new List<int>() { 0, period };
        //                while (points.Count() < numberOfPoints)
        //                {
        //                    int point = _rnd.Next(1, period - 1);
        //                    if (points.IndexOf(point) == -1)
        //                    {
        //                        points.Add(point);
        //                    }
        //                }
        //                points.Sort();

        //                individ._schedule = new List<double[]>();
        //                for (int i = 1; i < numberOfPoints; i++)
        //                {
        //                    TechnologicalSectionIndividual current = _rnd.NextDouble() < 0.5 ? this : pair;
        //                    for(int j = points[i - 1]; j < points[i]; j++)
        //                    {
        //                        individ._schedule.Add(current._schedule[j].Clone() as double[]);
        //                    }
        //                }

        //                individ.CalcFitness();

        //                return individ;
        //            }

        //            public List<DiscreteSchedule> GetSchedules()
        //            {
        //                List<DiscreteSchedule> schedules = _targetRegime.Select(x => new DiscreteSchedule(_schedule.Count())).ToList();
        //                for(int i = 0; i < _schedule.Count(); i++)
        //                {
        //                    for (int j = 0; j < _schedule[i].Count(); j++)
        //                    {
        //                        schedules[j][i] = _schedule[i][j];
        //                    }
        //                }
        //                return schedules;
        //            }
        //        }

        //        private class GeneticMachine
        //        {

        //            #region Параметры генетического алгоритма

        //            public struct GeneticParameters
        //            {
        //                public int populationNumber;
        //                public int childrenNumber;

        //                public double mutateProbability;

        //                public int eliteParentsNextGenerationNumber;
        //                public int childrenNextGenerationNumber;
        //                public int mutatsNextGenerationNumber;
        //            }
        //            private GeneticParameters _parameters;

        //            #endregion

        //            private TechnologicalSectionIndividual _baseIndivid;
        //            private List<TechnologicalSectionIndividual> _currentPopulation;
        //            private int _age;
        //            private Random _rnd = new Random();
        //            private TechnologicalSectionIndividual _bestFitness;

        //            public TechnologicalSectionIndividual BestFitness
        //            {
        //                get
        //                {
        //                    return _bestFitness;
        //                }
        //            }

        //            public GeneticMachine(double[] targets, List<List<double[]>> avaliableRegimesOnIntervals, GeneticParameters parameters)
        //            {
        //                _parameters = parameters;
        //                _baseIndivid = new TechnologicalSectionIndividual(targets, avaliableRegimesOnIntervals);
        //                _currentPopulation = CreateRandomPopulation(_parameters.populationNumber);
        //                _bestFitness = GetBestFitness();
        //                _age = 0;
        //            }

        //            private List<TechnologicalSectionIndividual> CreateRandomPopulation(int number)
        //            {
        //                ConcurrentBag<TechnologicalSectionIndividual> randomPopulation = new ConcurrentBag<TechnologicalSectionIndividual>();
        //#if (PARALLEL)
        //                Parallel.For(0, number, (i) =>
        //                {
        //                    randomPopulation.Add(_baseIndivid.CreateRandomIndividual());
        //                });
        //#else
        //                for (int i = 0; i < number; i++)
        //                    randomPopulation.Add(_baseIndivid.CreateRandomIndividual());
        //#endif
        //                return randomPopulation.ToList();
        //            }

        //            private List<TechnologicalSectionIndividual> MutatePopulation(List<TechnologicalSectionIndividual> population)
        //            {
        //                List<TechnologicalSectionIndividual> populationForMutate = new List<TechnologicalSectionIndividual>();
        //                foreach(var individ in population)
        //                {
        //                    if (_parameters.mutateProbability < _rnd.NextDouble())
        //                    {
        //                        populationForMutate.Add(individ);
        //                    }
        //                }

        //                ConcurrentBag<TechnologicalSectionIndividual> mutantPopulation = new ConcurrentBag<TechnologicalSectionIndividual>();
        //#if (PARALLEL)
        //                Parallel.ForEach(populationForMutate, (individ) =>
        //                {
        //                    mutantPopulation.Add(individ.Mutate());
        //                });
        //#else
        //                foreach (var individ in populationForMutate)
        //                    mutantPopulation.Add(individ.Mutate());
        //#endif

        //                return mutantPopulation.ToList();
        //            }

        //            private List<TechnologicalSectionIndividual> Crossbreeding(List<TechnologicalSectionIndividual> parents)
        //            {
        //                List<Tuple<TechnologicalSectionIndividual, TechnologicalSectionIndividual>> pairs = new List<Tuple<TechnologicalSectionIndividual, TechnologicalSectionIndividual>>();
        //                List<double> probabilities = parents.Select(x => 1 / (0.001 + x.Fitness)).ToList();
        //                probabilities = probabilities.Select(x => x / probabilities.Sum()).ToList();
        //                while (pairs.Count() < _parameters.childrenNumber)
        //                {
        //                    int p1 = AlgorithmHelper.GetRouletIndex(probabilities),
        //                        p2 = AlgorithmHelper.GetRouletIndex(probabilities);
        //                    if (p1 != p2)
        //                    {
        //                       pairs.Add(new Tuple<TechnologicalSectionIndividual, TechnologicalSectionIndividual>(parents[p1], parents[p2]));
        //                    }
        //                }

        //                ConcurrentBag<TechnologicalSectionIndividual> children = new ConcurrentBag<TechnologicalSectionIndividual>();
        //#if(PARALLEL)
        //                Parallel.ForEach(pairs, (pair) =>
        //                {
        //                    children.Add(pair.Item1.Crossingover(pair.Item2));
        //                });
        //#else
        //                foreach(var pair in pairs)
        //                {
        //                    children.Add(pair.Item1.Crossingover(pair.Item2));
        //                }
        //#endif

        //                return children.ToList();
        //            }

        //            private List<TechnologicalSectionIndividual> Selection(List<TechnologicalSectionIndividual> parents, List<TechnologicalSectionIndividual> children, List<TechnologicalSectionIndividual> mutants)
        //            {
        //                parents.Sort((x, y) => (x.Fitness > y.Fitness) ? 1 : (x.Fitness < y.Fitness ? -1 : 0));
        //                children.Sort((x, y) => (x.Fitness > y.Fitness) ? 1 : (x.Fitness < y.Fitness ? -1 : 0));
        //                mutants.Sort((x, y) => (x.Fitness > y.Fitness) ? 1 : (x.Fitness < y.Fitness ? -1 : 0));

        //                List<TechnologicalSectionIndividual> nextGeneration = new List<TechnologicalSectionIndividual>();

        //                // Элитных родителей
        //                nextGeneration.AddRange(parents.GetRange(0, _parameters.eliteParentsNextGenerationNumber));

        //                // Детей
        //                int eliteChildren = _parameters.childrenNextGenerationNumber / 2,
        //                    randomChildren = _parameters.childrenNextGenerationNumber - eliteChildren;
        //                // Элитных детей
        //                nextGeneration.AddRange(children.GetRange(0, _parameters.childrenNextGenerationNumber / 2));
        //                // Рандомных детей
        //                for (int i = 0; i < randomChildren; i++)
        //                {
        //                    nextGeneration.Add(children[_rnd.Next(eliteChildren, children.Count())]);
        //                }

        //                // Мутантов
        //                if (_parameters.mutatsNextGenerationNumber >= mutants.Count())
        //                {
        //                    nextGeneration.AddRange(mutants);
        //                }
        //                else
        //                {
        //                    int eliteMutants = _parameters.mutatsNextGenerationNumber / 2,
        //                        randomMutants = _parameters.mutatsNextGenerationNumber - eliteMutants;
        //                    // Элитных мутантов
        //                    nextGeneration.AddRange(mutants.GetRange(0, eliteMutants));
        //                    // Рандомных мутантов
        //                    for (int i = 0; i < randomMutants; i++)
        //                    {
        //                        nextGeneration.Add(mutants[_rnd.Next(eliteMutants, mutants.Count())]);
        //                    }
        //                }

        //                nextGeneration.AddRange(CreateRandomPopulation(_parameters.populationNumber - nextGeneration.Count()));

        //                return nextGeneration;
        //            }

        //            private TechnologicalSectionIndividual GetBestFitness()
        //            {
        //                double best = _currentPopulation.Min(x => x.Fitness);
        //                return _currentPopulation.First(x => x.Fitness == best);
        //            }

        //            public void Step()
        //            {
        //                List<TechnologicalSectionIndividual> parents = _currentPopulation;
        //                List<TechnologicalSectionIndividual> children = Crossbreeding(parents);
        //                List<TechnologicalSectionIndividual> mutants = MutatePopulation(parents.Concat(children).ToList());
        //                _currentPopulation = Selection(parents, children, mutants);
        //                _bestFitness = GetBestFitness();
        //                _age++;
        //            }

        //            public TechnologicalSectionIndividual Run(int maxAge)
        //            {
        //                for (int i = 0; i < maxAge; i++)
        //                    Step();

        //                return _bestFitness;
        //            }
        //        }

        //        #endregion

        #endregion
    }
}
