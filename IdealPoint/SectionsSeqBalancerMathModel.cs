#define PARALLEL

using GeneticSharp.Domain.Chromosomes;
using GeneticSharp.Domain.Fitnesses;
using GeneticSharp.Domain.Mutations;
using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using GeneticSharp.Domain.Crossovers;
using System.Collections.Concurrent;
using System.Threading.Tasks;
using System.Threading;

namespace Algorithms
{
    public class SectionsSeqBalancerMathModel
    {
        #region Поля

        List<double> _reservoirVolumes;
        List<double> _oilStartVolumes;
        List<ISection> _sections;
        int _period;
        int _resCount;

        List<List<Tuple<List<double[]>, List<int>>>> _tempSolutions;
        List<TargetVolumes> _tempTargetVolumes;
        List<List<double>> _tempPumpSchedules;

        #endregion

        #region Свойства



        #endregion

        #region Кострукторы

        public SectionsSeqBalancerMathModel(List<double> reservoirVolumes, List<double> oilStartVolumes, List<ISection> sections, int period)
        {
            if (reservoirVolumes == null)
                throw new Exception();
            if (reservoirVolumes.Count() == 0)
                throw new Exception();
            if (reservoirVolumes.Any(x => x < 0))
                throw new Exception();
            if (oilStartVolumes == null)
                throw new Exception();
            if (oilStartVolumes.Count() != reservoirVolumes.Count())
                throw new Exception();
            if (reservoirVolumes.Zip(oilStartVolumes, (x,y) => y < 0 || y > x).Any(x => x))
                throw new Exception();
            if (sections == null)
                throw new Exception();
            if (sections.Count() != reservoirVolumes.Count() + 1)
                throw new Exception();
            // Секции должны быть по обоим сторонам резервуара, пусть и null
            // Но в других местах null не может быть
            for (int i = 1; i < sections.Count() - 1; i++)
                if (sections[i] == null)
                    throw new Exception();
            if (sections.Any(x => x != null && x.Period != period))
                throw new Exception();


            _reservoirVolumes = reservoirVolumes.ToList();
            _sections = sections.ToList();
            _oilStartVolumes = oilStartVolumes.ToList();
            _period = period;
            _resCount = reservoirVolumes.Count();
        }

        #endregion

        #region Вспомогательные структуры и классы
        
        struct TempBalanceStruct
        {
            public List<double[]> initialSchedule;
            public List<double> startReservoirSchedule, endReservoirSchedule;
            public List<double> startReservoirPumpsSchedule, endReservoirPumpsSchedule;
            public List<int> startReservoirCrashIndexes, endReservoirCrashIndexes;
            public double startReservoirPrevOverfillVolume, endReservoirPrevOverfillVolume;
            public double startReservoirStartVolume, endReservoirStartVolume;
            public double startReservoirVolume, endReservoirVolume;
            public int sectionNumber;
        }

        struct TempBagStruct
        {
            public int sectionNumber;
            public double startReservoirOverfillVolume;
            public double endReservoirOverfillVolume;
            public int changesCount;
            public double changesVolume;
            public List<int> block1;
            public List<int> block2;
        }

        #endregion

        #region Методы

        public void SetTempParams(List<TargetVolumes> tempTargetVolumes, List<List<double>> tempPumpSchedules)
        {
            if (tempPumpSchedules == null)
                throw new Exception();
            if (tempPumpSchedules.Count() != _resCount)
                throw new Exception();
            if (tempPumpSchedules.Any(x => x == null))
                throw new Exception();
            if (tempPumpSchedules.Any(x => x.Count() != _period))
                throw new Exception();
            if (tempTargetVolumes == null)
                throw new Exception();
            else if (tempTargetVolumes.Count() != _sections.Count())
                throw new Exception();
            for (int i = 0; i < _sections.Count(); i++)
                if (_sections[i] != null && tempTargetVolumes[i] == null)
                    throw new Exception();
                else if (tempTargetVolumes[i] != null)
                    tempTargetVolumes[i].Check(_sections[i].Dimension);

            _tempPumpSchedules = tempPumpSchedules.Select(x => x.ToList()).ToList();
            _tempTargetVolumes = tempTargetVolumes.Select(x => x == null ? null : new TargetVolumes(x)).ToList();
            _tempSolutions = null;
        }

        private void CreateTempSolution()
        {
            _tempSolutions = new List<List<Tuple<List<double[]>, List<int>>>>();
            for (int i = 0; i < _sections.Count(); i++)
            {
                var section = _sections[i];
                if (section != null)
                {
                    var volume = _tempTargetVolumes[i];
                    section.CalcDefaultIntervalsParameters(volume);
                    _tempSolutions.Add(section.GetSolution(volume));
                }
                else
                {
                    _tempSolutions.Add(null);
                }
            }
        }
                
        public List<List<double[]>> Balance(Action<string> informationAction = null)
        {
            if (informationAction == null)
                informationAction = (str) => { };

            List<List<double[]>> initialSchedules = null;
            List<double> prevOverfillVolumes = new List<double>();

            List<List<double[]>> ConvertResult()
            {
                if (_sections[0] == null)
                    initialSchedules.RemoveAt(0);
                else if (_sections.Last() == null)
                    initialSchedules.RemoveAt(initialSchedules.Count() - 1);
                return initialSchedules;
            }
            void CalcInitialSchedules()
            {
                initialSchedules = _tempSolutions.Select((x, i) =>
                {
                    if (_sections[i] == null)
                        return AlgorithmHelper.CreateListOfArrays(_period, 1, 0.0);
                    else
                        return _sections[i].GetFullSchedule(CreateInitialSchedule(x));
                }).ToList();
            }
            void CalcOverfill()
            {
                prevOverfillVolumes = new List<double>();
                for (int i = 0; i < initialSchedules.Count() - 1; i++)
                {
                    var reservoirSchedule = GetReservoirSchedule(initialSchedules[i], initialSchedules[i + 1], _tempPumpSchedules[i], _oilStartVolumes[i]);
                    prevOverfillVolumes.Add(GetOverfillVolume(reservoirSchedule, 0, _reservoirVolumes[i]));
                }
            }

            // Сначала без анализа и перестановок
            CreateTempSolution();
            CalcInitialSchedules();
            CalcOverfill();
            if (prevOverfillVolumes.All(x => x == 0.0))
            {
                return ConvertResult();
            }

            int maxIter = 1;
            for (int i = 0; i < maxIter; i++)
            {
                // Перестановки
                CreateTempSolution();
                CalcInitialSchedules();
                initialSchedules = Permute(initialSchedules, 24, 100, informationAction);

                CalcOverfill();
                if (prevOverfillVolumes.All(x => x == 0.0))
                {
                    return ConvertResult();
                }

                // Анализ и изменение объемов
                if (!Analyse(initialSchedules))
                    break;

                CreateTempSolution();
                CalcInitialSchedules();
                CalcOverfill();
                if (prevOverfillVolumes.All(x => x == 0.0))
                {
                    return ConvertResult();
                }
            }

            if (initialSchedules == null)
                return null;
            else
                return ConvertResult();
        }
        
        private List<List<double[]>> Permute(List<List<double[]>> initialSchedules, int blockSize, int maxIter, Action<string> informationAction)
        {
            initialSchedules = initialSchedules.Select(x => x.Select(y => y.ToList().ToArray()).ToList()).ToList();
            List<List<double>> reservoirInitialSchedules = new List<List<double>>();
            List<List<int>> reservoirCrashIndexes = new List<List<int>>();
            List<double> prevOverfillVolumes = new List<double>();
            for (int i = 0; i < initialSchedules.Count() - 1; i++)
            {
                reservoirInitialSchedules.Add(GetReservoirSchedule(initialSchedules[i], initialSchedules[i + 1], _tempPumpSchedules[i], _oilStartVolumes[i]));
                prevOverfillVolumes.Add(GetOverfillVolume(reservoirInitialSchedules.Last(), 0, _reservoirVolumes[i]));
                reservoirCrashIndexes.Add(GetCrashIndexes(reservoirInitialSchedules.Last(), 0, _reservoirVolumes[i]));
            }

            int curIter = 0;
            TempBalanceStruct GetInitStruct(int sectionNumber)
            {
                var tempStruct = new TempBalanceStruct()
                {
                    sectionNumber = sectionNumber,
                    initialSchedule = initialSchedules[sectionNumber],
                    startReservoirCrashIndexes = new List<int>(),
                    endReservoirCrashIndexes = new List<int>()
                };
                if (sectionNumber != 0)
                {
                    tempStruct.startReservoirSchedule = reservoirInitialSchedules[sectionNumber - 1];
                    tempStruct.startReservoirCrashIndexes = reservoirCrashIndexes[sectionNumber - 1];
                    tempStruct.startReservoirPrevOverfillVolume = prevOverfillVolumes[sectionNumber - 1];
                    tempStruct.startReservoirStartVolume = _oilStartVolumes[sectionNumber - 1];
                    tempStruct.startReservoirPumpsSchedule = _tempPumpSchedules[sectionNumber - 1];
                    tempStruct.startReservoirVolume = _reservoirVolumes[sectionNumber - 1];
                }
                if (sectionNumber != _sections.Count() - 1)
                {
                    tempStruct.endReservoirSchedule = reservoirInitialSchedules[sectionNumber];
                    tempStruct.endReservoirCrashIndexes = reservoirCrashIndexes[sectionNumber];
                    tempStruct.endReservoirPrevOverfillVolume = prevOverfillVolumes[sectionNumber];
                    tempStruct.endReservoirStartVolume = _oilStartVolumes[sectionNumber];
                    tempStruct.endReservoirPumpsSchedule = _tempPumpSchedules[sectionNumber];
                    tempStruct.endReservoirVolume = _reservoirVolumes[sectionNumber];
                }
                return tempStruct;
            }
            //bool[][] crashIndexesBool = reservoirCrashIndexes.Select((x, i) => new bool[_period].Select((y,j) => x.Contains(j)).ToArray()).ToArray();
            while (maxIter > curIter++)
            {
                ConcurrentBag<TempBagStruct> variantBag = new ConcurrentBag<TempBagStruct>();
                for (int sectionNumber = 0; sectionNumber < _sections.Count(); sectionNumber++)
                {
                    informationAction($"Итерация {curIter}, переполненность {prevOverfillVolumes.Sum()}, расчет секции {sectionNumber}, найдено вариантов {variantBag.Count()}");
                    var section = _sections[sectionNumber];
                    if (section == null)
                        continue;

                    TempBalanceStruct tempStruct = GetInitStruct(sectionNumber);

                    foreach (var intervals in section.ControlAvaliableIntervals)
                    {
                        int intervalCount = intervals.Count();
                        // Получим все индексы, которые будем просматривать
                        List<int[]> checkIndexes =
                            Vector
                            .EnumerableRange(0, intervalCount * intervalCount)
                            .AsParallel()
                            .Select(idx =>
                            {
                                int i = idx / intervalCount, j = idx % intervalCount;
                                int i1 = intervals[i];
                                if (intervalCount - i > 2 * blockSize && j >= i + blockSize && j + blockSize < intervalCount)
                                {
                                    int i2 = intervals[j + blockSize];
                                    if (tempStruct.startReservoirCrashIndexes.Any(x => x <= i2 && x >= i1) || tempStruct.endReservoirCrashIndexes.Any(x => x <= i2 && x >= i1))
                                        return new int[] { i, j };
                                    else
                                        return null;
                                }
                                else
                                    return null;
                            })
                            .Where(x => x != null)
                            .ToList();
                        
                        Parallel.ForEach(checkIndexes, (indexes) =>
                        {
                            var block1 = intervals.GetRange(indexes[0], blockSize).ToList();
                            var block2 = intervals.GetRange(indexes[1], blockSize).ToList();

                            var newSchedule = SwapBlocks(tempStruct.initialSchedule, block1, block2);

                            double currentStartReservoirVolume = tempStruct.startReservoirStartVolume,
                                currentEndReservoirVolume = tempStruct.endReservoirStartVolume;
                            double currentStartReservoirOverfillVolume = 0.0,
                                currentEndReservoirOverfillVolume = 0.0;
                            double currentChangesVolume = 0.0;
                            int currentChangesCount = 0;

                            if (tempStruct.startReservoirSchedule != null)
                            {
                                currentStartReservoirVolume += -newSchedule[0].First() + tempStruct.startReservoirPumpsSchedule[0] + initialSchedules[sectionNumber - 1][0].Last();
                                if (currentStartReservoirVolume > tempStruct.startReservoirVolume)
                                {
                                    currentStartReservoirOverfillVolume += currentStartReservoirVolume - tempStruct.startReservoirVolume;
                                    /*if (!crashIndexesBool[sectionNumber - 1][0])
                                        return;*/
                                }
                                else if (currentStartReservoirVolume < 0)
                                {
                                    currentStartReservoirOverfillVolume -= currentStartReservoirVolume;
                                    /*if (!crashIndexesBool[sectionNumber - 1][0])
                                        return;*/
                                }
                            }
                            if (tempStruct.endReservoirSchedule != null)
                            {
                                currentEndReservoirVolume += newSchedule[0].Last() + tempStruct.endReservoirPumpsSchedule[0] - initialSchedules[sectionNumber + 1][0].First();
                                if (currentEndReservoirVolume > tempStruct.endReservoirVolume)
                                {
                                    currentEndReservoirOverfillVolume += currentEndReservoirVolume - tempStruct.endReservoirVolume;
                                    /*if (!crashIndexesBool[sectionNumber][0])
                                        return;*/
                                }
                                else if (currentEndReservoirVolume < 0)
                                {
                                    currentEndReservoirOverfillVolume -= currentEndReservoirVolume;
                                    /*if (!crashIndexesBool[sectionNumber][0])
                                        return;*/
                                }
                            }

                            for (int i = 1; i < _period; i++)
                            {
                                if (tempStruct.startReservoirSchedule != null)
                                {
                                    currentStartReservoirVolume += -newSchedule[i].First() + tempStruct.startReservoirPumpsSchedule[i] + initialSchedules[sectionNumber - 1][i].Last();
                                    if (currentStartReservoirVolume > tempStruct.startReservoirVolume)
                                    {
                                        currentStartReservoirOverfillVolume += currentStartReservoirVolume - tempStruct.startReservoirVolume;
                                        /*if (!crashIndexesBool[sectionNumber - 1][i])
                                            return;*/
                                    }
                                    else if (currentStartReservoirVolume < 0)
                                    {
                                        currentStartReservoirOverfillVolume -= currentStartReservoirVolume;
                                        /*if (!crashIndexesBool[sectionNumber - 1][i])
                                            return;*/
                                    }
                                }
                                if (tempStruct.endReservoirSchedule != null)
                                {
                                    currentEndReservoirVolume += newSchedule[i].Last() + tempStruct.endReservoirPumpsSchedule[i] - initialSchedules[sectionNumber + 1][i].First();
                                    if (currentEndReservoirVolume > tempStruct.endReservoirVolume)
                                    {
                                        currentEndReservoirOverfillVolume += currentEndReservoirVolume - tempStruct.endReservoirVolume;
                                        /*if (!crashIndexesBool[sectionNumber][i])
                                            return;*/
                                    }
                                    else if (currentEndReservoirVolume < 0)
                                    {
                                        currentEndReservoirOverfillVolume -= currentEndReservoirVolume;
                                        /*if (!crashIndexesBool[sectionNumber][i])
                                            return;*/
                                    }
                                }

                                double a = newSchedule[i - 1][0], b = newSchedule[i][0];
                                if (a != b)
                                {
                                    currentChangesCount++;
                                    currentChangesVolume += Math.Abs(a - b);
                                }
                            }

                            if ((currentStartReservoirOverfillVolume < tempStruct.startReservoirPrevOverfillVolume && !(currentEndReservoirOverfillVolume > tempStruct.endReservoirPrevOverfillVolume))
                                || (currentEndReservoirOverfillVolume < tempStruct.endReservoirPrevOverfillVolume && !(currentStartReservoirOverfillVolume > tempStruct.startReservoirPrevOverfillVolume)))
                                variantBag.Add(new TempBagStruct
                                {
                                    block1 = block1,
                                    block2 = block2,
                                    changesCount = currentChangesCount,
                                    changesVolume = currentChangesVolume,
                                    endReservoirOverfillVolume = currentEndReservoirOverfillVolume,
                                    startReservoirOverfillVolume = currentStartReservoirOverfillVolume,
                                    sectionNumber = sectionNumber
                                });
                        });
                    };
                }

                if (variantBag.Count() == 0)
                    break;

                List<TempBagStruct> variantsList = variantBag.ToList();

                // Отбираем по минимуму парелива
                List<double> zeroes = variantBag.Select(x =>
                {
                    var cur = prevOverfillVolumes.ToList();
                    if (x.sectionNumber != 0)
                        cur[x.sectionNumber - 1] = x.startReservoirOverfillVolume;
                    if (x.sectionNumber != _sections.Count() - 1)
                        cur[x.sectionNumber] = x.endReservoirOverfillVolume;
                    return cur.Sum();
                }).ToList();
                double minZeroes = zeroes.Min();
                variantsList = variantsList.Where((x, i) => zeroes[i] == minZeroes).ToList();

                // Отбираем по минимуму количества переключений
                int minChanges = variantsList.Min(x => x.changesCount);
                variantsList = variantsList.Where(x => x.changesCount == minChanges).ToList();

                // Отбираем по минимуму объема переключений
                double minChangeVolumes = variantsList.Min(x => x.changesVolume);
                variantsList = variantsList.Where(x => x.changesVolume == minChangeVolumes).ToList();

                var selected = variantsList[0];

                int sn = selected.sectionNumber;
                initialSchedules[sn] = SwapBlocks(initialSchedules[sn], selected.block1, selected.block2);
                if (sn != 0)
                {
                    reservoirInitialSchedules[sn - 1] = GetReservoirSchedule(initialSchedules[sn - 1], initialSchedules[sn], _tempPumpSchedules[sn - 1], _oilStartVolumes[sn - 1]);
                    reservoirCrashIndexes[sn - 1] = GetCrashIndexes(reservoirInitialSchedules[sn - 1], 0, _reservoirVolumes[sn - 1]);
                    prevOverfillVolumes[sn - 1] = selected.startReservoirOverfillVolume;
                }
                if (sn < _sections.Count() - 1)
                {
                    reservoirInitialSchedules[sn] = GetReservoirSchedule(initialSchedules[sn], initialSchedules[sn + 1], _tempPumpSchedules[sn], _oilStartVolumes[sn]);
                    reservoirCrashIndexes[sn] = GetCrashIndexes(reservoirInitialSchedules[sn], 0, _reservoirVolumes[sn]);
                    prevOverfillVolumes[sn] = selected.endReservoirOverfillVolume;
                }

                informationAction($"Итерация {curIter}, переполненность {prevOverfillVolumes.Sum()}");
                if (prevOverfillVolumes.All(x => x == 0.0))
                {
                    return initialSchedules;
                }
                //crashIndexesBool = reservoirCrashIndexes.Select((x, i) => new bool[_period].Select((y, j) => x.Contains(j)).ToArray()).ToArray();
            }

            return initialSchedules;
        }

        [Flags]
        enum Decision
        {
            NO_DECISION,
            DECREASE_FLOW_DIFFERENCE_INPUT_RESERVOIR,
            DECREASE_FLOW_DIFFERENCE_OUTPUT_RESERVOIR,
            DECREASE_LEVEL_INPUT_RESERVOIR,
            INCREASE_LEVEL_OUTPUT_RESERVOIR
        }

        private bool Analyse(List<List<double[]>> initialSchedules)
        {
            List<List<double>> reservoirInitialSchedules = new List<List<double>>();
            List<List<int>> reservoirCrashIndexes = new List<List<int>>();
            List<double> prevOverfillVolumes = new List<double>();
            for (int i = 0; i < initialSchedules.Count() - 1; i++)
            {
                reservoirInitialSchedules.Add(GetReservoirSchedule(initialSchedules[i], initialSchedules[i + 1], _tempPumpSchedules[i], _oilStartVolumes[i]));
                prevOverfillVolumes.Add(GetOverfillVolume(reservoirInitialSchedules.Last(), 0, _reservoirVolumes[i]));
                reservoirCrashIndexes.Add(GetCrashIndexes(reservoirInitialSchedules.Last(), 0, _reservoirVolumes[i]));
            }

            Decision decision = Decision.NO_DECISION;
            for (int i = 0; i < _sections.Count(); i++)
            {
                var section = _sections[i];
                if (section == null)
                    continue;

                for (int j = 0; j < section.RepairsIntervals.Count(); j++)
                {
                    var repair = section.RepairsIntervals[j];
                    int repairLen = repair.Count();
                    if (i > 0)
                    {
                        double resVolume = _reservoirVolumes[i - 1];
                        // Смотрим изменения во входном резервуаре на данном участке
                        double startVolume = reservoirInitialSchedules[i - 1][repair.First()],
                            endVolume = reservoirInitialSchedules[i - 1][repair.Last() + 1];

                        if (startVolume  < 0 || startVolume > resVolume || endVolume > 0 || endVolume > resVolume)
                        {
                            double changeVolume = endVolume - startVolume;

                            if (changeVolume < 0)
                            {
                                // Необъяснимо но факт
                                return false;
                            }
                            else if (changeVolume > resVolume)
                            {
                                if (_sections[i - 1] == null)
                                {
                                    // Неконтролируемый скачек больше объема резервуара
                                    return false;
                                }
                                else
                                {
                                    // Нужно уменьшать разницу между расходами
                                    decision = Decision.DECREASE_FLOW_DIFFERENCE_INPUT_RESERVOIR;
                                }
                            }
                            else if (startVolume > 0 && startVolume < resVolume && endVolume > resVolume)
                            {
                                // Нужно перед ремонтом опустошить резервуар 
                                decision = Decision.DECREASE_LEVEL_INPUT_RESERVOIR;
                                if (_sections[i - 1] == null)
                                {
                                    // за счет увеличения расхода в текущей секции
                                }
                                else
                                {
                                    // за счет уменьшения расхода на предыдущей секции
                                }
                            }
                        }
                    }

                    if (i < _resCount)
                    {
                        double resVolume = _reservoirVolumes[i];
                        // Смотрим изменения в выходном резервуаре на данном участке
                        double startVolume = reservoirInitialSchedules[i][repair.First()],
                            endVolume = reservoirInitialSchedules[i][repair.Last() + 1];

                        if (startVolume < 0 || startVolume > resVolume || endVolume > 0 || endVolume > resVolume)
                        {
                            double changeVolume = startVolume - endVolume;

                            if (changeVolume < 0)
                            {
                                // Необъяснимо но факт
                                return false;
                            }
                            else if (changeVolume > resVolume)
                            {
                                if (_sections[i + 1] == null)
                                {
                                    // Неконтролируемый скачек больше объема резервуара
                                    return false;
                                }
                                else
                                {
                                    // Нужно уменьшать разницу между расходами
                                    decision = decision | Decision.DECREASE_FLOW_DIFFERENCE_OUTPUT_RESERVOIR;
                                }
                            }
                            else if (startVolume > 0 && startVolume < resVolume && endVolume < 0)
                            {
                                // Нужно перед ремонтом немного наполнить резервуар
                                decision = decision | Decision.INCREASE_LEVEL_OUTPUT_RESERVOIR;
                                if (_sections[i + 1] == null)
                                {
                                    // за счет увеличения расхода в текущей секции
                                }
                                else
                                {
                                    // за счет уменьшения расхода на следующей секции
                                }
                            }
                        }
                    }
                }

                if (Decision.)
            }

            CreateTempSolution();
            return true;
        }

        private static List<double[]> CreateInitialSchedule(List<Tuple<List<double[]>, List<int>>> initialSolution)
        {
            int period = initialSolution.Sum(x => x.Item2.Count());
            var result = (new double[period][]).ToList();
            foreach (var tuple in initialSolution)
            {
                tuple.Item2.Sort();
                tuple.Item1.Sort((el1, el2) => el1[0] > el2[0] ? 1 : (el1[0] < el2[0] ? -1 : 0));
                for (int i = 0; i < tuple.Item1.Count(); i++)
                {
                    result[tuple.Item2[i]] = tuple.Item1[i];
                }
            }
            return result;
        }

        private static List<double[]> SwapBlocks(List<double[]> scheudule, List<int> block1, List<int> block2)
        {
            int blockSize = block1.Count();
            List<double[]> result = scheudule.ToList();
            for (int k = 0; k < blockSize; k++)
            {
                result[block1[k]] = scheudule[block2[k]];
                result[block2[k]] = scheudule[block1[k]];
            }
            return result;
        }

        private static List<double> GetReservoirSchedule(List<double[]> inputSchedule, List<double[]> outputSchedule, List<double> pumpsSchedule, double startVolume)
        {
            List<double> result = new List<double>() { startVolume };
            for (int i = 0; i < inputSchedule.Count(); i++)
                result.Add(result.Last() + inputSchedule[i].Last() - outputSchedule[i].First() + pumpsSchedule[i]);
            return result;
        }

        private static List<int> GetCrashIndexes(List<double> reservoirSchedule, double minVolume, double maxVolume)
        {
            var result = new List<int>();
            for (int i = 1; i < reservoirSchedule.Count(); i++)
            {
                var x = reservoirSchedule[i];
                if (x < minVolume || x > maxVolume)
                {
                    result.Add(i);
                }
            }
            return result;
        }

        private static double GetOverfillVolume(List<double> reservoirSchedule, double minVolume, double maxVolume)
        {
            double result = 0.0;
            for (int i = 1; i < reservoirSchedule.Count(); i++)
            {
                var x = reservoirSchedule[i];
                if (x < minVolume)
                {
                    result += minVolume - x;
                }
                else if (x > maxVolume)
                {
                    result += x - maxVolume;
                }
            }
            return result;
        }

        private static List<double[]> DecreaseRegime(ISection section, List<double[]> schedule, double volume, List<int> indexes, bool allowZero, bool input = true)
        {
            List<double[]> newSchedule = AlgorithmHelper.GetIndexes(schedule, indexes);
            double newVolume = newSchedule.Sum(x => input ? x.First() : x.Last());
            double oldSectionVolume = newVolume;
            int len = indexes.Count();
            while (true)
            {
                int zeroCounter = 0;
                bool end = false;
                for (int i = 0; i < len; i++)
                {
                    var idx = indexes[i];
                    var lowerRegime = section.GetLowerRegime(idx, newSchedule[i], !input);
                    if (lowerRegime != null)
                    {
                        if (input)
                            newVolume -= newSchedule[i].First() - lowerRegime.First();
                        else
                            newVolume -= newSchedule[i].Last() - lowerRegime.Last();

                        if (oldSectionVolume - newVolume > volume || newVolume < 0)
                        {
                            end = true;
                            break;
                        }
                        newSchedule[i] = lowerRegime;
                    }
                    else
                    {
                        zeroCounter++;
                    }
                }

                if (end)
                    break;

                if (zeroCounter == len)
                {
                    if (allowZero)
                        for (int i = len - 1; i >= 0; i--)
                        {
                            if (input)
                                newVolume -= newSchedule[i].First();
                            else
                                newVolume -= newSchedule[i].Last();

                            if (oldSectionVolume - newVolume > volume || newVolume < 0)
                                break;

                            newSchedule[i] = new double[section.Dimension];
                        }
                    break;
                }
            }

            return newSchedule;
        }

        #endregion
    }
}
