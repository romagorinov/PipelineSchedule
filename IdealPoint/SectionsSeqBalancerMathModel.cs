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
        
        //public class ReservoirMathModelChromosome: IChromosome
        //{
        //    public ReservoirMathModel _model;
        //    public List<int> _sequence;
        //    public List<double[]> _inputSchedule;
        //    public List<double[]> _outputSchedule;

        //    public double? Fitness
        //    {
        //        get;
        //        set;
        //    }

        //    public int Length
        //    {
        //        get;
        //        set;
        //    }

        //    public ReservoirMathModelChromosome(ReservoirMathModel model)
        //    {
        //        _model = model;
        //        Length = _model._period * 2;
        //        Initialize();
        //        Decode();
        //    }

        //    public ReservoirMathModelChromosome(ReservoirMathModelChromosome chromosome)
        //    {
        //        _model = chromosome._model;
        //        Length = _model._period * 2;
        //        _sequence = chromosome._sequence.Select(x => x).ToList();
        //        _inputSchedule = chromosome._inputSchedule.Select(x => x.Select(y => y).ToArray()).ToList();
        //        _outputSchedule = chromosome._outputSchedule.Select(x => x.Select(y => y).ToArray()).ToList();
        //        Fitness = chromosome.Fitness;
        //    }

        //    public IChromosome Clone()
        //    {
        //        return new ReservoirMathModelChromosome(this);
        //    }

        //    public int CompareTo(IChromosome other)
        //    {
        //        if (other == null)
        //        {
        //            return -1;
        //        }

        //        var otherFitness = other.Fitness;

        //        if (Fitness == otherFitness)
        //        {
        //            return 0;
        //        }
                
        //        return Fitness > otherFitness ? 1 : -1;
        //    }

        //    public IChromosome CreateNew()
        //    {
        //        return new ReservoirMathModelChromosome(_model);
        //    }
            
        //    private void Initialize()
        //    {
        //        _sequence = new List<int>();

        //        Func<List<int>, List<int>> initializer = (initialCount) =>
        //        {
        //            List<int> result = new List<int>();
        //            if (initialCount == null)
        //            {
        //                for (int i = 0; i < _model._period; i++)
        //                    result.Add(0);
        //            }
        //            else
        //            {
        //                for (int i = 0; i < initialCount.Count(); i++)
        //                {
        //                    int currentCount = initialCount[i];
        //                    if (currentCount == 1)
        //                        result.Add(0);
        //                    else
        //                        while (currentCount > 0)
        //                        {
        //                            result.Add(AlgorithmHelper.RandomInstance.Next(0, currentCount));
        //                            currentCount--;
        //                        }
        //                }
        //            }
        //            return result;
        //        };
                
        //        _sequence.AddRange(initializer(_model._tempInputCounts));
        //        _sequence.AddRange(initializer(_model._tempOutputCounts));
        //    }

        //    public void Decode()
        //    {
        //        int period = _model._period;
        //        _inputSchedule = new List<double[]>();
        //        _outputSchedule = new List<double[]>();

        //        Func<List<int>, List<Tuple<List<double[]>, List<int>>>, List<double[]>> decoder = (seq, initialSolution) =>
        //        {
        //            var result = new double[_model._period][].ToList();
        //            if (initialSolution == null)
        //            {
        //                return AlgorithmHelper.CreateListOfArrays(_model._period, 1, 0.0);
        //            }
        //            int seqIdx = 0;
        //            foreach (var el in initialSolution)
        //            {
        //                var indexes = el.Item2;
        //                var values = el.Item1.Select(x => x).ToList();
        //                for (int i = 0; i < indexes.Count(); i++)
        //                {
        //                    var val = values[seq[seqIdx]];
        //                    values.RemoveAt(seq[seqIdx]);
        //                    result[indexes[i]] = val;
        //                    seqIdx++;
        //                }
        //            }
        //            return result;
        //        };

        //        _inputSchedule = decoder(_sequence.GetRange(0, _model._period).ToList(), _model._tempInputSolution);
        //        _outputSchedule = decoder(_sequence.GetRange(_model._period, _model._period).ToList(), _model._tempOutputSolution);
        //    }

        //    public Gene GenerateGene(int geneIndex)
        //    {
        //        throw new NotImplementedException();
        //    }

        //    public void ReplaceGene(int index, Gene gene)
        //    {
        //        throw new NotImplementedException();
        //    }

        //    public void ReplaceGenes(int startIndex, Gene[] genes)
        //    {
        //        throw new NotImplementedException();
        //    }

        //    public void Resize(int newLength)
        //    {
        //        throw new NotImplementedException();
        //    }

        //    public Gene GetGene(int index)
        //    {
        //        throw new NotImplementedException();
        //    }

        //    public Gene[] GetGenes()
        //    {
        //        return _sequence.Select(x => new Gene(x)).ToArray();
        //    }
        //}

        //public class ReservoirMathModelMutation : IMutation
        //{
        //    public bool IsOrdered => true;

        //    public void Mutate(IChromosome chromosome, float probability)
        //    {
        //        var c = chromosome as ReservoirMathModelChromosome;

        //        Func<List<int>, List<int>, List<int>> mutator = (seq, initialCount) =>
        //        {
        //            var result = seq.Select(x => x).ToList();
        //            if (initialCount == null)
        //            {
        //                return seq;
        //            }
        //            else
        //            {
        //                int counter = 0;
        //                for (int i = 0; i < initialCount.Count(); i++)
        //                {
        //                    int currentCount = initialCount[i];
        //                    if (currentCount == 1)
        //                    {
        //                        counter++;
        //                    }
        //                    else
        //                        while (currentCount > 0)
        //                        {
        //                            if (AlgorithmHelper.RandomInstance.NextDouble() < probability)
        //                            {
        //                                seq[counter] = AlgorithmHelper.RandomInstance.Next(0, currentCount);
        //                            }
        //                            currentCount--;
        //                            counter++;
        //                        }
        //                }
        //                return result;
        //            }
        //        };

        //        var newSeq = new List<int>();
        //        int period = c._model._period;
        //        newSeq.AddRange(mutator(c._sequence.GetRange(0, period).ToList(), c._model._tempInputCounts));
        //        newSeq.AddRange(mutator(c._sequence.GetRange(period, period).ToList(), c._model._tempOutputCounts));
        //        c._sequence = newSeq;
        //        c.Decode();
        //    }
        //}

        //public class ReservoirMathModelFitness : IFitness
        //{
        //    public double Evaluate(IChromosome chromosome)
        //    {
        //        ReservoirMathModelChromosome c = chromosome as ReservoirMathModelChromosome;

        //        var reservoirSchedule = c._model.GetReservoirSchedule(c._inputSchedule, c._outputSchedule);

        //        var overvolume = reservoirSchedule.Sum(x =>
        //        {
        //            if (x < 0)
        //                return -x;
        //            else if (x > c._model._reservoirVolume)
        //                return c._model._reservoirVolume - x;
        //            else
        //                return 0.0;
        //        });

        //        var changesCount = 0.0;
        //        Func<List<double[]>, int> changesCounter = (schedule) =>
        //        {
        //            int result = 0;
        //            for (int i = 1; i < schedule.Count(); i++)
        //            {
        //                if (schedule[i][0] != schedule[i - 1][0])
        //                    result++;
        //            }
        //            return result;
        //        };
        //        changesCount += changesCounter(c._inputSchedule) + changesCounter(c._outputSchedule);

        //        return 1 / (changesCount + 0.01) /*+ overvolume*/;
        //    }
        //}

        //public class ReservoirMathModelCrossOver : ICrossover
        //{
        //    public int ParentsNumber => 2;

        //    public int ChildrenNumber => 2;

        //    public int MinChromosomeLength => 1488;

        //    public bool IsOrdered => true;

        //    public IList<IChromosome> Cross(IList<IChromosome> parents)
        //    {
        //        var p1 = parents[0] as ReservoirMathModelChromosome;
        //        var p2 = parents[1] as ReservoirMathModelChromosome;

        //        var model = p1._model;

        //        var c1 = new ReservoirMathModelChromosome(p1);
        //        var c2 = new ReservoirMathModelChromosome(p2);

        //        int pointIdx = AlgorithmHelper.RandomInstance.Next(1, model._period * 2 - 1);

        //        for (int i = 0; i < pointIdx; i++)
        //        {
        //            c1._sequence[i] = p1._sequence[i];
        //            c2._sequence[i] = p2._sequence[i];
        //        }
        //        for (int i = pointIdx; i < model._period * 2; i++)
        //        {
        //            c1._sequence[i] = p2._sequence[i];
        //            c2._sequence[i] = p1._sequence[i];
        //        }
        //        c1.Decode();
        //        c2.Decode();
        //        return new List<IChromosome>() { c1, c2 };
        //    }
        //}
        
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
            /*if (tempSolutions == null)
                throw new Exception();
            if (tempSolutions.Count() != _sections.Count())
                throw new Exception();
            for (int i = 0; i < tempSolutions.Count(); i++)
            {
                var tempSolution = tempSolutions[i];
                if (tempSolution == null && _sections[i] != null)
                    throw new Exception();
                if (tempSolution == null)
                    continue;
                if (tempSolution.Count() == 0)
                    throw new Exception();
                if (tempSolution.Any(x => x == null))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item1 == null))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item1.Count() == 0))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item1.Any(y => y == null)))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item1.Any(y => y.Count() != _sections[i].Dimension)))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item1.Any(y => y.Any(z => z < 0))))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item2 == null))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item2.Count() != x.Item1.Count()))
                    throw new Exception();
                if (tempSolution.Any(x => x.Item2.Any(y => y < 0 || y > _period - 1)))
                    throw new Exception();
                if (tempSolution.SelectMany(x => x.Item2).GroupBy(x => x).Any(x => x.Count() > 1))
                        throw new Exception();
            }*/

            _tempPumpSchedules = tempPumpSchedules.Select(x => x.ToList()).ToList();
            _tempTargetVolumes = tempTargetVolumes.Select(x => x == null ? null : new TargetVolumes(x)).ToList();
            _tempSolutions = null;
        }
        
        private static List<double[]> CreateInitialSchedule(List<Tuple<List<double[]>, List<int>>> initialSolution)
        {
            int period = initialSolution.Sum(x => x.Item2.Count());
            var result = (new double[period][]).ToList();
            foreach(var tuple in initialSolution)
            {
                tuple.Item2.Sort();
                tuple.Item1.Sort((el1, el2) => el1[0] > el2[0] ? 1 : (el1[0] < el2[0] ? -1 : 0));
                for(int i = 0; i < tuple.Item1.Count(); i++)
                {
                    result[tuple.Item2[i]] = tuple.Item1[i];
                }
            }
            return result;
        }

        private static List<double[]> SwapBlocks(List<double[]> scheudule, List<int> block1, List<int> block2)
        {
            int blockSize = block1.Count();
            if (blockSize != block2.Count())
                throw new Exception();

            List<double[]> result = scheudule.Select(x => x).ToList();
            for (int k = 0; k < blockSize; k++)
            {
                result[block1[k]] = scheudule[block2[k]];
                result[block2[k]] = scheudule[block1[k]];
            }
            return result;
        }

        private List<double> GetReservoirSchedule(List<double[]> inputSchedule, List<double[]> outputSchedule, List<double> pumpsSchedule, double startVolume)
        {
            List<double> result = new List<double>() { startVolume };
            for (int i = 0; i < inputSchedule.Count(); i++)
                result.Add(result.Last() + inputSchedule[i].Last() - outputSchedule[i].First() + pumpsSchedule[i]);
            return result;
        }

        private List<int> GetCrashIndexes(List<double> reservoirSchedule, double minVolume, double maxVolume)
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

        private double GetOverfillVolume(List<double> reservoirSchedule, double minVolume, double maxVolume)
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
        
        public List<List<double[]>> Balance(Action<string> informationAction = null)
        {
            if (informationAction == null)
                informationAction = (str) => { };

            List<List<double[]>> initialSchedules = null;
            List<double> prevOverfillVolumes = new List<double>();

            Func<List<List<double[]>>> ConvertResult = () =>
            {
                if (_sections[0] == null)
                    initialSchedules.RemoveAt(0);
                else if (_sections.Last() == null)
                    initialSchedules.RemoveAt(initialSchedules.Count() - 1);
                return initialSchedules;
            };
            Action CalcInitialSchedules = () =>
            {
                initialSchedules = _tempSolutions.Select((x, i) =>
                {
                    if (_sections[i] == null)
                        return AlgorithmHelper.CreateListOfArrays(_period, 1, 0.0);
                    else
                        return _sections[i].GetFullSchedule(CreateInitialSchedule(x));
                }).ToList();
            };
            Action CalcOverfill = () =>
            {
                prevOverfillVolumes = new List<double>();
                for (int i = 0; i < initialSchedules.Count() - 1; i++)
                {
                    var reservoirSchedule = GetReservoirSchedule(initialSchedules[i], initialSchedules[i + 1], _tempPumpSchedules[i], _oilStartVolumes[i]);
                    prevOverfillVolumes.Add(GetOverfillVolume(reservoirSchedule, 0, _reservoirVolumes[i]));
                }
            };

            // Сначала без анализа
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
                initialSchedules = Permute(initialSchedules, 24, 100, informationAction);
                CalcOverfill();
                if (prevOverfillVolumes.All(x => x == 0.0))
                {
                    return ConvertResult();
                }
                bool breaker = !Analyse(initialSchedules);
                if (breaker)
                    break;
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
            Func<int, TempBalanceStruct> GetInitStruct = (sectionNumber) =>
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
            };
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

            for (int i = 0; i < _sections.Count(); i++)
            {
                var section = _sections[i];
                if (section == null)
                    continue;

                var badStartRepairs = new List<List<int>>();
                var badEndRepairs = new List<List<int>>();
                for (int j = 0; j < section.RepairsIntervals.Count(); j++)
                {
                    var repair = section.RepairsIntervals[j];
                    int repairLen = repair.Count();
                    int test = 0;
                    if (i > 0)
                    {
                        // Смотрим изменения во входном резервуаре на данном участке
                        double startVolume = reservoirInitialSchedules[i - 1][repair.First()],
                            endVolume = reservoirInitialSchedules[i - 1][repair.Last() + 1];
                        // Неконтролируемый скачек больше объема резервуара
                        if (_sections[i - 1] == null && (endVolume - startVolume) > _reservoirVolumes[i - 1])
                            return false;

                        if (startVolume < 0 && endVolume > _reservoirVolumes[i - 1])
                        {
                            // Нужно на время ремонта уменьшить расход на предыдущей секции
                            test = 0;
                        }
                        else if (startVolume > 0 && endVolume > _reservoirVolumes[i - 1])
                        {
                            // Нужно немного опорожнить резервуар
                            test = 0;
                        }
                        else if (startVolume < 0 && endVolume < _reservoirVolumes[i - 1])
                        {
                            // Это странно....
                            // Нужно немного наполнить резервуар
                            test = 0;
                        }
                        else
                        {
                            // Этот ремонт не влияет ни на что
                        }
                    }

                    if (i < _resCount)
                    {
                        // Смотрим изменения в выходном резервуаре на данном участке
                        double startVolume = reservoirInitialSchedules[i][repair.First()],
                            endVolume = reservoirInitialSchedules[i][repair.Last() + 1];
                        // Неконтролируемый скачек больше объема резервуара
                        if (_sections[i + 1] == null && (startVolume - endVolume) > _reservoirVolumes[i])
                            return false;

                        if (startVolume > _reservoirVolumes[i] && endVolume < 0)
                        {
                            // Нужно на время ремонта уменьшить расход на следующей секции
                            test = 0;
                        }
                        else if (startVolume < _reservoirVolumes[i] && endVolume < 0)
                        {
                            // Нужно немного наполнить резервуар
                            test = 0;
                        }
                        else if (startVolume > _reservoirVolumes[i] && endVolume > 0)
                        {
                            // Это странно....
                            // Нужно немного опорожнить резервуар
                            test = 0;
                        }
                        else
                        {
                            // Этот ремонт не влияет ни на что
                        }

                    }
                }
            }

            CreateTempSolution();
            return true;
        }

        #endregion
    }
}
