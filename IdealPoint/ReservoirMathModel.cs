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
    public class ReservoirBalancerMathModel
    {
        #region Поля

        double _reservoirVolume;
        double _oilStartVolume;
        ISection _inputSection;
        ISection _outputSection;
        int _period;

        List<Tuple<List<double[]>, List<int>>> _tempInputSolution;
        List<Tuple<List<double[]>, List<int>>> _tempOutputSolution;
        List<double> _tempPumpSchedule;
        double _tempReservoirMaxVolume;
        double _tempReservoirMinVolume;

        #endregion

        #region Свойства



        #endregion

        #region Кострукторы

        public ReservoirBalancerMathModel(double reservoirVolume, double oilStartVolume, ISection inputSection, ISection outputSection, int period)
        {
            _reservoirVolume = reservoirVolume;
            _inputSection = inputSection;
            _outputSection = outputSection;
            _oilStartVolume = oilStartVolume;
            _period = period;
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

        public bool SetTempParams(List<Tuple<List<double[]>, List<int>>> tempInputSolution, List<Tuple<List<double[]>, List<int>>> tempOutputSolution, List<double> tempPumpSchedule,
            double tempReservoirMinVolume, double tempReservoirMaxVolume)
        {
            Func<List<Tuple<List<double[]>, List<int>>>, List<Tuple<List<double[]>, List<int>>>> converter = (solution) =>
            {
                var result = solution
                    .Select(x => new Tuple<List<double[]>, List<int>>(x.Item1.Select(y => y).ToList(), x.Item2.Select(y => y).ToList()))
                    .ToList();
                result.ForEach(x =>
                {
                    x.Item2.Sort();
                    x.Item1.Sort((el1, el2) => el1[0] > el2[0] ? 1 : (el1[0] < el2[0] ? -1 : 0));
                });
                return result;
            };

            if (_inputSection != null)
            {
                if (tempInputSolution == null)
                    throw new Exception();

                _tempInputSolution = converter(tempInputSolution);
            }

            if (_outputSection != null)
            {
                if (tempOutputSolution == null)
                    throw new Exception();

                _tempOutputSolution = converter(tempOutputSolution);
            }
            
            _tempPumpSchedule = tempPumpSchedule == null ?  AlgorithmHelper.CreateListOfElements(_period, 0.0) : tempPumpSchedule;
            _tempReservoirMaxVolume = tempReservoirMaxVolume;
            _tempReservoirMinVolume = tempReservoirMinVolume;

            return true;
        }

        List<double> GetReservoirSchedule(List<double[]> inputSchedule, List<double[]> outputSchedule)
        {
            List<double> result = new List<double>(_period) { _oilStartVolume };
            
            if (_inputSection == null)
                inputSchedule = AlgorithmHelper.CreateListOfArrays(_period, 1, 0.0);
            else if (_outputSection == null)
                outputSchedule = AlgorithmHelper.CreateListOfArrays(_period, 1, 0.0);

            for (int i = 0; i < _period; i++)
            {
                double inputVolume = inputSchedule[i].Last(),
                    outVolume = outputSchedule[i].First();
                result.Add(result.Last() + inputVolume - outVolume + _tempPumpSchedule[i]);
            }
            return result;
        }

        List<int> GetCrashIndexes(List<double> reservoirSchedule)
        {
            var result = new List<int>();
            for (int i = 1; i < _period + 1; i++)
            {
                var x = reservoirSchedule[i];
                if (x < 0 || x > _reservoirVolume)
                {
                    result.Add(i);
                }
            }
            return result;
        }
        
        int GetSumChanges(List<double[]> schedule)
        {
            var result = 0;
            for(int i = 1; i < _period; i++)
            {
                if (schedule[i][0] != schedule[i - 1][0])
                    result++;
            }
            return result;
        }

        public static List<double[]> CreateInitialSchedule(List<Tuple<List<double[]>, List<int>>> initialSolution)
        {
            if (initialSolution == null)
                return null;

            int period = initialSolution.Sum(x => x.Item2.Count());
            var result = (new double[period][]).ToList();
            foreach(var tuple in initialSolution)
            {
                for(int i = 0; i < tuple.Item1.Count(); i++)
                {
                    result[tuple.Item2[i]] = tuple.Item1[i];
                }
            }
            return result;
        }

        public Tuple<List<double[]>, List<double[]>> Balance()
        {
            if (_tempInputSolution == null && _tempOutputSolution == null)
                return null;

            var inputInitialSchedule = CreateInitialSchedule(_tempInputSolution);
            if (_inputSection != null)
                _inputSection.GetFullSchedule(inputInitialSchedule);
            var outputInitialSchedule = CreateInitialSchedule(_tempOutputSolution);
            var reservoirInitialSchedule = GetReservoirSchedule(inputInitialSchedule, outputInitialSchedule);
            var crashIndexes = GetCrashIndexes(reservoirInitialSchedule);
            if (crashIndexes.Count() == 0)
                return new Tuple<List<double[]>, List<double[]>>(inputInitialSchedule, outputInitialSchedule);

            int currentSumChanges = GetSumChanges(outputInitialSchedule);
            int blockSize = 24;
            int maxIter = 10;
            while (maxIter-- > 0)
            {
                ConcurrentBag<Tuple<List<int>, List<int>>> variantList = new ConcurrentBag<Tuple<List<int>, List<int>>>();
                foreach (var intervals in _outputSection.ControlAvaliableIntervals)
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
                            if (intervalCount - i > 2 * blockSize && j >= i + 24 && j + 24 < intervalCount)
                            {
                                int i2 = intervals[j + blockSize];
                                if (crashIndexes.Any(x => x <= i2 && x >= i1))
                                    return new int[] { i, j };
                                else
                                    return null;
                            }
                            else
                                return null;
                        })
                        .Where(x => x != null)
                        .ToList();

                    bool hasSolution = false;
                    Parallel.ForEach(checkIndexes, (indexes) =>
                    {
                        var block1 = intervals.GetRange(indexes[0], blockSize).ToList();
                        var block2 = intervals.GetRange(indexes[1], blockSize).ToList();

                        List<double[]> newSchedule = outputInitialSchedule.Select(x => x).ToList();
                        for (int k = 0; k < blockSize; k++)
                        {
                            newSchedule[block1[k]] = outputInitialSchedule[block2[k]];
                            newSchedule[block2[k]] = outputInitialSchedule[block1[k]];
                        }

                        List<double> newReservoirSchedule = GetReservoirSchedule(inputInitialSchedule, newSchedule);
                        List<int> newCrashIndexes = GetCrashIndexes(newReservoirSchedule);
                        if (newCrashIndexes.Count() == 0)
                            hasSolution = true;

                        if (newCrashIndexes.Count() < crashIndexes.Count())
                            if (hasSolution )
                            variantList.Add(new Tuple<List<int>, List<int>>(block1, block2));
                    });
                };

                if (variantList.Count() == 0)
                    return new Tuple<List<double[]>, List<double[]>>(inputInitialSchedule, outputInitialSchedule);

                var selected = variantList.ToList()[AlgorithmHelper.RandomInstance.Next(0, variantList.Count())];
                List<double[]> newSchedule1 = outputInitialSchedule.Select(x => x).ToList();
                for (int k = 0; k < blockSize; k++)
                {
                    newSchedule1[selected.Item1[k]] = outputInitialSchedule[selected.Item2[k]];
                    newSchedule1[selected.Item2[k]] = outputInitialSchedule[selected.Item1[k]];
                }
                outputInitialSchedule = newSchedule1;
                reservoirInitialSchedule = GetReservoirSchedule(inputInitialSchedule, outputInitialSchedule);
                crashIndexes = GetCrashIndexes(reservoirInitialSchedule);
                if (crashIndexes.Count() == 0)
                    return new Tuple<List<double[]>, List<double[]>>(inputInitialSchedule, outputInitialSchedule);
            }

            return new Tuple<List<double[]>, List<double[]>>(inputInitialSchedule, outputInitialSchedule);
        }

        #endregion
    }
}
