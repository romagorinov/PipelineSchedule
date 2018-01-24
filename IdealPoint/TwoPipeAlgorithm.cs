using GeneticSharp.Domain.Chromosomes;
using GeneticSharp.Domain.Fitnesses;
using GeneticSharp.Domain.Mutations;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class TwoPipeAlgorithm
    {
        #region Параметры

        private int _period;

        private List<double> _regimesIn;
        private List<double> _maxFlowIn;
        private double _maxVolumeIn;
        private SinglePipeAlgorithm _inAlgo;

        private List<double> _regimesOut;
        private List<double> _maxFlowOut;
        private double _maxVolumeOut;
        private SinglePipeAlgorithm _outAlgo;

        private double _bufferVolume;

        private class RegimeCombination
        {
            public double InVolume
            {
                get;
                private set;
            }

            public double OutVolume
            {
                get;
                private set;
            }

            public double Difference
            {
                get;
                private set;
            }

            public RegimeCombination(double inVolume, double outVolume)
            {
                InVolume = inVolume;
                OutVolume = outVolume;
                Difference = inVolume - outVolume;
            }
        }

        private List<RegimeCombination> _regimeCombinations;
        private List<Tuple<RegimeCombination, RegimeCombination, double>> _differences;
        
        private List<List<RegimeCombination>> _avaliableRegimeCombinations;

        #endregion

        #region Настройки алгоритма

        public static int GENERATION_NUMBER = 300;
        public static int MINIMUM_WORK_TIME = 1;
        public static int NUMBER_OF_TRIES = 90;

        #endregion

        #region Конструкторы

        public TwoPipeAlgorithm(List<double> regimesIn, List<double> maxFlowIn, List<double> regimesOut, List<double> maxFlowOut, double bufferVolume)
        {
            _regimesIn = regimesIn ?? throw new ArgumentNullException();
            _regimesIn.Sort();
            _maxFlowIn = maxFlowIn ?? throw new ArgumentNullException();

            _period = maxFlowIn.Count();
            _maxVolumeIn = _maxFlowIn.Sum(x => _regimesIn.Last(y => y <= x));

            _regimesOut = regimesOut ?? throw new ArgumentNullException();
            _regimesOut.Sort();
            _maxFlowOut = maxFlowOut ?? throw new ArgumentNullException();
            _maxVolumeOut = _maxFlowOut.Sum(x => _regimesOut.Last(y => y <= x));

            _bufferVolume = bufferVolume;

            _regimeCombinations = AlgorithmHelper.CartesianProduct(new List<List<double>>() { _regimesIn, _regimesOut }).Select(x => new RegimeCombination(x.First(), x.Last())).ToList();

            _avaliableRegimeCombinations = _maxFlowIn.Zip(_maxFlowOut, (mf1, mf2) => _regimeCombinations.Where(r => r.InVolume <= mf1 && r.OutVolume <= mf2).ToList()).ToList();

            _differences = new List<Tuple<RegimeCombination, RegimeCombination, double>>();

            foreach(var r1 in _regimeCombinations)
            {
                foreach(var r2 in _regimeCombinations)
                {
                    _differences.Add(new Tuple<RegimeCombination, RegimeCombination, double>(r1, r2, 
                        (Math.Abs(r1.InVolume - r2.InVolume) / _regimesIn.Max()  + Math.Abs(r1.OutVolume - r2.OutVolume) / _regimesOut.Max()) / 2));
                }
            }

            _inAlgo = new SinglePipeAlgorithm(_regimesIn, _maxFlowIn);
            _outAlgo = new SinglePipeAlgorithm(_regimesOut, _maxFlowOut);
        }

        #endregion

        #region Методы
        
        public enum LEADING
        {
            IN_LEADING,
            OUT_LEADING,
            NO_LEADING
        }

        public struct AlgorithmParametes
        {
            public double inTarget;
            public double outTarget;

            public double minLevel;
            public double maxLevel;
            public double startLevel;

            public List<double> pumpsSchedule;

            public List<double> inFixValues;
            public List<double> outFixValues;

            public void Check(TwoPipeAlgorithm algo)
            {
                if (inTarget < 0 || outTarget < 0) throw new ArgumentException();
                if (minLevel >= maxLevel || minLevel >= algo._bufferVolume || maxLevel <= 0 || startLevel > maxLevel || startLevel < minLevel) throw new ArgumentException();

                if (minLevel < 0) minLevel = 0;
                if (maxLevel > algo._bufferVolume) maxLevel = algo._bufferVolume;

                if (pumpsSchedule == null) pumpsSchedule = AlgorithmHelper.CreateListOfElements(algo._period, 0.0);
                if (inFixValues == null) inFixValues = AlgorithmHelper.CreateListOfElements(algo._period, -1.0);
                if (outFixValues == null) outFixValues = AlgorithmHelper.CreateListOfElements(algo._period, -1.0);

                if (pumpsSchedule.Count() != algo._period || inFixValues.Count() != algo._period || outFixValues.Count() != algo._period) throw new ArgumentException();
            }
        }
/*
        public Tuple<List<double>, List<double>> BalanceSchedules(AlgorithmParametes parameters)
        {
            parameters.Check(this);

            Random rnd = AlgorithmHelper.RandomInstance;

            List<double> inSchedule = _inAlgo.FormSchedule(parameters.inTarget, 0, parameters.inFixValues),
                outSchedule = _outAlgo.FormSchedule(parameters.outTarget, 0, parameters.outFixValues);

            List<double> reservoirSchedule = AlgorithmHelper.GetIntegralDifference(parameters.startLevel, inSchedule, outSchedule, parameters.pumpsSchedule);
            
            if (reservoirSchedule.Min() >= parameters.minLevel && reservoirSchedule.Max() <= parameters.maxLevel)
            {
                return new Tuple<List<double>, List<double>>(inSchedule, outSchedule);
            }

            List<List<RegimeCombination>> avaliableRegimeCombinations = _avaliableRegimeCombinations.Select((combinations, i) =>
            {
                double inFix = parameters.inFixValues[i], outFix = parameters.outFixValues[i];
                if (inFix < 0 && outFix < 0)
                {
                    return combinations;
                }
                else if (inFix >= 0 && outFix < 0)
                {
                    return combinations.Where(x => x.InVolume == inFix).ToList();
                }
                else if (inFix < 0 && outFix >= 0)
                {
                    return combinations.Where(x => x.OutVolume == outFix).ToList();
                }
                else
                {
                    return combinations.Where(x => x.InVolume == inFix && x.OutVolume == outFix).ToList();
                }
            }).ToList();
            
            ConcurrentBag<Tuple<List<double>, List<double>>> results = new ConcurrentBag<Tuple<List<double>, List<double>>>();
            Parallel.For(0, GENERATION_NUMBER, (iteration) =>
            {
                List<double> inLocalSchedule = inSchedule.Select(x => x).ToList(),
                    outLocalSchedule = outSchedule.Select(x => x).ToList();
                List<double> reservoirLocalSchedule;
                int iter = 0;
                while (iter++ < NUMBER_OF_TRIES)
                {
                    #region Поиск критической точки

                    reservoirLocalSchedule = AlgorithmHelper.GetIntegralDifference(parameters.startLevel, inLocalSchedule, outLocalSchedule, parameters.pumpsSchedule);
                    int crashIndex = reservoirLocalSchedule.FindIndex(x => x < parameters.minLevel || x > parameters.maxLevel);
                    if (crashIndex == -1)
                    {
                        break;
                    }
                    double sign = reservoirLocalSchedule[crashIndex] > parameters.maxLevel ? -1 : 1;
                    int crashInterval = crashIndex - 1;
                    double beforeCrashInVolume = AlgorithmHelper.GetSumOnInterval(inLocalSchedule, 0, crashInterval);

                    #endregion

                    #region Поиск режимов стабилизации

                    // Нужный знак
                    List<RegimeCombination> stabilizationAvaliableRegimes = avaliableRegimeCombinations[crashInterval].Where(x => x.Difference * sign >= 0).ToList();
                    // Работа не меньше, чем минимум
                    stabilizationAvaliableRegimes.Where(regime =>
                    {
                        double startVol = reservoirLocalSchedule[crashInterval]
                    });
                    // Если режимы стабилизации не найдены, то решение отсутствует
                    if (stabilizationAvaliableRegimes.Count() == 0)
                    {
                        break;
                    }

                    List<RegimeCombination> nonZeroRegimes = stabilizationAvaliableRegimes.Where(x => x.InVolume != 0 && x.OutVolume != 0).ToList();
                    if (nonZeroRegimes.Count() != 0)
                    {
                        stabilizationAvaliableRegimes = nonZeroRegimes;
                    }

                    #endregion

                    #region Выбор режимов стабилизации

                    // Чем ближе новая комбинация к текущей, тем лучше
                    RegimeCombination currentRegimeCombination = _regimeCombinations.First(x => x.InVolume == inSchedule[crashIndex] && x.OutVolume == outSchedule[crashIndex]);
                    List<double> probability = stabilizationAvaliableRegimes.Select(x => 1 - _differences.First(y => y.Item1 == currentRegimeCombination && y.Item2 == x).Item3).ToList();
                    probability = probability.Select(x => x / probability.Sum()).ToList();
                    RegimeCombination selectedRegimeCombination = stabilizationAvaliableRegimes[AlgorithmHelper.GetRouletIndex(probability)];

                    #endregion

                    #region Выбираем время

                    // Максимальный предел времени меньше crashIndex + maxStabilizationTime
                    int maxStabilizationTime = 1;
                    double reservoirVolume = reservoirSchedule[crashIndex] + selectedRegimeCombination.Difference;
                    for (int i = crashIndex + 1; i < period; i++)
                    {
                        reservoirVolume += selectedRegimeCombination.Difference;
                        if (_avaliableRegimeCombinations[i].IndexOf(selectedRegimeCombination) == -1 || reservoirVolume > maxLevel || reservoirVolume < minLevel)
                        {
                            break;
                        }
                        maxStabilizationTime++;
                    }

                    int stabilizationTime = rnd.Next(1, maxStabilizationTime + 1);

                    #endregion

                    #region Пересчитываем

                    for (int i = crashIndex; i < crashIndex + stabilizationTime; i++)
                    {
                        inSchedule[i] = selectedRegimeCombination.InVolume;
                        outSchedule[i] = selectedRegimeCombination.OutVolume;
                    }

                    if (crashIndex + stabilizationTime != period)
                    {
                        if (leading != LEADING.IN_LEADING)
                        {
                            inSchedule = new DiscreteSchedule(_inAlgo.FormScheduleFromPoint(inVolume, inSchedule.ToList(), crashIndex + stabilizationTime).ToArray());
                        }
                        if (leading != LEADING.OUT_LEADING)
                        {
                            outSchedule = new DiscreteSchedule(_outAlgo.FormScheduleFromPoint(outVolume, outSchedule.ToList(), crashIndex + stabilizationTime).ToArray());
                        }
                    }

                    #endregion
                }

                if (DiscreteSchedule.GetDiffereceIntegral(startLevel, inSchedule, outSchedule).FirstIndexOfOutOfRange(minLevel, maxLevel) == -1)
                {
                    res.Add(new DiscreteSchedule[] { inSchedule, outSchedule });
                }
            });
            results = res.ToList();

            double maxInDifference = inVolume > _inAlgo.MaxVolume / 2 ? inVolume : _inAlgo.MaxVolume - inVolume,
                maxOutDifference = outVolume > _outAlgo.MaxVolume / 2 ? outVolume : _outAlgo.MaxVolume - outVolume;
            List<double> volumeCriteria = results.Select((x) =>
            {
                return (Math.Abs(x[0].Sum - inVolume) / maxInDifference + Math.Abs(x[1].Sum - outVolume) / maxOutDifference) / 2;
            }).ToList();

            int maxChanges = period / 24 * 3;
            List<double> changesCriteria = results.Select((x) =>
            {
                double changes = 0;
                for (int i = 1; i < period; i++)
                {
                    if (x[0][i] != x[0][i - 1])
                        changes++;
                    if (x[1][i] != x[1][i - 1])
                        changes++;
                }
                changes /= maxChanges * 2;
                return changes;
            }).ToList();

            List<double> sumCriteria = volumeCriteria.Zip(changesCriteria, (x, y) => x + 0.2 * y).ToList();

            return results[sumCriteria.IndexOf(sumCriteria.Min())].Select(x => x.ToList()).ToArray();
        }
     */   
        #endregion

        #region Кастомизация GeneticSharp (пока не нужно)

        //public DiscreteSchedule[] GeneticBalanceSchedules(double minLevel, double maxLevel, double startLevel, DiscreteSchedule inBestSchedule, DiscreteSchedule outBestSchedule)
        //{
        //    DiscreteSchedule reservoirSchedule = DiscreteSchedule.GetDiffereceIntegral(startLevel, inBestSchedule, outBestSchedule);

        //    if (reservoirSchedule.Min >= minLevel && reservoirSchedule.Max <= maxLevel)
        //    {
        //        return new DiscreteSchedule[] { new DiscreteSchedule(inBestSchedule), new DiscreteSchedule(outBestSchedule) };
        //    }

        //    var selection = new GeneticSharp.Domain.Selections.EliteSelection();
        //    var crossover = new GeneticSharp.Domain.Crossovers.OnePointCrossover();
        //    var mutation = new TwoPipeMutation();
        //    var fitness = new TwoPipeFitness(
        //        inBestSchedule.Sum > _maxVolumeIn / 2 ? inBestSchedule.Sum : _maxVolumeIn - inBestSchedule.Sum,
        //        outBestSchedule.Sum > _maxVolumeOut / 2 ? outBestSchedule.Sum : _maxVolumeOut - outBestSchedule.Sum,
        //        startLevel, 0, maxLevel);
        //    var chromosome = new TwoPipeChromosome(
        //        _maxFlowIn.Select(x => _regimesIn.Where(y => y <= x).ToList()).ToList(),
        //        _maxFlowOut.Select(x => _regimesOut.Where(y => y <= x).ToList()).ToList(),
        //        inBestSchedule,
        //        outBestSchedule);
        //    var population = new GeneticSharp.Domain.Populations.Population(50, 100, chromosome);
        //    var ga = new GeneticSharp.Domain.GeneticAlgorithm(population, fitness, selection, crossover, mutation);
        //    ga.Termination = new GeneticSharp.Domain.Terminations.GenerationNumberTermination(1000);
        //    ga.MutationProbability = 0.3f;
        //    ga.Start();

        //    return (ga.BestChromosome as TwoPipeChromosome).GetSchedule();
        //}

        //public class TwoPipeFitness : IFitness
        //{
        //    private double _biggestInDifference;
        //    private double _biggestOutDifference;
        //    private double _startReservuirVolume;
        //    private double _minReservoirVolume;
        //    private double _maxReservoirVolume;

        //    public TwoPipeFitness(double biggestInDifference, double biggestOutDifference, double startReservuirVolume, double minReservoirVolume, double maxReservoirVolume)
        //    {
        //        _biggestInDifference = biggestInDifference;
        //        _biggestOutDifference = biggestOutDifference;
        //        _startReservuirVolume = startReservuirVolume;
        //        _minReservoirVolume = minReservoirVolume;
        //        _maxReservoirVolume = maxReservoirVolume;
        //    }

        //    public double Evaluate(IChromosome chromosome)
        //    {
        //        TwoPipeChromosome c = chromosome as TwoPipeChromosome;

        //        DiscreteSchedule[] schedules = c.GetSchedule();
        //        DiscreteSchedule inDifference = c.scheduleIn - schedules[0];
        //        DiscreteSchedule outDifference = c.scheduleOut - schedules[1];
                
        //        // Отличие от основных расписаний должно быть минимальным
        //        double differenceCriteria = - (inDifference.ToArray().Sum(x => Math.Abs(x)) / _biggestInDifference + outDifference.ToArray().Sum(x => Math.Abs(x)) / _biggestOutDifference) / 2;

        //        // Суммарный объем должен быть тем же
        //        double volumeCriteria = - (Math.Abs(inDifference.Sum) / _biggestInDifference + Math.Abs(outDifference.Sum) / _biggestOutDifference) / 2;

        //        // Стандартные штрафы за нули
        //        double zeroPenalty = - inDifference.ToArray().Zip(outDifference.ToArray(), (x, y) => (x == 0 ? 1 : 0) + (y == 0 ? 1 : 0)).Sum();

        //        // Большие пенальти за перелив!
        //        List<double> reservoir = DiscreteSchedule.GetDiffereceIntegral(_startReservuirVolume, schedules[0], schedules[1]).ToArray().ToList();
        //        double overvolumeCriteria = - reservoir.Sum((x) =>
        //        {
        //            if (x > _maxReservoirVolume)
        //                return 5;
        //            else if (x < _minReservoirVolume)
        //                return 5;

        //            return 0;
        //        });

        //        return volumeCriteria + differenceCriteria + zeroPenalty + overvolumeCriteria;
        //    }
        //}

        //public class TwoPipeChromosome : ChromosomeBase
        //{
        //    public List<List<double>> avaliableRegimesIn;
        //    public DiscreteSchedule scheduleIn;
        //    public List<List<double>> avaliableRegimesOut;
        //    public DiscreteSchedule scheduleOut;

        //    private static Random _rnd = new Random();

        //    public TwoPipeChromosome(List<List<double>> avaliableRegimesIn, List<List<double>> avaliableRegimesOut, DiscreteSchedule scheduleIn, DiscreteSchedule scheduleOut) : base(avaliableRegimesIn.Count())
        //    {
        //        this.avaliableRegimesIn = avaliableRegimesIn;
        //        this.avaliableRegimesOut = avaliableRegimesOut;
        //        this.scheduleIn = scheduleIn;
        //        this.scheduleOut = scheduleOut;
        //        CreateGenes();
        //    }

        //    public override Gene GenerateGene(int geneIndex)
        //    {
        //        double genIn = avaliableRegimesIn[geneIndex][AlgorithmHelper.GetRouletIndexCloserBetter(avaliableRegimesIn[geneIndex], scheduleIn[geneIndex])];
        //        double genOut = avaliableRegimesOut[geneIndex][AlgorithmHelper.GetRouletIndexCloserBetter(avaliableRegimesOut[geneIndex], scheduleOut[geneIndex])];
        //        return new Gene(new double[2] { genIn, genOut });
        //    }

        //    public override IChromosome CreateNew()
        //    {
        //        return new TwoPipeChromosome(avaliableRegimesIn, avaliableRegimesOut, scheduleIn, scheduleOut);
        //    }

        //    public DiscreteSchedule[] GetSchedule()
        //    {
        //        return new DiscreteSchedule[2] {
        //            new DiscreteSchedule(GetGenes().Select(x => (x.Value as double[])[0]).ToArray()),
        //            new DiscreteSchedule(GetGenes().Select(x => (x.Value as double[])[1]).ToArray())
        //        };
        //    }
        //}

        //public class TwoPipeMutation : MutationBase
        //{
        //    private static Random _rnd = new Random();

        //    protected override void PerformMutate(IChromosome chromosome, float probability)
        //    {
        //        if (_rnd.NextDouble() > probability)
        //        {
        //            return;
        //        }

        //        TwoPipeChromosome c = chromosome as TwoPipeChromosome;

        //        int startIndexIn = _rnd.Next(0, c.Length - 1), endIndexIn = _rnd.Next(startIndexIn + 1, c.Length),
        //            startIndexOut = _rnd.Next(0, c.Length - 1), endIndexOut = _rnd.Next(startIndexOut + 1, c.Length);
        //        double regimeIn = -1, regimeOut = -1;
        //        for (int i = startIndexIn; i <= endIndexIn; i++)
        //        {
        //            if (c.avaliableRegimesIn[i].IndexOf(regimeIn) == -1)
        //            {
        //                int maxIndex = c.avaliableRegimesIn[i].Count();
        //                regimeIn = c.avaliableRegimesIn[i][_rnd.Next(maxIndex)];
        //            }

        //            c.ReplaceGene(i, new Gene(new double[] { regimeIn, (c.GetGene(i).Value as double[])[1]}));
        //        }

        //        for (int i = startIndexOut; i <= endIndexOut; i++)
        //        {
        //            if (c.avaliableRegimesOut[i].IndexOf(regimeOut) == -1)
        //            {
        //                int maxIndex = c.avaliableRegimesOut[i].Count();
        //                regimeOut = c.avaliableRegimesOut[i][_rnd.Next(maxIndex)];
        //            }

        //            c.ReplaceGene(i, new Gene(new double[] { (c.GetGene(i).Value as double[])[0], regimeOut }));
        //        }
        //    }
        //}

        #endregion
    }
}
