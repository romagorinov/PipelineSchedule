using Accord.Math;
using GeneticSharp.Domain.Chromosomes;
using GeneticSharp.Domain.Fitnesses;
using GeneticSharp.Domain.Mutations;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class SinglePipeAlgorithm: ReservoirBalanceAlgorithm.TechnologicalSectionAlgorithm
    {
        #region Параметры

        private List<double> _regimes;
        private List<double> _normRegimes;
        private List<List<double>> _avaliableRegimesOnIntervals;
        private List<double> _upperBounds;
        private double _maxVolume;

        int _period;

        // Сортировка по ремонтным работам
        int[] _sortedIntervals;
        
        public double MaxVolume
        {
            get
            {
                return _maxVolume;
            }
        }   

        #endregion

        #region Конструкторы
        
        public SinglePipeAlgorithm(List<double> regimes, List<double> maxFlows)
        {
            if (regimes == null || maxFlows == null)
                throw new ArgumentNullException();

            _regimes = regimes.Select(x => x).ToList();
            _regimes.Sort();

            _avaliableRegimesOnIntervals = maxFlows.Select(maxFlow => regimes.Where(regime => regime <= maxFlow).ToList()).ToList();

            PreCalculations();
        }

        public SinglePipeAlgorithm(List<double> regimes, List<List<double>> avaliableRegimesOnIntervals)
        {
            if (regimes == null || avaliableRegimesOnIntervals == null)
                throw new ArgumentNullException();

            _regimes = regimes.Select(x => x).ToList();
            _regimes.Sort();

            _avaliableRegimesOnIntervals = avaliableRegimesOnIntervals.Select(x => x.Select(y => y).ToList()).ToList();

            PreCalculations();
        }

        private void PreCalculations()
        {
            for (int i = 0; i < _avaliableRegimesOnIntervals.Count(); i++)
            {
                _avaliableRegimesOnIntervals[i].Sort();
            }

            _upperBounds = _avaliableRegimesOnIntervals.Select(x => x.Last()).ToList();

            _maxVolume = _upperBounds.Sum();

            _period = _avaliableRegimesOnIntervals.Count();

            _sortedIntervals = (new DiscreteSchedule(_upperBounds.ToArray())).SortIndexesOnly();

            _normRegimes = _regimes.Select(x => x / _regimes.Max()).ToList();
        }

        #endregion

        #region Методы

        public List<double> FormSchedule(double volume, int point = 0, List<double> fixValues = null)
        {
            if (fixValues == null)
            {
                fixValues = _avaliableRegimesOnIntervals.Select(x => -1.0).ToList();
            }
            
            double sum = fixValues.Sum(x => x < 0 ? 0 : x);
            List<double> schedule = fixValues.Select(x => x < 0 ? 0 : x).ToList();

            // Начинаем заполнять по уровням
            foreach(var regime in _regimes)
            { 
                // Начинаем от point
                int leftPoint = point, rightPoint = point;

                // Пока не дошли до правой и левой границы
                while (leftPoint >= 0 || rightPoint < _period)
                {
                    // Если правая граница не выходит за пределы, и пропускная позволяет
                    if (rightPoint < _period && _avaliableRegimesOnIntervals[rightPoint].IndexOf(regime) != -1 && fixValues[rightPoint] < 0)
                    {
                        double newSum = sum - schedule[rightPoint] + regime;
                        if (newSum <= volume)
                        {
                            sum = newSum;
                            schedule[rightPoint] = regime;
                        }
                        else
                        {
                            return schedule;
                        }
                    }

                    if (leftPoint >= 0 && _avaliableRegimesOnIntervals[leftPoint].IndexOf(regime) != -1 && fixValues[leftPoint] == -1)
                    {
                        double newSum = sum - schedule[leftPoint] + regime;
                        if (newSum <= volume)
                        {
                            sum = newSum;
                            schedule[leftPoint] = regime;
                        }
                        else
                        {
                            return schedule;
                        }
                    }

                    leftPoint--;
                    rightPoint++;
                }
            }

            return schedule;
        }

        public List<double> FormScheduleFromPoint(double volume, List<double> baseSchedule, int startPoint, int point = 0)
        {
            return FormSchedule(volume, point, baseSchedule.ToArray().Select((x, i) => i < startPoint ? x : -1).ToList());
        }

        public List<double> FormScheduleContinuousRegimeField(double volume, List<double> fixValues = null, int neededPoint = -1)
        {
            if (fixValues == null)
            {
                fixValues = _avaliableRegimesOnIntervals.Select(x => -1.0).ToList();
            }

            if (neededPoint >= 0 && fixValues[neededPoint] >=0 )
            {
                return new List<double>() { fixValues[neededPoint] };
            }

            if (volume > _maxVolume)
            {
                volume = _maxVolume;
            }

            // Не стоит ограничивать, ибо могут качать на каком-нибудь спец. режиме
            //for (int i = 0; i < _period; i++)
            //{
            //    if (_upperBound[i] < fixValues[i])
            //        return null;
            //}

            double fixSum = fixValues.Sum(x => x < 0 ? 0 : x);
            if (fixSum >= volume)
            {
                if (neededPoint >= 0)
                {
                    return new List<double>() { 0 };
                }
                else
                {
                    return fixValues.Select(x => x < 0 ? 0 : x).ToList();
                }
            }

            double leftoverSum = volume - fixSum;
            int leftoverPeriod = _period - fixValues.Count(x => x > 0);
            double leftoverAverageRegime = leftoverSum / leftoverPeriod;
            List<double> schedule = fixValues.Select(x => x).ToList();
            foreach(var idx in _sortedIntervals)
            { 
                if (schedule[idx] < 0)
                {
                    double val = leftoverAverageRegime > _upperBounds[idx] ? _upperBounds[idx] : leftoverAverageRegime;

                    if (idx == neededPoint)
                    {
                        return new List<double>() { val };
                    }

                    schedule[idx] = val;
                    leftoverPeriod --;
                    leftoverSum -= val;
                    leftoverAverageRegime = leftoverSum / leftoverPeriod;
                }
            }

            return schedule;
        }
        
        public List<double> DinamicProgramming(double volume)
        {
            int volumeStepsCount = 1000,
                regimesCount = _regimes.Count(),
                gridSize = (volumeStepsCount + 1) * regimesCount;
            double volumeStepLength = volume / volumeStepsCount;

            double[] functionGrid = new double[gridSize];
            int[] distributionGrid = new int[gridSize];
            for (int i = 0; i < gridSize; i++)
            {
                int regimeIdx = i % regimesCount,
                    volumeIdx = i / regimesCount;
                double vol = (volumeIdx * volumeStepLength) - _regimes[regimeIdx];
                functionGrid[i] = _avaliableRegimesOnIntervals[0].Where(regime => regime <= vol).Max();
                distributionGrid[i] = _avaliableRegimesOnIntervals[0].IndexOf(functionGrid[i]);
            }

            Func<int, int, double> objective = (int interval, int regimeNumber) =>
            {
                return _avaliableRegimesOnIntervals[interval][regimeNumber];
            };

            //f1
            List<double[]> functions = new List<double[]>() { functionGrid };
            List<int[]> distributions = new List<int[]>() { distributionGrid };

            for (int interval = 1; interval < _period; interval++)
            {
                double[] prevFunctionGrid = functions.Last();
                functionGrid = new double[gridSize];
                distributionGrid = new int[gridSize];

                Parallel.For(0, gridSize, (i) =>
                {
                    int regimeIdx = i % regimesCount,
                        volumeIdx = i / regimesCount;
                    double vol = (volumeIdx * volumeStepLength) - _regimes[regimeIdx];
                    double maxValue = double.NegativeInfinity;
                    int maxDistribution = 0;

                    for (int j = 0; j < _regimes.Count(); j++)
                    {
                        if (_avaliableRegimesOnIntervals[interval].IndexOf(_regimes[j]) != -1)
                        {
                            double val = 1; 
                        }
                        else
                        {
                            continue;
                        }
                    }

                    functionGrid[i] = maxValue;
                    distributionGrid[i] = maxDistribution;
                });
                functions.Add(functionGrid);
                distributions.Add(distributionGrid);
            }

            int volumeStep = volumeStepsCount;
            double curRegime = _avaliableRegimesOnIntervals[_period - 1][distributions.Last()[volumeStep]],
                leftOverVolume = volume;
            var schedule = new List<double>() { curRegime };
            for (int i = _period - 2; i >= 0; i--)
            {
                leftOverVolume -= curRegime;
                if (leftOverVolume <= 0)
                {
                    leftOverVolume = 0;
                }
                volumeStep = (int)(leftOverVolume / volumeStepLength);
                curRegime = _avaliableRegimesOnIntervals[i][distributions[i][volumeStep]];
                schedule.Add(curRegime);
            }
            schedule.Reverse();
            return schedule;
        }

        #endregion

        #region Реализация интерфейса TechnologicalSectionAlgorithm

        public List<double[]> GetSchedule(double[] targets, List<double[]> fixValues = null)
        {
            if (fixValues == null)
            {
                fixValues = AlgorithmHelper.CreateListOfArrays(Period, Dimension, -1.0);
            }
            return FormSchedule(targets[0], 0, fixValues.Select(x => x[0]).ToList()).Select(x => new double[1] { x }).ToList();
        }

        public List<double[]> GetContinuousSchedule(double[] targets, List<double[]> fixValues = null)
        {
            if (fixValues == null)
            {
                fixValues = AlgorithmHelper.CreateListOfArrays(Period, Dimension, -1.0);
            }
            return FormScheduleContinuousRegimeField(targets[0], fixValues.Select(x => x[0]).ToList()).Select(x => new double[1] { x }).ToList();
        }

        public bool IsRegimeAvaliableOnInterval(double[] regime, int interval)
        {
            double r = regime[0];

            return _avaliableRegimesOnIntervals[interval].IndexOf(r) >= 0;
        }

        public List<double[]> Regimes
        {
            get
            {
                return _regimes.Select(x => new double[1] { x }).ToList();
            }
        }

        public List<double[]> NormRegimes
        {
            get
            {
                return _normRegimes.Select(x => new double[1] { x }).ToList();
            }
        }

        public int Period
        {
            get { return _period; }
        }

        public int Dimension
        {
            get { return 1; }
        }
        
        public List<double[]> UpperBounds
        {
            get
            {
                return _upperBounds.Select(x => new double[1] { x }).ToList();
            }
        }

        public double[] MaxFlows
        {
            get
            {
                return new double[] { _regimes.Max() };
            }
        }

        #endregion
    }

}
