﻿using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Algorithms
{
    public class SectionMathModel : ISection
    {
        #region Поля
        
        #region Не расчетные

        List<double> _regimes;
        List<double> _repairs;
        int _period;
        double _maxError;

        #endregion

        #region Расчетные

        List<List<double>> _avaliableRegimesOnIntervals;

        double _notRepair;

        #endregion

        #endregion

        #region Свойства 

        public int Dimension => 1;

        public double MaxError
        {
            get => _maxError;
            set
            {
                if (value <= 0)
                    _maxError = 0.001;
                else
                    _maxError = value;
            }
        }

        public List<List<int>> ControlAvaliableIntervals
        {
            get;
            private set;
        }

        public int Period => _period;

        public List<List<int>> RepairsIntervals
        {
            get;
            private set;
        }

        #endregion

        #region Конструкторы

        public SectionMathModel(List<double> regimes, List<double> maxFlows)
        {
            if (regimes == null)
                throw new Exception();
            if (maxFlows == null)
                throw new Exception();
            if (regimes.Count() == 0)
                throw new Exception();
            if (maxFlows.Count() == 0)
                throw new Exception();
            if (regimes.Any(x => x < 0))
                throw new Exception();
            if (maxFlows.Any(x => x < 0))
                throw new Exception();

            _regimes = regimes.Select(x => x).Distinct().ToList();
            _regimes.Sort();

            _repairs = maxFlows.Select(x => x).ToList();

            _notRepair = _repairs.Max();

            _period = maxFlows.Count();
            
            MaxError = 0.001;

            RepairsIntervals = new List<List<int>>();
            bool currentRepair = false;
            List<int> repairIndexes = null;
            for (int i = 0; i < Period; i++)
            {
                if (_repairs[i] != _notRepair)
                {
                    if (!currentRepair)
                    {
                        repairIndexes = new List<int>() { i };
                        currentRepair = true;
                    }
                    else
                        repairIndexes.Add(i);
                }
                else if (currentRepair == true)
                {
                    RepairsIntervals.Add(repairIndexes);
                    currentRepair = false;
                    repairIndexes = null;
                }
            }
            if (repairIndexes != null)
                RepairsIntervals.Add(repairIndexes);
        }

        #endregion

        #region Методы

        #region Проверки 

        void CheckVolumes(TargetVolumes volumes)
        {
            if (volumes == null)
                throw new Exception();
            volumes.Check(Dimension);
        }

        void CheckSchedule(List<double[]> schedule)
        {
            if (schedule == null)
                throw new Exception();
            if (schedule.Count() == 0)
                throw new Exception();
            if (schedule.Any(x => x == null))
                throw new Exception();
            if (schedule.Any(x => x.Count() != Dimension))
                throw new Exception();
            if (schedule.Any(x => x.Any(y => y < 0)))
                throw new Exception();
        }

        #endregion

        #region Функции для расчетов

        private List<Tuple<List<double>, List<int>>> DecomposeVolumes(TargetVolumes volumes)
        {
            List<Tuple<List<double>, List<int>>> result = new List<Tuple<List<double>, List<int>>>();

            if (GetContinuousSchedule(volumes) == null)
                return null;

            foreach(var tuple in volumes.targetVolumes)
            {
                var indexes = tuple.Item2;
                var volume = tuple.Item1;
                var schedule = GetDiscreteSchedule(volume[0], indexes);
                double maxDif = volume[0] * MaxError;
                if (Math.Abs(schedule.Sum() - volume[0]) > maxDif)
                    return null;

                var res = new List<Tuple<List<double>, List<int>>>() { new Tuple<List<double>, List<int>>(new List<double>(), new List<int>()) };

                for (int i = 0; i < indexes.Count(); i++)
                {
                    int idx = indexes[i];
                    if (_repairs[idx] != _notRepair || volumes.GetFix(i) != null)
                    {
                        res.Add(new Tuple<List<double>, List<int>>(new List<double>() { schedule[i] }, new List<int>() { idx }));
                    }
                    else
                    {
                        res[0].Item1.Add(schedule[i]);
                        res[0].Item2.Add(idx);
                    }
                }

                result.AddRange(res);
            }

            return result;
        }

        private List<double> GetDiscreteSchedule(double volume, List<int> indexes)
        {
            int period = indexes.Count();
            List<double> result = new double[period].ToList();
            double sum = 0.0;
            for(int i = 0; i < indexes.Count(); i++)
            {
                result[i] = _avaliableRegimesOnIntervals[indexes[i]][0];
                sum += result[i];
            }
            int[] currentRegimes = new int[period];

            bool end = false;
            while (!end)
            {
                end = true;
                for (int i = 0; i < period; i++)
                {
                    int idx = indexes[i];
                    int regimeIdx = currentRegimes[i];
                    List<double> avalRegimes = _avaliableRegimesOnIntervals[idx];

                    if (avalRegimes.Count() > regimeIdx + 1)
                    {
                        double newSum = sum - avalRegimes[regimeIdx] + avalRegimes[regimeIdx + 1];
                        if (newSum > volume)
                        {
                            return result;
                        }
                        else
                        {
                            sum = newSum;
                            currentRegimes[i]++;
                            result[i] = avalRegimes[regimeIdx + 1];
                            end = false;
                        }
                    }
                }
            }

            return result;
        }

        private List<double> GetContinuousSchedule(double volume, List<int> indexes)
        {
            int period = indexes.Count();
            List<double> result = AlgorithmHelper.CreateListOfElements(period, volume / period);
            List<double> upperBound = _avaliableRegimesOnIntervals.Where((x, i) => indexes.Contains(i)).Select(x => x.Max()).ToList();

            if (upperBound.Sum() < volume)
                return null;

            List<bool> canChange = AlgorithmHelper.CreateListOfElements(period, true);
            int k = 1000;
            while (k -- > 0)
            {
                double vol = 0.0;
                for (int i = 0; i < period; i++)
                {
                    if (result[i] > upperBound[i])
                    {
                        vol += result[i] - upperBound[i];
                        result[i] = upperBound[i];
                        canChange[i] = false;
                    }
                }

                if (vol == 0.0)
                    return result;

                vol /= canChange.Count(x => x);
                for (int i = 0; i < period; i++)
                {
                    if (canChange[i])
                    {
                        result[i] += vol;
                    }
                }
            }

            throw new Exception();
        }

        #endregion

        #region Реализация интерфейсов

        public List<Tuple<List<double[]>, List<int>>> GetSolution(TargetVolumes volumes)
        {
            CheckVolumes(volumes);
            var result = DecomposeVolumes(volumes);
            if (result == null)
                return null;
            else
                return result.Select(x => new Tuple<List<double[]>, List<int>>(x.Item1.Select(y => new double[] { y }).ToList(), x.Item2)).ToList();
        }

        public void CalcDefaultIntervalsParameters(TargetVolumes volumes)
        {
            CheckVolumes(volumes);

            _avaliableRegimesOnIntervals = new List<double>[_period].ToList();
            ControlAvaliableIntervals = new List<List<int>>();

            foreach(var tuple in volumes.targetVolumes)
            {
                List<int> indexes = tuple.Item2;
                double volume = tuple.Item1[0];

                if (volume == 0.0)
                    for (int i = 0; i < indexes.Count(); i++)
                        _avaliableRegimesOnIntervals[indexes[i]] = new List<double> { 0.0 };
                else
                    for (int i = 0; i < indexes.Count(); i++)
                    {
                        _avaliableRegimesOnIntervals[indexes[i]] = _regimes.Where(regime => regime <= _repairs[indexes[i]]).ToList();

                        if (_avaliableRegimesOnIntervals[indexes[i]].Count() == 0)
                            throw new Exception();

                        if (_avaliableRegimesOnIntervals[indexes[i]].Count() > 1)
                        {
                            _avaliableRegimesOnIntervals[indexes[i]].Remove(0.0);
                            _avaliableRegimesOnIntervals[indexes[i]].Sort();
                        }
                    }
                
                foreach(var idx in indexes)
                { 
                    double[] fixVal = volumes.GetFix(idx);
                    if (fixVal != null)
                        _avaliableRegimesOnIntervals[idx] = new List<double> { fixVal[0] };
                }

                ControlAvaliableIntervals.Add(indexes.Where(x => _repairs[x] == _notRepair && volumes.GetFix(x) == null).ToList());
                ControlAvaliableIntervals.Last().Sort();
            }
        }

        public List<double[]> GetFullSchedule(List<double[]> schedule)
        {
            return schedule;
        }

        public List<double[]> GetShortSchedule(List<double[]> schedule)
        {
            return schedule;
        }

        public List<double[]> GetContinuousSchedule(TargetVolumes volumes)
        {
            CheckVolumes(volumes);

            List<double[]> result = AlgorithmHelper.CreateListOfArrays(_period, 1, 0.0);

            foreach (var tuple in volumes.targetVolumes)
            {
                var indexes = tuple.Item2;
                var volume = tuple.Item1[0];

                var notFixIndexes = new List<int>();
                double notFixVolume = volume;
                for(int i = 0; i < indexes.Count(); i++)
                {
                    int idx = indexes[i];
                    double[] fix = volumes.GetFix(idx);
                    if (fix == null)
                    {
                        notFixIndexes.Add(idx);
                    }
                    else
                    {
                        notFixVolume -= fix[0];
                        result[idx] = fix;
                    }
                }
                var res = GetContinuousSchedule(notFixVolume, notFixIndexes);
                for (int i = 0; i < notFixIndexes.Count(); i++)
                {
                    result[notFixIndexes[i]] = new double[] { res[i] };
                }
            }

            return result;
        }

        public RepairMathModel GetRepair(int interval)
        {
            double rep = _repairs[interval];
            return new RepairMathModel(rep, new double[] { }, rep);
        }

        public double[] GetLowerRegime(int interval, double[] currentRegime, bool output = false)
        {
            double val = currentRegime[0];
            var avalRegimes = _avaliableRegimesOnIntervals[interval].Where(x => x < val).ToList();
            if (avalRegimes.Count() == 0)
                return null;
            return new double[] { avalRegimes[AlgorithmHelper.NearestByModul(avalRegimes, val)] };
        }

        #endregion

        #endregion
    }

}
