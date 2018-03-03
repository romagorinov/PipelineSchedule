using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Algorithms
{
    public class SectionMathModel : ISection
    {
        #region Поля
        
        #region Не расчетные

        private List<double> _regimes;
        private List<double> _repairs;
        int _period;

        #endregion

        #region Расчетные

        List<List<double>> _avaliableRegimesOnIntervals;

        double _notRepair;

        List<int> _reapirDayIntervals;

        #endregion

        #endregion

        #region Свойства 

        public int Dimension => 1;

        public double MaxError
        {
            get;
            set;
        }

        public List<List<int>> ControlAvaliableIntervals
        {
            get;
            private set;
        }

        #endregion

        #region Конструкторы

        public SectionMathModel(List<double> regimes, List<double> maxFlows)
        {
            if (regimes == null || maxFlows == null)
                throw new ArgumentNullException();

            _regimes = regimes.Select(x => x).ToList();
            _regimes.Sort();

            _repairs = maxFlows.Select(x => x).ToList();

            _notRepair = _repairs.Max();

            _period = maxFlows.Count();

            HashSet<int> repairDays = new HashSet<int>();
            for (int i =  0; i < _period; i++)
            {
                int day = i / 24;
                if (_repairs[i] != _notRepair)
                    repairDays.Add(day);
            }

            _reapirDayIntervals = new List<int>();
            foreach (var day in repairDays.ToList())
                for (int hour = 0; hour < 24; hour++)
                    _reapirDayIntervals.Add(day * 24 + hour);

            MaxError = 0.001;
        }

        #endregion

        #region Методы

        #region Функции для расчетов

        public List<Tuple<List<double>, List<int>>> DecomposeVolumes(List<Tuple<double, int[]>> volumes)
        {
            List<Tuple<List<double>, List<int>>> result = new List<Tuple<List<double>, List<int>>>();
            
            foreach(var tuple in volumes)
            {
                var indexes = tuple.Item2;
                var volume = tuple.Item1;
                var schedule = GetDiscreteSchedule(volume, indexes);
                double maxDif = volume * MaxError;
                if (Math.Abs(schedule.Sum() - volume) > maxDif)
                    return null;

                var res = new List<Tuple<List<double>, List<int>>>() { new Tuple<List<double>, List<int>>(new List<double>(), new List<int>()) };

                for (int i = 0; i < indexes.Count(); i++)
                {
                    int idx = indexes[i];
                    if (_repairs[idx] != _notRepair)
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

        public List<double> GetDiscreteSchedule(double volume, int[] indexes)
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

        #endregion

        #region Реализация интерфейсов

        public List<Tuple<List<double[]>, List<int>>> GetSchedule(List<Tuple<double[], int[]>> volumes)
        {
            var convert = volumes.Select(x => new Tuple<double, int[]>(x.Item1[0], x.Item2)).ToList();
            var result = DecomposeVolumes(convert);
            if (result == null)
                return null;
            else
                return result.Select(x => new Tuple<List<double[]>, List<int>>(x.Item1.Select(y => new double[] { y }).ToList(), x.Item2)).ToList();
        }

        public void CalcDefaultIntervalsParameters(List<Tuple<double[], int[]>> volumes)
        {
            var temp = volumes.SelectMany(x => x.Item2);
            if (temp.Count() != temp.Distinct().Count())
                throw new Exception();

            _avaliableRegimesOnIntervals = new List<double>[_period].ToList();
            ControlAvaliableIntervals = new List<List<int>>();

            foreach(var tuple in volumes)
            {
                int[] indexes = tuple.Item2;
                double volume = tuple.Item1[0];

                if (volume == 0.0)
                    for (int i = 0; i < indexes.Count(); i++)
                        _avaliableRegimesOnIntervals[i] = new List<double> { 0.0 };
                else
                    for (int i = 0; i < indexes.Count(); i++)
                    {
                        _avaliableRegimesOnIntervals[i] = _regimes.Where(regime => regime <= _repairs[i]).ToList();

                        if (_avaliableRegimesOnIntervals[i].Count() == 0)
                            throw new Exception();

                        if (_avaliableRegimesOnIntervals[i].Count() > 1)
                        {
                            _avaliableRegimesOnIntervals[i].Remove(0.0);
                            _avaliableRegimesOnIntervals[i].Sort();
                        }
                    }

                ControlAvaliableIntervals.Add(indexes.Where(x => _repairs[x] == _notRepair).ToList());
                ControlAvaliableIntervals.Last().Sort();
            }
        }

        public List<double[]> GetFullSchedule(List<double[]> schedule)
        {
            return schedule;
        }

        #endregion

        #endregion
    }

}
