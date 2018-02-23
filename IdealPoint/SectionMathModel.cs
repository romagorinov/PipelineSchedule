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
    public class SectionMathModel : ISection
    {
        #region Поля
        
        #region Параметры секции

        private List<double> _regimes;
        private List<List<double>> _avaliableRegimesOnIntervals;
        int _period;

        #endregion

        #endregion

        #region Свойства 

        public double MaxError
        {
            get;
            set;
        }

        #endregion

        #region Конструкторы

        public SectionMathModel(List<double> regimes, List<double> maxFlows)
        {
            if (regimes == null || maxFlows == null)
                throw new ArgumentNullException();

            _regimes = regimes.Select(x => x).ToList();
            _regimes.Sort();

            _avaliableRegimesOnIntervals = maxFlows.Select(maxFlow => regimes.Where(regime => regime <= maxFlow).ToList()).ToList();

            for (int i = 0; i < _avaliableRegimesOnIntervals.Count(); i++)
                _avaliableRegimesOnIntervals[i].Sort();

            _period = _avaliableRegimesOnIntervals.Count();

            MaxError = 0.001;
        }

        #endregion

        #region Методы

        public Tuple<List<double>, List<int>> DecomposeVolumes(List<Tuple<double, int[]>> volumes)
        {
            Tuple<List<double>, List<int>> result = new Tuple<List<double>, List<int>>(new List<double>(), new List<int>());
            
            foreach(var tuple in volumes)
            {
                var indexes = tuple.Item2;
                var volume = tuple.Item1;
                var schedule = GetDiscreteSchedule(volume, indexes);
                double maxDif = volume * MaxError;
                if (Math.Abs(schedule.Sum() - volume) > maxDif)
                    return null;

                result.Item1.AddRange(schedule);
                result.Item2.AddRange(indexes);
            }

            return result;
        }

        public List<double> GetDiscreteSchedule(double volume, int[] indexes)
        {
            int period = indexes.Count();
            List<double> result = new double[period].ToList();
            double sum = 0.0;
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

        public Tuple<List<double[]>, List<int>> GetSchedule(List<Tuple<double[], int[]>> volumes)
        {
            var convert = volumes.Select(x => new Tuple<double, int[]>(x.Item1[0], x.Item2)).ToList();
            var result = DecomposeVolumes(convert);
            return new Tuple<List<double[]>, List<int>>(result.Item1.Select(x => new double[] { x }).ToList(), result.Item2);
        }

        #endregion
    }

}
