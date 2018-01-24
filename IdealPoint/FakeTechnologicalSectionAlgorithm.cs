using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class FakeTechnologicalSectionAlgorithm : ReservoirBalanceAlgorithm.TechnologicalSectionAlgorithm
    {
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
            get;
            private set;
        }

        public int Dimension
        {
            get;
            private set;
        }

        public List<double[]> UpperBounds
        {
            get
            {
                return _schedule.Select(x => x.Select(y => y).ToArray()).ToList();
            }
        }

        public List<double[]> GetContinuousSchedule(double[] targets, List<double[]> fixValues = null)
        {
            return _schedule.Select(x => x.Select(y => y).ToArray()).ToList();
        }

        public List<double[]> GetSchedule(double[] targets, List<double[]> fixValues = null)
        {
            return _schedule.Select(x => x.Select(y => y).ToArray()).ToList();
        }

        public bool IsRegimeAvaliableOnInterval(double[] regime, int interval)
        {
            return _schedule[interval].SequenceEqual(regime);
        }

        public double[] MaxFlows
        {
            get
            {
                return AlgorithmHelper.GetMaxInListByComponents(_regimes);
            }
        }

        List<double[]> _regimes;
        List<double[]> _normRegimes;

        List<double[]> _schedule;

        public FakeTechnologicalSectionAlgorithm(List<double[]> schedule)
        {
            _schedule = schedule;
            _regimes = AlgorithmHelper.RemoveDuplicates(_schedule);
            _normRegimes = AlgorithmHelper.NormalizeByComponents(_regimes);

            Period = schedule.Count();
            Dimension = schedule[0].Count();
        }
    }
}
