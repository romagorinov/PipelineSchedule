using System;
using System.Collections.Generic;

namespace Algorithms
{
    public interface ISection
    {
        #region Свойства

        int Dimension
        {
            get;
        }

        int Period
        {
            get;
        }

        List<List<int>> ControlAvaliableIntervals
        {
            get;
        }

        List<List<int>> RepairsIntervals
        {
            get;
        }

        #endregion

        #region Методы

        List<Tuple<List<double[]>, List<int>>> GetSolution(TargetVolumes volumes);

        List<double[]> GetContinuousSchedule(TargetVolumes volumes);

        void CalcDefaultIntervalsParameters(TargetVolumes volumes);

        List<double[]> GetFullSchedule(List<double[]> schedule);

        List<double[]> GetShortSchedule(List<double[]> schedule);

        RepairMathModel GetRepair(int interval);

        double[] GetLowerRegime(int interval, double[] currentRegime, bool output = false);

        #endregion
    }
}
