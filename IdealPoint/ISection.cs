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

        #endregion

        #region Методы

        List<Tuple<List<double[]>, List<int>>> GetSchedule(List<Tuple<double[], int[]>> volumes);

        void CalcDefaultIntervalsParameters(List<Tuple<double[], int[]>> volumes);

        List<double[]> GetFullSchedule(List<double[]> schedule);

        List<double[]> GetShortSchedule(List<double[]> schedule);

        #endregion
    }
}
