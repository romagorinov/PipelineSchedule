using System;
using System.Collections.Generic;

namespace Algorithms
{
    public interface ISection
    {
        List<Tuple<List<double[]>, List<int>>> GetSchedule(List<Tuple<double[], int[]>> volumes);        

        int Dimension
        {
            get;
        }

        List<List<int>> ControlAvaliableIntervals
        {
            get;
        }

        void CalcDefaultIntervalsParameters(List<Tuple<double[], int[]>> volumes);

        List<double[]> GetFullSchedule(List<double[]> schedule);
    }
}
