using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public interface ISection
    {
        Tuple<List<double[]>, List<int>> GetSchedule(List<Tuple<double[], int[]>> volumes);
        List<double[]> AddOutputComponent(List<double[]> schedule);
    }
}
