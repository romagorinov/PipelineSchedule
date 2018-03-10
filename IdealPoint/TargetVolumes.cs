using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class TargetVolumes
    {
        public List<Tuple<double[], List<int>>> targetVolumes;
        public int period;
        public int dimension;

        public TargetVolumes(int period, int dimension)
        {
            if (period < 1)
                throw new Exception();
            if (dimension < 1)
                throw new Exception();

            this.period = period;
            this.dimension = dimension;
            targetVolumes = new List<Tuple<double[], List<int>>>();
        }

        public void AddVolume(double[] volume, List<int> indexes)
        {
            if (volume == null)
                throw new Exception();
            if (volume.Count() < dimension)
                throw new Exception();
            if (indexes == null)
                throw new Exception();
            if (indexes.Count() == 0)
                throw new Exception();
            if (targetVolumes.Any(x => x.Item2.Intersect(indexes).Count() > 0))
                throw new Exception();

            if (volume.Count() > dimension)
                volume = volume.ToList().GetRange(0, dimension).ToArray();

            targetVolumes.Add(new Tuple<double[], List<int>>(volume, indexes));
        }

        public void Check(int dim)
        {
            if (targetVolumes == null)
                throw new Exception();
            if (targetVolumes.Count() == 0)
                throw new Exception();
            if (targetVolumes.Any(x => x == null))
                throw new Exception();
            if (targetVolumes.Any(x => x.Item1 == null))
                throw new Exception();
            if (targetVolumes.Any(x => x.Item2 == null))
                throw new Exception();
            if (targetVolumes.Any(x => x.Item1.Count() != dim))
                throw new Exception();
            if (targetVolumes.Any(x => x.Item1.Any(y => y < 0)))
                throw new Exception();
            if (targetVolumes.Any(x => x.Item2.Any(y => y < 0 || y > period - 1)))
                throw new Exception();
            if (targetVolumes.SelectMany(x => x.Item2).GroupBy(x => x).Any(x => x.Count() > 1))
                throw new Exception();
        }
    }
}
