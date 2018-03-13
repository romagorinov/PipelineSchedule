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
        public Dictionary<int, double[]> fixValues;
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
            fixValues = new Dictionary<int, double[]>();
        }

        public TargetVolumes(TargetVolumes p)
        {
            period = p.period;
            dimension = p.dimension;
            targetVolumes = p.targetVolumes.Select(x => new Tuple<double[], List<int>>(x.Item1.ToArray(), x.Item2.ToList())).ToList();
            fixValues = p.fixValues.ToDictionary(kv => kv.Key, kv => kv.Value.ToArray());
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

            indexes.Sort();
            targetVolumes.Add(new Tuple<double[], List<int>>(volume, indexes));
        }

        public void AddFixValue(int interval, double[] value)
        {
            if (fixValues.ContainsKey(interval))
                fixValues[interval] = value;
            else
                fixValues.Add(interval, value);
        }

        public double[] GetFix(int interval)
        {
            if (fixValues.ContainsKey(interval))
                return fixValues[interval];
            else
                return null;
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
            if (fixValues == null)
                throw new Exception();
            if (fixValues.Any(x => x.Key < 0 || x.Key > period - 1))
                throw new Exception();
            if (fixValues.Any(x => x.Value == null))
                throw new Exception();
            if (fixValues.Any(x => x.Value.Count() != dim))
                throw new Exception();
            if (fixValues.Any(x => x.Value.Any(y => y < 0)))
                throw new Exception();

        }

        public int GetBiggestBatch(int startIndex, int endIndex)
        {
            if (startIndex < 0 || startIndex > period - 1 || endIndex < 1 || endIndex > period)
                throw new Exception();
            if (endIndex <= startIndex)
                throw new Exception();

            var indexesCount = targetVolumes.Select(x => x.Item2.Count(y => y >= startIndex && y < endIndex)).ToList();
            int max = indexesCount.Max();
            return indexesCount.IndexOf(max);
        }
        
        public bool SplitBatch(int batchIdx, int splitIdx, double[] beforeVolume)
        {
            if (batchIdx < 0 || batchIdx > targetVolumes.Count() - 1)
                throw new Exception();
            if (splitIdx < 1 || splitIdx > period - 1)
                throw new Exception();
            if (beforeVolume == null)
                throw new Exception();
            if (beforeVolume.Count() != dimension)
                throw new Exception();
            if (beforeVolume.Any(x => x < 0))
                throw new Exception();

            var volume = targetVolumes[batchIdx].Item1;
            var indexes = targetVolumes[batchIdx].Item2;

            int batchSplitPoint = 0;
            for (int i = 1; i < indexes.Count(); i++)
            {
                if (indexes[i] >= splitIdx)
                {
                    batchSplitPoint = indexes[i];
                    break;
                }
            }

            List<int> indexes1 = indexes.GetRange(0, batchSplitPoint),
                indexes2 = indexes.GetRange(batchSplitPoint, indexes.Count() - batchSplitPoint);

            double[] volume1 = beforeVolume;
            double[] volume2 = volume.Zip(beforeVolume, (x,y) => x - y).ToArray();

            targetVolumes.RemoveAt(batchIdx);
            targetVolumes.Add(new Tuple<double[], List<int>>(volume1, indexes1));
            targetVolumes.Add(new Tuple<double[], List<int>>(volume1, indexes1));

            return true;
        }
                
        public List<double[]> GetUniformVolumes()
        {
            List<double[]> result = new double[period][].ToList();

            foreach(var kv in fixValues)
            {
                result[kv.Key] = kv.Value;
            }

            foreach(var tuple in targetVolumes)
            {
                List<int> indexes = tuple.Item2.ToList();
                double[] volume = tuple.Item1.ToList().ToArray();
                foreach(var idx in indexes)
                {
                    if (fixValues.ContainsKey(idx))
                    {
                        indexes.Remove(idx);
                        volume = volume.Zip(fixValues[idx], (x, y) => x - y).ToArray();
                    }
                }
                volume = volume.Select(x => x / indexes.Count()).ToArray();
                foreach (var idx in indexes)
                {
                    result[idx] = volume.ToList().ToArray();
                }
            }

            return result;
        }
    }
}
