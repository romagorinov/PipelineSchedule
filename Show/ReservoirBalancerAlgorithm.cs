using Algorithms;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace Show
{
    public class ReservoirBalancerAlgorithm
    {
        public TechnologicalSectionAlgorithm tu1Algorithm;
        public TechnologicalSectionAlgorithm tu2Algorithm;
        public double rpVolume;
        public int period;

        public static double[] weights = { 0.4, 0.3, 0.3 };
        public static int Tmin = 1;
        public static int MaxIter = 10000;

        public ReservoirBalancerAlgorithm(TechnologicalSectionAlgorithm tuInAlgorithm, TechnologicalSectionAlgorithm tuOutAlgorithm, double rpVolume)
        {
            this.tu1Algorithm = tuInAlgorithm;
            this.tu2Algorithm = tuOutAlgorithm;
            this.rpVolume = rpVolume;
            period = tu1Algorithm.period;
        }
        
        public double[][] Balance(double Vplan1, double Vplan2, double[] fi1, double[] fi2, double Vstart, double[] z, Action<string> info)
        {
            double[] schedule1 = tu1Algorithm.FormSchedule(Vplan1, fi1),
                schedule2 = tu2Algorithm.FormSchedule(Vplan2, fi2);

            int idx = -1;
            double curVolume = Vstart;
            for (int i = 0; i < period; i++)
            {
                double addVolume = schedule1[i] - schedule2[i] + z[i];
                curVolume += addVolume;
                if (curVolume < 0 || curVolume > rpVolume)
                {
                    idx = i;
                    break;
                }
            }

            if (idx == -1)
                return new double[][] { schedule1, schedule2 };

            bool parallel = true;
            ConcurrentBag<double[][]> solutionsBag = new ConcurrentBag<double[][]>();
            int counter = 0,
                part = 0,
                partSize = 10;
            int threadId = 0;
            Parallel.Invoke(() =>
            {
                threadId = Thread.CurrentThread.ManagedThreadId;
            });
            if (parallel)
                Parallel.For(0, MaxIter, (it) =>
                {
                    var sol = GenerateSolution(Vplan1, Vplan2, fi1, fi2, Vstart, z);
                    if (sol != null)
                        solutionsBag.Add(sol);

                    Interlocked.Increment(ref counter);
                    if (counter > part * partSize && Thread.CurrentThread.ManagedThreadId == threadId)
                    {
                        info("Итерация " + counter.ToString() + ", сгенерированно решений " + solutionsBag.Count().ToString());
                        part++;
                    }
                });
            else
                for (int i = 0; i < MaxIter; i ++)
                {
                    var sol = GenerateSolution(Vplan1, Vplan2, fi1, fi2, Vstart, z);
                    if (sol != null)
                        solutionsBag.Add(sol);
                }

            if (solutionsBag.Count() == 0)
                return null;
            if (solutionsBag.Count() == 1)
                return solutionsBag.First();

            info("Cгенерированно решений " + solutionsBag.Count().ToString() + ", расчет критериев");
            var solutionList = solutionsBag.ToList();
            var criterias = solutionList.AsParallel().Select(sol => GetCriteria(sol, Vplan1, Vplan2)).ToList();
            var result = solutionList[criterias.IndexOf(criterias.Min())];

            return result;
        }

        public double[][] GenerateSolution(double Vplan1, double Vplan2, double[] fi1, double[] fi2, double Vstart, double[] z)
        {
            double[] fi1Cur = fi1.Select(x => x).ToArray(),
                fi2Cur = fi2.Select(x => x).ToArray();
            double[] schedule1 = tu1Algorithm.FormSchedule(Vplan1, fi1Cur),
                schedule2 = tu2Algorithm.FormSchedule(Vplan2, fi2Cur);

            int maxIter = 2 * tu1Algorithm.Nsut;
            while (maxIter-- > 0)
            {
                // Поиск индекса idx
                int idx = -1;
                double curVolume = Vstart;
                for (int i = 0; i < period; i++)
                {
                    double addVolume = schedule1[i] - schedule2[i] + z[i];
                    curVolume += addVolume;
                    if(curVolume < 0 || curVolume > rpVolume)
                    {
                        curVolume -= addVolume;
                        idx = i;
                        break;
                    }
                }

                if (idx == -1)
                    return new double[][] { schedule1, schedule2 };

                // Формируем множество пар
                double[] tu1U = fi1Cur[idx] >= 0 ? new double[] { fi1Cur[idx] } : tu1Algorithm.U[idx],
                    tu2U = fi2Cur[idx] >= 0 ? new double[] { fi2Cur[idx] } : tu2Algorithm.U[idx];
                var allPairs = AlgorithmHelper.CartesianProduct(new double[][] { tu1U, tu2U }).Select(x => x.ToArray());

                var avalPairs = allPairs.Select(pair =>
                {
                    double vol = curVolume;
                    for (int i = idx; i < period; i++)
                    {
                        if ((fi1Cur[i] >= 0 && fi1Cur[i] != pair[0]) || (fi2Cur[i] >= 0 && fi2Cur[i] != pair[1]))
                            return new { PAIR = pair, T = i - 1 };

                        if (!tu1Algorithm.U[i].Contains(pair[0]) || !tu2Algorithm.U[i].Contains(pair[1]))
                            return new { PAIR = pair, T = i - 1 };

                        double addVolume = pair[0] - pair[1] + z[i];
                        vol += addVolume;
                        if (vol < 0 || vol > rpVolume)
                            return new { PAIR = pair, T = i - 1 };
                    }
                    return new { PAIR = pair, T = period - 1 };
                })
                .Where(x => (x.T - idx + 1) >= Tmin || x.T == (period - 1)).ToArray();

                if (avalPairs.Count() == 0)
                    return null;

                int minZeros = avalPairs.Min(a => a.PAIR.Count(x => x == 0));
                avalPairs = avalPairs.Where(a => a.PAIR.Count(x => x == 0) == minZeros).ToArray();

                // Формируем функции предпочтения выбора пары
                double[] curPair = new double[] { schedule1[idx], schedule2[idx] };
                var pref = avalPairs.Select(a =>  1 / AlgorithmHelper.GetDistance(curPair, a.PAIR)).ToArray();

                // Вычисляем вероятности
                var prefSum = pref.Sum();
                var probabl = pref.Select(x => x / prefSum).ToList();

                // Выбираем пару 
                var selectedPair = avalPairs[AlgorithmHelper.GetRouletIndex(probabl)];

                // Выбираем, до какого момента работать
                int minT = Math.Min(idx + Tmin - 1, period - 1),
                    maxT = selectedPair.T;
                int selectedT = AlgorithmHelper.RandomInstance.Next(minT, maxT + 1);

                for (int i = idx; i <= selectedT; i++)
                {
                    schedule1[i] = selectedPair.PAIR[0];
                    schedule2[i] = selectedPair.PAIR[1];
                }

                for (int i = 0; i <= selectedT; i++)
                {
                    fi1Cur[i] = schedule1[i];
                    fi2Cur[i] = schedule2[i];
                }

                schedule1 = tu1Algorithm.FormSchedule(Vplan1, fi1Cur);
                schedule2 = tu2Algorithm.FormSchedule(Vplan2, fi2Cur);
            }

            return null;
        }

        public double GetCriteria(double[][] sol, double Vplan1, double Vplan2)
        {
            return tu1Algorithm.GetCriteria(sol[0], Vplan1, weights) + tu2Algorithm.GetCriteria(sol[1], Vplan2, weights);
        }

    }
}
