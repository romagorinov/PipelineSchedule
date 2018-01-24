using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Show
{
    public class TechnologicalSectionAlgorithm
    {
        public double[] G;
        public double[] gmax;
        public int[] Ucount;
        public double[][] U;
        public int period;
        public int Nsut;

        public TechnologicalSectionAlgorithm(double[] gmax, double[] G, int Nsut)
        {
            this.gmax = gmax.Select(x => x).ToArray();
            this.G = G.Select(x => x).ToArray();
            this.Nsut = Nsut;
            Array.Sort(this.G);
            period = this.gmax.Count();
            Ucount = new int[period];
            U = new double[period][];
            for (int i = 0; i < period; i++)
            {
                U[i] = this.G.Where(x => x <= this.gmax[i]).ToArray();
                Ucount[i] = U[i].Count();
            }
        }

        public double[] FormSchedule(double Vplan, double[] fi, int point = 0)
        {
            double[] schedule = fi.Select(x => x < 0 ? 0 : x).ToArray();
            double sum = schedule.Sum();

            // Начинаем заполнять по уровням
            for (int i = 1; i < G.Count(); i++)
            {
                // Начинаем от point
                int leftPoint = point, rightPoint = point;

                // Пока не дошли до правой и левой границы
                while (leftPoint >= 0 || rightPoint < period)
                {
                    // Если правая граница не выходит за пределы, и пропускная позволяет
                    if (rightPoint < period && Ucount[rightPoint] > i && fi[rightPoint] < 0)
                    {
                        double newSum = sum - schedule[rightPoint] + G[i];
                        if (newSum <= Vplan)
                        {
                            sum = newSum;
                            schedule[rightPoint] = G[i];
                        }
                        else
                        {
                            return schedule;
                        }
                    }

                    if (leftPoint >= 0 && Ucount[leftPoint] > i && fi[leftPoint] < 0)
                    {
                        double newSum = sum - schedule[leftPoint] + G[i];
                        if (newSum <= Vplan)
                        {
                            sum = newSum;
                            schedule[leftPoint] = G[i];
                        }
                        else
                        {
                            return schedule;
                        }
                    }

                    leftPoint--;
                    rightPoint++;
                }
            }

            return schedule;
        }

        private double GetVolumeCriteria(double[] x, double Vplan)
        {
            double criteria = Math.Pow((x.Sum() - Vplan), 2);
            double maxVal = Math.Max(
                Math.Pow(Vplan, 2),
                Math.Pow(Ucount.Sum(i => G[i - 1]) - Vplan, 2));
            return criteria / maxVal;
        }

        private double GetChangeCriteria(double[] x)
        {
            double criteria = 0;
            double maxVal = (2 * Nsut);
            for (int i = 1; i < x.Count(); i++)
            {
                if (x[i] != x[i - 1])
                    criteria++;
            }
            return criteria / maxVal;
        }

        private double GetUniformCriteria(double[] x, double Vplan)
        {
            double uniform = Vplan / period;
            double criteria = 0;
            double maxVal = 0;
            for (int i = 0; i < x.Count(); i++)
            {
                criteria += Math.Pow(x[i] - uniform, 2);
                maxVal += Math.Max(
                    Math.Pow(uniform, 2),
                    Math.Pow(G[Ucount[i] - 1] - uniform, 2));
            }
            return criteria / maxVal;
        }

        public double GetCriteria(double[] x, double Vplan, double[] weights)
        {
            return weights[0] * GetVolumeCriteria(x, Vplan)
                + weights[1] * GetChangeCriteria(x)
                + weights[2] * GetUniformCriteria(x, Vplan);
        }

        public double UpperBound
        {
            get => U.Sum(u => u.Last());
        }
    }
}
