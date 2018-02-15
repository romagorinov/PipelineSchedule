using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using System.Text;
using System.Threading.Tasks;
using Accord.Math.Optimization;

namespace Algorithms
{
    public class RegimeMathModel
    {
        const double EPS = 1E-8;

        public double Gin;

        public double[] GpumpMax;
        public double[] GpumpMin;

        public double[] pumpSign;

        public int PumpsCount
        {
            get;
            set;
        }

        public RegimeMathModel(double Gin, double[] GpumpMax, double[] GpumpMin, double[] pumpSign)
        {
            this.Gin = Gin;
            this.GpumpMax = GpumpMax;
            this.GpumpMin = GpumpMin;
            this.pumpSign = pumpSign;
            PumpsCount = pumpSign.Count();
        }

        public Tuple<double[][], double[]> GetSystemOfInequalities(RepairMathModel repair)
        {
            int rowsCount = PumpsCount * 2 + 1;
            double[][] A = new double[rowsCount][];
            double[] b = new double[rowsCount];

            A[rowsCount - 1] = new double[PumpsCount];
            for (int i = 0; i < PumpsCount; i++)
            {
                A[i * 2] = new double[PumpsCount];
                A[i * 2 + 1] = new double[PumpsCount];
                A[i * 2][i] = 1.0;
                A[i * 2 + 1][i] = -1.0;
                b[i * 2] = Math.Min(GpumpMax[i], repair.MaxPumps[i]);
                b[i * 2 + 1] = - GpumpMin[i];

                A[rowsCount - 1][i] = pumpSign[i];
            }
            b[rowsCount - 1] = repair.MaxOutput - Gin;

            return new Tuple<double[][], double[]>(A, b);
        }

        public bool CanUse(RepairMathModel repair)
        {
            return repair.MaxInput >= Gin
                && GpumpMin.Zip(repair.MaxPumps, (m, u) => m <= u).All(x => x)
                && GpumpMin.Zip(pumpSign, (x, y) => x * y).Sum() + Gin <= repair.MaxOutput;
        }
        
        public bool CanPump(double[] Gpumps, RepairMathModel repair)
        {
            Tuple<double[][], double[]> ineq = GetSystemOfInequalities(repair);
            // Смотрим, можно или нет качать
            for (int i = 0; i < ineq.Item1.Count(); i++)
            {
                var multipleResult = ineq.Item1[i].Zip(Gpumps, (x, y) => x * y).Sum();
                if (multipleResult - ineq.Item2[i] > EPS)
                    return false;
            }

            return true;
        }

        public bool MeetsZeroConditions(bool inZero, bool[] pumpsZero)
        {
            if (!inZero && Gin == 0)
                return false;

            if (GpumpMax.Zip(pumpsZero, (x, y) => !y && x == 0).Any(x => x))
                return false;

            return true;
        }
        
        public double[] GetNearestPumpsPoit(double[] pumpsVolume, RepairMathModel repair, double[] pumpsPriority = null)
        {
            if (pumpsPriority == null)
                pumpsPriority = pumpSign.Select(x => 1.0).ToArray();

            var system = GetSystemOfInequalities(repair);
            double[,] A = system.Item1.ToMatrix().Multiply(-1.0);
            double[] b = system.Item2.Multiply(-1.0);

            double[,] Q = new double[PumpsCount, PumpsCount];
            double[] d = new double[PumpsCount];

            for (int i =  0; i < PumpsCount; i++)
            {
                Q[i, i] = 2.0 * pumpsPriority[i] * pumpsPriority[i];
                d[i] = pumpsPriority[i] * pumpsPriority[i] * -2.0 * pumpsVolume[i];
            }

            return AlgorithmHelper.SolveQP(Q, d, A, b, 0);
        }

        /*public static bool operator == (RegimeMathModel r1, RegimeMathModel r2)
        {
            return r1.Gin == r2.Gin
                && r1.GpumpMax.SequenceEqual(r2.GpumpMax)
                && r1.GpumpMin.SequenceEqual(r2.GpumpMin)
                && r1.pumpSign.SequenceEqual(r2.pumpSign);
        }

        public static bool operator != (RegimeMathModel r1, RegimeMathModel r2)
        {
            return r1.Gin != r2.Gin
                || !r1.GpumpMax.SequenceEqual(r2.GpumpMax)
                || !r1.GpumpMin.SequenceEqual(r2.GpumpMin)
                || !r1.pumpSign.SequenceEqual(r2.pumpSign);
        }*/
    }
}
