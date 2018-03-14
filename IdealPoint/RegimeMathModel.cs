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
        #region Поля

        const double EPS = 1E-8;

        public readonly double Gin;

        public readonly double[] GpumpMax;

        public readonly double[] GpumpMin;

        public readonly double[] pumpSign;

        private readonly int _pumpsCount;

        #endregion

        #region Свойства

        public int PumpsCount => _pumpsCount;

        #endregion

        #region Конструкторы

        public RegimeMathModel(double Gin, double[] GpumpMax, double[] GpumpMin, double[] pumpSign)
        {
            if (pumpSign == null)
                throw new Exception();
            if (pumpSign.Count() == 0)
                throw new Exception();
            if (pumpSign.Any(x => x != -1 && x != 1))
                throw new Exception();
            if (Gin < 0)
                throw new Exception();
            if (GpumpMax == null)
                throw new Exception();
            if (GpumpMax.Count() != pumpSign.Count())
                throw new Exception();
            if (GpumpMax.Any(x => x < 0))
                throw new Exception();
            if (GpumpMin == null)
                throw new Exception();
            if (GpumpMin.Count() != pumpSign.Count())
                throw new Exception();
            if (GpumpMin.Any(x => x < 0))
                throw new Exception();
            if (GpumpMax.Zip(GpumpMin, (max, min) => max < min).Any(x => x))
                throw new Exception();
            
            this.Gin = Gin;
            this.GpumpMax = GpumpMax.ToList().ToArray();
            this.GpumpMin = GpumpMin.ToList().ToArray();
            this.pumpSign = pumpSign.ToList().ToArray();
            _pumpsCount = pumpSign.Count();
        }

        #endregion

        #region Методы

        private void CheckRepair(RepairMathModel repair)
        {
            if (repair.MaxPumps.Count() != PumpsCount)
                throw new Exception();
        }

        private void CheckPumps(double[] Gpumps)
        {
            if (Gpumps == null)
                throw new Exception();
            if (Gpumps.Count() != PumpsCount)
                throw new Exception();
            if (Gpumps.Any(x => x < 0))
                throw new Exception();
        }

        /// <summary>
        /// Система "МЕНЬШЕ ИЛИ РАВНО"
        /// </summary>
        /// <param name="repair"></param>
        /// <returns></returns>
        public Tuple<double[][], double[]> GetSystemOfInequalities(RepairMathModel repair)
        {
            CheckRepair(repair);

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
            CheckRepair(repair);

            return repair.MaxInput >= Gin
                && GpumpMin.Zip(repair.MaxPumps, (m, u) => m <= u).All(x => x)
                && GpumpMin.Zip(pumpSign, (x, y) => x * y).Sum() + Gin <= repair.MaxOutput;
        }
        
        public bool CanPump(double[] Gpumps, RepairMathModel repair)
        {
            CheckRepair(repair);
            CheckPumps(Gpumps);

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
        
        public bool CanProvide(double[] G)
        {
            if (G == null)
                throw new Exception();
            if (G.Count() != _pumpsCount + 1)
                throw new Exception();
            if (G.Any(x => x < 0))
                throw new Exception();

            if (G[0] != Gin)
                return false;

            for (int i = 0; i < _pumpsCount; i++)
            {
                double providePump = G[i + 1];
                if (providePump < GpumpMin[i] || providePump > GpumpMax[i])
                    return false;
            }

            return true;
        }
        
        public double[] GetNearestPumpsPoit(double[] Gpumps, RepairMathModel repair, double[] pumpsPriority = null)
        {
            CheckRepair(repair);
            CheckPumps(Gpumps);
            if (pumpsPriority == null)
                pumpsPriority = pumpSign.Select(x => 1.0).ToArray();
            else
            {
                if (pumpsPriority.Count() != PumpsCount)
                    throw new Exception();
                if (pumpsPriority.Any(x => x < 0))
                    throw new Exception();
            }

            var system = GetSystemOfInequalities(repair);
            double[,] A = system.Item1.ToMatrix().Multiply(-1.0);
            double[] b = system.Item2.Multiply(-1.0);

            double[,] Q = new double[PumpsCount, PumpsCount];
            double[] d = new double[PumpsCount];

            for (int i =  0; i < PumpsCount; i++)
            {
                Q[i, i] = 2.0 * pumpsPriority[i] * pumpsPriority[i];
                d[i] = pumpsPriority[i] * pumpsPriority[i] * -2.0 * Gpumps[i];
            }

            return AlgorithmHelper.SolveQP(Q, d, A, b, 0);
        }
        
        public double GetOutputElement(double[] Gpumps)
        {
            CheckPumps(Gpumps);

            return Gin + Gpumps.Zip(pumpSign, (x, y) => x * y).Sum();
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

        #endregion
    }
}
