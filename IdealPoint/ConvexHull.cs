using System;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using System.Threading.Tasks;
using Accord.Math.Optimization;

namespace Algorithms
{    
    public class ConvexHull
    {
        public const double EPS = 1E-8;

        protected List<double[]> _basePoints;
        protected double[] _max, _min;
        protected bool[] _notEqMask;

        protected List<double[]> _realPoints;
        protected int _realDimension;
        protected double[] _realMax, _realMin;

        protected MIConvexHull.ConvexHull<MIConvexHull.DefaultVertex, MIConvexHull.DefaultConvexFace<MIConvexHull.DefaultVertex>> _convexHull;

        protected Tuple<double[,], double[]> _systemOfRealInequalities;

        public int Dimension
        {
            get;
            protected set;
        }

        public List<double[]> ConvexHullPoints
        {
            get
            {
                if (_realDimension == 0)
                    return new List<double[]> { _max.Select(x => x).ToArray() };
                else if (_realDimension == 1)
                    return new List<double[]> { _min.Select(x => x).ToArray(), _max.Select(x => x).ToArray() };
                else
                    return _convexHull.Points.Select(p => ConvertFromReal(p.Position)).ToList();
            }
        }

        public double[,] Amatrix
        {
            get;
            private set;
        }

        public double[] Bvector
        {
            get;
            private set;
        }

        public int EqNumber
        {
            get;
            private set;
        }

        public double[] Max => _max;

        public double[] Min => _min;

        public ConvexHull(List<double[]> basePoints)
        {
            int dimension = basePoints[0].Count();
            if (dimension == 0) throw new ArgumentException();
            if (basePoints.Any(x => x.Count() != dimension)) throw new ArgumentException();
            
            Dimension = dimension;
            _basePoints = AlgorithmHelper.RemoveDuplicates(basePoints);
            _max = AlgorithmHelper.GetMaxInListByComponents(_basePoints);
            _min = AlgorithmHelper.GetMinInListByComponents(_basePoints);
            _notEqMask = _max.Zip(_min, (x, y) => x != y).ToArray();
            _realDimension = _notEqMask.Count(x => x);

            if (_realDimension > 0)
            {
                _realPoints = AlgorithmHelper.RemoveDuplicates(AlgorithmHelper.MaskListOfArrays(_basePoints, _notEqMask));
                _realMax = AlgorithmHelper.GetArrayByMask(_max, _notEqMask);
                _realMin = AlgorithmHelper.GetArrayByMask(_min, _notEqMask);

                if (_realDimension > 1)
                    _convexHull = MIConvexHull.ConvexHull.Create(_realPoints);
            }

            CalculateInequalities();
        }

        private void CalculateInequalities()
        {
            EqNumber = Dimension - _realDimension;
            List<double[]> AEqList = new List<double[]>();
            List<double> BEqList = new List<double>();
            for (int i = 0; i < Dimension; i++)
            {
                if (!_notEqMask[i])
                {
                    double[] a = new double[Dimension];
                    a[i] = 1.0;
                    double b = _max[i];
                    AEqList.Add(a);
                    BEqList.Add(b);
                }
            }

            double[,] AIneq = new double[0,0];
            double[] BIneq = new double[0];
            if (_realDimension == 1)
            {
                AIneq = new double[2, 1];
                BIneq = new double[2];
                AIneq[0,0] = 1.0;
                AIneq[0,1] = -1.0;
                BIneq[0] = _realMin[0];
                BIneq[1] = -_realMax[0];
            }
            else if (_realDimension > 1)
            {
                var faces = _convexHull.Faces.ToList();
                int fCount = faces.Count();
                // Получаем свободный член уравнений поверхностей
                double[] facesD = faces.Select(face =>
                {
                    return face.Normal.Zip(face.Vertices[0].Position, (x, y) => x * y).Sum();
                }).ToArray();
                AIneq = new double[fCount, _realDimension];
                BIneq = new double[fCount];
                for (int i = 0; i < fCount; i++)
                {
                    BIneq[i] = facesD[i];
                    for (int j = 0; j < _realDimension; j++)
                    {
                        AIneq[i, j] = faces[i].Normal[j];
                    }
                }
            }
            _systemOfRealInequalities = new Tuple<double[,], double[]>(AIneq, BIneq);

            int constrNumber = BEqList.Count() + BIneq.Count();
            double[,] Aall = new double[constrNumber, Dimension];
            double[] Ball = new double[constrNumber];
            for (int i = 0; i < EqNumber; i++)
            {
                for (int j = 0; j < Dimension; j++)
                {
                    Aall[i, j] = AEqList[i][j];
                    Ball[i] = BEqList[i];
                }
            }

            for (int i = EqNumber; i < constrNumber; i++)
            {
                int idx = 0;
                for (int j = 0; j < Dimension; j++)
                {
                    if (_notEqMask[j])
                    {
                        Aall[i,j] = AIneq[i - EqNumber, idx];
                        idx++;
                    }
                    else
                        Aall[i,j] = _max[j];
                }
                Ball[i] = BIneq[i];
            }

            Amatrix = Aall;
            Bvector = Ball;
        }

        /// <summary>
        /// Создает ConvexHull из системы линейных уравнений A*x <= b
        /// </summary>
        /// <param name="A">Матрица A</param>
        /// <param name="b">вектор b</param>
        /// <returns></returns>
        public static ConvexHull CreateFromHPolytope(double[][] A, double[] b)
        {
            int[] rowsNumbers = A.Select((x, i) => i).ToArray();
            int varsCount = A[0].Count();
            List<double[]> vertices = new List<double[]>();
            foreach(var combination in Combinatorics.Combinations(rowsNumbers, varsCount))
            {
                double[,] Acomb = new double[varsCount, varsCount];
                double[] bcomb = new double[varsCount];
                for (int i = 0; i < varsCount; i++)
                {
                    for (int j = 0; j < varsCount; j++)
                    {
                        Acomb[i,j] = A[combination[i]][j];
                    }
                    bcomb[i] = b[combination[i]];
                }

                if (!Acomb.IsSingular())
                {
                    var potentialV = Acomb.Solve(bcomb);
                    if (A.Dot(potentialV).Subtract(b).All(x => x <= EPS))
                        vertices.Add(potentialV);
                }
            }

            if (vertices.Count() > 0)
                return new ConvexHull(vertices);
            else return null;
        }
        
        //public List<double[]> GetGridPoints(int[] stepsCount)
        //{
        //    if (_realDimension == 0)
        //        return new List<double[]> { _basePoints[0] };

        //    List<double[]> grid;
        //    int[] realStepsCount = AlgorithmHelper.GetArrayByMask(stepsCount, _notEqMask);
        //    if (_realDimension == 1)
        //        grid = GetSingleDimensionGrid(realStepsCount);
        //    else
        //        grid = GetMultiDimensionGrid(realStepsCount);

        //    if (_realDimension != Dimension)
        //    {
        //        grid = grid.Select(realPoint => ConvertFromReal(realPoint)).ToList();
        //    }

        //    return grid;
        //}

        public bool IsPointInConvexHull(double[] point)
        {
            if (point.Count() != Dimension) throw new Exception();

            if (_realDimension == 0)
                return point.SequenceEqual(_max);

            if (_realDimension != Dimension)
            {
                for (int i = 0; i < Dimension; i++)
                {
                    if (!_notEqMask[i] && point[i] != _max[i])
                        return false;
                }
            }

            if (_realDimension == 1)
            {
                int idx = _notEqMask.IndexOf(true);
                return point[idx] <= _max[idx] && point[idx] >= _min[idx];
            }

            return IsRealPointInConvexHull(AlgorithmHelper.GetArrayByMask(point, _notEqMask));
        }
        
        public double[] FindNearestInteriorPoint(double[] point, bool[] notUpper, double[] weights = null)
        {
            if (weights == null)
                weights = _max.Select(x => 1.0).ToArray();

            if (point.Count() != Dimension) throw new Exception();

            if (_realDimension == 0)
                return _max;

            var nearestPoint = point.Select(x => x).ToArray();
            if (_realDimension != Dimension)
            {
                for (int i = 0; i < Dimension; i++)
                {
                    if (!_notEqMask[i])
                        nearestPoint[i] = _max[i];
                }
            }


            if (_realDimension == 1)
            {

                int idx = _notEqMask.IndexOf(true);
                if (nearestPoint[idx] < _min[idx])
                {
                    if (notUpper[idx])
                        return null;
                    else
                        nearestPoint[idx] = _min[idx];
                }
                else if (nearestPoint[idx] > _max[idx])
                {
                    nearestPoint[idx] = _max[idx];
                }

                return nearestPoint;
            }

            var nearestRealPoint = NearestRealPoint(AlgorithmHelper.GetArrayByMask(nearestPoint, _notEqMask), AlgorithmHelper.GetArrayByMask(weights, _notEqMask), AlgorithmHelper.GetArrayByMask(notUpper, _notEqMask));
            return ConvertFromReal(nearestRealPoint);
        }

        //private List<double[]> GetSingleDimensionGrid(int[] stepsCount)
        //{
        //    var result = new List<double[]>();
        //    double stepLen = (_realMax[0] - _realMin[0]) / (stepsCount[0] - 1);
        //    for (int i = 0; i < stepsCount[0]; i++)
        //    {
        //        result.Add(new double[] { stepLen * i });
        //    }
        //    return result;
        //}

        //private List<double[]> GetMultiDimensionGrid(int[] stepsCount)
        //{
        //    var result = new List<double[]>();
        //    double[] stepLen = _realMax.Zip(_realMin, (x, y) => x - y).Zip(stepsCount, (x,y) => x / (y - 1)).ToArray();

        //    int[][] truthTable = Combinatorics.TruthTable(stepsCount);
        //    List<double[]> pointsInHull = (new double[truthTable.Count()][]).ToList();
        //    Parallel.For(0, pointsInHull.Count(), (idx) =>
        //    { 
        //        double[] point = stepLen.Zip(truthTable[idx], (x,y) => x * y).ToArray();
                
        //        if (IsRealPointInConvexHull(point))
        //        {
        //            pointsInHull[idx] = point;
        //        }
        //    });

        //    pointsInHull = pointsInHull.Where(x => x != null).ToList();
            
        //    return pointsInHull;
        //}

        private bool IsRealPointInConvexHull(double[] realPoint)
        {
            return !_systemOfRealInequalities.Item1.Dot(realPoint).Subtract(_systemOfRealInequalities.Item2).Any(x => x > EPS);
        }
        
        /// <summary>
        /// Получает точку на границе выпуклой оболочки, ближайшую к заданной точке, решая задачу квадратичного программирования
        /// </summary>
        /// <param name="realPoint">Точка, не принадлежащая к выпуклой оболочке</param>
        /// <returns></returns>
        private double[] NearestRealPoint(double[] realPoint, double[] weights, bool[] notUpper)
        {
            var faces = _convexHull.Faces.ToList();
            int fCount = faces.Count();
            // Получаем свободный член уравнений поверхностей
            double[] facesD = faces.Select(face =>
            {
                return face.Normal.Zip(face.Vertices[0].Position, (x, y) => x * y).Sum();
            }).ToArray();

            double[,] Q = new double[_realDimension, _realDimension];
            double[] d = new double[_realDimension];
            for (int i = 0; i < _realDimension; i++)
            {
                Q[i, i] = weights[i] * weights[i] * 2.0;
                d[i] = weights[i] * weights[i] * (-2.0 * realPoint[i]);
            }

            bool hasUpperBound = notUpper.Any(x => x);
            int constrCount = hasUpperBound ? fCount + _realDimension : fCount;
            double[,] A = new double[constrCount, _realDimension];
            double[] b = new double[constrCount];
            for (int i = 0; i < fCount; i++)
            {
                b[i] = -facesD[i];
                for (int j = 0; j < _realDimension; j++)
                {
                    A[i, j] = -faces[i].Normal[j];
                }
            }
            if (hasUpperBound)
                for (int i = fCount; i < fCount + _realDimension; i++)
                {
                    int dim = i - fCount;
                    b[i] = notUpper[dim] ? -realPoint[dim] : -_realMax[dim] * 0.1;
                    A[i, dim] = -1.0;
                }

            return AlgorithmHelper.SolveQP(Q, d, A, b, 0);
        }

        private double[] ConvertFromReal(double[] realPoint)
        {
            if (realPoint == null)
                return null;

            int idx = 0;
            double[] point = new double[Dimension];
            for (int i = 0; i < Dimension; i++)
            {
                if (_notEqMask[i])
                {
                    point[i] = realPoint[idx];
                    idx++;
                }
                else
                    point[i] = _max[i];
            }
            return point;
        }
        
    }
}
