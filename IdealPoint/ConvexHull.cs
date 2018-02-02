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
        public const double EPS = 1E-10;

        protected List<double[]> _basePoints;
        protected double[] _max, _min;
        protected bool[] _notEqMask;

        protected List<double[]> _realPoints;
        protected int _realDimension;
        protected double[] _realMax, _realMin;

        protected MIConvexHull.ConvexHull<MIConvexHull.DefaultVertex, MIConvexHull.DefaultConvexFace<MIConvexHull.DefaultVertex>> _convexHull;

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
                    return new List<double[]> { _max };
                else if (_realDimension == 1)
                    return new List<double[]> { _min, _max };
                else
                    return _convexHull.Points.Select(p => ConvertFromReal(p.Position.Select(x => x).ToArray())).ToList();
            }
        }

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
                {
                    _convexHull = MIConvexHull.ConvexHull.Create(_realPoints);
                }
            }
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
                    if (A.Dot(potentialV).Subtract(b).All(x => x <= 0))
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
        
        public double[] FindNearestInteriorPoint(double[] point, bool notUpper, double[] weights = null)
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
                    if (notUpper)
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

            var nearestRealPoint = NearestRealPoint(AlgorithmHelper.GetArrayByMask(nearestPoint, _notEqMask), AlgorithmHelper.GetArrayByMask(weights, _notEqMask), notUpper);
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
            foreach (var face in _convexHull.Faces)
            {
                // Вектор от тестируемой точки
                double[] vec = AlgorithmHelper.GetDifOfVectors(face.Vertices[0].Position, realPoint);

                // Проверяем сонаправленность с нормалью
                if (vec.Dot(face.Normal) < -EPS)
                {
                    return false;
                }
            }

            return true;
        }
        
        /// <summary>
        /// Получает точку на границе выпуклой оболочки, ближайшую к заданной точке, решая задачу квадратичного программирования
        /// </summary>
        /// <param name="realPoint">Точка, не принадлежащая к выпуклой оболочке</param>
        /// <returns></returns>
        private double[] NearestRealPoint(double[] realPoint, double[] weights, bool notUpper)
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

            int constrCount = notUpper ? fCount + _realDimension : fCount;
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
            if (notUpper)
                for (int i = fCount; i < fCount + _realDimension; i++)
                {
                    b[i] = -realPoint[i - fCount];
                    A[i, i - fCount] = -1.0;
                }

            return AlgorithmHelper.SolveQP(Q, d, A, b, 0);
        }

        private double[] ConvertFromReal(double[] realPoint)
        {
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
