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
        
        public double[] FindNearestInteriorPoint(double[] point)
        {
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
                    nearestPoint[idx] = _min[idx];
                }
                else if (nearestPoint[idx] > _max[idx])
                {
                    nearestPoint[idx] = _max[idx];
                }

                return nearestPoint;
            }

            var nearestRealPoint = NearestRealPoint(AlgorithmHelper.GetArrayByMask(nearestPoint, _notEqMask));
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
        private double[] NearestRealPoint(double[] realPoint)
        {
            var faces = _convexHull.Faces.ToList();

            // Получаем свободный член уравнений поверхностей
            double[] facesD = faces.Select(face =>
            {
                return face.Normal.Zip(face.Vertices[0].Position, (x, y) => x * y).Sum();
            }).ToArray();

            double[,] Q = new double[_realDimension, _realDimension];
            double[] d = new double[_realDimension];
            for (int i =  0; i < _realDimension; i++)
            {
                for (int j = 0; j < _realDimension; j++)
                    Q[i, j] = i != j ? 0.0 : 2.0;

                d[i] = -2.0 * realPoint[i];
            }

            double[,] A = new double[faces.Count(), _realDimension];
            double[] b = new double[faces.Count()];
            for (int i = 0; i < faces.Count(); i++)
            {
                b[i] = - facesD[i];
                for (int j = 0; j < _realDimension; j++)
                {
                    A[i, j] = - faces[i].Normal[j];
                }
            }

            var solver = new GoldfarbIdnani(new QuadraticObjectiveFunction(Q, d), A, b, 0);
            solver.Minimize();
            if (solver.Status != GoldfarbIdnaniStatus.Success) throw new Exception();
            return solver.Solution;
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
