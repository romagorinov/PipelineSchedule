using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

using AH = Algorithms.AlgorithmHelper;

namespace Algorithms
{
    public class ReservoirBalanceAlgorithm
    {
        public interface TechnologicalSectionAlgorithm
        {
            List<double[]> GetSchedule(double[] targets, List<double[]> fixValues = null);

            List<double[]> GetContinuousSchedule(double[] targets, List<double[]> fixValues = null);

            bool IsRegimeAvaliableOnInterval(double[] regime, int interval);

            List<double[]> Regimes
            {
                get;
            }

            List<double[]> NormRegimes
            {
                get;
            }

            int Period
            {
                get;
            }

            int Dimension
            {
                get;
            }

            List<double[]> UpperBounds
            {
                get;
            }

            double[] MaxFlows
            {
                get;
            }
        }

        #region Параметры объектов

        private TechnologicalSectionAlgorithm _inAlgorithm;
        private TechnologicalSectionAlgorithm _outAlgorithm;
        private Tuple<int, int> _connection;
        private Tuple<bool[], bool[]> _constraints;
        private double _bufferVolume;
        private int _period;
        private int _inDimension;
        private int _outDimension;

        #endregion

        #region Параметры алгоритма
        
        public int GENERATION_NUMBER = 2000;
        public int MINIMUM_WORK_TIME = 1;

        #endregion

        #region Дополнительные расчеты

        private List<Tuple<double[], double[]>> _regimesPairs;
        private List<Tuple<double[], double[]>> _normRegimesPairs;
        private List<List<double>> _pairsDifferences;
        private List<double> _reservoirDifferences;
        private List<List<Tuple<double[], double[]>>> _avaliablePairsOnIntervals;

        #endregion

        #region Конструкторы

        public ReservoirBalanceAlgorithm(TechnologicalSectionAlgorithm inAlgorithm, TechnologicalSectionAlgorithm outAlgorithm, double bufferVolume, Tuple<int, int> connection, Tuple<bool[], bool[]> constraints)
        {
            _inAlgorithm = inAlgorithm;
            _outAlgorithm = outAlgorithm;

            _period = _inAlgorithm.Period;
            _inDimension = _inAlgorithm.Dimension;
            _outDimension = _outAlgorithm.Dimension;

            _connection = connection;
            _constraints = constraints;

            _bufferVolume = bufferVolume;

            _regimesPairs = AH.CartesianProduct(new List<List<double[]>>() { _inAlgorithm.Regimes, _outAlgorithm.Regimes }).Select(x => new Tuple<double[], double[]>(x.First(), x.Last())).ToList();
            _normRegimesPairs = AH.CartesianProduct(new List<List<double[]>>() { _inAlgorithm.NormRegimes, _outAlgorithm.NormRegimes }).Select(x => new Tuple<double[], double[]>(x.First(), x.Last())).ToList();

            _reservoirDifferences = _regimesPairs.Select(x => x.Item1[_connection.Item1] - x.Item2[_connection.Item2]).ToList();

            _pairsDifferences = new List<List<double>>();
            int num = _regimesPairs.Count();
            for (int i = 0; i < num; i++)
            {
                List<double> l = new List<double>();
                _pairsDifferences.Add(l);
                for (int j = 0; j < num; j++)
                {
                    l.Add(GetRegimesDifference(_normRegimesPairs[i], _normRegimesPairs[j]));
                }
            }

            _avaliablePairsOnIntervals = new List<List<Tuple<double[], double[]>>>();
            for (int i = 0; i < _period; i++)
            {
                _avaliablePairsOnIntervals.Add(_regimesPairs.Where(x => IsRegimeAvaliable(x, i)).ToList());
            }
        }

        #endregion

        #region Начальные параметры расчета
        
        public enum AlgorithmType
        {
            GRAPH,
            DINAMIC
        }

        public struct InitialValues
        {
            public double[] inTargets;
            public double[] outTargets;
            public double startVolume;
            public double minVolume;
            public double maxVolume;
            public List<double> pumpsSchedule;
            
            public AlgorithmType algorithmType;

            public GreedyNode initRootNode;

            public int maxFragmentation;
        }

        #endregion

        #region Методы

        public Tuple<List<double[]>, List<double[]>> FormSchedule(ref InitialValues initVals, Action<string> action = null)
        {
            List<double[]> inFixValues = AH.CreateListOfArrays(_period, _inDimension, -1.0),
                outFixValues = AH.CreateListOfArrays(_period, _outDimension, -1.0),
                inSchedule = _inAlgorithm.GetSchedule(initVals.inTargets, inFixValues),
                outSchedule = _outAlgorithm.GetSchedule(initVals.outTargets, outFixValues);

            DiscreteSchedule reservoirSchedule = GetReservoirSchedule(inSchedule, outSchedule, initVals.startVolume, initVals.pumpsSchedule);
            int idx = reservoirSchedule.FirstIndexOfOutOfRange(initVals.minVolume, initVals.maxVolume);
            if (idx == -1)
            {
                return new Tuple<List<double[]>, List<double[]>>(inSchedule, outSchedule);
            }

            List<Tuple<List<double[]>, List<double[]>>> schedules;
            switch (initVals.algorithmType)
            {
                case AlgorithmType.GRAPH:
                    schedules = GreedyRandomizeSearch(ref initVals, inSchedule, outSchedule, action);
                    break;
                default:
                    schedules = new List<Tuple<List<double[]>, List<double[]>>>() { DinamicProgramming(initVals, action) };
                    break;
            }

            action("Выбираем оптимальное на основе критериев");

            if (schedules.Count() == 0)
            {
                return null;
            }
            else if (schedules.Count() == 1)
            {
                return schedules[0];
            }
            else
            {
                List<double> criterias = GetCriteria(initVals.inTargets, initVals.outTargets, schedules);
                var best = schedules[criterias.IndexOf(criterias.Min())];
                return best;
            }
        }

        #region Жадный рандомизированный поиск

        private List<Tuple<List<double[]>, List<double[]>>> GreedyRandomizeSearch(ref InitialValues initVals, List<double[]> inSchedule, List<double[]> outSchedule, Action<string> action)
        {
            ConcurrentBag<Tuple<List<double[]>, List<double[]>>> schedulesConcurrent = new ConcurrentBag<Tuple<List<double[]>, List<double[]>>>();
            int counter = 0,
                part = 0,
                partSize = 10;
            int threadId = 0;
            Parallel.Invoke(() =>
            {
                threadId = Thread.CurrentThread.ManagedThreadId;
            });
            GreedyNode rootNode;
            if (initVals.initRootNode == null)
            {
                rootNode = new GreedyNode(inSchedule, outSchedule, this, initVals);
                rootNode.CalculateEdges();
            }
            else
            {
                rootNode = initVals.initRootNode;
            }

            Parallel.For(0, GENERATION_NUMBER, (generation, state) =>
            {
                if (rootNode.edgesList.Count() == 0)
                {
                    state.Break();
                }

                GreedyNode currentNode = rootNode;

                while (currentNode.crashIndex != -1)
                {
                    if (!currentNode.IsCalculated)
                        currentNode.CalculateEdges();

                    GreedyEdge edge = currentNode.GetRandomEdge();

                    if (edge == null)
                        break;
                    else if (!edge.IsCalculated)
                        edge.CalculateNode();

                    currentNode = edge.end;
                }

                if (currentNode.crashIndex == -1)
                {
                    schedulesConcurrent.Add(new Tuple<List<double[]>, List<double[]>>(currentNode.inSchedule, currentNode.outSchedule));
                }

                if (currentNode.inputEdge != null)
                    currentNode.inputEdge.start.RemoveEdge(currentNode.inputEdge);
                
                Interlocked.Increment(ref counter);
                if (counter > part * partSize && Thread.CurrentThread.ManagedThreadId == threadId)
                {
                    action("Просмотрено листьев " + counter.ToString() + ", найдено решений " + schedulesConcurrent.Count().ToString());
                    part++;
                }
            });
            initVals.initRootNode = rootNode;

            return schedulesConcurrent.ToList();
        }

        public class GreedyNode
        {
            public bool IsCalculated
            {
                get;
                private set;
            }
            public bool IsFullExplored
            {
                get;
                private set;
            }

            public ReservoirBalanceAlgorithm algorithm;
            public InitialValues initVals;

            public List<double[]> inSchedule;
            public List<double[]> outSchedule;
            public DiscreteSchedule reservoirSchedule;
            public int crashIndex;
            public int sign;

            public List<GreedyEdge> edgesList = new List<GreedyEdge>();
            public GreedyEdge inputEdge = null;
            public int gridSize = 5;

            public List<Tuple<double[], double[]>> avaliablePairs;
            public List<double> pairProbabilities;
            public List<int> edgesCountForPair = new List<int>();

            private object nodeLocker = new object();

            public GreedyNode(List<double[]> inSchedule, List<double[]> outSchedule, ReservoirBalanceAlgorithm algorithm, InitialValues initVals)
            {
                this.inSchedule = inSchedule;
                this.outSchedule = outSchedule;
                this.initVals = initVals;
                this.algorithm = algorithm;

                reservoirSchedule = algorithm.GetReservoirSchedule(inSchedule, outSchedule, initVals.startVolume, initVals.pumpsSchedule);
                crashIndex = reservoirSchedule.FirstIndexOfOutOfRange(initVals.minVolume, initVals.maxVolume);

                if (crashIndex != -1)
                {
                    sign = reservoirSchedule[crashIndex] > initVals.maxVolume ? -1 : 1;
                }

                IsCalculated = false;
                IsFullExplored = false;
            }

            public void CalculateEdges()
            {
                lock (nodeLocker)
                {
                    if (IsCalculated)
                        return;

                    // Отбираем режимы
                    double[] beforeCrashInVolumes = AH.GetSumOnInterval(inSchedule, 0, crashIndex - 1),
                        beforeCrashOutVolumes = AH.GetSumOnInterval(outSchedule, 0, crashIndex - 1);
                    avaliablePairs = algorithm._avaliablePairsOnIntervals[crashIndex - 1];
                    List<int> timesOfAvaliability = avaliablePairs.Select(pair =>
                    {
                        int start = crashIndex - 1;
                        double resVol = reservoirSchedule[start];
                        double[] curInVolumes = beforeCrashInVolumes.Select(x => x).ToArray(),
                            curOutVolumes = beforeCrashOutVolumes.Select(x => x).ToArray();
                        double difference = algorithm._reservoirDifferences[algorithm._regimesPairs.IndexOf(pair)];

                        for (int i = start; i < algorithm._period; i++)
                        {
                            if (algorithm._avaliablePairsOnIntervals[i].IndexOf(pair) == -1)
                                return i;

                            resVol += (difference + initVals.pumpsSchedule[i]);
                            if (resVol < initVals.minVolume || resVol > initVals.maxVolume)
                                return i;

                            curInVolumes = curInVolumes.Zip(pair.Item1, (x, y) => x + y).ToArray();
                            for (int j = 0; j < algorithm._inDimension; j++)
                            {
                                if (algorithm._constraints.Item1[j] && initVals.inTargets[j] < curInVolumes[j])
                                    return i;
                            }

                            curOutVolumes = curOutVolumes.Zip(pair.Item2, (x, y) => x + y).ToArray();
                            for (int j = 0; j < algorithm._outDimension; j++)
                            {
                                if (algorithm._constraints.Item2[j] && initVals.outTargets[j] < curOutVolumes[j])
                                    return i;
                            }
                        }

                        return algorithm._period;
                    }).ToList();
                    avaliablePairs = avaliablePairs.Where((x, i) => timesOfAvaliability[i] == algorithm._period || timesOfAvaliability[i] - crashIndex + 1 >= algorithm.MINIMUM_WORK_TIME).ToList();
                    timesOfAvaliability = timesOfAvaliability.Where(x => x == algorithm._period || x - crashIndex + 1 >= algorithm.MINIMUM_WORK_TIME).ToList();
                    
                    if (avaliablePairs.Count() > 0)
                    {                        
                        // Вычисляем вероятности
                        Tuple<double[], double[]> currentPair = algorithm._regimesPairs.First(x => x.Item1.SequenceEqual(inSchedule[crashIndex - 1]) && x.Item2.SequenceEqual(outSchedule[crashIndex - 1]));
                        int currentIdx = algorithm._regimesPairs.IndexOf(currentPair);
                        pairProbabilities = avaliablePairs.Select((x) =>
                        {
                            int r1 = algorithm._regimesPairs.IndexOf(x);
                            return 1 / (algorithm._pairsDifferences[currentIdx][r1] + 0.1);
                        }).ToList();
                        pairProbabilities = pairProbabilities.Select(x => x / pairProbabilities.Sum()).ToList();

                        // Делаем временную сетку
                        for (int i = 0; i < avaliablePairs.Count(); i++)
                        {
                            int minT = Math.Min(algorithm._period, crashIndex - 1 + algorithm.MINIMUM_WORK_TIME),
                                maxT = timesOfAvaliability[i];
                            Tuple<double[], double[]> p = avaliablePairs[i];
                            List<int> times = AH.GetGridOnInterval(minT, maxT, gridSize);
                            if (times.Count() > 1)
                            {
                                times.RemoveAt(0);
                            }
                            edgesCountForPair.Add(times.Count());
                            times.ForEach(x => edgesList.Add(new GreedyEdge(this, avaliablePairs[i], x)));
                        }
                    }

                    reservoirSchedule = null;
                    IsCalculated = true;
                }
            }

            public GreedyEdge GetRandomEdge()
            {
                lock (nodeLocker)
                {
                    if (edgesList.Count() == 0)
                        return null;

                    int pairIndex = AH.GetRouletIndex(pairProbabilities);
                    int timeIndex = AH.RandomInstance.Next(0, edgesCountForPair[pairIndex]);

                    return edgesList.Where(x => x.regimesPair == avaliablePairs[pairIndex]).ToList()[timeIndex];
                }
            }
            
            public void RemoveEdge(GreedyEdge edge)
            {
                lock(nodeLocker)
                {
                    int edgeIdx = edgesList.IndexOf(edge);
                    if (edgeIdx == -1)
                        return;

                    edgesList.RemoveAt(edgeIdx);
                    int probabilityIdx = avaliablePairs.IndexOf(edge.regimesPair);
                    edgesCountForPair[probabilityIdx]--;

                    bool recalcProbs = false;
                    if (edgesCountForPair[probabilityIdx] == 0)
                    {
                        pairProbabilities.RemoveAt(probabilityIdx);
                        edgesCountForPair.RemoveAt(probabilityIdx);
                        avaliablePairs.RemoveAt(probabilityIdx);
                        recalcProbs = true;
                    }

                    if (edgesList.Count() == 0)
                    {
                        IsFullExplored = true;
                        if (inputEdge != null)
                            inputEdge.start.RemoveEdge(inputEdge);
                    }
                    else if (recalcProbs)
                    {
                        pairProbabilities = pairProbabilities.Select(x => x / pairProbabilities.Sum()).ToList();
                    }
                }
            }
        }

        public class GreedyEdge
        {
            public bool IsCalculated
            {
                get;
                private set;
            }

            public Tuple<double[], double[]> regimesPair;
            public int workTime;

            public GreedyNode start;
            public GreedyNode end;

            public object edgeLocker = new object();

            public GreedyEdge(GreedyNode start, Tuple<double[], double[]> regimesPair, int workTime)
            {
                this.start = start;
                this.end = null;
                this.regimesPair = regimesPair;
                this.workTime = workTime;
                IsCalculated = false;
            }

            public void CalculateNode()
            {
                lock(edgeLocker)
                {
                    if (IsCalculated) return;

                    List<double[]> inSchedule, outSchedule,
                        inFixValues = AH.CreateListOfArrays(start.algorithm._period, start.algorithm._inDimension, -1.0),
                        outFixValues = AH.CreateListOfArrays(start.algorithm._period, start.algorithm._outDimension, -1.0);

                    for (int i = 0; i < start.crashIndex - 1; i++)
                    {
                        inFixValues[i] = start.inSchedule[i].Select(x => x).ToArray();
                        outFixValues[i] = start.outSchedule[i].Select(x => x).ToArray();
                    }
                    for (int i = start.crashIndex - 1; i < workTime; i++)
                    {
                        inFixValues[i] = regimesPair.Item1.Select(x => x).ToArray();
                        outFixValues[i] = regimesPair.Item2.Select(x => x).ToArray();
                    }

                    inSchedule = start.algorithm._inAlgorithm.GetSchedule(start.initVals.inTargets, inFixValues);
                    outSchedule = start.algorithm._outAlgorithm.GetSchedule(start.initVals.outTargets, outFixValues);

                    end = new GreedyNode(inSchedule, outSchedule, start.algorithm, start.initVals);
                    end.inputEdge = this;

                    IsCalculated = true;
                }
            }
        }

        #endregion
        
        #region Динамическое программирование
        
        Tuple<List<double[]>, List<double[]>> DinamicProgramming(InitialValues initVals, Action<string> action)
        {
            // Сформируем развернутую сетку
            return null;
        }

        #endregion

        #region Вспомогательные методы

        private List<double> GetCriteria(double[] inTargets, double[] outTargets, List<Tuple<List<double[]>, List<double[]>>> schedules)
        {
            List<double> changesCriteria = AH.CreateListOfElements(schedules.Count(), 0.0),
                differenceCriteria = AH.CreateListOfElements(schedules.Count(), 0.0);

            Parallel.For(0, schedules.Count(), (scheduleNumber) =>
            {
                var schedule = schedules[scheduleNumber];

                int cangeCount = 0;
                List<double[]> s1 = schedule.Item1, s2 = schedule.Item2;
                for (int i = 1; i < _period; i++)
                {
                    if (!s1[i].SequenceEqual(s1[i - 1]))
                        cangeCount++;

                    if (!s2[i].SequenceEqual(s2[i - 1]))
                        cangeCount++;
                }

                changesCriteria[scheduleNumber] = cangeCount;

                double[] sum = AH.GetSumOnInterval(schedule.Item1, 0, _period).Concat(AH.GetSumOnInterval(schedule.Item2, 0, _period).ToArray()).ToArray();
                double[] dif = AH.GetDifOfVectors(sum, inTargets.Concat(outTargets).ToArray());
                double[] maxDif = AH.GetSumOnInterval(_inAlgorithm.UpperBounds, 0, _period).Concat(AH.GetSumOnInterval(_outAlgorithm.UpperBounds, 0, _period)).ToArray();
                dif = dif.Zip(maxDif, (x, y) => x / y).ToArray();
                differenceCriteria[scheduleNumber] = dif.Sum() / dif.Count();
            });
            
            changesCriteria = changesCriteria.Select(x => x / (_period * 2)).ToList();

            //List<double> differenceCriteria = sumDif.Select(x => x.Sum() / x.Count()).ToList();

            //// Суммы перекачанной нефти
            //List<Tuple<double[], double[]>> sum = schedules.Select((x) => 
            //    new Tuple<double[], double[]>( AH.GetSumOnInterval(x.Item1, 0, _period), AH.GetSumOnInterval(x.Item2, 0, _period))).ToList();

            //// Отличие перекачки от целей (по модулю)
            //List<double[]> sumDif = sum.Select((x) =>
            //{
            //    double[] inDif = AH.GetDifOfVectors(x.Item1, inTargets),
            //        outDif = AH.GetDifOfVectors(x.Item2, outTargets);
            //    return inDif.Concat(outDif).Select(Math.Abs).ToArray();
            //}).ToList();

            //// Нормализуем отличия
            //// sumDif = AH.NormalizeByComponents(sumDif);
            //double[] maxSum = AH.GetSumOnInterval(_inAlgorithm.UpperBounds, 0, _period).Concat(AH.GetSumOnInterval(_outAlgorithm.UpperBounds, 0, _period)).ToArray();
            //sumDif = sumDif.Select(x => x.Zip(maxSum, (y, yMax) => y / yMax).ToArray()).ToList();

            //// КРИТЕРИЙ: отличие от целей
            //List<double> differenceCriteria = sumDif.Select(x => x.Sum() / x.Count()).ToList();
            ////differenceCriteria = AH.NormalizeList(differenceCriteria);

            //// Сумма отличий на каждом интервале от среднего по модулю
            //List<double[]> difModulSum = schedules.Select((x,i) =>
            //{
            //    // Получим средние расписания
            //    List<double[]> inAvgSchedule = _inAlgorithm.GetContinuousSchedule(sum[i].Item1),
            //        outAvgSchedule = _outAlgorithm.GetContinuousSchedule(sum[i].Item2);

            //    // Вычислим отклонения
            //    List<double[]> inDeviationsByModul = AH.GetDeviations(x.Item1, inAvgSchedule).Select(y => y.Select(Math.Abs).ToArray()).ToList(),
            //        outDeviationsByModul = AH.GetDeviations(x.Item2, outAvgSchedule).Select(y => y.Select(Math.Abs).ToArray()).ToList();

            //    return AH.GetSumOnInterval(inDeviationsByModul, 0, _period).Concat(AH.GetSumOnInterval(outDeviationsByModul, 0, _period)).ToArray();
            //}).ToList();

            //// Нормализуем
            ////difModulSum = AH.NormalizeByComponents(difModulSum);
            //difModulSum = difModulSum.Select(x => x.Zip(maxSum, (y, yMax) => y / yMax).ToArray()).ToList();

            //// КРИТЕРИЙ: Неравномерность
            //List<double> uniformityCriteria = difModulSum.Select(x => x.Sum() / x.Count()).ToList();
            ////uniformityCriteria = AH.NormalizeList(uniformityCriteria);

            return schedules.Select((x,i) => 0.0 * changesCriteria[i] + 1.0 * differenceCriteria[i]).ToList();
        }

        private DiscreteSchedule GetReservoirSchedule(List<double[]> inSchedule, List<double[]> outSchedule, double startVolume, List<double> pumpsSchedule)
        {
            DiscreteSchedule result = new DiscreteSchedule(_period + 1, 0);
            result[0] = startVolume;
            for (int i = 1; i < _period + 1; i++)
            {
                result[i] = result[i-1] + inSchedule[i - 1][_connection.Item1] - outSchedule[i - 1][_connection.Item2] + pumpsSchedule[i - 1];
            }
            return result;
        }

        private bool IsRegimeAvaliable(Tuple<double[], double[]> pair, int interval)
        {
            return _inAlgorithm.IsRegimeAvaliableOnInterval(pair.Item1, interval) && _outAlgorithm.IsRegimeAvaliableOnInterval(pair.Item2, interval);
        }

        private delegate int GetTimeOfAvaliability(Tuple<double[], double[]> pair, int start);
        
        private double GetRegimesDifference(Tuple<double[], double[]> p1, Tuple<double[], double[]> p2)
        {
            return AH.GetDistance(p1.Item1.Concat(p1.Item2).ToArray(), p2.Item1.Concat(p2.Item2).ToArray());
        }

        #endregion

        #endregion
    }
}
