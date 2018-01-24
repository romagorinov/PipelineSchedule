using Algorithms;
using Pipeline.PipelineObjects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AntColony
{
    public class Ant
    {
        private List<Node> _way = new List<Node>();
        private List<Edge> _edge_way = new List<Edge>();
        private Dictionary<object, double> _length;
        private Random _rnd = new Random();
        private Dictionary<object, double> _target;
        private Dictionary<string, int> _regime_change_count = new Dictionary<string, int>();
        private int _sumChangeCount = 0;

        public static double MAX_VOLUME_ERROR = 0.5;
        public static int REGIME_CHANGE_MAXIMUM_COUNT = 3;

        public static double ALPHA = 2;
        public static double BETA = 2.5;

        public List<Edge> Way
        {
            get
            {
                return _edge_way;
            }
        }
                
        public Ant(Node startNode, Dictionary<object, double> initialLength, Dictionary<object, double> target)
        {
            _target = target;
            _length = initialLength;
            _way.Add(startNode);
            foreach(var r in startNode.State.Regimes)
            {
                _regime_change_count.Add(r.TechnologicalSectionName, 0);
            }
        }

        public List<Node> Run()
        {
            while (_way.Last().Outputs.Count() != 0)
            {
                Edge edge = SelectEdge();
                foreach(var len in edge.Transition.PipelineObjectChanges)
                {
                    _length[len.Key] += len.Value;
                }

                foreach (var s in edge.Transition.RegimeChanges)
                {
                    _regime_change_count[s]++;
                    _sumChangeCount++;
                }

                if ((_way.Count() - 1) % 24 == 0)
                {
                    foreach (var x in _regime_change_count.Keys.ToList())
                    {
                        _regime_change_count[x] = 0;
                    }
                }

                _way.Add(edge.Target);
                _edge_way.Add(edge);
            }
            return _way;
        }
        
        private Edge SelectEdge()
        {
            Node currentNode = _way.Last();

            // Проверяем ограничения задачи
            List<Edge> avaliableEdges = new List<Edge>();
            foreach (var edge in currentNode.Outputs)
            {
                // Проверяем все резервуары на переливы недоливы и что не превысили максимальное количество нефти по поставщикам и потребителям
                bool badPoints = edge.Transition.PipelineObjectChanges.Any(x =>
                {
                    if (typeof(Reservoir) == x.Key.GetType())
                    {
                        double reservoirLevel = _length[x.Key] + x.Value;

                        if (reservoirLevel < 0 || reservoirLevel > (x.Key as Reservoir).Volume)
                        {
                            return true;
                        }
                    }

                    if ((typeof(IOObject) == x.Key.GetType()) && (Math.Abs(_length[x.Key] + x.Value) > Math.Abs(_target[x.Key]) + MAX_VOLUME_ERROR))
                    {
                        return true;
                    }

                    return false;
                });

                if (!badPoints)
                    avaliableEdges.Add(edge);
            }

            // Из оставшихся переходов выбираем на базе уровня феромона и эвристики
            List<double> probability = new List<double>();
            double sum = 0;
            foreach(var edge in avaliableEdges)
            {
                double pheromone_alpha = Math.Pow(edge.Pheromone, ALPHA);
                                
                double heuristic_beta = Math.Pow(GetEdgeHeuristic(edge), BETA);

                double pr = pheromone_alpha * heuristic_beta;
                probability.Add(pr);
                sum += pr;
            }

            // Вероятности выбора
            probability = probability.Select(x => x / sum).ToList();

            return avaliableEdges[AlgorithmHelper.GetRouletIndex(probability)];
        }

        private double GetEdgeHeuristic(Edge edge)
        {
            // Ребро приклекательно, если режим не изменяется
            double sameness = 1 - edge.Transition.DiametralDifference;

            // Ребро привлекательно, если режим близок к среднему
            double averageSameness = 1 - _way.First().Outputs.First(x => x.Target.State == edge.Target.State).Transition.DiametralDifference;

            double heuristic =  0.4 * sameness +  0.7 * averageSameness;

            int bad = 0;
            foreach(var r in edge.Target.State.Regimes)
            {
                foreach(var p in  r.Q)
                {
                    object o = null;
                    if (typeof(IOObject) == p.Key.Target.GetType())
                    {
                        o = p.Key.Target;
                    }
                    else if (typeof(IOObject) == p.Key.Source.GetType())
                    {
                        o = p.Key.Source;
                    }

                    if (o != null && p.Value == 0 && Math.Abs(_length[o]) <  MAX_VOLUME_ERROR + Math.Abs(_target[o]))
                    {
                        bad++;
                    }
                }
            }

            return heuristic / Math.Exp(bad);
        }

        public double GetVolumeCriteria()
        { 
            double sumVolumeDifference = 0;
            foreach(var o in _target)
            {
                sumVolumeDifference += Math.Abs(o.Value - _length[o.Key]);
            }

            return sumVolumeDifference;
        }

        public int GetChangesCriteria()
        {
            return _sumChangeCount;
        }

        public Dictionary<Pipe, DiscreteSchedule> GetSchedule()
        {
            Dictionary<Pipe, DiscreteSchedule> result = new Dictionary<Pipe, DiscreteSchedule>();
            foreach(var r in _way.First().State.Regimes)
            {
                foreach(var p in r.Q)
                {
                    result.Add(p.Key, new DiscreteSchedule(_way.Count() - 1));
                }
            }

            for (int i = 1; i < _way.Count(); i++)
            {
                foreach (var r in _way[i].State.Regimes)
                {
                    foreach (var p in r.Q)
                    {
                        result[p.Key][i - 1] = p.Value;
                    }
                }
            }

            return result;
        }
    }
}
