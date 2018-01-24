using Algorithms;
using Pipeline;
using Pipeline.Needles;
using Pipeline.PipelineObjects;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AntColony
{
    public class Runner
    {
        public static double TEST_VOLUME = 1500.0;

        public static double EP = 0.6;
        public static double EP_1 = 1 - EP;
        public static double Q = 1;
        public static int ANTS_COUNT = 50;
        public static int STEPS_COUNT = 500;
        public static int LAYERS_COUNT = 744;

        public Node startNode;
        
        public PipelineSystem system = new PipelineSystem();

        public List<Edge> allEdges = new List<Edge>();


        private Dictionary<object, double> GetTargets()
        {
            return new Dictionary<object, double>() {
                { system.GetPointByName("input"), -TEST_VOLUME }, { system.GetPointByName("output"), TEST_VOLUME }
            };
        }

        private Dictionary<object, double> GetInitialStates()
        {
            Dictionary<object, double> initialState = new Dictionary<object, double>();
            foreach (var o in system.points)
            {
                initialState.Add(o, 0);
            }
            return initialState;
        }

        public Runner()
        {
            // Фэйковая нода со среним расходом по кажой трубе
            startNode = new Node( new PipelineSystem.SystemState(system,
                new List<Regime>() {
                    new Regime("1", "fake-1", 10, new Dictionary<Pipe, double>() { { system.pipes.First(x => x.Name == "inputPipe"), TEST_VOLUME / LAYERS_COUNT } } ),
                    new Regime("2", "fake-2", 10, new Dictionary<Pipe, double>() { { system.pipes.First(x => x.Name == "outputPipe"), TEST_VOLUME / LAYERS_COUNT } } )
                }
            ));

            List<Node> oldLayer = system.States.Select(x => new Node(x)).ToList();

            foreach (var n in oldLayer)
            {
                Edge edge = new Edge(startNode, n, new PipelineSystem.SystemTransition(system, startNode.State, n.State));
                allEdges.Add(edge);
                startNode.Outputs.Add(edge);
            }
            
            for (int i = 1; i < LAYERS_COUNT; i++)
            {
                List<Node> newLayer = system.States.Select(x => new Node(x)).ToList();
                foreach (var oldNode in oldLayer)
                {
                    foreach(var newNode in newLayer)
                    {
                        Edge edge = new Edge(oldNode, newNode, system.Transitions.First(x => x.StateFrom == oldNode.State && x.StateTo == newNode.State));
                        allEdges.Add(edge);
                        oldNode.Outputs.Add(edge);
                    }
                }

                oldLayer = newLayer;
            }
        }

        public Ant RunAnt(Node startNode)
        { 
            Ant ant = new Ant(startNode, GetInitialStates(), GetTargets());

            ant.Run();

            return ant;
        }

        public Tuple<DiscreteSchedule, DiscreteSchedule> Run()
        {
            Ant bestAnt = null;
            double bestCriteria = double.MaxValue;
            for (int age = 0; age < STEPS_COUNT; age++)
            {
                ConcurrentBag<Ant> concurrent_ants = new ConcurrentBag<Ant>();
                // Пустили муравьев
                /*for (int i = 0; i< ANTS_COUNT; i++)
                {
                    concurrent_ants.Add(RunAnt(startNode));
                }*/
                Parallel.For(0, ANTS_COUNT, iter =>
                {
                    concurrent_ants.Add(RunAnt(startNode));
                });

                List<Ant> ants = concurrent_ants.ToList();
                List<double> changeCriteria = new List<double>(), volumeCriteria = new List<double>();
                foreach (var a in ants)
                {
                    changeCriteria.Add(a.GetChangesCriteria());
                    volumeCriteria.Add(a.GetVolumeCriteria());
                }
                // Получили значения критерия для каждого пути
                double maxRegimeChanges = Ant.REGIME_CHANGE_MAXIMUM_COUNT * LAYERS_COUNT / 24;
                List<double> sumCriteria = changeCriteria.Zip(volumeCriteria, (x, y) => 0.1 * x / maxRegimeChanges + y / (2* TEST_VOLUME) ).ToList();

                // Испарили феромоны
                foreach(var e in allEdges)
                {
                    e.Pheromone *= EP_1;
                }

                for (int i = 0; i < ants.Count(); i++)
                {
                    Ant a = ants[i];
                    foreach(var e in a.Way)
                    {
                        e.Pheromone += Q / (sumCriteria[i] + 1);
                    }
                }

                double currentCriteria = sumCriteria.Min();
                if (currentCriteria < bestCriteria)
                {
                    bestCriteria = currentCriteria;
                    bestAnt = ants[sumCriteria.IndexOf(currentCriteria)];
                }

                ants.Clear();
            }
            
            Dictionary<Pipe, DiscreteSchedule> s = bestAnt.GetSchedule();
            return new Tuple<DiscreteSchedule, DiscreteSchedule>(s.First(x => x.Key.Name == "inputPipe").Value, s.First(x => x.Key.Name == "outputPipe").Value);
        }
    }
}
