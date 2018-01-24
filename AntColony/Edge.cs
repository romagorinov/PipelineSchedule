using Pipeline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AntColony
{
    public class Edge
    {
        public static double MAX_PHEROMONE_LEVEL = double.MaxValue;

        private double _pheromone = 0.01;

        public Node Source
        {
            get;
            private set;
        }

        public Node Target
        {
            get;
            private set;
        }

        public double Pheromone
        {
            get
            {
                return _pheromone;
            }
            set
            {
                if (value > MAX_PHEROMONE_LEVEL)
                {
                    _pheromone = MAX_PHEROMONE_LEVEL;
                }
                else
                {
                    _pheromone = value;
                }
            }
        }
        
        public PipelineSystem.SystemTransition Transition
        {
            get;
            private set;
        }

        public Edge(Node source, Node target, PipelineSystem.SystemTransition transition)
        {
            Source = source;
            Target = target;

            Transition = transition;
        }
    }
}
