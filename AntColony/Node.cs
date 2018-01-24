using Pipeline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AntColony
{
    public class Node
    {
        public List<Edge> Outputs
        {
            get;
            private set;
        }

        public PipelineSystem.SystemState State
        {
            get;
            private set;
        }

        public Node(PipelineSystem.SystemState state)
        {
            Outputs = new List<Edge>();
            State = state;
        }
    }
}
