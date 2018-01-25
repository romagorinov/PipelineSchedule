using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pipeline
{
    namespace PipelineObjects
    {
        public class SubSection: NamedObject
        {
            public NamedObject Source
            {
                get;
                private set;
            }

            public NamedObject Target
            {
                get;
                private set;
            }
            
            public double MaxFlows
            {
                get;
                private set;
            }

            public string TechnologicalSectionName
            {
                get;
                private set;
            }

            public string Name
            {
                get;
                private set;
            }

            public SubSection(string name, NamedObject source = null, NamedObject target = null, string tsName = null, double maxFlow = 10000.0)
            {
                Name = name;
                Source = source;
                Target = target;
                TechnologicalSectionName = tsName;
                MaxFlows = maxFlow;
            }
        }
    }
}
