using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pipeline
{
    namespace PipelineObjects
    {
        public class Pipe: NamedObject
        {
            public NamedObject Source
            {
                get;
                set;
            }

            public NamedObject Target
            {
                get;
                set;
            }

            public List<double> MaxFlows
            {
                get;
                set;
            }

            public string TechnologicalSectionName
            {
                get;
                set;
            }

            public string Name
            {
                get;
                private set;
            }

            public Pipe(string name)
            {
                Name = name;
            }
        }
    }
}
