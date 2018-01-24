using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pipeline
{
    namespace PipelineObjects
    {
        public class IOObject: NamedObject
        {
            public string Name
            {
                get;
                private set;
            }

            public IOObject(string name)
            {
                Name = name;
            }
        }
    }
}
