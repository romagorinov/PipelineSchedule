using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pipeline
{
    namespace PipelineObjects
    {
        public class Reservoir: NamedObject
        {
            public string Name
            {
                get;
                private set;
            }

            #region Технические ограничения

            public double Volume
            {
                get;
                private set;
            }

            #endregion
            
            #region Конструкторы

            public Reservoir(string name, double volume)
            {
                Name = name;
                Volume = volume;
            }

            #endregion

        }
    }
}
