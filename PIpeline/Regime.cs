using Pipeline.PipelineObjects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pipeline
{
    namespace Needles
    {
        public class Regime
        {
            #region Параметры режима

            public double W
            {
                get;
                private set;
            }

            public Dictionary<Pipe, double> Q
            {
                get;
                private set;
            }

            public double? this[string pipeName]
            {
                get
                {
                    if (Q.Any(x => x.Key.Name == pipeName))
                        return Q.First(x => x.Key.Name == pipeName).Value;
                    else
                        return null;
                }
            }

            public double? this[Pipe pipe]
            {
                get
                {
                    if (Q.ContainsKey(pipe))
                        return Q[pipe];
                    else
                        return null;
                }
            }

            public string Name
            {
                get;
                private set;
            }

            public string TechnologicalSectionName
            {
                get;
                private set;
            }

            #endregion

            #region Конструкторы

            public Regime(string tsName, string name, double w, Dictionary<Pipe, double> q)
            {
                TechnologicalSectionName = tsName;
                Name = name;
                W = w;
                Q = q;
            }

            #endregion
        }
    }
}
