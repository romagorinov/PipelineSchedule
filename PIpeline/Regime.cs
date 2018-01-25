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
            #region Свойства

            public double W
            {
                get;
                private set;
            }

            public Tuple<double, double[][]> G
            {
                get;
                set;
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

            public Regime(string tsName, string name, double w, Tuple<double, double[][]> g)
            {
                TechnologicalSectionName = tsName;
                Name = name;
                W = w;
                G = g;
            }

            #endregion
        }
    }
}
