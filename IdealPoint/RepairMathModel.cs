using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class RepairMathModel
    {
        double _maxInput;
        double[] _maxPumps;
        double _maxOutput;

        public double MaxInput => _maxInput;
        public double[] MaxPumps => _maxPumps;
        public double MaxOutput => _maxOutput;

        public RepairMathModel(double maxInput, double[] maxPumps, double maxOutput)
        {
            _maxInput = maxInput;
            _maxPumps = maxPumps;
            _maxOutput = maxOutput;
        }

        public static bool operator == (RepairMathModel r1, RepairMathModel r2)
        {
            return r1._maxInput == r2._maxInput
                && r1._maxPumps.SequenceEqual(r2._maxPumps)
                && r1._maxOutput == r2._maxOutput;
        }

        public static bool operator != (RepairMathModel r1, RepairMathModel r2)
        {
            return r1._maxInput != r2._maxInput
                || !r1._maxPumps.SequenceEqual(r2._maxPumps)
                || r1._maxOutput != r2._maxOutput;
        }
    }
}
