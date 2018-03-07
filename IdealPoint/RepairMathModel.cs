using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class RepairMathModel
    {
        #region Поля

        readonly double _maxInput;
        readonly double[] _maxPumps;
        readonly double _maxOutput;

        #endregion

        #region Свойства

        public double MaxInput => _maxInput;
        public double[] MaxPumps => _maxPumps;
        public double MaxOutput => _maxOutput;

        #endregion

        #region Конструкторы

        public RepairMathModel(double maxInput, double[] maxPumps, double maxOutput)
        {
            if (maxInput < 0)
                throw new Exception();
            if (maxPumps == null)
                throw new Exception();
            if (maxPumps.Count() == 0)
                throw new Exception();
            if (maxOutput < 0)
                throw new Exception();

            _maxInput = maxInput;
            _maxPumps = maxPumps.ToList().ToArray();
            _maxOutput = maxOutput;
        }

        #endregion

        #region Методы

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

        #endregion
    }
}
