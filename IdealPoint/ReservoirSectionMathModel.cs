using Algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class ReservoirSectionMathModel
    {
        #region Поля

        double _reservoirVolume;
        ISection _inputSection;
        ISection _outputSection;
        double _oilStartVolume;
        int _period;

        List<double[]> _tempInputSchedule;
        List<double[]> _tempOutputSchedule;
        List<double> _tempPumpSchedule;

        #endregion

        #region Свойства



        #endregion

        #region Кострукторы

        public ReservoirSectionMathModel(double reservoirVolume, double oilStartVolume, ISection inputSection, ISection outputSection, int period)
        {
            _reservoirVolume = reservoirVolume;
            _inputSection = inputSection;
            _outputSection = outputSection;
            _oilStartVolume = oilStartVolume;
            _period = period;
        }

        #endregion

        #region Методы

        public bool SetTempParams(List<double[]> tempInputSchedule, List<double[]> tempOutputSchedule, List<double> tempPumpSchedule)
        {
            _tempInputSchedule = tempInputSchedule;
            _tempOutputSchedule = tempOutputSchedule;
            _tempPumpSchedule = tempPumpSchedule;

            return true;
        }

        List<double> GetReservoirSchedule()
        {
            List<double> result = new List<double>(_period) { _oilStartVolume };
            for (int i = 0; i < _period; i++)
            {
                result.Add(result.Last() + _tempInputSchedule[i].Last() - _tempOutputSchedule[i].First() + _tempPumpSchedule[i]);
            }
            return result;
        }

        List<double> 

        #endregion
    }
}
