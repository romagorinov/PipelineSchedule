using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    public class Result
    {

        #region Ошибки

        public enum Error
        {
            NO_ERRORS
        }

        Error errors;

        #endregion
            
        #region Конструкторы

        public Result()
        {
            errors = Error.NO_ERRORS;
        }

        #endregion
    }
}
