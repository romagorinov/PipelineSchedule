using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Algorithms
{
    /// <summary>
    /// Класс дискретного расписания с одинаковой длинной интервалов
    /// </summary>
    public class DiscreteSchedule
    {
        private double[] _schedule;

        public int Size
        {
            get { return _schedule.Length; }
        }

        public double this[int key]
        {
            get { return _schedule[key]; }
            set { _schedule[key] = value; }
        }

        public double Sum
        {
            get { return _schedule.Sum(); }
        }

        public double Max
        {
            get { return _schedule.Max(); }
        }

        public double Min
        {
            get { return _schedule.Min(); }
        }

        private DiscreteSchedule() { }

        public DiscreteSchedule(int size, double defaultValues = 0)
        {
            if (size < 1)
            {
                throw new ArgumentOutOfRangeException("Размер расписания не должен быть меньше 1");
            }
            _schedule = Enumerable.Repeat(defaultValues, size).ToArray<double>();
        }

        public DiscreteSchedule(DiscreteSchedule p)
        {
            if (p == null)
            {
                throw new ArgumentNullException();
            }

            _schedule = p._schedule.Clone() as double[];
        }

        public DiscreteSchedule(double[] array)
        {
            if (array == null)
            {
                throw new ArgumentNullException();
            }

            _schedule = array.Clone() as double[];
        }

        /// <summary>
        /// Разбивает каждый интервал расписания на k элементов
        /// </summary>
        /// <param name="k">Количество элементов, на которые разбивается каждый элемент расписания</param>
        /// <returns>Дискретизированное расписание</returns>
        public DiscreteSchedule Split(int k)
        {
            if (k < 2)
            {
                throw new ArgumentOutOfRangeException("Нельзя разделить меньше чем на 2 части");
            }

            DiscreteSchedule result = new DiscreteSchedule(Size * k);

            for (int i = 0; i < Size; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    result._schedule[i * k + j] = _schedule[i] / k;
                }
            }

            return result;
        }

        /// <summary>
        /// Объединяет каждые k элементов в расписании
        /// </summary>
        /// <param name="k">Количество объединияемых подряд элементов</param>
        /// <returns>Новое расписание с объединенными интервалами</returns>
        public DiscreteSchedule UnSplit(int k)
        {
            if (k < 2)
            {
                throw new ArgumentOutOfRangeException("Нельзя объединять больше чем 2 части");
            }
            if (Size % k != 0)
            {
                throw new ArgumentException("Размер расписания не кратен величине объединения");
            }

            DiscreteSchedule result = new DiscreteSchedule(Size / k);

            for (int i = 0; i < Size; i++)
            {
                result._schedule[i / k] += _schedule[i];
            }

            return result;
        }

        /// <summary>
        /// Линейная комбинация расписаний где каждый элемент нового расписания равен p1+k*p2
        /// </summary>
        /// <param name="p1">Первое расписание</param>
        /// <param name="p2">Второе расписание</param>
        /// <param name="k">Коэффициент</param>
        /// <returns>Линейная комбинация расписаний</returns>
        public static DiscreteSchedule Combine(DiscreteSchedule p1, DiscreteSchedule p2, int k)
        {
            if (p1 == null || p2 == null)
            {
                throw new ArgumentNullException();
            }

            if (p1.Size != p2.Size)
            {
                throw new ArgumentException("Не совпадают размеры расписаний");
            }

            DiscreteSchedule result = new DiscreteSchedule();

            result._schedule = p1._schedule.Zip(p2._schedule, (x, y) => x + k * y).ToArray();

            return result;
        }

        public static DiscreteSchedule operator -(DiscreteSchedule p1, DiscreteSchedule p2)
        {
            return Combine(p1, p2, -1);
        }

        public static DiscreteSchedule operator +(DiscreteSchedule p1, DiscreteSchedule p2)
        {
            return Combine(p1, p2, 1);
        }

        public bool HasNegative()
        {
            return _schedule.Min() < 0;
        }

        public double[] ToArray()
        {
            return _schedule.ToArray();
        }

        public List<double> ToList()
        {
            return _schedule.ToList();
        }

        /// <summary>
        /// Проверяем, правильно ли заданы начальный и конечный индекс для данного расписания
        /// </summary>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        public void CheckInterval(int startIndex, ref int endIndex)
        {
            if (startIndex < 0 || startIndex > Size - 1)
            {
                throw new ArgumentOutOfRangeException("Неверный начальный индекс");
            }

            if (endIndex < 0)
            {
                endIndex = Size;
            }

            if (endIndex <= startIndex)
            {
                throw new ArgumentOutOfRangeException("Конечный индекс должен быть больше начального");
            }
        }

        /// <summary>
        /// Сумма значений интервала
        /// </summary>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        /// <returns></returns>
        public double SumInterval(int startIndex = 0, int endIndex = -1)
        {
            CheckInterval(startIndex, ref endIndex);

            double result = 0;
            for (int i = startIndex; i < endIndex; i++)
            {
                result += _schedule[i];
            }
            return result;
        }

        /// <summary>
        /// Заполнить интервал определенным значением
        /// </summary>
        /// <param name="value">Значение, которым заполняется интервал</param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        public void FillInterval(double value, int startIndex = 0, int endIndex = -1)
        {
            CheckInterval(startIndex, ref endIndex);

            for (int i = startIndex; i < endIndex; i++)
            {
                _schedule[i] = value;
            }
        }

        /// <summary>
        /// Сортировка
        /// </summary>
        /// <param name="ASC">true - по возрастанию, false - по убыванию</param>
        public void Sort(bool ASC = true)
        {
            Array.Sort(_schedule);
            if (!ASC)
            {
                Array.Reverse(_schedule);
            }
        }

        public int[] SortIndexesOnly(bool ASC = true)
        {
            int[] indexes = new int[Size];

            for (int i = 0; i < Size; i++)
            {
                indexes[i] = i;
            }
            Array.Sort(_schedule.Clone() as double[], indexes);

            return indexes;
        }

        /// <summary>
        /// Получает уникальные значения расписания и подсчитывает количество каждого из значений
        /// </summary>
        /// <returns></returns>
        public Dictionary<double, int> GetUniqueWithCount()
        {
            return _schedule.GroupBy(x => x).ToDictionary(x => x.Key, x => x.Count());
        }

        /// <summary>
        /// Перекрытие поменшему значению
        /// </summary>
        /// <param name="p2"></param>
        /// <returns></returns>
        public DiscreteSchedule Overlap(DiscreteSchedule p2)
        {
            DiscreteSchedule result = new DiscreteSchedule(Size);

            for (int i = 0; i < Size; i++)
            {
                result[i] = _schedule[i] < p2._schedule[i] ? _schedule[i] : p2._schedule[i];
            }

            return result;
        }

        public DiscreteSchedule GetIntegral()
        {
            DiscreteSchedule result = new DiscreteSchedule(this);

            for (int i = 1; i < result.Size; i++)
            {
                result[i] = result[i] + result[i - 1];
            }

            return result;
        }

        public DiscreteSchedule GetReverseIntegral()
        {
            DiscreteSchedule result = new DiscreteSchedule(this);

            for (int i = result.Size - 2; i >= 0; i--)
            {
                result[i] = result[i] + result[i + 1];
            }

            return result;

        }

        public DiscreteSchedule GetReverse()
        {
            DiscreteSchedule result = new DiscreteSchedule(this);
            result._schedule.Reverse();
            return result;
        }

        public DiscreteSchedule GetPart(int startIndex, int endIndex = -1)
        {
            CheckInterval(startIndex, ref endIndex);

            DiscreteSchedule result = new DiscreteSchedule(endIndex - startIndex);
            for (int i = startIndex; i < endIndex; i++)
            {
                result[i - startIndex] = _schedule[i];
            }
            return result;
        }

        public static DiscreteSchedule GetDiffereceIntegral(double startValue, DiscreteSchedule p1, DiscreteSchedule p2)
        {
            if (p1 == null | p2 == null)
            {
                throw new ArgumentNullException();
            }
            if (p1.Size != p2.Size)
            {
                throw new ArgumentException();
            }

            DiscreteSchedule result = new DiscreteSchedule(p1.Size + 1);
            result[0] = startValue;
            for (int i = 1; i < result.Size; i++)
            {
                result[i] = result[i - 1] + p1[i - 1] - p2[i - 1];
            }

            return result;
        }

        public static DiscreteSchedule GetDiffereceIntegral(double startValue, List<double> p1, List<double> p2)
        {
            return DiscreteSchedule.GetDiffereceIntegral(startValue, new DiscreteSchedule(p1.ToArray()), new DiscreteSchedule(p2.ToArray()));
        }

        public bool LowerThanPart(double value, int startIndex, int endIndex = -1)
        {
            CheckInterval(startIndex, ref endIndex);

            for (int i = startIndex; i < endIndex; i++)
            {
                if (value > _schedule[i])
                {
                    return false;
                }
            }

            return true;
        }

        public int FirstIndexOfOutOfRange(double minValue, double maxValue)
        {
            for(int i = 0; i < Size; i++)
            {
                if (_schedule[i] < minValue || _schedule[i] > maxValue)
                {
                    return i;
                }
            }

            return -1;
        }
    }
}    