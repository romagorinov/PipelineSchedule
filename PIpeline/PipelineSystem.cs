using Algorithms;
using Pipeline.Needles;
using Pipeline.PipelineObjects;
using Spire.Xls;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using System.Threading;
using System.Windows;

namespace Pipeline
{
    public class PipelineSystem
    {
        #region Структура

        public List<Pipe> pipes = new List<Pipe>();
        public List<NamedObject> points = new List<NamedObject>();
        public List<Regime> regimes = new List<Regime>();

        #endregion

        #region Конструкторы

        public PipelineSystem()
        {
            // Создаем систему
            CreateTransoilNorth();
        }

        #endregion

        #region Методы
        
        public void CreateTransoilNorth()
        {
            Days = 31;
            Period = Days * 24;
            double bigMaxflowValue = 1000;

            Reservoir rp1 = new Reservoir("РП1", 47.4),
                rp2 = new Reservoir("РП2", 114.3),
                rp3 = new Reservoir("РП3", 47.4);

            IOObject mainInput = new IOObject("Прием от НГДУ (1-5)"),
                mainOutput = new IOObject("Сдача в соседний ОСТ"),
                input1 = new IOObject("Подкачка 1"),
                input2 = new IOObject("Подкачка 2"),
                input3 = new IOObject("ГНПС 2"),
                output1 = new IOObject("У НПЗ"),
                output2 = new IOObject("ж/д Светлый");

            Pipe pipe1 = new Pipe("ТУ0") { TechnologicalSectionName = "ТУ0", Source = mainInput, Target = rp1, MaxFlows = AlgorithmHelper.CreateListOfElements(Period, 51.1 / 24) },
                pipe2_in = new Pipe("ТУ1 вход") { TechnologicalSectionName = "ТУ1", Source = rp1 },
                pipe2_1 = new Pipe("ТУ1 подкачка 1") { TechnologicalSectionName = "ТУ1", Source = input1, MaxFlows = AlgorithmHelper.CreateListOfElements(Period, bigMaxflowValue) },
                pipe2_2 = new Pipe("ТУ1 подкачка 2") { TechnologicalSectionName = "ТУ1", Source = input2, MaxFlows = AlgorithmHelper.CreateListOfElements(Period, bigMaxflowValue) },
                pipe2_out = new Pipe("ТУ1 выход") { TechnologicalSectionName = "ТУ1", Target = rp2 },
                pipe3 = new Pipe("ТУ2") { TechnologicalSectionName = "ТУ2", Source = rp2, Target = rp3 },
                pipe4_in = new Pipe("ТУ3 вход") { TechnologicalSectionName = "ТУ3", Source = rp3},
                pipe4_1 = new Pipe("ТУ3 откачка ж/д Светлый") { TechnologicalSectionName = "ТУ3", Target = output2, MaxFlows = AlgorithmHelper.CreateListOfElements(Period, bigMaxflowValue) },
                pipe4_out = new Pipe("ТУ3 вход") { TechnologicalSectionName = "ТУ3", Target = mainOutput },
                pipe5 = new Pipe("Труба откачка У НПЗ") { TechnologicalSectionName = "Откачка У НПЗ", Source = rp2, Target = output1, MaxFlows = AlgorithmHelper.CreateListOfElements(Period, 21.6 / 24) },
                pipe6 = new Pipe("Труба подкачка ГНПС 2") { TechnologicalSectionName = "Подкачка ГНПС 2", Source = input3, Target = rp2, MaxFlows = AlgorithmHelper.CreateListOfElements(Period, 10.0 / 24) };
            
            //////////////////////
            // Ремонты
            //////////////////////
            DiscreteSchedule pipe2MaxFlow = new DiscreteSchedule(Period, bigMaxflowValue);
            pipe2MaxFlow.FillInterval(61.5 / 24, 10 * 24 + 12, 10 * 24 + 13);
            pipe2MaxFlow.FillInterval(46.1 / 24, 11 * 24 + 12, 11 * 24 + 16);
            pipe2MaxFlow.FillInterval(61.5 / 24, 12 * 24 + 12, 12 * 24 + 14);
            pipe2MaxFlow.FillInterval(53.5 / 24, 17 * 24 + 12, 17 * 24 + 16);
            pipe2MaxFlow.FillInterval(46.1 / 24, 18 * 24 + 12, 18 * 24 + 14);
            pipe2MaxFlow.FillInterval(46.1 / 24, 26 * 24 + 12, 26 * 24 + 14);
            pipe2_out.MaxFlows = pipe2MaxFlow.ToList();
            pipe2_in.MaxFlows = pipe2MaxFlow.ToList();
            DiscreteSchedule pipe3MaxFlow = new DiscreteSchedule(Period, bigMaxflowValue);
            pipe3MaxFlow.FillInterval(56.1 / 24, 10 * 24 + 12, 10 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 11 * 24 + 12, 11 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 12 * 24 + 12, 12 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 13 * 24 + 12, 13 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 24 * 24 + 12, 24 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 25 * 24 + 12, 25 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 26 * 24 + 12, 26 * 24 + 14);
            pipe3MaxFlow.FillInterval(56.1 / 24, 27 * 24 + 12, 27 * 24 + 14);
            pipe3.MaxFlows = pipe3MaxFlow.ToList();
            DiscreteSchedule pipe4MaxFlow = new DiscreteSchedule(Period, bigMaxflowValue);
            pipe4MaxFlow.FillInterval(47.3 / 24, 17 * 24 + 12, 17 * 24 + 14);
            pipe4MaxFlow.FillInterval(47.3 / 24, 18 * 24 + 12, 18 * 24 + 14);
            pipe4MaxFlow.FillInterval(32.2 / 24, 19 * 24 + 12, 19 * 24 + 16);
            pipe4MaxFlow.FillInterval(47.3 / 24, 20 * 24 + 12, 20 * 24 + 14);
            pipe4MaxFlow.FillInterval(47.3 / 24, 24 * 24 + 12, 24 * 24 + 14);
            pipe4MaxFlow.FillInterval(47.3 / 24, 25 * 24 + 12, 25 * 24 + 14);
            pipe4MaxFlow.FillInterval(33 / 24, 26 * 24 + 12, 26 * 24 + 16);
            pipe4MaxFlow.FillInterval(47.3 / 24, 27 * 24 + 12, 27 * 24 + 14);
            pipe4_out.MaxFlows = pipe4MaxFlow.ToList();
            pipe4_in.MaxFlows = pipe4MaxFlow.ToList();

            points = new List<NamedObject>() { rp1, rp2, rp3, mainInput, mainOutput, input1, input2, input3, output1, output2 };
            pipes = new List<Pipe>() { pipe1, pipe2_in, pipe2_1, pipe2_2, pipe2_out, pipe3, pipe4_in, pipe4_1, pipe4_out, pipe5, pipe6 };

            /////////////////////
            // Режимы
            /////////////////////
            // Первый ТУ
            regimes.Add(new Regime("ТУ1", "1.2-0-0-0-0-0", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 14.4 },
                { pipe2_1, 6.552 },
                { pipe2_2, 0.696 },
                { pipe2_out, 21.648 }
            }));
            regimes.Add(new Regime("ТУ1", "10.22-22-2-22-22-2", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 51.84 },
                { pipe2_1, 6.66 },
                { pipe2_2, 0.516 },
                { pipe2_out, 59.016 }
            }));
            regimes.Add(new Regime("ТУ1", "12.22-22-22-22-22-12", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 55.2 },
                { pipe2_1, 5.292 },
                { pipe2_2, 1.032 },
                { pipe2_out, 61.524 }
            }));
            regimes.Add(new Regime("ТУ1", "13.522-22-22-22-22-22", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 54.576 },
                { pipe2_1, 7.512 },
                { pipe2_2, 0.432 },
                { pipe2_out, 62.52 }
            }));
            regimes.Add(new Regime("ТУ1", "14.222-22-22-222-22-22", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 58.416 },
                { pipe2_1, 5.112 },
                { pipe2_2, 1.152 },
                { pipe2_out, 64.68 }
            }));
            regimes.Add(new Regime("ТУ1", "16.222-222-22-222-222-22", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 61.68 },
                { pipe2_1, 5.016 },
                { pipe2_2, 1.344 },
                { pipe2_out, 68.04 }
            }));
            regimes.Add(new Regime("ТУ1", "17.222-222-222-222-222-22", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 62.5 },
                { pipe2_1, 5.4 },
                { pipe2_2, 1.4 },
                { pipe2_out, 69.32640234375 }
            }));
            regimes.Add(new Regime("ТУ1", "2.2-0-0-2-0-0", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 23.4 },
                { pipe2_1, 6.3 },
                { pipe2_2, 0.7 },
                { pipe2_out, 30.444 }
            }));
            regimes.Add(new Regime("ТУ1", "3.2-2-0-2-0-0", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 24.8 },
                { pipe2_1, 7.2 },
                { pipe2_2, 0.7 },
                { pipe2_out, 32.748 }
            }));
            regimes.Add(new Regime("ТУ1", "4.2-2-0-2-2-0", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 33.4 },
                { pipe2_1, 6.2 },
                { pipe2_2, 0.8 },
                { pipe2_out, 40.4172 }
            }));
            regimes.Add(new Regime("ТУ1", "6.2-2-2-2-2-1", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 38.8 },
                { pipe2_1, 6.6 },
                { pipe2_2, 0.7 },
                { pipe2_out, 46.146 }
            }));
            regimes.Add(new Regime("ТУ1", "7.22-2-2-2-2-2", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 43.4 },
                { pipe2_1, 6.6 },
                { pipe2_2, 0.7 },
                { pipe2_out, 50.728 }
            }));
            regimes.Add(new Regime("ТУ1", "8.22-2-2-22-2-2", 10, new Dictionary<Pipe, double>()
            {
                { pipe2_in, 46.62 },
                { pipe2_1, 6.168 },
                { pipe2_2, 0.684 },
                { pipe2_out, 53.46857 }
            }));
            // Второй ТУ
            regimes.Add(new Regime("ТУ2", "1.1-0-0-0", 10, new Dictionary<Pipe, double>() {
                { pipe3, 14.76 }
            }));
            regimes.Add(new Regime("ТУ2", "10.33-11-511-511", 10, new Dictionary<Pipe, double>() {
                { pipe3, 62.904 }
            }));
            regimes.Add(new Regime("ТУ2", "2.3-1-0-0", 10, new Dictionary<Pipe, double>() {
                { pipe3, 24 }
            }));
            regimes.Add(new Regime("ТУ2", "4.1-5-5-5", 10, new Dictionary<Pipe, double>() {
                { pipe3, 37.4 }
            }));
            regimes.Add(new Regime("ТУ2", "4.3-1-1-1", 10, new Dictionary<Pipe, double>() {
                { pipe3, 40.6 }
            }));
            regimes.Add(new Regime("ТУ2", "5.13-1-1-1", 10, new Dictionary<Pipe, double>() {
                { pipe3, 44.3 }
            }));
            regimes.Add(new Regime("ТУ2", "6.13-1-11-1", 10, new Dictionary<Pipe, double>() {
                { pipe3, 48.428841796875 }
            }));
            regimes.Add(new Regime("ТУ2", "8.31-11-51-11", 10, new Dictionary<Pipe, double>() {
                { pipe3, 56.064 }
            }));
            regimes.Add(new Regime("ТУ2", "8.31-51-11-51", 10, new Dictionary<Pipe, double>() {
                { pipe3, 54.1 }
            }));
            // Третий ТУ
            regimes.Add(new Regime("ТУ3", "10.511-11-511-11 с отбором", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 51.4 },
                { pipe4_1, 2.8 },
                { pipe4_out, 48.6 }
            }));
            regimes.Add(new Regime("ТУ3", "12.511-511-511-151", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 54.48 },
                { pipe4_1, 0 },
                { pipe4_out, 54.48 }
            }));
            regimes.Add(new Regime("ТУ3", "2.1-1-0-0", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 18 },
                { pipe4_1, 0 },
                { pipe4_out, 18 }
            }));
            regimes.Add(new Regime("ТУ3", "3.1-1-1-0", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 26.124 },
                { pipe4_1, 0 },
                { pipe4_out, 26.124 }
            }));
            regimes.Add(new Regime("ТУ3", "4.1-1-1-7", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 33 },
                { pipe4_1, 0 },
                { pipe4_out, 33 }
            }));
            regimes.Add(new Regime("ТУ3", "6.11-1-11-1", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 41.6 },
                { pipe4_1, 0 },
                { pipe4_out, 41.6 }
            }));
            regimes.Add(new Regime("ТУ3", "6.11-1-11-1 с отбором", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 42.2 },
                { pipe4_1, 2.2 },
                { pipe4_out, 40 }
            }));
            regimes.Add(new Regime("ТУ3", "8.11-11-11-15", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 47.3 },
                { pipe4_1, 0 },
                { pipe4_out, 47.3 }
            }));
            regimes.Add(new Regime("ТУ3", "8.11-11-11-51 с отбором", 10, new Dictionary<Pipe, double>() {
                { pipe4_in, 48.5 },
                { pipe4_1, 2.9 },
                { pipe4_out, 45.6 }
            }));
            // Уровняем сумму вход-выход
            regimes.ForEach(regime =>
            {
                switch (regime.TechnologicalSectionName)
                {
                    case "ТУ1":
                        regime.Q[pipe2_out] = regime[pipe2_in].Value + regime[pipe2_1].Value + regime[pipe2_2].Value;
                        break;
                    case "ТУ3":
                        regime.Q[pipe4_out] = regime[pipe4_in].Value - regime[pipe4_1].Value;
                        break;
                }
            });
            // Сделаем значения "в час"
            regimes.ForEach(regime =>
            {
                foreach (var key in regime.Q.Keys.ToList())
                {
                    regime.Q[key] /= 24.0;
                }
            });


            //string path = Directory.GetParent(Directory.GetCurrentDirectory()).Parent.FullName + "\\Режимы.xlsx";
            //if (!File.Exists(path))
            //{
            //    return;
            //}
            //Workbook book = new Workbook();
            //book.LoadFromFile(path);
            //Worksheet tu = book.Worksheets["ТУ Уса-Ухта"];
            //for(int i = 0; i < tu.Rows.Count(); i++)
            //{
            //    var row = tu.Rows[i];
            //    if (/*row.Cells[1].Value == "Вход в РП Ухта-1"*/ row.Cells[0].Value == "Январь")
            //    {
            //        //i++;
            //        row = tu.Rows[i];
            //        double val1 = Double.Parse(row.Cells[4].Value) / 24,
            //            val2 = Double.Parse(row.Cells[3].Value) / 24,
            //            val3 = Double.Parse(row.Cells[2].Value) / 24,
            //            val4 = val1 + val2 + val3;
            //        regimes.Add(new Regime("ТУ1", tu.Rows[i - 2].Cells[0].Value, 10, new Dictionary<Pipe, double>() {
            //            { pipe2_in, val1 },
            //            { pipe2_1, val2 },
            //            { pipe2_2, val3 },
            //            { pipe2_out, val4 }
            //        }));
            //    }
            //    /*if (row.Cells[0].Value == "Январь")
            //    {
            //        regimes.Add(new Regime("ТУ1", tu.Rows[i - 2].Cells[0].Value, 10, new Dictionary<Pipe, double>() {
            //            { pipe2_in, Double.Parse(row.Cells[1].Value) / 24 },
            //            { pipe2_1, Double.Parse(row.Cells[3].Value) / 24 },
            //            { pipe2_2, Double.Parse(row.Cells[2].Value) / 24 },
            //            { pipe2_out, Double.Parse(row.Cells[4].Value) / 24 }
            //        }));
            //    }*/
            //}
            //tu = book.Worksheets["ТУ Ухта-Приводино"];
            //for (int i = 0; i < tu.Rows.Count(); i++)
            //{
            //    var row = tu.Rows[i];
            //    if (row.Cells[0].Value == "Январь")
            //    {
            //        regimes.Add(new Regime("ТУ2", tu.Rows[i - 2].Cells[0].Value, 10, new Dictionary<Pipe, double>() {
            //            { pipe3, Double.Parse(row.Cells[1].Value) / 24 }
            //        }));
            //    }
            //}
            //tu = book.Worksheets["ТУ Приводино-Ярославль"];
            //for (int i = 0; i < tu.Rows.Count(); i++)
            //{
            //    var row = tu.Rows[i];
            //    if (row.Cells[0].Value == "Январь")
            //    {
            //        double val1 = Double.Parse(row.Cells[3].Value) / 24,
            //            val2 = Double.Parse(row.Cells[2].Value) / 24,
            //            val3 = val1 - val2;
            //        regimes.Add(new Regime("ТУ3", tu.Rows[i - 2].Cells[0].Value, 10, new Dictionary<Pipe, double>() {
            //            { pipe4_in, val1 },
            //            { pipe4_1, val2 },
            //            { pipe4_out, val3 }
            //        }));
            //    }
            //}
        }


        public NamedObject GetPoint(string name)
        {
            return points.FirstOrDefault(x => x.Name == name);
        }

        public Pipe GetPipe(string name)
        {
            return pipes.FirstOrDefault(x => x.Name == name);
        }

        public List<string> GetTechnologicalSectionNames()
        {
            return pipes.Select(p => p.TechnologicalSectionName).Distinct().ToList();
        }

        public List<Pipe> GetTechnologicalSectionPipes(string name)
        {
            return pipes.Where(p => p.TechnologicalSectionName == name).ToList();
        }

        public List<double[]> GetTechnologicalSectionMaxFlows(string name)
        {
            var tsPipes = GetTechnologicalSectionPipes(name);

            if (pipes.Count() == 0)
                return new List<double[]>();

            return AlgorithmHelper.CombineIntoList(pipes.Select(x => x.MaxFlows).ToArray());
        }

        public List<double[]> GetTechnologicalSectionRegimes(string name)
        {
            List<Pipe> tsPipes = GetTechnologicalSectionPipes(name);

            if (tsPipes.Count() == 0)
                return new List<double[]>();

            List<Regime> tsRegimes = regimes.Where(r => r.TechnologicalSectionName == name).ToList();

            List<double[]> result = new List<double[]>();
            foreach(var regime in tsRegimes)
            {
                result.Add(tsPipes.Select(x => regime[x].Value).ToArray());
            }

            return AlgorithmHelper.RemoveDuplicates(result);
        }

        public List<double[]> GetTechnologicalSectionRegimes(string name, bool[] mask)
        {
            if (mask.All(x => !x)) throw new ArgumentException();

            var tsRegimes = GetTechnologicalSectionRegimes(name);

            if (tsRegimes.Count() == 0)
                return new List<double[]>();

            return AlgorithmHelper.RemoveDuplicates(AlgorithmHelper.MaskListOfArrays(tsRegimes, mask));
        }

        public List<List<double[]>> GetStates()
        {
            var ts = GetTechnologicalSectionNames();
            var regimes = ts.Select(x => GetTechnologicalSectionRegimes(x)).Where(x => x.Count() != 0).ToList();
            var states = AlgorithmHelper.CartesianProduct(regimes).Select(x => x.ToList()).ToList();
            return states;
        }

        public Dictionary<Reservoir, List<double[]>> GetTankersMasks()
        {
            var tankers = points.Where(x => x.GetType() == typeof(Reservoir)).Select(x => x as Reservoir).ToList();
            var tankerDependences = tankers.Select(tanker =>
            {
                return pipes.Where(pipe => pipe.Target == tanker || pipe.Source == tanker)
                    .ToDictionary(pipe => pipe, pipe => pipe.Target == tanker ? 1 : -1);
            }).ToList();
            var masks = new List<List<double[]>>();
            var tsNames = GetTechnologicalSectionNames();
            var tsPipes = tsNames.Select(x => GetTechnologicalSectionPipes(x)).ToList();
            foreach (var tankerDependence in tankerDependences)
            {
                var mask = tsPipes.Select(x => x.Select(y => 0.0).ToArray()).ToList();

                foreach (var pipe in tankerDependence)
                {
                    string tsName = pipe.Key.TechnologicalSectionName;
                    int tsNumber = tsNames.IndexOf(tsName);
                    int pipeNumber = tsPipes[tsNumber].IndexOf(pipe.Key);
                    mask[tsNumber][pipeNumber] = pipe.Value;
                }
                masks.Add(mask);
            }

            return masks.Select((x,i) => new { IDX = i, EL = x}).ToDictionary(el => tankers[el.IDX], el => el.EL);
        }
        
        public void Algorithm()
        {
            var targets =
                 new List<double[]>() {
                    new double[] { 1584.27 },
                    new double[] { 1586.6, 221.6, 9.8, 1818.0},
                    new double[] { 1677.0 },
                    new double[] { 1677.0, 15.0, 1662.0},
                    new double[] { 269.0 },
                    new double[] { 162.847 } };
            var tankersStartVolume = new double[] { 12.8, 30.6, 17.7 };

            var targetsVector = AlgorithmHelper.ListToVector(targets);

            // Проверка на объем РП
            var tankers = points.Where(x => x.GetType() == typeof(Reservoir)).Select(x => x as Reservoir).ToList();
            var tankersDependence = GetTankersMasks();
            for (int i = 0; i < tankers.Count(); i++)
            {
                var dependence = AlgorithmHelper.ListToVector(tankersDependence[tankers[i]]);
                var difVol = dependence.Zip(targetsVector, (x, y) => x * y).Sum();
                var endVol = tankersStartVolume[i] + difVol;
                if (endVol > tankers[i].Volume || endVol < 0)
                {
                    MessageBox.Show($"Ошибка исходных данных. К концу месяца объем нефти в {tankers[i].Name} будет {endVol} тыс.т.");
                    return;
                }
            }

            var tsNames = GetTechnologicalSectionNames();
            var tsMathModels = tsNames.Select(tsName => new SectionMathModel(GetTechnologicalSectionRegimes(tsName), ));


            var tsPipes = tsNames.Select(x => GetTechnologicalSectionPipes(x)).ToList();
            var tsMaxFlows = tsNames.Select(x => GetTechnologicalSectionMaxFlows(x)).ToList();

            var tsSurplusMask = new List<bool[]>
            {
                null,
                new bool[] { true, true, true, false },
                new bool[] { true },
                new bool[] { true, true, false},
                null,
                null
            };

            var tsMathModels = tsRegimes.Zip(tsSurplusMask, (x,y) => 
                y == null ? null : new SectionMathModel(AlgorithmHelper.RemoveDuplicates(AlgorithmHelper.MaskListOfArrays(x, y)))).ToList();

            var minVol = targets[1].Select(x => x * 0.99).ToArray();
            minVol[0] = targets[1][0] * 0.99;
            var maxVol = targets[1].Select(x => x).ToArray();
            maxVol[0] = targets[1][0] * 1.01;
            var schedule = tsMathModels[1].GetOptimalDecomposition(targets[1], minVol, maxVol, Period, new double[] { 1, 2, 1.5 });

            var vol = AlgorithmHelper.GetSumOnInterval(schedule.Select(x => x.Item1.Select(y => y * x.Item2).ToArray()).ToList(), 0, schedule.Count());
        }

        #endregion

        #region Свойства

        public int Period
        {
            get;
            private set;
        }

        public int Days
        {
            get;
            private set;
        }

        #endregion
    }
}
