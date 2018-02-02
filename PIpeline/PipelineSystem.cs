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

        public List<SubSection> pipes = new List<SubSection>();
        public List<NamedObject> points = new List<NamedObject>();
        public List<Regime> regimes = new List<Regime>();
        public Dictionary<string, List<Tuple<double, double[]>>> sectionRepairs = new Dictionary<string, List<Tuple<double, double[]>>>();

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

            SubSection pipe1 = new SubSection("ТУ0", mainInput, rp1, "ТУ0", 51.1 / 24),
                pipe2_in = new SubSection("ТУ1 вход", rp1, null, "ТУ1", 69.5 / 24),
                pipe2_1 = new SubSection("ТУ1 подкачка 1", input1, null, "ТУ1", 250.0 / (31 * 24)),
                pipe2_2 = new SubSection("ТУ1 подкачка 2", input2, null, "ТУ1", 50.0 / (31 * 24)),
                pipe2_out = new SubSection("ТУ1 выход", null, rp2, "ТУ1", 69.5 / 24),
                pipe3 = new SubSection("ТУ2", rp2, rp3, "ТУ2", 64.3 / 24),
                pipe4_in = new SubSection("ТУ3 вход", rp3, null, "ТУ3", 56.4 / 24),
                pipe4_1 = new SubSection("ТУ3 откачка ж/д Светлый", null, output2, "ТУ3", 10.8 / 24),
                pipe4_out = new SubSection("ТУ3 вход", null, mainInput, "ТУ3", 54.5 / 24),
                pipe5 = new SubSection("Труба откачка У НПЗ", rp2, output1, "Труба откачка У НПЗ", 21.6 / 24),
                pipe6 = new SubSection("Труба подкачка ГНПС 2", input3, rp2, "Труба подкачка ГНПС 2", 9.72 / 24);

            points = new List<NamedObject>() { rp1, rp2, rp3, mainInput, mainOutput, input1, input2, input3, output1, output2 };
            pipes = new List<SubSection>() { pipe1, pipe2_in, pipe2_1, pipe2_2, pipe2_out, pipe3, pipe4_in, pipe4_1, pipe4_out, pipe5, pipe6 };

            //////////////////////
            // Ремонты
            //////////////////////
            DiscreteSchedule tu1MainRepairs = new DiscreteSchedule(Period, bigMaxflowValue);
            tu1MainRepairs.FillInterval(61.5 / 24, 10 * 24 + 12, 10 * 24 + 13);
            tu1MainRepairs.FillInterval(46.1 / 24, 11 * 24 + 12, 11 * 24 + 16);
            //tu1MainRepairs.FillInterval(0, 11 * 24 + 12, 11 * 24 + 24);
            tu1MainRepairs.FillInterval(61.5 / 24, 12 * 24 + 12, 12 * 24 + 14);
            tu1MainRepairs.FillInterval(53.5 / 24, 17 * 24 + 12, 17 * 24 + 16);
            tu1MainRepairs.FillInterval(46.1 / 24, 18 * 24 + 12, 18 * 24 + 14);
            tu1MainRepairs.FillInterval(46.1 / 24, 26 * 24 + 12, 26 * 24 + 14);
            var tu1Rep = tu1MainRepairs.ToList().Select(x => new Tuple<double, double[]>(x, new double[] { bigMaxflowValue, bigMaxflowValue })).ToList();
            DiscreteSchedule tu2MainRepairs = new DiscreteSchedule(Period, bigMaxflowValue);
            tu2MainRepairs.FillInterval(56.1 / 24, 10 * 24 + 12, 10 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 11 * 24 + 12, 11 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 12 * 24 + 12, 12 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 13 * 24 + 12, 13 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 24 * 24 + 12, 24 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 25 * 24 + 12, 25 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 26 * 24 + 12, 26 * 24 + 14);
            tu2MainRepairs.FillInterval(56.1 / 24, 27 * 24 + 12, 27 * 24 + 14);
            var tu2Rep = tu2MainRepairs.ToList().Select(x => new Tuple<double, double[]>(x, new double[] { })).ToList();
            DiscreteSchedule tu3MainRepairs = new DiscreteSchedule(Period, bigMaxflowValue);
            tu3MainRepairs.FillInterval(47.3 / 24, 17 * 24 + 12, 17 * 24 + 14);
            tu3MainRepairs.FillInterval(47.3 / 24, 18 * 24 + 12, 18 * 24 + 14);
            tu3MainRepairs.FillInterval(32.2 / 24, 19 * 24 + 12, 19 * 24 + 16);
            tu3MainRepairs.FillInterval(47.3 / 24, 20 * 24 + 12, 20 * 24 + 14);
            tu3MainRepairs.FillInterval(47.3 / 24, 24 * 24 + 12, 24 * 24 + 14);
            tu3MainRepairs.FillInterval(47.3 / 24, 25 * 24 + 12, 25 * 24 + 14);
            tu3MainRepairs.FillInterval(33 / 24, 26 * 24 + 12, 26 * 24 + 16);
            tu3MainRepairs.FillInterval(47.3 / 24, 27 * 24 + 12, 27 * 24 + 14);
            var tu3Rep = tu3MainRepairs.ToList().Select(x => new Tuple<double, double[]>(x, new double[] { bigMaxflowValue })).ToList();
            sectionRepairs = new Dictionary<string, List<Tuple<double, double[]>>>()
            {
                {"ТУ1",  tu1Rep},
                {"ТУ2", tu2Rep },
                {"ТУ3", tu3Rep },
                {"ТУ0", AlgorithmHelper.CreateListOfElements(Period, bigMaxflowValue).Select(x => new Tuple<double, double[]>(bigMaxflowValue,  new double[] { })).ToList() },
                {"Труба подкачка ГНПС 2", AlgorithmHelper.CreateListOfElements(Period, bigMaxflowValue).Select(x => new Tuple<double, double[]>(bigMaxflowValue,  new double[] { })).ToList() },
                {"Труба откачка У НПЗ", AlgorithmHelper.CreateListOfElements(Period, bigMaxflowValue).Select(x => new Tuple<double, double[]>(bigMaxflowValue,  new double[] { })).ToList() }
            };

            /////////////////////
            // Режимы
            /////////////////////
            // Первый ТУ
            regimes.Add(new Regime("ТУ1", "0", 10, new Tuple<double, double[][]>(0, new double[][] { new double[] { 0, 0 }, new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ1", "1.2-0-0-0-0-0", 10, new Tuple<double, double[][]>(14.4, new double[][] { new double[] { 0, 6.552 }, new double[] {0, 0.696 }})));
            regimes.Add(new Regime("ТУ1", "10.22-22-2-22-22-2", 10, new Tuple<double, double[][]>(51.84, new double[][] { new double[] { 0, 6.66 }, new double[] { 0, 0.516 } })));
            regimes.Add(new Regime("ТУ1", "12.22-22-22-22-22-12", 10, new Tuple<double, double[][]>(55.2, new double[][] { new double[] { 0, 5.292 }, new double[] { 0, 1.032 } })));
            regimes.Add(new Regime("ТУ1", "13.522-22-22-22-22-22", 10, new Tuple<double, double[][]>(54.576, new double[][] { new double[] { 0, 7.512 }, new double[] { 0, 0.432 } })));
            regimes.Add(new Regime("ТУ1", "14.222-22-22-222-22-22", 10, new Tuple<double, double[][]>(58.416, new double[][] { new double[] { 0, 5.112 }, new double[] { 0, 1.152 } })));
            regimes.Add(new Regime("ТУ1", "16.222-222-22-222-222-22", 10, new Tuple<double, double[][]>(61.68, new double[][] { new double[] { 0, 5.016 }, new double[] { 0, 1.344 } })));
            regimes.Add(new Regime("ТУ1", "17.222-222-222-222-222-22", 10, new Tuple<double, double[][]>(62.5, new double[][] { new double[] { 0, 5.4 }, new double[] { 0, 1.4 } })));
            regimes.Add(new Regime("ТУ1", "2.2-0-0-2-0-0", 10, new Tuple<double, double[][]>(23.4, new double[][] { new double[] { 0, 6.3 }, new double[] { 0, 0.7 } })));
            regimes.Add(new Regime("ТУ1", "3.2-2-0-2-0-0", 10, new Tuple<double, double[][]>(24.8, new double[][] { new double[] { 0, 7.2 }, new double[] { 0, 0.7 } })));
            regimes.Add(new Regime("ТУ1", "4.2-2-0-2-2-0", 10, new Tuple<double, double[][]>(33.4, new double[][] { new double[] { 0, 6.2 }, new double[] { 0, 0.8 } })));
            regimes.Add(new Regime("ТУ1", "6.2-2-2-2-2-1", 10, new Tuple<double, double[][]>(38.8, new double[][] { new double[] { 0, 6.6 }, new double[] { 0, 0.7 } })));
            regimes.Add(new Regime("ТУ1", "7.22-2-2-2-2-2", 10, new Tuple<double, double[][]>(43.4, new double[][] { new double[] { 0, 6.6 }, new double[] { 0, 0.7 } })));
            regimes.Add(new Regime("ТУ1", "8.22-2-2-22-2-2", 10, new Tuple<double, double[][]>(46.62, new double[][] { new double[] { 0, 6.168 }, new double[] { 0, 0.684 } })));

            // Второй ТУ
            regimes.Add(new Regime("ТУ2", "0", 10, new Tuple<double, double[][]>(0, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "1.1-0-0-0", 10, new Tuple<double, double[][]>(14.76, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "10.33-11-511-511", 10, new Tuple<double, double[][]>(62.904, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "2.3-1-0-0", 10, new Tuple<double, double[][]>(24, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "4.1-5-5-5", 10, new Tuple<double, double[][]>(37.4, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "4.3-1-1-1", 10, new Tuple<double, double[][]>(40.6, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "5.13-1-1-1", 10, new Tuple<double, double[][]>(44.3, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "6.13-1-11-1", 10, new Tuple<double, double[][]>(48.428841796875, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "8.31-11-51-11", 10, new Tuple<double, double[][]>(56.064, new double[][] { })));
            regimes.Add(new Regime("ТУ2", "8.31-51-11-51", 10, new Tuple<double, double[][]>(54.1, new double[][] { })));

            // Третий ТУ
            regimes.Add(new Regime("ТУ3", "0", 10, new Tuple<double, double[][]>(0, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "10.511-11-511-11 с отбором", 10, new Tuple<double, double[][]>(51.4, new double[][] { new double[] { 0, 2.8 } })));
            regimes.Add(new Regime("ТУ3", "12.511-511-511-151", 10, new Tuple<double, double[][]>(54.48, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "2.1-1-0-0", 10, new Tuple<double, double[][]>(18, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "3.1-1-1-0", 10, new Tuple<double, double[][]>(26.124, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "4.1-1-1-7", 10, new Tuple<double, double[][]>(33, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "6.11-1-11-1", 10, new Tuple<double, double[][]>(41.6, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "6.11-1-11-1 с отбором", 10, new Tuple<double, double[][]>(42.2, new double[][] { new double[] { 0, 2.2 } })));
            regimes.Add(new Regime("ТУ3", "8.11-11-11-15", 10, new Tuple<double, double[][]>(47.3, new double[][] { new double[] { 0, 0 } })));
            regimes.Add(new Regime("ТУ3", "8.11-11-11-51 с отбором", 10, new Tuple<double, double[][]>(48.5, new double[][] { new double[] { 0, 2.9 } })));
            
            // Сделаем значения "в час"
            regimes.ForEach(regime =>
            {
                regime.G = new Tuple<double, double[][]>(regime.G.Item1 / 24, regime.G.Item2.Select(x => new double[] {
                    //x[0] / 24
                    (x[1] / 24) * 0.8
                    , (x[1] / 24) * 1.2 }).ToArray());
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

        public SubSection GetPipe(string name)
        {
            return pipes.FirstOrDefault(x => x.Name == name);
        }

        public List<string> GetTechnologicalSectionNames()
        {
            return pipes.Select(p => p.TechnologicalSectionName).Distinct().ToList();
        }

        public List<SubSection> GetTechnologicalSectionPipes(string name)
        {
            return pipes.Where(p => p.TechnologicalSectionName == name).ToList();
        }

        public List<Tuple<double, double[], double>> GetTechnologicalSectionMaxFlows(string name)
        {
            var tsPipes = GetTechnologicalSectionPipes(name);
            List<Tuple<double, double[], double>> maxFlows = new List<Tuple<double, double[], double>>();

            if (pipes.Count() == 0)
                return maxFlows;

            var tsRepairs = sectionRepairs[name];
            for (int i = 0; i < Period; i++)
            {
                var pipesMaxFlows = tsPipes.Select(x => x.MaxFlows).ToList();
                var repair = tsRepairs[i];
                if (tsPipes.Count() == 1)
                {
                    var inputMax = Math.Min(pipesMaxFlows.First(), repair.Item1);
                    maxFlows.Add(new Tuple<double, double[], double>(inputMax, new double[] { }, inputMax));
                }
                else
                {
                    var inputMax = Math.Min(pipesMaxFlows.First(), repair.Item1);
                    var outputMax = Math.Min(pipesMaxFlows.Last(), repair.Item1);
                    var pumpsMax = pipesMaxFlows.GetRange(1, tsPipes.Count() - 2).Zip(repair.Item2, (x, y) => y < x ? y : x).ToArray();
                    maxFlows.Add(new Tuple<double, double[], double>(inputMax, pumpsMax, outputMax));
                }
            }

            return maxFlows;
        }

        public List<Regime> GetTechnologicalSectionRegimes(string name)
        {
            return regimes.Where(regime => regime.TechnologicalSectionName == name).ToList();
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
        
        public SectionMathModel CreateSectionMathModel(string name)
        {
            var regimes = GetTechnologicalSectionRegimes(name);

            if (regimes.Count() == 0)
                return null;

            var maxFlows = GetTechnologicalSectionMaxFlows(name);
            var tsPipes = GetTechnologicalSectionPipes(name);
            List<SubSection> tsPumps = new List<SubSection>();
            if (tsPipes.Count() > 1)
                tsPumps = tsPipes.GetRange(1, tsPipes.Count() - 2);

            return new SectionMathModel(regimes.Select(regime => regime.G).ToList(), tsPumps.Select(pump => pump.Target == null ? +1.0 : -1.0).ToArray(), maxFlows);
        }

        public List<double[]> Algorithm()
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
                    return null;
                }
            }

            // Проверка, можно ли перекачать по секцийм такие объемы


            var tu1MathModel = CreateSectionMathModel("ТУ1");
            return tu1MathModel.AddOutputElement(tu1MathModel.GetRegimesDecomposition(targets[1][0], targets[1].ToList().GetRange(1, 2).ToArray(), 0, Period - 1, new double[] { 1, 1, 1 }).Select(x => SectionMathModel.Convert(x)).ToList());

           // return tu1MathModel.AddOutputElement(tu1MathModel.GetContinuousSchedule(targets[1][0], targets[1].ToList().GetRange(1, 2).ToArray(), 0, Period - 1, tu1MathModel.GetNormQuadraticWeights(new double[] { 1, 1, 1 })));
            /*var tsNames = GetTechnologicalSectionNames();
            var tsMathModels = GetTechnologicalSectionNames().Select(name =>
            {
                var regimes = GetTechnologicalSectionRegimes(name);

                if (regimes.Count() == 0)
                    return null;

                var maxFlows = GetTechnologicalSectionMaxFlows(name);
                var tsPipes = GetTechnologicalSectionPipes(name);
                List<SubSection> tsPumps = new List<SubSection>();
                if (tsPipes.Count() > 1)
                    tsPumps = tsPipes.GetRange(1, tsPipes.Count() - 2);
                
                return new SectionMathModel(regimes.Select(regime => regime.G).ToList(), tsPumps.Select(pump => pump.Target == null ? +1.0 : -1.0).ToArray(), maxFlows);
            }).ToList();*/


            /*var tsPipes = tsNames.Select(x => GetTechnologicalSectionPipes(x)).ToList();
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

            var vol = AlgorithmHelper.GetSumOnInterval(schedule.Select(x => x.Item1.Select(y => y * x.Item2).ToArray()).ToList(), 0, schedule.Count());*/
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

        #region Для графика, потом удалить

        public List<double[]> GetTechnologicalSectionInputOutputMax(string name)
        {
            var maxFlows = GetTechnologicalSectionMaxFlows(name);
            var arrMaxFlows = new List<double[]>();
            for (int i = 0; i < Period; i++)
            {
                var maxFlow = new List<double>();
                maxFlow.Add(maxFlows[i].Item1);
                maxFlow.AddRange(maxFlows[i].Item2);
                maxFlow.Add(maxFlows[i].Item3);
                arrMaxFlows.Add(maxFlow.ToArray());
            }

            // Ограничиваем режимами
            var regimes = GetTechnologicalSectionRegimes(name);
            for(int i = 0; i < Period; i++)
            {
                var maxRegime = regimes.Max(x => x.G.Item1 > arrMaxFlows[i][0] ? double.NegativeInfinity : x.G.Item1);
                arrMaxFlows[i][0] = maxRegime;
            }

            return arrMaxFlows;
        }

        #endregion
    }
}
