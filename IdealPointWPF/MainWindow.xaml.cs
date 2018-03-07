using OxyPlot;
using Pipeline;
using Pipeline.PipelineObjects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using OxyPlot.Series;
using System.Threading;
using Accord.Math;

namespace IdealPointWPF
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            var a = DataContext as MainViewModel;
            a.Win = this;
            model = a;
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            new Thread(() => 
            {
                model.MainAlgo();
            }).Start();
        }

        public MainViewModel model
        {
            get;
            set;
        }
    }

    public class MainViewModel
    {
        private void GetBasePlots(PipelineSystem p)
        {
            LineSeries
                tu1InMaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu1Pump1MaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu1Pump2MaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu1OutMaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu2MaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu3InMaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu3PumpMaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                tu3OutMaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
                zeros = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                };

            Func<LineSeries, LineSeries> CopySerie = (source) =>
            {
                return new LineSeries()
                {
                    ItemsSource = source.ItemsSource as List<DataPoint>,
                    StrokeThickness = source.StrokeThickness,
                    Color = source.Color,
                    LineStyle = source.LineStyle
                };
            };

            List<double[]>
                tu1MaxFlowsList = p.GetTechnologicalSectionInputOutputMax("ТУ1"),
                tu2MaxFlowsList = p.GetTechnologicalSectionInputOutputMax("ТУ2"),
                tu3MaxFlowsList = p.GetTechnologicalSectionInputOutputMax("ТУ3");
            
            zeros.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, 0));
                total.Add(new DataPoint(i + 1, 0));
                return total;
            });

            InSchedule = new PlotModel();
            TU1InputSchedule = new PlotModel();
            TU1Pump1Schedule = new PlotModel();
            TU1Pump2Schedule = new PlotModel();
            TU1OutputSchedule = new PlotModel();
            tu1InMaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[0]));
                total.Add(new DataPoint(i + 1, current[0]));
                return total;
            });
            tu1Pump1MaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[1]));
                total.Add(new DataPoint(i + 1, current[1]));
                return total;
            });
            tu1Pump2MaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[2]));
                total.Add(new DataPoint(i + 1, current[2]));
                return total;
            });
            tu1OutMaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[3]));
                total.Add(new DataPoint(i + 1, current[3]));
                return total;
            });
            TU1InputSchedule.Series.Add(tu1InMaxFlows);
            TU1Pump1Schedule.Series.Add(tu1Pump1MaxFlows);
            TU1Pump2Schedule.Series.Add(tu1Pump2MaxFlows);
            TU1OutputSchedule.Series.Add(tu1OutMaxFlows);
            TU1InputSchedule.Series.Add(CopySerie(zeros));
            TU1Pump1Schedule.Series.Add(CopySerie(zeros));
            TU1Pump2Schedule.Series.Add(CopySerie(zeros));
            TU1OutputSchedule.Series.Add(CopySerie(zeros));

            TU2Schedule = new PlotModel();
            tu2MaxFlows.ItemsSource = tu2MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[0]));
                total.Add(new DataPoint(i + 1, current[0]));
                return total;
            });
            TU2Schedule.Series.Add(tu2MaxFlows);
            TU2Schedule.Series.Add(CopySerie(zeros));


            TU3InputSchedule = new PlotModel();
            TU3PumpSchedule = new PlotModel();
            TU3OutputSchedule = new PlotModel();
            tu3InMaxFlows.ItemsSource = tu3MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[0]));
                total.Add(new DataPoint(i + 1, current[0]));
                return total;
            });
            tu3PumpMaxFlows.ItemsSource = tu3MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[1]));
                total.Add(new DataPoint(i + 1, current[1]));
                return total;
            });
            tu3OutMaxFlows.ItemsSource = tu3MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, current[2]));
                total.Add(new DataPoint(i + 1, current[2]));
                return total;
            });
            TU3InputSchedule.Series.Add(tu3InMaxFlows);
            TU3PumpSchedule.Series.Add(tu3PumpMaxFlows);
            TU3OutputSchedule.Series.Add(tu3OutMaxFlows);
            TU3InputSchedule.Series.Add(CopySerie(zeros));
            TU3PumpSchedule.Series.Add(CopySerie(zeros));
            TU3OutputSchedule.Series.Add(CopySerie(zeros));

            Reservoir tank1 = (p.GetPoint("РП1") as Reservoir);
            Tank1 = new PlotModel();
            Tank1.Series.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, tank1.Volume),
                    new DataPoint(p.Period, tank1.Volume)
                },
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.LongDash
            });
            Tank1.Series.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, 0),
                    new DataPoint(p.Period, 0)
                },
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.LongDash
            });
            Reservoir tank2 = (p.GetPoint("РП2") as Reservoir);
            Tank2 = new PlotModel();
            Tank2.Series.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, tank2.Volume),
                    new DataPoint(p.Period, tank2.Volume)
                },
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.LongDash
            });
            Tank2.Series.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, 0),
                    new DataPoint(p.Period, 0)
                },
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.LongDash
            });
            Reservoir tank3 = (p.GetPoint("РП3") as Reservoir);
            Tank3 = new PlotModel();
            Tank3.Series.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, tank3.Volume),
                    new DataPoint(p.Period, tank3.Volume)
                },
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.LongDash
            });
            Tank3.Series.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, 0),
                    new DataPoint(p.Period, 0)
                },
                Color = OxyColors.Red,
                StrokeThickness = 2,
                LineStyle = LineStyle.LongDash
            });

            UpdateDatacontext();
        }

        private List<DataPoint> ConvertToDatapoints(double[] arr)
        {
            List<DataPoint> result = new List<DataPoint>();
            for (int i = 0; i < arr.Count(); i++)
            {
                result.Add(new DataPoint(i + 0.001, arr[i]));
                result.Add(new DataPoint(i + 1, arr[i]));
            }
            return result;
        }

        private void SetPlots(double[] p1Schedule, double[] p2Schedule, double[] p3Schedule, double[] p4Schedule, 
            double[] p5Schedule, double[] p6Schedule, double[] p7Schedule, double[] p8Schedule, double[] p9Schedule,
            double tank1StartVolume, double tank2StartVolume, double tank3StartVolume,
            double[] rp1PumpSchedule, double[] rp2PumpSchedule, double[] rp3PumpSchedule)
        {
            int period = p1Schedule.Count();
            List<DataPoint> 
                tank1Schedule = new List<DataPoint>() { new DataPoint(0, tank1StartVolume) },
                tank2Schedule = new List<DataPoint>() { new DataPoint(0, tank2StartVolume) },
                tank3Schedule = new List<DataPoint>() { new DataPoint(0, tank3StartVolume) };
            for (int i = 0; i < period; i++)
            {
                tank1Schedule.Add(new DataPoint(i + 1, tank1Schedule.Last().Y + p1Schedule[i] - p2Schedule[i] + rp1PumpSchedule[i]));
                tank2Schedule.Add(new DataPoint(i + 1, tank2Schedule.Last().Y + p5Schedule[i] - p6Schedule[i] + rp2PumpSchedule[i]));
                tank3Schedule.Add(new DataPoint(i + 1, tank3Schedule.Last().Y + p6Schedule[i] - p7Schedule[i] + rp3PumpSchedule[i]));
            }

            InSchedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p1Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            InSchedule.Title = "Q = " + p1Schedule.Sum().ToString();
            TU1InputSchedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p2Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU1InputSchedule.Title = "Q = " + p2Schedule.Sum().ToString();
            TU1Pump1Schedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p3Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU1Pump1Schedule.Title = "Q = " + p3Schedule.Sum().ToString();
            TU1Pump2Schedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p4Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU1Pump2Schedule.Title = "Q = " + p4Schedule.Sum().ToString();
            TU1OutputSchedule.Series.Add(new LineSeries()
             {
                 ItemsSource = ConvertToDatapoints(p5Schedule),
                 Color = OxyColors.Blue,
                 StrokeThickness = 2,
                 LineStyle = LineStyle.Solid
            });
            TU1OutputSchedule.Title = "Q = " + p5Schedule.Sum().ToString();
            TU2Schedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p6Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU2Schedule.Title = "Q = " + p6Schedule.Sum().ToString();
            TU3InputSchedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p7Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU3InputSchedule.Title = "Q = " + p7Schedule.Sum().ToString();
            TU3PumpSchedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p8Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU3PumpSchedule.Title = "Q = " + p8Schedule.Sum().ToString();
            TU3OutputSchedule.Series.Add(new LineSeries()
            {
                ItemsSource = ConvertToDatapoints(p9Schedule),
                Color = OxyColors.Blue,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid
            });
            TU3OutputSchedule.Title = "Q = " + p9Schedule.Sum().ToString();
            Tank1.Series.Add(new LineSeries()
             {
                 ItemsSource = tank1Schedule,
                 Color = OxyColors.Blue,
                 StrokeThickness = 2,
                 LineStyle = LineStyle.Solid
             });
             Tank2.Series.Add(new LineSeries()
             {
                 ItemsSource = tank2Schedule,
                 Color = OxyColors.Blue,
                 StrokeThickness = 2,
                 LineStyle = LineStyle.Solid
             });
             Tank3.Series.Add(new LineSeries()
             {
                 ItemsSource = tank3Schedule,
                 Color = OxyColors.Blue,
                 StrokeThickness = 2,
                 LineStyle = LineStyle.Solid
             });

            UpdateDatacontext();
        }

        public void MainAlgo()
        {
            PipelineSystem p = new PipelineSystem();
            GetBasePlots(p);

            var tu3Interval = Accord.Math.Vector.EnumerableRange(12 * 24, 19 * 24).ToArray();
            var interval = Accord.Math.Vector.EnumerableRange(0, p.Period).ToArray();
            PipelineSystem.PumpsParameters targets = new PipelineSystem.PumpsParameters(new List<string> { "ТУ0", "ТУ1", "ТУ2", "ТУ3", "Труба откачка У НПЗ", "Труба подкачка ГНПС 2" });
            targets.SetUniformity("ТУ0", true);
            targets.SetUniformity("Труба откачка У НПЗ", true);
            targets.SetUniformity("Труба подкачка ГНПС 2", true);
            targets.AddBatch("ТУ0", new double[] { 1584.27 }, interval);
            targets.AddBatch("ТУ1", new double[] { 1586.6, 221.6, 15, 1586.6 + 221.6 + 14 }, interval);
            targets.AddBatch("ТУ2", new double[] { 1650.0 }, interval);
            targets.AddBatch("ТУ3", new double[] { 350, 15.0, 350 - 15.0 }, tu3Interval);
            targets.AddBatch("ТУ3", new double[] { 1640 - 350, 0, 1640 - 350 }, interval.Except(tu3Interval).ToArray());
            targets.AddBatch("Труба откачка У НПЗ", new double[] { 269.0 }, interval);
            targets.AddBatch("Труба подкачка ГНПС 2", new double[] { 162.847 }, interval);

            var reservoir2Pumps = Algorithms.AlgorithmHelper.CreateListOfElements(p.Period, (-targets.batches["Труба откачка У НПЗ"][0].Item1[0] + targets.batches["Труба подкачка ГНПС 2"][0].Item1[0]) / p.Period).ToArray();

            var tankersStartVolume = new double[] { 12.8, 30.6, 17.7 };

            DateTime time = DateTime.Now;
            Dictionary<string, List<double[]>> tuSchedules = p.Algorithm(targets, tankersStartVolume, (str) => SetWinTitle(str));
            //MessageBox.Show((DateTime.Now - time).TotalSeconds.ToString());

            if (tuSchedules == null)
                return;

            double[] zeros = new double[p.Period];
            double[] inputSchedule = Algorithms.AlgorithmHelper.CreateListOfElements(p.Period, targets.batches["ТУ0"][0].Item1[0] / p.Period).ToArray();

            SetPlots(
                inputSchedule,
                tuSchedules["ТУ1"].Select(x => x[0]).ToArray(),
                tuSchedules["ТУ1"].Select(x => x[1]).ToArray(),
                tuSchedules["ТУ1"].Select(x => x[2]).ToArray(),
                tuSchedules["ТУ1"].Select(x => x[3]).ToArray(),
                tuSchedules["ТУ2"].Select(x => x[0]).ToArray(),
                tuSchedules["ТУ3"].Select(x => x[0]).ToArray(),
                tuSchedules["ТУ3"].Select(x => x[1]).ToArray(),
                tuSchedules["ТУ3"].Select(x => x[2]).ToArray(),
                tankersStartVolume[0],
                tankersStartVolume[1],
                tankersStartVolume[2],
                zeros,
                reservoir2Pumps,
                zeros);
            
        }

        public MainViewModel()
        {
        }
       
        public PlotModel InSchedule
        {
            get;
            private set;
        }
        
        public PlotModel TU1InputSchedule
        {
            get;
            private set;
        }
        
        private PlotModel TU1InputScheduleCopy
        {
            get
            {
                if (TU1InputSchedule == null)
                {
                    return null;
                }
                PlotModel m = new PlotModel();
                foreach (LineSeries el in TU1InputSchedule.Series)
                {
                    m.Series.Add(new LineSeries()
                    {
                        ItemsSource = el.ItemsSource as List<DataPoint>,
                        StrokeThickness = el.StrokeThickness,
                        Color = el.Color,
                        LineStyle = el.LineStyle
                    });
                }
                m.Title = TU1InputSchedule.Title;
                return m;
            }
        }
        
        public PlotModel TU1InputSchedule1
        {
            get
            {
                return TU1InputScheduleCopy;
            }
        }
        
        public PlotModel TU1InputSchedule2
        {
            get
            {
                return TU1InputScheduleCopy;
            }
        }

        public PlotModel TU1Pump1Schedule
        {
            get;
            private set;
        }

        public PlotModel TU1Pump2Schedule
        {
            get;
            private set;
        }

        public PlotModel TU1OutputSchedule
        {
            get;
            private set;
        }

        private PlotModel TU1OutputScheduleCopy
        {
            get
            {
                if (TU1OutputSchedule == null)
                {
                    return null;
                }
                PlotModel m = new PlotModel();
                foreach (LineSeries el in TU1OutputSchedule.Series)
                {
                    m.Series.Add(new LineSeries()
                    {
                        ItemsSource = el.ItemsSource as List<DataPoint>,
                        StrokeThickness = el.StrokeThickness,
                        Color = el.Color,
                        LineStyle = el.LineStyle
                    });
                }
                m.Title = TU1OutputSchedule.Title;
                return m;
            }
        }
        
        public PlotModel TU1OutputSchedule1
        {
            get
            {
                return TU1OutputScheduleCopy;
            }
        }
        
        public PlotModel TU1OutputSchedule2
        {
            get
            {
                return TU1OutputScheduleCopy;
            }
        }

        public PlotModel TU2Schedule
        {
            get;
            private set;
        }

        public PlotModel TU2ScheduleCopy
        {
            get
            {
                if (TU2Schedule == null)
                {
                    return null;
                }
                PlotModel m = new PlotModel();
                foreach (LineSeries el in TU2Schedule.Series)
                {
                    m.Series.Add(new LineSeries()
                    {
                        ItemsSource = el.ItemsSource as List<DataPoint>,
                        StrokeThickness = el.StrokeThickness,
                        Color = el.Color,
                        LineStyle = el.LineStyle
                    });
                }
                m.Title = TU2Schedule.Title;
                return m;
            }
        }
        
        public PlotModel TU2Schedule1
        {
            get
            {
                return TU2ScheduleCopy;
            }
        }
        
        public PlotModel TU2Schedule2
        {
            get
            {
                return TU2ScheduleCopy;
            }
        }

        public PlotModel TU3InputSchedule
        {
            get;
            private set;
        }
        
        public PlotModel TU3InputScheduleCopy
        {
            get
            {
                if (TU3InputSchedule == null)
                {
                    return null;
                }
                PlotModel m = new PlotModel();
                foreach (LineSeries el in TU3InputSchedule.Series)
                {
                    m.Series.Add(new LineSeries()
                    {
                        ItemsSource = el.ItemsSource as List<DataPoint>,
                        StrokeThickness = el.StrokeThickness,
                        Color = el.Color,
                        LineStyle = el.LineStyle
                    });
                }
                m.Title = TU3InputSchedule.Title;
                return m;
            }
        }
        
        public PlotModel TU3InputSchedule1
        {
            get
            {
                return TU3InputScheduleCopy;
            }
        }
        
        public PlotModel TU3InputSchedule2
        {
            get
            {
                return TU3InputScheduleCopy;
            }
        }

        public PlotModel TU3PumpSchedule
        {
            get;
            private set;
        }
        
        public PlotModel TU3OutputSchedule
        {
            get;
            private set;
        }

        public PlotModel Tank1
        {
            get;
            private set;
        }

        public PlotModel Tank2
        {
            get;
            private set;
        }

        public PlotModel Tank3
        {
            get;
            private set;
        }

        public MainWindow Win
        {
            get;
            set;
        }

        public void SetWinTitle(string str)
        {
            Win.Dispatcher.Invoke(() =>
            {
                Win.Title = str;
            });
        }

        public void UpdateDatacontext()
        {
            Win.Dispatcher.Invoke(() =>
            {
                Win.DataContext = null;
                Win.DataContext = this;
            });
        }
    }
}
