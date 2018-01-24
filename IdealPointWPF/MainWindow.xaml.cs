using Algorithms;
using OxyPlot;
using Pipeline;
using Pipeline.PipelineObjects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using OxyPlot.Series;
using System.Threading;
using Show;

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
                inMaxFlows = new LineSeries()
                {
                    Color = OxyColors.Red,
                    StrokeThickness = 2,
                    LineStyle = LineStyle.LongDash
                },
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
                };

            List<LineSeries>
                inRegimes = new List <LineSeries>(),
                tu1InRegimes = new List<LineSeries>(),
                tu1Pump1Regimes = new List<LineSeries>(),
                tu1Pump2Regimes = new List<LineSeries>(),
                tu1OutRegimes = new List<LineSeries>(),
                tu2Regimes = new List<LineSeries>(),
                tu3InRegimes = new List<LineSeries>(),
                tu3PumpRegimes = new List<LineSeries>(),
                tu3OutRegimes = new List<LineSeries>();


            List<double[]>
                inMaxFlowsList = p.GetTechnologicalSectionMaxFlows("ТУ0"),
                inRegimesList = p.GetTechnologicalSectionRegimes("ТУ0"),
                tu1MaxFlowsList = p.GetTechnologicalSectionMaxFlows("ТУ1"),
                tu1RegimesList = p.GetTechnologicalSectionRegimes("ТУ1"),
                tu2MaxFlowsList = p.GetTechnologicalSectionMaxFlows("ТУ2"),
                tu2RegimesList = p.GetTechnologicalSectionRegimes("ТУ2"),
                tu3MaxFlowsList = p.GetTechnologicalSectionMaxFlows("ТУ3"),
                tu3RegimesList = p.GetTechnologicalSectionRegimes("ТУ3");

            InSchedule = new PlotModel();
            double[] inRegimesMax = AlgorithmHelper.GetMaxInListByComponents(inRegimesList);
            inMaxFlows.ItemsSource = inMaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(inRegimesMax[0], current[0])));
                total.Add(new DataPoint(i + 1, Math.Min(inRegimesMax[0], current[0])));
                return total;
            });
            InSchedule.Series.Add(inMaxFlows);
            inRegimesList.ForEach(x => inRegimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[0]),
                    new DataPoint(p.Period, x[0])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            inRegimes.ForEach(x => InSchedule.Series.Add(x));

            TU1InputSchedule = new PlotModel();
            TU1Pump1Schedule = new PlotModel();
            TU1Pump2Schedule = new PlotModel();
            TU1OutputSchedule = new PlotModel();
            double[] tu1RegimesMax = AlgorithmHelper.GetMaxInListByComponents(tu1RegimesList);
            tu1InMaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu1RegimesMax[0], current[0])));
                total.Add(new DataPoint(i + 1, Math.Min(tu1RegimesMax[0], current[0])));
                return total;
            });
            tu1Pump1MaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu1RegimesMax[1], current[1])));
                total.Add(new DataPoint(i + 1, Math.Min(tu1RegimesMax[1], current[1])));
                return total;
            });
            tu1Pump2MaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu1RegimesMax[2], current[2])));
                total.Add(new DataPoint(i + 1, Math.Min(tu1RegimesMax[2], current[2])));
                return total;
            });
            tu1OutMaxFlows.ItemsSource = tu1MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu1RegimesMax[3], current[3])));
                total.Add(new DataPoint(i + 1, Math.Min(tu1RegimesMax[3], current[3])));
                return total;
            });
            TU1InputSchedule.Series.Add(tu1InMaxFlows);
            TU1Pump1Schedule.Series.Add(tu1Pump1MaxFlows);
            TU1Pump2Schedule.Series.Add(tu1Pump2MaxFlows);
            TU1OutputSchedule.Series.Add(tu1OutMaxFlows);
            tu1RegimesList.ForEach(x => tu1InRegimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[0]),
                    new DataPoint(p.Period, x[0])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu1InRegimes.ForEach(x => TU1InputSchedule.Series.Add(x));
            tu1RegimesList.ForEach(x => tu1Pump1Regimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[1]),
                    new DataPoint(p.Period, x[1])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu1Pump1Regimes.ForEach(x => TU1Pump1Schedule.Series.Add(x));
            tu1RegimesList.ForEach(x => tu1Pump2Regimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[2]),
                    new DataPoint(p.Period, x[2])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu1Pump2Regimes.ForEach(x => TU1Pump2Schedule.Series.Add(x));
            tu1RegimesList.ForEach(x => tu1OutRegimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[3]),
                    new DataPoint(p.Period, x[3])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu1OutRegimes.ForEach(x => TU1OutputSchedule.Series.Add(x));
            
            TU2Schedule = new PlotModel();
            double[] tu2RegimesMax = AlgorithmHelper.GetMaxInListByComponents(tu2RegimesList);
            tu2MaxFlows.ItemsSource = tu2MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu2RegimesMax[0], current[0])));
                total.Add(new DataPoint(i + 1, Math.Min(tu2RegimesMax[0], current[0])));
                return total;
            });
            TU2Schedule.Series.Add(tu2MaxFlows);
            tu2RegimesList.ForEach(x => tu2Regimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[0]),
                    new DataPoint(p.Period, x[0])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu2Regimes.ForEach(x => TU2Schedule.Series.Add(x));


            TU3InputSchedule = new PlotModel();
            TU3PumpSchedule = new PlotModel();
            TU3OutputSchedule = new PlotModel();
            double[] tu3RegimesMax = AlgorithmHelper.GetMaxInListByComponents(tu3RegimesList);
            tu3InMaxFlows.ItemsSource = tu3MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu3RegimesMax[0], current[0])));
                total.Add(new DataPoint(i + 1, Math.Min(tu3RegimesMax[0], current[0])));
                return total;
            });
            tu3PumpMaxFlows.ItemsSource = tu3MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu3RegimesMax[1], current[1])));
                total.Add(new DataPoint(i + 1, Math.Min(tu3RegimesMax[1], current[1])));
                return total;
            });
            tu3OutMaxFlows.ItemsSource = tu3MaxFlowsList.Aggregate(new List<DataPoint>(), (total, current) =>
            {
                int i = total.Count() / 2;
                total.Add(new DataPoint(i + 0.001, Math.Min(tu3RegimesMax[2], current[2])));
                total.Add(new DataPoint(i + 1, Math.Min(tu3RegimesMax[2], current[2])));
                return total;
            });
            TU3InputSchedule.Series.Add(tu3InMaxFlows);
            TU3PumpSchedule.Series.Add(tu3PumpMaxFlows);
            TU3OutputSchedule.Series.Add(tu3OutMaxFlows);
            tu3RegimesList.ForEach(x => tu3InRegimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[0]),
                    new DataPoint(p.Period, x[0])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu3InRegimes.ForEach(x => TU3InputSchedule.Series.Add(x));
            tu3RegimesList.ForEach(x => tu3PumpRegimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[1]),
                    new DataPoint(p.Period, x[1])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu3PumpRegimes.ForEach(x => TU3PumpSchedule.Series.Add(x));
            tu3RegimesList.ForEach(x => tu3OutRegimes.Add(new LineSeries()
            {
                ItemsSource = new List<DataPoint>()
                {
                    new DataPoint(0, x[2]),
                    new DataPoint(p.Period, x[2])
                },
                StrokeThickness = 0.5,
                Color = OxyColors.Green,
                LineStyle = LineStyle.Dash
            }));
            tu3OutRegimes.ForEach(x => TU3OutputSchedule.Series.Add(x));

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
        }

        public void ShowAlgo()
        {
            // Формирование системы
            SetWinTitle("Формируем систему");
            PipelineSystem p = new PipelineSystem();

            // Объем, принимаемый от соседнего ОСТ
            double p1Volume = 1584.27;
            // Объем, откачки1 на ТУ1
            double p3Volume = 221.6;
            // Объем, откачки2 на ТУ1
            double p4Volume = 9.8;
            // Объем перекачки в конце ТУ1
            double p5Volume = 1818;
            // Объем перекачки в начале ТУ1
            double p2Volume = p5Volume - p3Volume - p4Volume;
            // Объем перекачки ТУ2
            double p6Volume = 1677;
            // Объем откачки ТУ3
            double p8Volume = 15;
            // Объем перекачки в конце ТУ3
            double p9Volume = 1662;
            // Объем перекачки в начале ТУ3
            double p7Volume = p9Volume - p8Volume;

            // объемы РП
            double 
                rp1Volume = 47.4,
                rp2Volume = 114.3,
                rp3Volume = 47.4;

            // Подкачки к РП
            double 
                rp1PumpVolume = 0,
                rp2PumpVolume = 162.847 - 269,
                rp3PumpVolume = 0;
            double[]
                rp1PumpSchedule = AlgorithmHelper.CreateListOfElements(p.Period, rp1PumpVolume / p.Period).ToArray(),
                rp2PumpSchedule = AlgorithmHelper.CreateListOfElements(p.Period, rp2PumpVolume / p.Period).ToArray(),
                rp3PumpSchedule = AlgorithmHelper.CreateListOfElements(p.Period, rp3PumpVolume / p.Period).ToArray();

            // Приурачиваем путевые откачки и подкачки к РП
            rp1PumpSchedule = rp1PumpSchedule.Select(x => x + p3Volume / p.Period + p4Volume / p.Period).ToArray();
            rp3PumpSchedule = rp3PumpSchedule.Select(x => x - p8Volume / p.Period).ToArray();

            // Начальные объемы в РП
            double
                rp1StartVolume = 3,//12.8,
                rp2StartVolume = 1,//30.6,
                rp3StartVolume = 1;//17.7;

            // Проверяем данные
            //РП1
            double rp1End = rp1StartVolume + p1Volume + p3Volume + p4Volume - p5Volume + rp1PumpVolume;
            if (rp1End < 0 || rp1End > rp1Volume)
            {
                SetWinTitle("Не соблюден баланс в РП1, решение задачи невозможно");
                return;
            }
            //РП2
            double rp2End = rp2StartVolume + p5Volume - p6Volume + rp2PumpVolume;
            if (rp2End < 0 || rp2End > rp2Volume)
            {
                SetWinTitle("Не соблюден баланс в РП2, решение задачи невозможно");
                return;
            }
            //РП3
            double rp3End = rp3StartVolume + p6Volume - p8Volume - p9Volume + rp3PumpVolume;
            if (rp3End < 0 || rp3End > rp3Volume)
            {
                SetWinTitle("Не соблюден баланс в РП3, решение задачи невозможно");
                return;
            }

            // Алгоритм для первой трубы
            var p1MaxFlow = AlgorithmHelper.CreateListOfElements(p.Period, 100.0).ToArray();
            var p1Regimes = new double[] { p1Volume / p.Period };
            var p1Algo = new TechnologicalSectionAlgorithm(p1MaxFlow, p1Regimes, 31);

            // Алгоритм для пятой трубы
            var p5MaxFlow = p.GetTechnologicalSectionMaxFlows("ТУ1").Select(x => x[3]).ToArray();
            var p5Regimes = p.GetTechnologicalSectionRegimes("ТУ1").Select(x => x[3]).Distinct().ToArray();
            var p5Algo = new TechnologicalSectionAlgorithm(p5MaxFlow, p5Regimes, 31);

            // Алгоритм для шестой трубы
            var p6MaxFlow = p.GetTechnologicalSectionMaxFlows("ТУ2").Select(x => x[0]).ToArray();
            var p6Regimes = p.GetTechnologicalSectionRegimes("ТУ2").Select(x => x[0]).ToArray();
            var p6Algo = new TechnologicalSectionAlgorithm(p6MaxFlow, p6Regimes, 31);

            // Алгоритм для девятой трубы
            var p9MaxFlow = p.GetTechnologicalSectionMaxFlows("ТУ3").Select(x => x[2]).ToArray();
            var p9Regimes = p.GetTechnologicalSectionRegimes("ТУ3").Select(x => x[2]).Distinct().ToArray();
            var p9Algo = new TechnologicalSectionAlgorithm(p9MaxFlow, p9Regimes, 31);

            // Отсутствие фиксированных расходов
            double[] noFix = AlgorithmHelper.CreateListOfElements(p.Period, -1.0).ToArray();
            double[] p1Schedule = AlgorithmHelper.CreateListOfElements(p.Period, p1Volume / p.Period).ToArray();

            // Балансируем РП1
            SetWinTitle("Баланс РП1.");
            var rp1BalanceAlgorithm = new ReservoirBalancerAlgorithm(p1Algo, p5Algo, rp1Volume);
            var rp1Balance = rp1BalanceAlgorithm.Balance(p1Volume, p5Volume, p1Schedule, noFix, rp1StartVolume, rp1PumpSchedule,
                (string info) => SetWinTitle("Баланс РП1. " + info));
            if (rp1Balance == null)
            {
                SetWinTitle("Не удалось сбалансировать РП1. Решение не найдено.");
                return;
            }
            double[]
                p2Schedule = rp1Balance[1],
                p3Schedule = AlgorithmHelper.CreateListOfElements(p.Period, p3Volume / p.Period).ToArray(),
                p4Schedule = AlgorithmHelper.CreateListOfElements(p.Period, p4Volume / p.Period).ToArray(),
                p5Schedule = rp1Balance[1];

            // Балансируем РП2
            SetWinTitle("Баланс РП2.");
            var rp2BalanceAlgorithm = new ReservoirBalancerAlgorithm(p5Algo, p6Algo, rp2Volume);
            var rp2Balance = rp2BalanceAlgorithm.Balance(p5Volume, p6Volume, p5Schedule, noFix, rp2StartVolume, rp2PumpSchedule,
                (string info) => SetWinTitle("Баланс РП2. " + info));
            if (rp2Balance == null)
            {
                SetWinTitle("Не удалось сбалансировать РП2. Решение не найдено.");
                return;
            }
            double[] p6Schedule = rp2Balance[1];

            // Балансируем РП3
            SetWinTitle("Баланс РП3.");
            var rp3BalanceAlgorithm = new ReservoirBalancerAlgorithm(p6Algo, p9Algo, rp3Volume);
            var rp3Balance = rp3BalanceAlgorithm.Balance(p6Volume, p9Volume, p6Schedule, noFix, rp3StartVolume, rp3PumpSchedule,
                (string info) => SetWinTitle("Баланс РП3. " + info));
            if (rp3Balance == null)
            {
                SetWinTitle("Не удалось сбалансировать РП3. Решение не найдено.");
                return;
            }
            double[] 
                p7Schedule = rp3Balance[1],
                p8Schedule = AlgorithmHelper.CreateListOfElements(p.Period, p8Volume / p.Period).ToArray(),
                p9Schedule = rp3Balance[1];

            //проверка баланса
            double rp1Cur = rp1StartVolume,
                rp2Cur = rp2StartVolume,
                rp3Cur = rp3StartVolume;
            for (int i = 0; i < p.Period; i++)
            {
                rp1Cur += p1Schedule[i] - p2Schedule[i] + rp1PumpSchedule[i];
                rp2Cur += p5Schedule[i] - p6Schedule[i] + rp2PumpSchedule[i];
                rp3Cur += p6Schedule[i] - p7Schedule[i] + rp3PumpSchedule[i];

                if (rp1Cur > rp1Volume || rp1Cur < 0)
                    throw new Exception();
                if (rp2Cur > rp2Volume || rp2Cur < 0)
                    throw new Exception();
                if (rp3Cur > rp3Volume || rp3Cur < 0)
                    throw new Exception();
            }

            SetWinTitle("Решение найдено.");
            GetBasePlots(p);
            SetPlots(
                rp1Balance[0], //p1
                p2Schedule, //p2
                p3Schedule, //p3
                p4Schedule, //p4
                p5Schedule, //p5
                p6Schedule, //p6
                p7Schedule, //p7
                p8Schedule, //p8
                p9Schedule, //p9
                rp1StartVolume,
                rp2StartVolume,
                rp3StartVolume,
                rp1PumpSchedule,
                rp2PumpSchedule,
                rp3PumpSchedule);
            UpdateDatacontext();
        }

        public void MainAlgo()
        {
            PipelineSystem p = new PipelineSystem();

            p.Algorithm();
            
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
