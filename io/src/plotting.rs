// SPDX-FileCopyrightText: 2026 <Oliver A. Bramley; Technical University of Denmark>
//
// SPDX-License-Identifier: AGPL-3.0-only

//! Plot consolidated Monte Carlo CSV results with Plotters.
//!
//! [`SimulationResults::from_csv`] loads the `average_fill.csv` and
//! `average_event.csv` files produced after a simulation. Individual methods
//! create filling, temperature, and event-frequency plots, while
//! [`plot_default_results`] creates a useful standard set in one call.

use crate::outputs::{AverageEventRow, ContinuousValueRow};
use crate::errors::PlotError;
use common::numeric::{Float, TimeFloat};
use common::constants::time::TimeUnit;
use common::constants::temperature::TemperatureUnit;
use plotters::coord::Shift;
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;
use std::fmt;
use std::fs;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Expand a finite data extent so Plotters always receives a non-empty range.
fn padded_range(values: impl IntoIterator<Item = Float>) -> Range<Float> {
    let mut values = values.into_iter();
    let first = values.next().expect("validated result data is non-empty");
    let (mut minimum, mut maximum) = (first, first);
    for value in values {
        minimum = minimum.min(value);
        maximum = maximum.max(value);
    }
    let span = maximum - minimum;
    let padding = if span > 0.0 {
        span * 0.05
    } else {
        minimum.abs().max(1.0) * 0.05
    };
    (minimum - padding)..(maximum + padding)
}

/// Central filling statistic drawn as a line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillAverage {
    /// Arithmetic mean filling fraction.
    Mean,
    /// Median filling fraction.
    Median,
}
impl FillAverage{
    pub fn get_label(&self) -> &str{
        match self {
            FillAverage::Mean => "Mean filling",
            FillAverage::Median => "Median filling",
        }
    
    }
}
/// Optional uncertainty region for a filling plot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillBand {
    /// Draw only the central line.
    None,
    /// Shade one population standard deviation below and above the mean.
    StandardDeviation,
    /// Shade the interquartile range from the 25th to 75th percentile.
    InterquartileRange,
}
impl FillBand {
    pub fn get_label(&self) -> &str{
        match self {
            FillBand::StandardDeviation => "Mean ± standard deviation",
            FillBand::InterquartileRange => "25th–75th percentile",
            FillBand::None => unreachable!(),
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FillStatistic {
    /// Arithmetic mean filling fraction.
    average: FillAverage,
    /// Median filling fraction.
    band: FillBand,
}

impl FillStatistic{

    pub fn new(mean:bool, fill:&str) -> Self{
        let average = if mean{
            FillAverage::Mean
        }else {
            FillAverage::Median
        };

        let band = if fill == "sd" {
            FillBand::StandardDeviation
        } else if fill == "iqr" {
            FillBand::InterquartileRange
        } else {
            FillBand::None
        }; 
        Self { average, band }
    }
    pub fn check_data(&self, fill_rows: &Vec<ContinuousValueRow>, path: PathBuf) -> Result<(),PlotError> {
        if self.band == FillBand::InterquartileRange
            && fill_rows.iter().any(|row| {
                !row.fill_quantile_0_25.is_finite() || !row.fill_quantile_0_75.is_finite()
            })
        {
            return Err(PlotError::InvalidData {
                path: path,
                message: "interquartile plots require fill_quantile_0_25 and fill_quantile_0_75 columns; regenerate this CSV with the current version".into(),
            });
        } else {
            Ok(())
        }
    }

    pub fn get_fill_average(&self,) -> 
        Result< impl Fn(&ContinuousValueRow) -> f64,  PlotError> 
    {
        Ok(move |row: &ContinuousValueRow| match self.average  {
            FillAverage::Mean => row.fill,
            FillAverage::Median => row.fill_median,
        })
       
        
    }
    pub fn get_fill_band(&self,) -> Result< impl Fn(&ContinuousValueRow) 
    -> (f64,f64), PlotError> {
        Ok(move |row: &ContinuousValueRow| match self.band {
            FillBand::None => {
                let bcentral =  match self.average {
                                FillAverage::Mean => row.fill,
                                FillAverage::Median => row.fill_median,
                            };   
                (bcentral, bcentral) 
            }
            
            FillBand::StandardDeviation => (
                row.fill - row.fill_standard_deviation,
                row.fill + row.fill_standard_deviation,
            ),
            FillBand::InterquartileRange => (row.fill_quantile_0_25, row.fill_quantile_0_75),
        })
        
    }
    pub fn get_y_range(&self,fill_rows: &Vec<ContinuousValueRow>,) 
    -> Result< Range< f64 >, PlotError> {
        let average = self.get_fill_average()?;
        let band = self.get_fill_band()?;

        let mut y_values = vec![0.0, 1.0];
        for row in fill_rows {
            let (lower, upper) = band(row);
            y_values.extend([average(row), lower, upper]);
        }

        Ok(padded_range(y_values))
       
    } 

}

/// Horizontal coordinate used for a filling plot.
#[derive(Debug, Clone, Copy, PartialEq,)]
pub enum Axis {
    /// Plot elapsed/profile time in seconds.
    Time {unit: TimeUnit},
    /// Plot against temperature.
    Temperature {unit: TemperatureUnit},
    /// Plot trap filling ratio
    Fill{statistics: FillStatistic },
    /// Plot events/luminescence glow curve
    Event,
}

impl Axis{

    pub fn new(name: &str, unit: &str, mean: bool, fill: &str) -> Result<Self, PlotError> {
        match name {
            "Time" => {
                let unit = TimeUnit::from_str(unit)
                    .map_err(|e| PlotError::Setup{ source: "TimeUnit".into(),
                                                           message: e})?;
                Ok(Self::Time { unit })
            },  
            "Temperature" => {
                let unit = TemperatureUnit::from_str(unit)
                .map_err(|e| PlotError::Setup{ source: "TemperatureUnit".into(),
                                                           message: e})?;
                Ok(Self::Temperature { unit })
            },
            "Fill" => Ok(Self::Fill { statistics: FillStatistic::new(mean,fill)}),
            "Event" => Ok(Self::Event),
            _ => Err(PlotError::Setup { 
                                source: "YAxis".into(), 
                                message: format!("Invalid type {name}") })

        }
    } 

    pub fn get_label(&self) -> String {
        match self {
            Axis::Time{ unit} => {
                format!("Time ({unit})") 
            }
            Axis::Temperature { unit} => {
                format!("Temperature ({unit})")
            }
            Axis::Event => { format!("Events / bin width / repetition (s⁻¹)")},
            Axis::Fill{..} => { format!("n/N filling")},

        }
    }

}

/// Raster dimensions shared by plot-producing methods.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlotOptions {
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,

    pub font:  &'static str,

    pub caption_font_size: u32,

    pub margin: u32
}

impl PlotOptions {

    pub fn validate_dimensions(&self) -> Result<(), PlotError> {
        if self.width == 0 || self.height == 0 {
            return Err(PlotError::Setup { 
                                source: "PlotOptions".into(), 
                                message:"plot width and height must both be greater than zero".into() });
        }
        Ok(())
    }

    pub fn get_caption_options(&self) -> (&'static str, u32){
        (self.font, self.caption_font_size)
    }

}

impl Default for PlotOptions {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 800,
            font: "sans-serif".into(),
            caption_font_size: 32,
            margin: 20
        }
    }
}

pub struct PlotWindow{

    pub output: PathBuf,

    pub x_axis: Axis,

    pub y_axis: Axis,

    pub options: PlotOptions,

}

impl PlotWindow{

    pub fn new(output: PathBuf, x_axis: &str, x_unit: &str, y_axis: &str, y_unit: &str, mean:bool, fill:&str) -> Result<Self, PlotError> {

        let x_axis = Axis::new(x_axis,x_unit, mean, fill)?;

        let y_axis = Axis::new(y_axis, y_unit, mean, fill)?; 

        let options = PlotOptions::default(); 
        options.validate_dimensions()?;
        Ok( Self {
                output,
                x_axis,
                y_axis,
                options
            })
    } 
    pub fn get_x_filter(&self) -> Result< impl Fn(&ContinuousValueRow) -> f64,  PlotError>  {
        let x = match self.x_axis {
            Axis::Time { .. } => {
                Box::new(|row: &ContinuousValueRow| row.time)
                as Box<dyn Fn(&ContinuousValueRow) -> f64>
            }
            Axis::Temperature { .. } => {
                Box::new(|row: &ContinuousValueRow| row.temperature)
            }
        _ => {
            return Err(PlotError::Draw {
                path: self.output.clone(),
                message: "Attempted to access Fill or Events as an X-axis".into(),
            });
            }
        };
        Ok(x)
    }

    pub fn get_continuous_x_range(&self,fill_rows: &Vec<ContinuousValueRow>,)
        -> Result< Range< f64 >, PlotError> 
    {
        let x=self.get_x_filter()?;

        Ok(padded_range(fill_rows.iter().map(x)))

    }

    pub fn get_y_filter(&self) -> Result< impl Fn(&ContinuousValueRow) -> f64,  PlotError>  {
        
        let y = match self.y_axis {
            Axis::Time { .. } => {
                Box::new(|row: &ContinuousValueRow| row.time)
                as Box<dyn Fn(&ContinuousValueRow) -> f64>
            
            },
            Axis::Temperature { .. } => {
                Box::new(|row: &ContinuousValueRow| row.temperature)
            },
            Axis::Fill { statistics } => {
                match statistics.average {
                    FillAverage::Mean => { 
                        Box::new( |row: &ContinuousValueRow| row.fill)
                    },
                    FillAverage::Median => { 
                        Box::new( |row: &ContinuousValueRow| row.fill_median)
                        as Box<dyn Fn(&ContinuousValueRow) -> f64> 
                    },
                }
            },
            _ => return Err(PlotError::Draw {
                path: self.output.clone(),
                message: "Attempted to access Events as an Y-axis".into(),
            })

        };
        Ok(y)

    }


    pub fn get_continuous_y_range(&self,fill_rows: &Vec<ContinuousValueRow>,)
        -> Result< Range< f64 >, PlotError> 
    {
        let y = self.get_y_filter()?;
        match self.y_axis {
            Axis::Time { .. } => {
                
                return Ok(padded_range(fill_rows.iter().map(y)));
            }
            Axis::Temperature { .. } => {
                let y = |row: &ContinuousValueRow| row.temperature;
                return Ok(padded_range(fill_rows.iter().map(y)));
            },
            Axis::Fill { statistics } => {
                let range = statistics.get_y_range(fill_rows)?;
                return Ok(range);
            }
            _ => return Err(PlotError::Draw {
                path: self.output.clone(),
                message: "Not yet written".into(),
            }),

        }
        
    }

    pub fn plot_setup(&self, caption: &str, x_range: Range<Float>, y_range: Range<Float>) -> Result<
        (
            DrawingArea<BitMapBackend<'_>, Shift>,
            ChartContext<
                '_,
                BitMapBackend<'_>,
                Cartesian2d<RangedCoordf64, RangedCoordf64>,
            > 
        ), PlotError>{
        
        let root = BitMapBackend::new(
            &self.output, 
            ( self.options.width, self.options.height))
            .into_drawing_area();

        root.fill(&WHITE)
            .map_err(
                    |error| PlotError::Draw {
                            path: self.output.clone(),
                            message: format!("{error:?}"),
                        }
            )?;

        let mut chart = ChartBuilder::on(&root)
                    .caption(caption, self.options.get_caption_options())
                    .margin(self.options.margin)
                    .x_label_area_size(55)
                    .y_label_area_size(85)
                    .build_cartesian_2d(x_range, y_range)
                    .map_err(
                              | error | PlotError::Draw {
                                    path: self.output.clone(),
                                    message: format!("{error:?}"),
                                }
                    )?;
        chart
                .configure_mesh()
                .x_desc(self.x_axis.get_label())
                .y_desc(self.y_axis.get_label())
                .draw()
                 .map_err(
                        |error| PlotError::Draw {
                                path: self.output.clone(),
                                message: format!("{error:?}"),
                            }
                        )?;
        Ok((root, chart))
    }


     /// Plot mean or median filling against time or temperature.
    pub fn plot_fill(
        &self,
        fill_rows: &Vec<ContinuousValueRow>,
        caption: &str,
    ) -> Result<(), PlotError> {
        
        if let Axis::Fill { statistics } = &self.y_axis {
            statistics.check_data(fill_rows, self.output.clone())?;
        }
       
        let x_range = self.get_continuous_x_range(fill_rows)?;
        let y_range = self.get_continuous_y_range(fill_rows)?; 

        let (root, mut chart) = self.plot_setup(caption, x_range, y_range)?;
        
        

        let x = self.get_x_filter()?;
        if let Axis::Fill { statistics } = &self.y_axis {
            if statistics.band != FillBand::None {
                let bounds = statistics.get_fill_band()?;
                
                let mut polygon = fill_rows
                    .iter()
                    .map(|row| (x(row), bounds(row).1))
                    .collect::<Vec<_>>();
                
                polygon.extend(fill_rows
                        .iter()
                        .rev()
                        .map(|row| (x(row), bounds(row).0)),
                );
                
                let band_label = statistics.band.get_label();
                
                chart.draw_series(
                                std::iter::once(Polygon::new(
                                            polygon,
                                             BLUE.mix(0.18).filled(),
                                                )
                                            )
                                )
                                .map_err(|error| 
                                    PlotError::Draw {
                                            path: self.output.clone(),
                                            message: format!("{error:?}"),
                                            }
                                        )?
                                .label(band_label)
                                .legend(|(x, y)| {
                                        Rectangle::new(
                                            [(x, y - 5), (x + 20, y + 5)], 
                                            BLUE.mix(0.18).filled()
                                        )
                                    }
                                );
            }
            let central = statistics.get_fill_average()?;
            
            chart.draw_series(LineSeries::new(
                fill_rows.iter().map(|row| (x(row), central(row))),
                BLUE.stroke_width(3),
                            )
                        )
                        .map_err(|error| 
                            PlotError::Draw {
                                path: self.output.clone(),
                                message: format!("{error:?}"),
                            }
                        )?
                        .label(statistics.average.get_label())
                        .legend(|(x, y)| 
                            PathElement::new(
                                [(x, y), (x + 20, y)], 
                                BLUE.stroke_width(3))
                            );
        
        
            chart
                .configure_series_labels()
                .background_style(WHITE.mix(0.85))
                .border_style(BLACK)
                .draw()
                .map_err(|error| PlotError::Draw {
                    path: self.output.clone(),
                    message: format!("{error:?}")}
                )?;

            root.present()
                .map_err(|error| PlotError::Draw {
                    path: self.output.clone(),
                    message: format!("{error:?}")}
                )?;
       
    }
    
    Ok(())
    }


      /// Plot temperature in kelvin against time in seconds.
    pub fn plot_temperature_vs_time(
        &self,
        fill_rows: &Vec<ContinuousValueRow>,
        caption: &str,
    ) -> Result<(), PlotError> {

        let x_range = self.get_continuous_x_range(fill_rows)?;
        let y_range = self.get_continuous_y_range(fill_rows)?; 
       
        let (root, mut chart) = self.plot_setup(caption, x_range, y_range)?;
        let x = self.get_x_filter()?;
        let y = self.get_y_filter()?;
        chart.draw_series(LineSeries::new(
                fill_rows.iter().map(|row| (x(row), y(row))),
                RED.stroke_width(3),
            )).map_err(|error| 
                                    PlotError::Draw {
                                            path: self.output.clone(),
                                            message: format!("{error:?}"),
                                            }
                                        )?;
        root.present()
                .map_err(|error| PlotError::Draw {
                    path: self.output.clone(),
                    message: format!("{error:?}")}
                )?;
        Ok(())
    }

    // /// Plot time in seconds against temperature in kelvin.
    // ///
    // /// This is the axis-reversed form of [`Self::plot_temperature_vs_time`].
    // pub fn plot_time_vs_temperature(
    //     &self,
    //     output: impl AsRef<Path>,
    //     options: PlotOptions,
    // ) -> Result<(), PlotError> {
    //     let output = output.as_ref();
    //     validate_dimensions(output, options)?;
    //     let x_range = padded_range(self.fill_rows.iter().map(|row| row.temperature));
    //     let y_range = padded_range(self.fill_rows.iter().map(|row| row.time));
    //     let root = BitMapBackend::new(output, (options.width, options.height)).into_drawing_area();
    //     draw_result(output, root.fill(&WHITE))?;
    //     let mut chart = draw_result(
    //         output,
    //         ChartBuilder::on(&root)
    //             .caption("Time vs temperature", ("sans-serif", 32))
    //             .margin(20)
    //             .x_label_area_size(55)
    //             .y_label_area_size(70)
    //             .build_cartesian_2d(x_range, y_range),
    //     )?;
    //     draw_result(
    //         output,
    //         chart
    //             .configure_mesh()
    //             .x_desc("Temperature (K)")
    //             .y_desc("Time (s)")
    //             .draw(),
    //     )?;
    //     draw_result(
    //         output,
    //         chart.draw_series(LineSeries::new(
    //             self.fill_rows.iter().map(|row| (row.temperature, row.time)),
    //             RED.stroke_width(3),
    //         )),
    //     )?;
    //     draw_result(output, root.present())
    // }

    // /// Plot selected event frequencies against time.
    // ///
    // /// Pass `None` to use the CSV's original bins or `Some(width)` to smooth
    // /// the data with a larger bin width in seconds. CSV timestamps are treated
    // /// as right bin edges; plotted timestamps are the midpoint between adjacent
    // /// edges. Each averaged count is divided by that original or new bin width
    // /// during drawing.
    // pub fn plot_events_vs_time(
    //     &self,
    //     output: impl AsRef<Path>,
    //     series: &[EventSeries],
    //     new_bin_width: Option<TimeFloat>,
    //     options: PlotOptions,
    // ) -> Result<(), PlotError> {
    //     let output = output.as_ref();
    //     validate_dimensions(output, options)?;
    //     if series.is_empty() {
    //         return Err(PlotError::InvalidData {
    //             path: self.event_path.clone(),
    //             message: "at least one event series must be selected".into(),
    //         });
    //     }
    //     let rebinned;
    //     let rows = if let Some(width) = new_bin_width {
    //         rebinned = self.rebin_events(width)?;
    //         rebinned.as_slice()
    //     } else {
    //         self.event_rows.as_slice()
    //     };
    //     let bins = event_plot_bins(rows);
    //     if bins.is_empty() {
    //         return Err(PlotError::InvalidData {
    //             path: self.event_path.clone(),
    //             message: "at least one completed event bin is required for plotting".into(),
    //         });
    //     }
    //     let x_range = padded_range(bins.iter().map(|bin| bin.centre));
    //     let maximum = bins
    //         .iter()
    //         .flat_map(|bin| {
    //             series
    //                 .iter()
    //                 .map(move |column| event_frequency(bin, *column))
    //         })
    //         .fold(0.0_f64, Float::max);
    //     let y_range = 0.0..if maximum > 0.0 { maximum * 1.08 } else { 1.0 };

    //     let root = BitMapBackend::new(output, (options.width, options.height)).into_drawing_area();
    //     draw_result(output, root.fill(&WHITE))?;
    //     let mut chart = draw_result(
    //         output,
    //         ChartBuilder::on(&root)
    //             .caption("Event frequency vs time", ("sans-serif", 32))
    //             .margin(20)
    //             .x_label_area_size(55)
    //             .y_label_area_size(85)
    //             .build_cartesian_2d(x_range, y_range),
    //     )?;
    //     draw_result(
    //         output,
    //         chart
    //             .configure_mesh()
    //             .x_desc("Time (s)")
    //             .y_desc("Events / bin width / repetition (s⁻¹)")
    //             .draw(),
    //     )?;

    //     for (index, column) in series.iter().copied().enumerate() {
    //         let color = Palette99::pick(index).to_rgba();
    //         draw_result(
    //             output,
    //             // Draw adjacent pairs independently. Plotters' bitmap backend
    //             // can generate incorrect polygon joins for a large, dense
    //             // polyline, which appeared as negative spikes even though all
    //             // input frequencies were non-negative.
    //             chart.draw_series(bins.windows(2).map(|pair| {
    //                 PathElement::new(
    //                     [
    //                         (pair[0].centre, event_frequency(&pair[0], column)),
    //                         (pair[1].centre, event_frequency(&pair[1], column)),
    //                     ],
    //                     color.stroke_width(1),
    //                 )
    //             })),
    //         )?
    //         .label(column.label())
    //         .legend(move |(x, y)| PathElement::new([(x, y), (x + 20, y)], color.stroke_width(3)));
    //     }
    //     draw_result(
    //         output,
    //         chart
    //             .configure_series_labels()
    //             .background_style(WHITE.mix(0.85))
    //             .border_style(BLACK)
    //             .draw(),
    //     )?;
    //     draw_result(output, root.present())
    // }



}


/// Event-count column that can be included in an event-frequency plot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventSeries {
    /// Ground-state localised recombination.
    LocalisedRecombinationGround,
    /// Excited-state localised recombination.
    LocalisedRecombinationExcited,
    /// Ground-state delocalised recombination.
    DelocalisedRecombinationGround,
    /// Excited-state delocalised recombination.
    DelocalisedRecombinationExcited,
    /// Ground-state localised retrapping.
    LocalisedRetrappingGround,
    /// Excited-state localised retrapping.
    LocalisedRetrappingExcited,
    /// Ground-state delocalised retrapping.
    DelocalisedRetrappingGround,
    /// Excited-state delocalised retrapping.
    DelocalisedRetrappingExcited,
    /// All events originating in ground states.
    Ground,
    /// All events originating in excited states.
    Excited,
    /// All recombination events.
    Recombination,
    /// All retrapping events.
    Retrapping,
    /// Irradiation-driven filling events.
    Filling,
}

impl EventSeries {
    /// Every event series in CSV column order.
    pub const ALL: [Self; 13] = [
        Self::LocalisedRecombinationGround,
        Self::LocalisedRecombinationExcited,
        Self::DelocalisedRecombinationGround,
        Self::DelocalisedRecombinationExcited,
        Self::LocalisedRetrappingGround,
        Self::LocalisedRetrappingExcited,
        Self::DelocalisedRetrappingGround,
        Self::DelocalisedRetrappingExcited,
        Self::Ground,
        Self::Excited,
        Self::Recombination,
        Self::Retrapping,
        Self::Filling,
    ];

    /// Human-readable label used in the plot legend.
    pub fn label(self) -> &'static str {
        match self {
            Self::LocalisedRecombinationGround => "Localised recombination (ground)",
            Self::LocalisedRecombinationExcited => "Localised recombination (excited)",
            Self::DelocalisedRecombinationGround => "Delocalised recombination (ground)",
            Self::DelocalisedRecombinationExcited => "Delocalised recombination (excited)",
            Self::LocalisedRetrappingGround => "Localised retrapping (ground)",
            Self::LocalisedRetrappingExcited => "Localised retrapping (excited)",
            Self::DelocalisedRetrappingGround => "Delocalised retrapping (ground)",
            Self::DelocalisedRetrappingExcited => "Delocalised retrapping (excited)",
            Self::Ground => "Ground-state events",
            Self::Excited => "Excited-state events",
            Self::Recombination => "Recombination",
            Self::Retrapping => "Retrapping",
            Self::Filling => "Filling",
        }
    }

    /// Extract this series from an event row.
    fn value(self, row: &AverageEventRow) -> Float {
        match self {
            Self::LocalisedRecombinationGround => row.localised_recombination_ground_count,
            Self::LocalisedRecombinationExcited => row.localised_recombination_excited_count,
            Self::DelocalisedRecombinationGround => row.delocalised_recombination_ground_count,
            Self::DelocalisedRecombinationExcited => row.delocalised_recombination_excited_count,
            Self::LocalisedRetrappingGround => row.localised_retrapping_ground_count,
            Self::LocalisedRetrappingExcited => row.localised_retrapping_excited_count,
            Self::DelocalisedRetrappingGround => row.delocalised_retrapping_ground_count,
            Self::DelocalisedRetrappingExcited => row.delocalised_retrapping_excited_count,
            Self::Ground => row.ground_count,
            Self::Excited => row.excited_count,
            Self::Recombination => row.recombination_count,
            Self::Retrapping => row.retrapping_count,
            Self::Filling => row.filling_count,
        }
    }
}



/// Paths written by [`plot_default_results`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneratedPlots {
    /// Mean filling against time, with a standard-deviation band.
    pub mean_fill_vs_time: PathBuf,
    /// Mean filling against temperature, with a standard-deviation band.
    pub mean_fill_vs_temperature: PathBuf,
    /// Median filling against time, with an interquartile band.
    pub median_fill_vs_time: PathBuf,
    /// Median filling against temperature, with an interquartile band.
    pub median_fill_vs_temperature: PathBuf,
    /// Temperature against time.
    pub temperature_vs_time: PathBuf,
    /// Aggregate event frequencies against time.
    pub events_vs_time: PathBuf,
}

/// Consolidated fill and event results loaded from the two simulation CSVs.
#[derive(Debug, Clone)]
pub struct SimulationResults {
    fill_path: PathBuf,
    event_path: PathBuf,
    fill_rows: Vec<ContinuousValueRow>,
    event_rows: Vec<AverageEventRow>,
}

/// One completed event bin represented by its centre, width, and frequency row.
#[derive(Debug, Clone, Copy)]
struct EventPlotBin<'a> {
    /// Midpoint between the bin's left and right edges.
    centre: TimeFloat,
    /// Difference between the bin's right and left edges.
    width: TimeFloat,
    /// Average counts recorded at the bin's right edge in the CSV.
    row: &'a AverageEventRow,
}

impl SimulationResults {
    /// Load and validate one fill CSV and one event CSV.
    pub fn from_csv(
        fill_csv: impl AsRef<Path>,
        event_csv: impl AsRef<Path>,
    ) -> Result<Self, PlotError> {
        let fill_path = fill_csv.as_ref().to_path_buf();
        let event_path = event_csv.as_ref().to_path_buf();
        let fill_rows = read_csv::<ContinuousValueRow>(&fill_path)?;
        let event_rows = read_csv::<AverageEventRow>(&event_path)?;
        validate_fill_rows(&fill_path, &fill_rows)?;
        validate_event_rows(&event_path, &event_rows)?;

        Ok(Self {
            fill_path,
            event_path,
            fill_rows,
            event_rows,
        })
    }

    /// Borrow the loaded continuous fill rows.
    pub fn fill_rows(&self) -> &[ContinuousValueRow] {
        &self.fill_rows
    }

    /// Borrow the loaded averaged event-count rows.
    pub fn event_rows(&self) -> &[AverageEventRow] {
        &self.event_rows
    }

    /// Increase the event-bin width and return counts in the wider bins.
    ///
    /// If old and new boundaries do not align, a source bin's count is split
    /// in proportion to its overlap with each new bin. This preserves the
    /// integrated count, including in a final partial bin. Plotting subsequently
    /// divides each result by its new, actual width.
    pub fn rebin_events(
        &self,
        new_bin_width: TimeFloat,
    ) -> Result<Vec<AverageEventRow>, PlotError> {
        rebin_event_rows(&self.event_path, &self.event_rows, new_bin_width)
    }

}







/// Load both result CSVs and create six standard PNG plots.
///
/// `event_bin_width` can request smoothing of the event plot. The four fill
/// plots cover both available central statistics and both horizontal axes.
pub fn plot_default_results(
    fill_csv: impl AsRef<Path>,
    event_csv: impl AsRef<Path>,
    output_directory: impl AsRef<Path>,
    event_bin_width: Option<TimeFloat>,
) -> Result<GeneratedPlots, PlotError> {
    let results = SimulationResults::from_csv(fill_csv, event_csv)?;
    let output_directory = output_directory.as_ref();
    fs::create_dir_all(output_directory).map_err(|source| PlotError::CreateDirectory {
        path: output_directory.to_path_buf(),
        source,
    })?;

    let to_plot = PlotWindow::new(
        output_directory.join("mean_fill_vs_time.png"),
        "Time",
        "second",
        "Fill",
        "None",
        true,
        "sd"
    )?;
    to_plot.plot_fill(&results.fill_rows, "Mean fill vs Time")?;
    let to_plot = PlotWindow::new(
        output_directory.join("mean_fill_vs_temperature.png"),
        "Temperature",
        "Kelvin",
        "Fill",
        "None",
        true,
        "sd"
    )?;
    to_plot.plot_fill(&results.fill_rows, "Mean fill vs Temperature")?;

   let to_plot = PlotWindow::new(
        output_directory.join("median_fill_vs_time.png"),
        "Time",
        "second",
        "Fill",
        "None",
        false,
        "sd"
    )?;
    to_plot.plot_fill(&results.fill_rows, "Median fill vs Time")?;

    let to_plot = PlotWindow::new(
        output_directory.join("median_fill_vs_temperature.png"),
        "Temperature",
        "Kelvin",
        "Fill",
        "None",
        false,
        "sd"
    )?;
    to_plot.plot_fill(&results.fill_rows, "Meanian fill vs Temperature")?;

    let paths = GeneratedPlots {
        mean_fill_vs_time: output_directory.join("mean_fill_vs_time.png"),
        mean_fill_vs_temperature: output_directory.join("mean_fill_vs_temperature.png"),
        median_fill_vs_time: output_directory.join("median_fill_vs_time.png"),
        median_fill_vs_temperature: output_directory.join("median_fill_vs_temperature.png"),
        temperature_vs_time: output_directory.join("temperature_vs_time.png"),
        events_vs_time: output_directory.join("events_vs_time.png"),
    };
    
    let to_plot = PlotWindow::new(
        output_directory.join("temperature_vs_time.png"),
        "Time",
        "second",
        "Temperature",
        "Kelvin",
        false,
        "sd"
    )?;
    to_plot.plot_temperature_vs_time(&results.fill_rows, "Time vs Temperature")?;
     let to_plot = PlotWindow::new(
        output_directory.join("time_vs_temperature.png"),
        "Temperature",
        "Kelvin",
        "Time",
        "second",
        false,
        "sd"
    )?;
    to_plot.plot_temperature_vs_time(&results.fill_rows, "Temperature vs Time ")?;
    // results.plot_temperature_vs_time(&paths.temperature_vs_time, options)?;
    // results.plot_events_vs_time(
    //     &paths.events_vs_time,
    //     &[
    //         EventSeries::Recombination,
    //         EventSeries::Retrapping,
    //         EventSeries::Filling,
    //     ],
    //     event_bin_width,
    //     options,
    // )?;
    Ok(paths)
}

/// Deserialize every record in a CSV file, trimming accidental header spaces.
fn read_csv<T: for<'de> serde::Deserialize<'de>>(path: &Path) -> Result<Vec<T>, PlotError> {
    let mut reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(path)
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })?;
    reader
        .deserialize()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|source| PlotError::Csv {
            path: path.to_path_buf(),
            source,
        })
}

/// Validate all mandatory fill fields without rejecting absent legacy quartiles.
fn validate_fill_rows(path: &Path, rows: &[ContinuousValueRow]) -> Result<(), PlotError> {
    if rows.is_empty() {
        return invalid(path, "the file contains no fill rows");
    }
    for (index, row) in rows.iter().enumerate() {
        let mandatory = [
            row.time,
            row.temperature,
            row.fill,
            row.fill_standard_deviation,
            row.fill_median,
            row.fill_quantile_0_1,
            row.fill_quantile_0_9,
        ];
        if mandatory.iter().any(|value| !value.is_finite()) {
            return invalid(
                path,
                format!("row {} contains a non-finite value", index + 2),
            );
        }
        if row.fill_standard_deviation < 0.0 {
            return invalid(
                path,
                format!("row {} has a negative standard deviation", index + 2),
            );
        }
        let quartiles = (row.fill_quantile_0_25, row.fill_quantile_0_75);
        if quartiles.0.is_finite() != quartiles.1.is_finite()
            || (quartiles.0.is_finite() && quartiles.0 > quartiles.1)
        {
            return invalid(
                path,
                format!("row {} has invalid quartile limits", index + 2),
            );
        }
    }
    Ok(())
}

/// Validate monotonic event-bin endpoints and finite, non-negative frequencies.
fn validate_event_rows(path: &Path, rows: &[AverageEventRow]) -> Result<(), PlotError> {
    if rows.is_empty() {
        return invalid(path, "the file contains no event rows");
    }
    for (index, row) in rows.iter().enumerate() {
        if !row.time.is_finite() {
            return invalid(
                path,
                format!("row {} contains a non-finite time", index + 2),
            );
        }
        if index > 0 && row.time <= rows[index - 1].time {
            return invalid(path, "event-bin times must be strictly increasing");
        }
        if EventSeries::ALL
            .iter()
            .map(|column| column.value(row))
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return invalid(
                path,
                format!("row {} contains an invalid event count", index + 2),
            );
        }
    }
    Ok(())
}

/// Convert right-edge CSV rows into bins located at their temporal midpoints.
fn event_plot_bins(rows: &[AverageEventRow]) -> Vec<EventPlotBin<'_>> {
    rows.windows(2)
        .map(|window| {
            let left = window[0].time;
            let right = window[1].time;
            EventPlotBin {
                centre: left + (right - left) / 2.0,
                width: right - left,
                row: &window[1],
            }
        })
        .collect()
}

/// Convert an averaged bin count to an event frequency.
fn event_frequency(bin: &EventPlotBin<'_>, series: EventSeries) -> Float {
    debug_assert!(bin.width > 0.0);
    series.value(bin.row) / bin.width
}

/// Rebin average counts while conserving their sum across destination bins.
fn rebin_event_rows(
    path: &Path,
    rows: &[AverageEventRow],
    new_bin_width: TimeFloat,
) -> Result<Vec<AverageEventRow>, PlotError> {
    if !new_bin_width.is_finite() || new_bin_width <= 0.0 {
        return invalid(path, "the new event-bin width must be finite and positive");
    }
    if rows.len() < 2 {
        return invalid(path, "at least two event rows are required for rebinning");
    }
    let source_width = rows
        .windows(2)
        .map(|window| window[1].time - window[0].time)
        .fold(0.0_f64, Float::max);
    let tolerance = source_width.abs().max(new_bin_width.abs()) * 1e-12;
    if new_bin_width + tolerance < source_width {
        return invalid(
            path,
            format!(
                "the new event-bin width ({new_bin_width}) is smaller than the source width ({source_width})"
            ),
        );
    }

    let origin = rows[0].time;
    let final_time = rows.last().expect("at least two rows were checked").time;
    let mut rebinned = vec![zero_event_row(origin)];
    let mut bin_start = origin;
    while bin_start < final_time {
        let bin_end = (bin_start + new_bin_width).min(final_time);
        let mut rebinned_counts = [0.0; 13];

        let mut index = rows.partition_point(|row| row.time <= bin_start).max(1);
        while index < rows.len() {
            let segment_start = rows[index - 1].time;
            let segment_end = rows[index].time;
            let overlap_start = segment_start.max(bin_start);
            let overlap_end = segment_end.min(bin_end);
            if overlap_end > overlap_start {
                let overlap = overlap_end - overlap_start;
                let source_bin_width = segment_end - segment_start;
                let overlap_fraction = overlap / source_bin_width;
                for (value, column) in rebinned_counts.iter_mut().zip(EventSeries::ALL) {
                    *value += column.value(&rows[index]) * overlap_fraction;
                }
            }
            if segment_end >= bin_end {
                break;
            }
            index += 1;
        }
        rebinned.push(event_row_from_values(bin_end, rebinned_counts));
        bin_start = bin_end;
    }
    Ok(rebinned)
}

/// Construct an all-zero event row at an arbitrary starting time.
fn zero_event_row(time: TimeFloat) -> AverageEventRow {
    event_row_from_values(time, [0.0; 13])
}

/// Convert values in [`EventSeries::ALL`] order back into a typed row.
fn event_row_from_values(time: TimeFloat, values: [Float; 13]) -> AverageEventRow {
    AverageEventRow {
        time,
        localised_recombination_ground_count: values[0],
        localised_recombination_excited_count: values[1],
        delocalised_recombination_ground_count: values[2],
        delocalised_recombination_excited_count: values[3],
        localised_retrapping_ground_count: values[4],
        localised_retrapping_excited_count: values[5],
        delocalised_retrapping_ground_count: values[6],
        delocalised_retrapping_excited_count: values[7],
        ground_count: values[8],
        excited_count: values[9],
        recombination_count: values[10],
        retrapping_count: values[11],
        filling_count: values[12],
    }
}


/// Convert a validation failure into the common error type.
fn invalid<T>(path: &Path, message: impl Into<String>) -> Result<T, PlotError> {
    Err(PlotError::InvalidData {
        path: path.to_path_buf(),
        message: message.into(),
    })
}




/// Erase Plotters' backend-specific error type while retaining its message.
fn draw_result<T, E: fmt::Debug>(path: &Path, result: Result<T, E>) -> Result<T, PlotError> {
    result.map_err(|error| PlotError::Draw {
        path: path.to_path_buf(),
        message: format!("{error:?}"),
    })
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::fs;
//     use std::sync::atomic::{AtomicU64, Ordering};
//     use std::time::{SystemTime, UNIX_EPOCH};

//     static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

//     fn temporary_directory() -> PathBuf {
//         let unique = SystemTime::now()
//             .duration_since(UNIX_EPOCH)
//             .unwrap()
//             .as_nanos();
//         let sequence = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
//         let path = std::env::temp_dir().join(format!(
//             "mcrustlum_plotting_{}_{}_{}",
//             std::process::id(),
//             unique,
//             sequence
//         ));
//         fs::create_dir(&path).unwrap();
//         path
//     }

//     fn test_results() -> (PathBuf, SimulationResults) {
//         let directory = temporary_directory();
//         let fill = directory.join("average_fill.csv");
//         let events = directory.join("average_event.csv");
//         fs::write(
//             &fill,
//             concat!(
//                 "time,temperature,fill,fill_standard_deviation,fill_median,fill_quantile_0_1,fill_quantile_0_9,fill_quantile_0_25,fill_quantile_0_75\n",
//                 "0,300,0.2,0.02,0.19,0.15,0.25,0.17,0.22\n",
//                 "0.1,310,0.4,0.03,0.39,0.34,0.46,0.36,0.42\n",
//                 "0.2,320,0.6,0.04,0.59,0.52,0.68,0.55,0.63\n",
//                 "0.3,330,0.8,0.05,0.79,0.70,0.88,0.74,0.84\n",
//             ),
//         )
//         .unwrap();
//         fs::write(
//             &events,
//             concat!(
//                 "time,localised_recombination_ground,localised_recombination_excited,delocalised_recombination_ground,delocalised_recombination_excited,localised_retrapping_ground,localised_retrapping_excited,delocalised_retrapping_ground,delocalised_retrapping_excited,ground,excited,recombination,retrapping,filling_count\n",
//                 "0,0,0,0,0,0,0,0,0,0,0,0,0,0\n",
//                 "0.1,1,0,0,0,0,0,0,0,1,0,1,0,2\n",
//                 "0.2,3,0,0,0,0,0,0,0,3,0,3,0,4\n",
//                 "0.3,5,0,0,0,0,0,0,0,5,0,5,0,6\n",
//             ),
//         )
//         .unwrap();
//         let results = SimulationResults::from_csv(&fill, &events).unwrap();
//         (directory, results)
//     }

//     #[test]
//     fn reads_both_result_files() {
//         let (directory, results) = test_results();
//         assert_eq!(results.fill_rows().len(), 4);
//         assert_eq!(results.event_rows().len(), 4);
//         fs::remove_dir_all(directory).unwrap();
//     }

//     #[test]
//     fn wider_bins_conserve_counts_and_keep_a_partial_tail() {
//         let (directory, results) = test_results();
//         let rows = results.rebin_events(0.2).unwrap();
//         assert_eq!(rows.len(), 3);
//         assert!((rows[1].recombination_count - 4.0).abs() < 1e-12);
//         assert!((rows[1].filling_count - 6.0).abs() < 1e-12);
//         assert!((rows[2].recombination_count - 5.0).abs() < 1e-12);
//         assert!((rows[2].time - 0.3).abs() < 1e-12);
//         let bins = event_plot_bins(&rows);
//         assert!((bins[0].centre - 0.1).abs() < 1e-12);
//         assert!((bins[0].width - 0.2).abs() < 1e-12);
//         assert!((bins[1].centre - 0.25).abs() < 1e-12);
//         assert!((bins[1].width - 0.1).abs() < 1e-12);
//         assert!((event_frequency(&bins[0], EventSeries::Recombination) - 20.0).abs() < 1e-12);
//         assert!((event_frequency(&bins[1], EventSeries::Recombination) - 50.0).abs() < 1e-12);
//         fs::remove_dir_all(directory).unwrap();
//     }

//     #[test]
//     fn rebinning_splits_counts_at_unaligned_boundaries() {
//         let (directory, results) = test_results();
//         let rows = results.rebin_events(0.15).unwrap();
//         assert_eq!(rows.len(), 3);
//         assert!((rows[1].recombination_count - 2.5).abs() < 1e-12);
//         assert!((rows[2].recombination_count - 6.5).abs() < 1e-12);
//         assert!((rows[1].recombination_count + rows[2].recombination_count - 9.0).abs() < 1e-12);
//         fs::remove_dir_all(directory).unwrap();
//     }

//     #[test]
//     fn event_bins_are_plotted_at_their_centres() {
//         let (directory, results) = test_results();
//         let bins = event_plot_bins(results.event_rows());
//         assert_eq!(bins.len(), 3);
//         assert!((bins[0].centre - 0.05).abs() < 1e-12);
//         assert!((bins[1].centre - 0.15).abs() < 1e-12);
//         assert!((bins[2].centre - 0.25).abs() < 1e-12);
//         assert!((bins[0].width - 0.1).abs() < 1e-12);
//         assert!((event_frequency(&bins[0], EventSeries::Recombination) - 10.0).abs() < 1e-12);
//         fs::remove_dir_all(directory).unwrap();
//     }

//     #[test]
//     fn narrower_or_non_positive_bins_are_rejected() {
//         let (directory, results) = test_results();
//         assert!(results.rebin_events(0.05).is_err());
//         assert!(results.rebin_events(0.0).is_err());
//         fs::remove_dir_all(directory).unwrap();
//     }

//     #[test]
//     fn creates_the_default_png_set() {
//         let (directory, _) = test_results();
//         let output = directory.join("plots");
//         let plots = plot_default_results(
//             directory.join("average_fill.csv"),
//             directory.join("average_event.csv"),
//             &output,
//             Some(0.2),
//         )
//         .unwrap();
//         for path in [
//             plots.mean_fill_vs_time,
//             plots.mean_fill_vs_temperature,
//             plots.median_fill_vs_time,
//             plots.median_fill_vs_temperature,
//             plots.temperature_vs_time,
//             plots.events_vs_time,
//         ] {
//             assert!(fs::metadata(path).unwrap().len() > 0);
//         }
//         fs::remove_dir_all(directory).unwrap();
//     }
// }
