Contains code and example data sets related to chapter 4 of my dissertation: Behavior & Habitat Use of Three Sportfish in a Restored Seascape
Example datasets are from red drum tracking data collected at Northeast Point aux Pins (NEPaP), off the coast of Bayou la Batre, Alabama, USA

Run the code in the following order:
(1) DataPrep_SplitTracks&Interpolate.R
  This code requires Theo Michelot's split tracks function: splittracks_Michelot22.R
  And goes with the example dataset: habitats_SO368b.28093_NEPaP.csv
  A series of fine-scale acoustic telemetry positions triangulated for a 368 mm red drum

(2) HMM_NEPaP_RedDrum.R
  Fits HMMs with various numbers of behavior states and underlying data distributions to interpolated tracks of red drum positions
  Goes with the example data set: RedDrumHMMData.csv
  The split and interpolated tracking data from all red drum tracked at NEPaP

(3) SO_NEPaP_GLMM_habdielbehav.R
  Fits a multinomial mixed effects model to determine the effect habitat type and diel period on red drum behavior, as inferred using a HMM
  Goes with the example dataset: habitats_NEPaP_SO_HMM3wc_diel.csv
  In which triangulated red drum positions have been assigned a behavior state, a habitat type, and diel period
