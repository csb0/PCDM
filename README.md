## Pupil common drive model (PCDM) toolbox ########################
A pupillometry toolbox for MATLAB.

## Authors ########################################################
Charlie S. Burlingham &lt;<cs.burlingham@gmail.com>&gt; and Saghar Mirbagheri &lt;<sagharm@uw.edu>&gt;

## Use ############################################################
To use the toolbox, navigate to PCDM directory in MATLAB, run the command addpath(genpath(pwd)). Then, create an input data struct “in” containing your eye data, task event timing, and trial types (see DataAnalysis.m for details on formatting) and run fitModel. The toolbox takes in the raw gaze position and pupil area time series. The code was tested on data collected with an EyeLink eye tracker, but should work with data from other systems. The function fitModel.m returns a struct “f” containing the model’s parameter estimates, prediction of the task-evoked pupil response, and goodness of fit (R^2), as well as plots of the model fits and parameter estimates for each run/block of data.

If you are using pupil data that already has blinks removed, set “op.blinkInterpolation” to 0 in the preamble of the function dataAnalysis.m

For estimating gain correctly, 5 minutes of eye data suffices. If you’d like to make a relatively noiseless estimate of the generator function at a very high sampling rate (e.g., 500 or 1000 Hz), it is important to use a large amount of data (around 200 minutes of eye data). If you want to estimate the generator function at a lower sampling rate, simply downsample your input data, and you will require much less data for good estimation of the generator function. This is common practice in analyzing saccade rates and will often be desirable because the underlying neural process is very unlikely to occur at more than 200 Hz. So don’t downsample your inpur data below 200 Hz, to be conservative.

## Variants/Extensions ############################################

### Locking the pupil response to task events other than trial onset
If you’d like to examine the pupil response to other task events (e.g., button press, cues), simply create a column vector with the start times of these events (in units of samples) called “eventTimes” and replace each instance of “in.startInds{kk}(:,1)” on lines 100 and 123 of dataAnalysis.m with it. If you are analyzing multiple trial types, make the variable a cell array with one cell per trial type containing the vector of event times. This will estimate an event-locked pupil response and saccade rate function.

If you’d like to do this with multiple events simultaneously and take into account their shared variance, you must create a block design matrix with a block for each event separately for deconvolving the saccade rate function and pupil response (see Knapen et al 2016, *PLOS One* for example).

### Estimating multiple gains per trial
If you have some reason to think the gain changes within a trial, you can estimate multiple gains per trial by …… [add]

## Citing #########################################################

To cite PCDM, please reference the following:
* Burlingham CS*, Mirbagheri S*, Heeger DJ (2022). A unified model of the task-evoked pupil response. Science Advances. * Equal Authors

## References #####################################################

* Engbert R, Kliegl R. (2003) Microsaccades uncover the orientation of covert attention. *Vision Research*, 43: 1035-1045.
* Engbert R, Mergenthaler K (2006) Microsaccades are triggered by low retinal image slip. *Proceedings of the National Academy of Sciences of the United States of America*, 103: 7192-7197.
* Gardner J, Larsson J. (2012). MGL Toolbox. In: Gardner Research Unit | mgl:overview (Online). http://gru.brain.riken.jp/doku.php/mgl/overview
* Mathôt, S (2013). A simple way to reconstruct pupil size during eye blinks. *FigShare*. 10.6084/m9.figshare.688001
* Denison RN, Parker, JA, Carrasco, M (2020). Modeling pupil responses to rapid sequential events. *Behav Res Methods.* doi: 10.3758/s13428-020-01368-6.

## License ########################################################

This README file is part of the PCDM library.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
Affero General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>
