# eyetrack
Perform dual-angle video-oculography calibration for magnetic eye tracking

Hannah Payne <hpayne@stanford.edu>
Raymond Lab
2017


Matlab code is included to perform three functions:
1. Collect dual-angle video-oculography data from 2 cameras
2. Analyze the video-oculography data to extract pupil and CR positions
3. Calibrate simultaneously recorded (in Spike2, CED, Cambridge) magnetic eye tracking data


Example data is included. Instructions for analyzing example data:
1. In Matlab, run runEYETRACK.m (Add it to the path when prompted. It will add subfolders to the path automatically.)
2. Select the second option in the pop-up: analyze video data.
3. Select the folder you'd like to analyze (the 20140821_H2_calibL folder)
4. It will prompt to input video analysis settings -- just leave at default
5. It will analyze the video and then apply the calibration
6. The final calibration file, *.mat, contains a number of results. The most importantare scaleCh1 and scaleCh2. One of these will be 0. 
The other is the calibration factor to apply to that magnetic sensor channel, for future experiments with the video cameras removed.





