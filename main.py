'''
READ ME
    Run this file ONLY after running loading_matlab_file.py first!!
    Put in the paths to the matlab file and the place to store the data in loading_matlab_file.py first!!

DESCRIPTION
    This script puts out all the data needed for the BSc Thesis of Nubia Middelkoop.
    This code uses 3 classes from the file ClassMosquito.py:
        Track (to load single tracks), Trials (to load single trials) and Dataset (to load the whole dataset).

PARAMETERS
    The class Dataset does not take any parameters
    The class Trial takes only one integer as an input: The trial nuber
    The class Track takes two integers parameters: the trial number and the track number
        Some tracks don't exist in this dataset,then an error message arises.

'''

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('Hi reader, this is main!\n Run this script AFTER loading the data into the correct folders/directory with the loading_matlab_file.py.\n Close the plots to continue in script\n')


import ClassMosquito


# Run loading_matlab_file.py first!


# Initializing the dataset
dataset = ClassMosquito.Dataset()
# Initializing the highlighted trial
trial = ClassMosquito.Trial(1)
# Initializing the highlighted track
track = ClassMosquito.Track(1, 2)


#  Basic dataset analysis

#  --> Amount and percentage of landing mosquitoes and captured mosquitoes
landing_count = dataset.countLandingPoints()
capturing_count = dataset.countCapturingPoints()
total_tracks = sum(dataset.getNumTracksPerTrial())
print(f'\nTotal tracks = {total_tracks}.')
print(f'In the whole dataset {landing_count} mosquitoes are landing. ({round(landing_count/total_tracks*100, 2)} % of the tracks end in landing)')
print(f'From this testing: In the whole dataset {capturing_count} mosquitoes are captured. ({round(capturing_count/total_tracks*100, 2)} % of the tracks end in capture)')
print(f'From the capture rate of the article: In the whole dataset 1335 mosquitoes are captured. ({round(1335/total_tracks*100, 2)} % of the tracks end in capture)')


# --> Normal distribution of number of tracks per trial
dataset.testNormalDistribution(dataset.getNumTracksPerTrial())
print(f'Confidence interval of the number of tracks per trial: {dataset.testConfidenceInterval(dataset.getNumTracksPerTrial())}\n')

# --> average duration of all the tracks
print(f'The average duration all the tracks is {round(dataset.getAvarageLengthTracks(), 2)} seconds.\n')
# --> Durations of all trials visualized
dataset.plotDuration()
# --> start and end times visualized in a violin plot
dataset.plotViolinStartEndTimes()

# Single track plotting
track.plotTrack()  # 3D
track.plotTheta2DTrack()  # 2D

# single trial plotting
trial.plotTrial()  # 3D
trial.plot2DTrial()  # 2D


# Boxplot statistics to see significant differences  (significance in terminal)

# --> Amount of captures vs landings per trial
dataset.plotBoxplotCatchesVsLandings()
# --> Caches per short-range cue conditions
dataset.plotBoxplotCatchesConditions()
# --> Landings per short-range cue conditions
dataset.plotBoxplotLandingCondition()


# Resting time statistics

# --> Average resting time
print(f'\nAverage resting time of the whole dataset = {dataset.getAvarageRestingTime()}')
# --> Median reating time
print(f'The median resting time of the whole dataset = {round(dataset.getMedianRestingTime()[1], 2)}')
print(f'With the first quarter being  = {round(dataset.getMedianRestingTime()[0], 2)} and the third being {round(dataset.getMedianRestingTime()[2], 2)}.')
# --> Longest resting time
print(f'The longest resting time of the whole dataset is = {round(dataset.getLongestRestingTime(), 2)}')

# Statistics above 0.5 seconds resting time
print('\nExcluding all the resting times below 0.5 gives the following average and median: ')
# --> Average resting time (above 0.5)
print(f'Average resting time (above 0.5 s resting times) = {dataset.getAvarageRestingTime(0.5)}')
# --> Median resting time (above 0.5)
print(f'The median resting time (above 0.5 s resting times) = {round(dataset.getMedianRestingTime(0.5)[1], 2)}')
print(f'With the first quarter being  = {round(dataset.getMedianRestingTime(0.5)[0], 2)} and the third being {round(dataset.getMedianRestingTime(0.5)[2], 2)}.\n')

# Resting time histogram (general)
dataset.plotHistogramRestingTime()
# --> Resting time histogram between 0 - 1 sec
dataset.plotHistogramRestingTimeZoomedIn(0, 1)
# --> Resting time histogram between 0 - 2 sec
dataset.plotHistogramRestingTimeZoomedIn(0, 2)
# --> Resting time histogram between 0 - 10 sec
dataset.plotHistogramRestingTimeZoomedIn(0, 10)


# Resting times amount below 0.5 sec
print(f'The amount of very small resting times (below 0.5 seconds) = {round(dataset.countSmallRestingTimes(0.5), 2)}, this is {round(dataset.countSmallRestingTimes(0.5)/dataset.countTotalRestingTimes(), 2)*100} %')
# Resting times amount below 1 sec
print(f'The amount of very small resting times (below 1 second) = {round(dataset.countSmallRestingTimes(1), 2)}, this is {round(dataset.countSmallRestingTimes(1)/dataset.countTotalRestingTimes()*100, 2)} %')

# Resting time histograms per conditions
dataset.plotHistogramRestingTimeCondition('without')
dataset.plotHistogramRestingTimeCondition('with_heat')
dataset.plotHistogramRestingTimeCondition('with_heat_water')
# Resting time violin plot with all conditions visualized
dataset.plotViolinplotRestingTimeCondition()


# Heatmaps

# Landing points
dataset.plotHeatmapLandingPointNormilized()

# Resting points (all)
dataset.plotHeatmapRestingPoints(title = 'Density Resting Points all.')
# --> Resting points (< 0.5s)
dataset.plotHeatmapRestingPoints(upper_time_boundary=0.5, title='Density Resting Points below 0.5 s.')
# --> Resting points (> 0.5s)
dataset.plotHeatmapRestingPoints(lower_time_boundary=0.5, title='Density Resting Points above 0.5 s.')
# --> Resting points (< 1s)
dataset.plotHeatmapRestingPoints(upper_time_boundary=1, title='Density Resting Points below 1 s.')
# --> Resting points (> 1s)
dataset.plotHeatmapRestingPoints(lower_time_boundary=1, title='Density Resting Points above 1 s.')

# Resting times (all)
dataset.plotHeatmapRestingTimes()
# --> Resting times (> 0.5s)
dataset.plotHeatmapRestingTimes(lower_time_boundary=0.5, title = 'bigger then 0.5s')
# --> Resting times (< 0.5s)
dataset.plotHeatmapRestingTimes(upper_time_boundary=0.5, title = 'smaller than 0.5s')
# --> Resting times (> 1s)
dataset.plotHeatmapRestingTimes(lower_time_boundary=1, title = 'bigger then 1s')
# --> Resting times (< 1s)
dataset.plotHeatmapRestingTimes(upper_time_boundary=1, title = 'smaller then 1s')

# --> Resting times per condition
dataset.plotHeatmapRestingTimePerCondition()


# Landing Again visualisation
trial.plotTrackLandingAgainTrial(8)
# Capture after landing visualisation
trial.plotTrackLandingToCaptureTrial(1)

# Sensitivity / association boxplot
dataset.plotBoxplotRadiusAssociations()  # Landing -- take-off pairing radius
dataset.plotBoxplotCountLandingTakeOffAccosiation()  # Landing -- take-off pairing counting

# Heatmap  Landing again
dataset.plotHeatmapLandingAgainProbability()
# Heatmap capture after landing
dataset.plotHeatmapLandingToCaptureProbability()


# Probability capture/land again after landing
percentage_land_again = dataset.calculatingPercentagesLandingAgain()
percentage_land_to_capture = dataset.calculatingPercentagesLandingToCapture()
print(f'\nPercentage of take offs that land again is {round(percentage_land_again, 2)} %\nPercentage of take offs that lead to capture {round(percentage_land_to_capture, 2)} %')

