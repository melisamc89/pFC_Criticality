'''
Created on mrt 18 11:41

author: Melisa

This script converts npy files from calcium analysis to .mat data.
It converts the calcium traces files as well as the timeline files.

'''

from scipy.io import savemat
import numpy as np
import os
import pickle

data_dir = '/home/melisa/Documents/criticality/data/calcium_activity_day_wise/'
new_data_dir = '/home/melisa/Documents/criticality/data/calcium_activity/'

time_dir = '/home/melisa/Documents/criticality/data/timeline/'

mouse = [32363,32364,32365]
mouse = [56165, 56166]
trials = [1,6,11,16,21]
sessions = [[1,2],[1,2],[2,3]]
sessions = [[1,2,4],[1,3]]


#mouse_32363_session_1_trial_1_v1.4.20.3.0.1.1.0

i = 0
for mouse_id in mouse:
    print(sessions[i])
    for session in sessions[i]:
        print(session)
        for days in trials:
            file_name = 'mouse_' + f'{mouse_id}' + '_session_' + f'{session}' + '_trial_' + f'{days}' + '_v1.4.20.3.0.1.1.0'
            file_dir = data_dir + file_name + '.npy'
            data = np.load(file_dir)
            print(data.shape)

            time_file_session_1 = 'mouse_' + f'{mouse_id}' + '_session_' + f'{session}' + '_trial_1_v' + '1.4.' + f'{1}' + \
                                  '.' + f'{0}' + '_10.pkl'
            timeline_file = open(time_dir + time_file_session_1, 'rb')
            timeline_info = pickle.load(timeline_file)
            timeline_1 = np.zeros(len(timeline_info) + 1)
            for i in range(len(timeline_info)):
                timeline_1[i] = timeline_info[i][1]
            timeline_1[len(timeline_info)] = data.shape[1]
            savemat(new_data_dir + file_name + '.mat' , {'calcium_trace': data, 'timeline':timeline_1})
    i=i+1
