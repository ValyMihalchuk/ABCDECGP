#####################################
cgp_parameters_file = "parameters.txt"
train_big_file='train_big.data'
train_tiny_file='train_tiny.data'
test_big_file='test_big.data'
test_tiny_file='test_tiny.data'

cgp_exe_name = 'D:/test_pr/Project1/Release/project1.exe'
#####################################


import sys
import os
import subprocess
import configparser
import numpy as np


text_file = open("D:/ABCDECGP/Connection/parameters_pattern.txt", "r")
data = text_file.read()
text_file.close()
def InsertParamsInCGP(parameters_file_cgp, deep_ini_file='../ABCDE/deepdeep.ini',  section_name = 'test' ):
    config = configparser.ConfigParser()
    config.read(deep_ini_file)
    params = dict(config.items('test'))
    params=list(params.values())


    text_file = open(parameters_file_cgp, "w")
    new_data = data.format(*params)
    text_file.write(new_data)
    text_file.close()

def split_test_train(current_db, FilenameForTrain, FilenameForTest):
    X= [x[1:-1] for x in current_db]
    y = [x[-1] for x in current_db]

    X = np.array(X)
    y = np.array(y)
    
    from sklearn.model_selection import KFold 
    test_fold = KFold(n_splits=7, shuffle = True)
    test_fold.get_n_splits(X)
    train_index, test_index = list(test_fold.split(X))[0]

    test_X = X[test_index]
    test_y = y[test_index]
    
    
    inputs = len(test_X[0])
    outputs = 1

    with open(FilenameForTest, 'w') as csvfile:
        csvfile.write(str(inputs) + "," + str(outputs)+","+ str(len(test_X)) +",\n")
        for xx, yy in zip(test_X, test_y):
            zz = [*xx, yy]
            line = str(zz)[1:-1]+",\n"
            csvfile.write(line)
            
    train_X = X[train_index]
    train_y = y[train_index]

    with open(FilenameForTrain, 'w') as csvfile:
        csvfile.write(str(inputs) + "," + str(outputs)+","+ str(len(train_X)) +",\n")
        for xx, yy in zip(train_X, train_y):
            zz = [*xx, yy]
            line = str(zz)[1:-1]+",\n"
            csvfile.write(line)


arg = sys.argv[1]
deep_ini_file = arg[len("--default-name="):]
InsertParamsInCGP(deep_ini_file[deep_ini_file.index("."):]+cgp_parameters_file, deep_ini_file)


current_db_big = np.load('D:/ABCDECGP/Connection/start_seedling_complete_maturing_big_all.npy').astype(float)
current_db_tiny = np.load('D:/ABCDECGP/Connection/start_seedling_complete_maturing_tiny_all.npy').astype(float)
split_test_train(current_db_big, train_big_file, test_big_file)
split_test_train(current_db_tiny, train_tiny_file, test_tiny_file)


subprocess.call(cgp_exe_name + ' 1 '+ cgp_parameters_file +' '+train_big_file+' '+train_tiny_file, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
proc = subprocess.Popen([cgp_exe_name, '3', cgp_parameters_file, train_big_file,train_tiny_file, test_big_file, test_tiny_file], stdout=subprocess.PIPE)

val = proc.communicate()[0].decode("utf-8")

import numpy as np
errors = [float(s.split(",")[1]) for s in val.split("\n") if len(s.split(",")) == 3]
error = np.mean(errors)
print(":"+str(error))
#split_val = [float(s) for s in val.split()]
#error = ...
#print(":"+error)