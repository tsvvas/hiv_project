import pandas as pd
import os
from os.path import isfile, join
from joblib import dump
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report


def preprocessing_data(path):
    """
    Reading all csv_files from {your path} dir
    
    Args:
        path: str, path to dir with csv_files

    Returns:
        feature_matrix: array, array with features needed to train classificator
        labels_categorized: array, array with labels for objects
    """
    # reading files' paths from csv_data folder
    csv_files = [f for f in os.listdir(path) if isfile(join(path, f))]

    # changing paths
    for i in range(len(csv_files)):
        csv_files[i] = f'{path}/{csv_files[i]}'

    # making dataframe
    proteins_data = pd.DataFrame()

    # reading rest of the files
    for file in csv_files:
        res_data = pd.read_csv(file, index_col=0)

        # adding data into dataframe
        proteins_data = pd.concat([proteins_data, res_data], ignore_index=True)


    # Getting labels and feature_matrix from Proteins data
    labels_categorized = proteins_data[proteins_data.columns[0]].values
    feature_matrix = proteins_data[proteins_data.columns[2:]].values

    # Initializing labels
    labels = []

    if 'human_proteome' not in labels_categorized:
        raise ValueError('human_proteome is not in labels_categorized, be careful')

    # categorizing labels
    for x in labels_categorized:
        if x == 'human_proteome':
            labels.append(1)
        else:
            labels.append(0)
    
    # fixing feature_matrix
    if type(feature_matrix[0][0]) == str:
        n, m = feature_matrix.shape
        for i in range(n):
            for j in range(m):
                tmp = float(feature_matrix[i][j].replace('[','').replace(']',''))
                feature_matrix[i][j] = tmp
    
    return feature_matrix, labels


def making_classificator(feature_matrix, labels, name, n_estimators=350, max_depth=15, bootstrap='True',
                         max_features='sqrt', n_jobs=-1):
    """
        Making Random Forest classificator and storing it at 'saved_sklean_models' dir
        If you want to change classificator parameters, please see scikit-learn documentation
        Note: class of human proteins is 1, class of others' proteins is 0

    Args:
        feature_matrix: array, array with features needed to train classificator
        labels: array, array with labels connected to feature_matrix to train classificator
        name: str, name for the classificator to save
        n_estimators: int, default = 300
        max_depth: int, default = 15
        bootstrap: str, default = 'True'
        max_features: str, default = 'sqrt'
        n_jobs: int, default = -1
    """

    # Making 2 Train-Test split
    train_feature_matrix, test_feature_matrix, \
    train_labels, test_labels = train_test_split(feature_matrix, labels, test_size=0.15, random_state=42)

    # Fitting our classificator
    forest = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, bootstrap=bootstrap,
                                    max_features=max_features, n_jobs=n_jobs)
    forest.fit(train_feature_matrix, train_labels)

    # showing score for user
    y_pred = forest.predict(test_feature_matrix)
    print(classification_report(test_labels, y_pred, labels=[0, 1]))

    # saving classificator in 'saved_sklearn_models' folder
    print(r'Saving classificator in ../saved_sklearn_models/ folder')

    # making dir
    if not os.path.isdir('saved_sklearn_models'):
        os.mkdir('saved_sklearn_models')
        print('saved_sklearn_models directory is created')

    # saving model
    dump(forest, f'saved_sklearn_models/{name}.joblib')

    return
