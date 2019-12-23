import pandas as pd
import os
from os.path import isfile, join
from joblib import dump

# reading all files we have

def creating_RF_classificator(n_estimators=300, max_depth=15, bootstrap='True', max_features='sqrt', n_jobs=-1):
    """
        Making Random Forest classificator and storing it at 'saved_sklean_models' dir
        If you want to change classificator parameters, please see scikit-learn documentation
    Args:
        n_estimators: int, default = 300
        max_depth: int, default = 15
        bootstrap: str, default = 'True'
        max_features: str, default = 'sqrt'
        n_jobs: int, default = -1

    Returns:

    """
    # reading files from csv_data folder
    csv_files = [f for f in os.listdir('csv_data') if isfile(join('csv_data', f))]
    for i in range(len(csv_files)):
        csv_files[i] = 'csv_data/' + csv_files[i]

    # reading first file
    Proteins_data = pd.read_csv(csv_files.pop(0), index_col=0)

    # reading rest of the files
    for x in csv_files:
        res_data = pd.read_csv(x,  index_col=0)
        Proteins_data = pd.concat([Proteins_data, res_data], ignore_index=True)


    # Be careful with your data

    labels_categorized = Proteins_data[Proteins_data.columns[0]].values
    feature_matrix = Proteins_data[Proteins_data.columns[2:403]].values
    del Proteins_data


    # Making labels not categorized (slightly)

    labels = []

    for x in labels_categorized:
        if x == 'human_proteome':
            labels.append(1)
        else:
            labels.append(0)


    # Importing SKlearn magic

    from sklearn.model_selection import train_test_split
    from sklearn.metrics import f1_score, accuracy_score
    from sklearn import preprocessing
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import classification_report

    # Making 2 Train-Test split

    train_feature_matrix, test_feature_matrix,\
    train_labels, test_labels = train_test_split(feature_matrix, labels, test_size=0.2, random_state=42)

    # Fitting our classificator

    forest = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, bootstrap=bootstrap,
                                    max_features=max_features, n_jobs=n_jobs)
    forest.fit(train_feature_matrix, train_labels)

    # showing score for user

    y_pred = forest.predict(test_feature_matrix)
    print(classification_report(test_labels, y_pred, labels=[0, 1]))

    # saving classificator in 'saved_sklearn_models' folder
    print(r'Saving classificator in ../saved_sklearn_models folder')

    if not os.path.isdir('saved_sklearn_models'):
        os.mkdir('saved_sklearn_models')

    dump(forest, '../saved_sklearn_models/random_forest.joblib')
    del forest