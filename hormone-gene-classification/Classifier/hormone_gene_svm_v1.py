#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
import fasttext
import csv
import numpy as np
import random
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
import math
import matplotlib.pyplot as plt
from matplotlib import pyplot
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

model = fasttext.load_model("BioWordVec_PubMed_MIMICIII_d200.bin")


# In[2]:


def cosine_similarity(a,b):
    return np.inner(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def train_test_split(dict_hormone):
    training_pairs = []
    testing_pairs = []
    pair_cnt = 0
    train = 0
    test = 0
    for hormone in dict_hormone.keys():
        i = 0
        random.shuffle(dict_hormone[hormone])
        cnt = math.floor(0.7*len(dict_hormone[hormone]))
        for gene in dict_hormone[hormone]:
            if i < cnt:
                pair_cnt += 1
                train += 1
                training_pairs.append((hormone,gene.lower()))
            else:
                pair_cnt += 1
                test += 1
                testing_pairs.append((hormone,gene.lower()))
            i += 1 
    return training_pairs,testing_pairs

def transform_X_values(pair_list):
    X_pairs = []
    for pair in pair_list:
        np1 = model.get_word_vector(pair[0])
        np2 = model.get_word_vector(pair[1])
        X_pairs.append(np.concatenate([np1,np2]))
    return np.array(X_pairs)

def get_best_parameters(X_train,y_train):
    C_range = np.logspace(-2, 10, 13)
    gamma_range = np.logspace(-9, 3, 13)
    param_grid = dict(gamma=gamma_range, C=C_range)
    cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
    grid = GridSearchCV(SVC(), param_grid=param_grid, n_jobs=-1, cv=cv)
    grid.fit(X_train,y_train)
    print("The best parameters are %s with a score of %0.2f" % (grid.best_params_, grid.best_score_))
    C = grid.best_params_['C']
    gamma = grid.best_params_['gamma']
    return C,gamma


# In[3]:


f, axes = plt.subplots(figsize=(7,7))
for file_index in range(5):
    with open('hormone_gene_graph'+str(file_index)+'.txt') as json_file:
        all_hormone_genes = json.load(json_file)
    
    pos_training_pairs, pos_testing_pairs = train_test_split(all_hormone_genes["pos"])
    neg_training_pairs, neg_testing_pairs = train_test_split(all_hormone_genes["neg"])

    all_training_pairs = pos_training_pairs + neg_training_pairs
    all_testing_pairs = pos_testing_pairs + neg_testing_pairs

    X_pos_training_pairs = transform_X_values(pos_training_pairs)
    X_neg_training_pairs = transform_X_values(neg_training_pairs)
    X_pos_testing_pairs = transform_X_values(pos_testing_pairs)
    X_neg_testing_pairs = transform_X_values(neg_testing_pairs)

    pos_y_train = np.ones((X_pos_training_pairs.shape[0],), dtype=int)
    neg_y_train = np.zeros((X_neg_training_pairs.shape[0],), dtype=int)
    y_train = np.concatenate([pos_y_train, neg_y_train])

    pos_y_test = np.ones((X_pos_testing_pairs.shape[0],), dtype=int)
    neg_y_test = np.zeros((X_neg_testing_pairs.shape[0],), dtype=int)
    y_test = np.concatenate([pos_y_test, neg_y_test])

    X_train = np.concatenate([X_pos_training_pairs, X_neg_training_pairs])
    X_test = np.concatenate([X_pos_testing_pairs, X_neg_testing_pairs])
    
    
    C, gamma = get_best_parameters(X_train,y_train)
    rbf_svclassifier = SVC(kernel='rbf',C=C,gamma=gamma,probability=True)
    rbf_svclassifier.fit(X_train, y_train)
    y_pred_test = rbf_svclassifier.predict(X_test)
    y_probs_test = rbf_svclassifier.predict_proba(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, y_probs_test[:,1])
    roc_lab = 'Dataset %d AUC=%.4f' % (file_index, roc_auc_score(y_test, y_probs_test[:,1]))
    axes.step(fpr, tpr, label=roc_lab)
    print("Testing results")
    print(confusion_matrix(y_test, y_pred_test))
    print(classification_report(y_test, y_pred_test))
          
plt.plot([0, 1], [0, 1], 'k--')
axes.set_xlabel('False Positive Rate')
axes.set_ylabel('True Positive Rate')
axes.set_title("Receiver Operating Characteristics Curve")
axes.legend(loc='lower right', fontsize='small')
f.savefig('hormone_gene_pred_unbal_roc.pdf')
pyplot.show()

