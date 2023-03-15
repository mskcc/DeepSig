from sklearn.model_selection import train_test_split
import tensorflow as tf
import pandas as pd
import numpy as np
import os
import pdb
import random

S = 'SBS.A'
mid = '1a'

mepochs = pd.read_csv('../../mrange2.txt',header=0, sep='\t')
epochs = mepochs['epochs'][mid]
tf.random.set_seed(135)
e0 = pd.read_csv('../data/binary.csv', header=0, sep=',')
e0.index = e0['Unnamed: 0'].values.tolist()
e = e0[S]

x0 = pd.read_csv('../data/catalog.csv',header=0, sep=',')
x0.index = x0['Unnamed: 0'].values.tolist()
x = x0.drop(columns='Unnamed: 0')

ex = pd.concat([e,x], axis=1)

train = ex

def format_output(data):
  y = data.pop(S)
  y = np.array(y)
  return y

train_Y = format_output(train)
train_Y = tf.convert_to_tensor(train_Y, dtype=tf.int64)

train_X = tf.convert_to_tensor(np.array(train),dtype=tf.float64)

def build_model():
  input_layer = tf.keras.Input(shape=(len(train.columns), ))
  L1 = tf.keras.layers.Dense(units='64',activation='sigmoid')(input_layer)
  L2 = tf.keras.layers.Dense(units='64',activation='sigmoid')(L1)
  y_output = tf.keras.layers.Dense(units='1', name=S)(L2)
  model = tf.keras.Model(inputs=input_layer, outputs=y_output)
  return model


model = build_model()
model.compile(optimizer='adam', loss=tf.keras.losses.BinaryCrossentropy(from_logits=True), 
              metrics=[tf.keras.metrics.BinaryAccuracy()])

history = model.fit(train_X, train_Y, batch_size=8, epochs=epochs)
z = history.history

model.save(S+'_'+mid)
