from __future__ import absolute_import
from __future__ import print_function
from one_hot_encode import do_all_the_shit, rev_one_hot, to_aa, index_to_aa, set_k, tf_maxpool, tf_maxpool_with_argmax, tf_preprocess_spectrum, ion_current_normalize, parse, fire_up_generator, get_dataset
import tensorflow as tf
import numpy as np
import random
import os
from tensorflow.keras.callbacks import CSVLogger

os.environ["CUDA_VISIBLE_DEVICES"]="-1"
tf.compat.v1.disable_eager_execution()

def euclidean_distance(vects):
    x, y = vects
    sum_square = tf.math.reduce_sum(tf.math.square(x - y), axis=1, keepdims=True)
    return tf.math.sqrt(tf.math.maximum(sum_square, tf.keras.backend.epsilon()))  # fuzz factor


def eucl_dist_output_shape(shapes):
    shape1, shape2 = shapes
    return (shape1[0], 1)


def contrastive_loss(y_true, y_pred):
    margin = 1
    square_pred = tf.math.square(y_pred)
    margin_square = tf.math.square(tf.math.maximum((margin - y_pred), 0))
    return tf.math.reduce_mean(y_true * square_pred + (1 - y_true) * margin_square)

def compute_accuracy(y_true, y_pred):
    '''Compute classification accuracy with a fixed threshold on distances.
    '''
    pred = y_pred.ravel() < 0.5  #returns a single array
    return np.mean(pred == y_true)


def accuracy(y_true, y_pred):
    '''Compute classification accuracy with a fixed threshold on distances.
    '''
    return tf.keras.backend.mean(tf.keras.backend.equal(y_true, tf.keras.backend.cast(y_pred < 0.5, y_true.dtype))) #cast into y_true type


class MyOwnFunction(tf.keras.callbacks.Callback):
  def on_epoch_end(self, epoch, logs=None):
    seq_model.save('/hpi/fs00/home/jiahao.hu/data/save_model/sequence_model')
    siam_model.save('/hpi/fs00/home/jiahao.hu/data/save_model/siamese_model')
    spec_model.save('/hpi/fs00/home/jiahao.hu/data/save_model/spectrum_model')


seq_shape=(1,24,22)
spec_shape=(3600,2)

seq_inputs = tf.keras.Input(shape=seq_shape)
x = tf.keras.layers.Reshape(target_shape=(24,22))(seq_inputs)
x = tf.keras.layers.Conv1D(filters=8,kernel_size=3,activation=tf.nn.relu, strides=1, padding='same', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=16, kernel_size=3, activation=tf.nn.relu, strides=1, padding='same', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='same')(x)
x = tf.keras.layers.Conv1D(filters=16, kernel_size=3, activation=tf.nn.relu, strides=1, padding='same', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=32, kernel_size=3, activation=tf.nn.relu, strides=1, padding='same', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='same')(x)
x = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation=tf.nn.relu, strides=1, padding='same', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=128, kernel_size=3, activation=tf.nn.relu, strides=1, padding='same', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='same')(x)
x = tf.keras.layers.Flatten()(x)
x = tf.keras.layers.Dense(128, activation=tf.nn.relu)(x)
seq_outputs = tf.keras.layers.Dense(128, activation=tf.nn.relu)(x)
seq_model = tf.keras.Model(inputs=seq_inputs, outputs=seq_outputs)


spec_inputs = tf.keras.Input(shape=spec_shape)
x = tf.keras.layers.Conv1D(filters=8, kernel_size=3, activation=tf.nn.relu, strides=2, padding='valid', kernel_initializer='glorot_uniform')(spec_inputs)
x = tf.keras.layers.Conv1D(filters=8, kernel_size=3, activation=tf.nn.relu, strides=2, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=16, kernel_size=3, activation=tf.nn.relu, strides=2, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')(x)
x = tf.keras.layers.Conv1D(filters=16, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=32, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=32, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')(x)
x = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=128, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')(x)
x = tf.keras.layers.Conv1D(filters=128, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=128, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=256, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')(x)
x = tf.keras.layers.Conv1D(filters=256, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Conv1D(filters=512, kernel_size=3, activation=tf.nn.relu, strides=1, padding='valid', kernel_initializer='glorot_uniform')(x)
x = tf.keras.layers.Flatten()(x)
x = tf.keras.layers.Dense(128, activation=tf.nn.relu)(x)
x = tf.keras.layers.Dropout(0.1)(x)
spec_outputs = tf.keras.layers.Dense(128, activation=tf.nn.relu)(x)
spec_model = tf.keras.Model(inputs=spec_inputs, outputs=spec_outputs)

inputs = tf.keras.Input(shape=(128,))
x = tf.keras.layers.Flatten()(inputs)
x = tf.keras.layers.Dense(128, activation=tf.nn.relu)(x)
x = tf.keras.layers.Dense(128, activation=tf.nn.relu)(x)
x = tf.keras.layers.Dense(64, activation=tf.nn.relu)(x)
x = tf.keras.layers.Dense(64, activation=tf.nn.relu)(x)
x = tf.keras.layers.Dense(32, activation=tf.nn.relu)(x)
outputs = tf.keras.layers.Dense(32)(x)
siam_model = tf.keras.Model(inputs=inputs, outputs=outputs)


#embeddings
embedding_spec = siam_model(spec_outputs)
embedding_seq = siam_model(seq_outputs)
#Compute the euclidean distance
distance = tf.keras.layers.Lambda(euclidean_distance,
                  output_shape=eucl_dist_output_shape)([embedding_spec, embedding_seq])

# define train sets
train_set = get_dataset(fire_up_generator('/hpi/fs00/home/jiahao.hu/data/files/training.hdf5',preshuffle=True),training=True)
#train_set = train_set.prefetch(tf.data.experimental.AUTOTUNE)
test_set = get_dataset(fire_up_generator('/hpi/fs00/home/jiahao.hu/data/files/test.hdf5'),training=False)


# save hstory

my_callbacks = [tf.keras.callbacks.ModelCheckpoint('/hpi/fs00/home/jiahao.hu/data/training/checkpoint_path.{epoch:02d}-{val_loss:.2f}.hdf5', monitor='val_loss', verbose=0, save_best_only=False, mode='auto'),
                tf.keras.callbacks.CSVLogger("/hpi/fs00/home/jiahao.hu/data/model_history_log.csv", append=True, separator=";"),
                MyOwnFunction()]


# define model (inputs,outputs)
model = tf.keras.Model([spec_inputs, seq_inputs], distance)

model.compile(loss=contrastive_loss, optimizer=tf.keras.optimizers.Adam(0.000001), metrics=[accuracy])
model.fit(train_set,
          epochs=450,
          steps_per_epoch = 2000,
          validation_data = test_set,
          validation_steps = 1000,
          initial_epoch = 0,
          callbacks= my_callbacks) 



