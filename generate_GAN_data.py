############################################################################################################################################
# 
# python generate_GAN_data.py --generator model_generator_100.h5 --constrainer model_constrainer_e_115.h5 --label labels_5872.pkl
#
############################################################################################################################################

import numpy as np
import tensorflow as tf
import pickle
import argparse

def create_GAN_waveforms(label,g_path,ce_path):
    model_path = "/gpfs/slac/staas/fs1/g/g.exo/shaolei/GAN_models/"
    label_path = "/gpfs/slac/staas/fs1/g/g.exo/shaolei/GAN_event_recon/label/"
    save_path = "/gpfs/slac/staas/fs1/g/g.exo/shaolei/GAN_event_recon/gan_wf/"
    ############# import models ###############################
    generator = tf.keras.models.load_model(model_path + g_path)
    constrainer_e = tf.keras.models.load_model(model_path + ce_path)
    with open(label_path + label,"rb") as file:
        Y = pickle.load(file,encoding='latin1')
    file.close()
    print("Y shape",Y.shape)

    Y = np.array([[i[0], i[1], i[2], i[3] / 10] for i in Y])
    size = Y.shape[0]
    seed = [Y,np.random.rand(size,100)]
    X_gen = generator.predict(seed)
    print("nEvent: ",X_gen.shape[0])
    '''
    with open(save_path + "GAN_generated_waveform_5872.pkl","wb") as f:
        pickle.dump(X_gen,f,protocol=2)
    '''
    X_gen = X_gen.reshape(X_gen.shape[0],74*350)
    print("Done creating wwavforms.")
    f = open(save_path + 'GAN_generated_waveform_5872(1).csv','w')
    for i in X_gen:
        row = ''
        for k in range(i.shape[0]):
            row += '%f'%i[k]
            if not k == i.shape[0]-1:
                row += ','
        row += '\n'
        f.write(row)
    f.close()
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=str)
    parser.add_argument("--constrainer",type=str)
    parser.add_argument("--label",type=str,default="labels_5872.pkl")
    args = parser.parse_args()
    generator = args.generator
    constrainer = args.constrainer
    label = args.label

    create_GAN_waveforms(label,generator,constrainer)
