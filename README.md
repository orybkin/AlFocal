# Algebraic estimation of Focal length 

These are the tools I developed while working on my Bachelor's thesis on Robust Focal Length Estimation. For more information about the project, see the thesis repo at https://github.com/orybkin/Bachelor-Thesis.

The project was developed in MatLab. It contains implementations of multiple methods for camera focal length estimation. It also contains the framework for evaluation of the methods and various plotting scripts. All figures in the thesis are reproducible via respective scripts.

## References
The project grew out of Pavel Trutman's work on applying deep learning to the focal length estimation here: https://github.com/PavelTrutman/FocalNet

It also makes heavy use of the toolkit for camera geometry privately developed by the group of Tomas Pajdla (under a free license).

## Folder structure ##
 - `deepag/`
   - `bin/`
     - `deepag\`: Main folder with the most important scripts.
       - `results\`: Folder containing results saved by the scripts. The results can be viewed by the `plot*` scripts.
       - `deepag\`: Scripts used to generate the training, validating and testing datasets from the raw data. The results are saved in `data/` folder as `features.mat`, which contains generated feature vectors and as `correspondences.mat`, which contains selected correspondences.
       - `cons_predictor.m`, `closed_form.m`: Scripts to compute error of the constant predictor and of the closed form linear regressor.
       - `k_nn.m`, `k_nn_GPU.m`: Scripts to compute the error of the k nearest neighbours regressor in the feature vector space. The version `_GPU` is designed to use on the GPU. 
       - `k_nn_pixels.m`, `k_nn_pixels_generateData.m`: Scripts to compute the error of the k nearest neighbours regressor in the image coordinates space. `_generateData` version is used to generate synthetic data for testing. The generated data are saved in `data/` folder as `correspondences_syntetic.mat`.
       - `nn*`: Scripts with neural network regressors. The scipts differ in the number of non-linear layers. The `_GPU` versions are used to be executed on the GPU.
       - `plotKNN.m`, `plotNN.m`: Scripts to plot and inspect the results from `k_nn*` and `nn*` scripts respectively. Very handy when the scripts are executed on the cluster on which you can not plot anything.
     - `lib`: Folder with libraries, like `MatconvNet` and utils functions.
   - `data/`: Folder containing all data files.
	 - `RS_paris_i.mat` - original data
	 - `features_c.mat` is corresponding to correspondences.mat
	 - `features.mat` doesn't have correspondences file
	 - `features_F*` - feature vector is fundamental matrix
	   - `*synth*` - no noise, synthentic data
	 - `correspondences_synteticKNN` - correspondences generated for efficient KNN
	 - `correspondences_syntetic` - synthetic correspondences
 - `export/`: Folder with exported images and other documents, usually used in the gDoc with partial results.
 - `logs/`: Folder with log files from batch jobs executed on the cluster
 - `coef_mat_two_focal.m`: Matlab function generating feature vector from seven correspondences as Zuzana sent it.
 - `*.pbs`: Script filed used to submit jobs to the cluster.
 - `Whiteboard-*.jpg`: Captured whiteboards after some discussions.

## Example figure

![alt text](https://github.com/orybkin/AlFocal/blob/master/results/opt_comp_f1.jpg "Figure")
