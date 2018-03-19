##########################################################
This is the readme file for Single Cell Migration Paper
Mar.15,2018
##########################################################

For mitochondria three-class classification (filamentous,
intermediate, dots), you may run the program randomfor.m, which 
uses the random forest method to do this
classification. The program named as nnetw.m is also doing
this three-class classification but uses a three-layer 
neural net. The program named as svmc.m is doing the three-
class classification using supporting vector machine.

For cell migration direction prediction, there are two 
relative programs. The first one is RFforCell.m, which
uses a random forest method for cell migration direction
prediction. The second one is npradaption.m, which uses
a three-layer neural net to predict migration direction
for each individual cell.

For cell migration speed (motile/non-motile) prediction,
you may run NNTrain2.m for 10/90 and 20/80 partition with 
other supporting data (MAT files) included in this folder.
It will do feature reduction and provide XLSX files like 
"Feature2080_3" at the end. This program is time-consuming
(needs ~10 hrs running on HPZ240 using Intel Xeon processor).
Bascially, NNTrain1.m does the same thing but with train/test
sequence in a different order. All values recorded in XLSX 
files are the misprediction rate. For the prediction accuracy,
we need to substract 1 with these values.





