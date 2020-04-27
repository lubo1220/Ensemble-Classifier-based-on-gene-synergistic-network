Writing in the front
GSE2034_ma2.mat and test_example.mat data sets necessary for the program are too large to be uploaded here. If you need them, please contact the author (516267582@qq.com).


//Execute the program in sequence
background_vector_calculation.m //The uncertainty coefficient of the gene was calculated to obtain the background distribution
synergistic_network.m //According to the three synergistic relationships, three gene synergistic networks are built
centroid_classifier.m //Center of mass classifier model function
classifier_construct.m //According to the communities, weak classifiers are established, and the performance of the classifiers are evaluated using cross validation, and then the well-performing classifier is integrated as the ensemble classifier
standard.m //Standardize all data sets
independ_verified.m //Independent validation, run standards.m before running, one independent test data set at a time
independ_verified_optimization.m //New method of independent validation, run standards.m before running, one independent test data set at a time
using_geneId_obtained.m //Obtain the community gene Id that performs well

//Description of data
test_example.mat //five data sets
Background_vector.mat //The background distribution of the uncertainty coefficient
Relation.mat //There are three types of significant pairs of synergies
Subnetwork.mat //dense subnetworks of three synergistic network
R.mat //the performance results of all weak classifiers
Using_geneId.mat //The gene Id of a useful subnetwork
GeneId_Symbol.mat //The id of the gene corresponds to the symbol

mRMR_0.9_compiled2 //A toolkit for calculating entropy