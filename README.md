# large_random_fields

developed by Alexander Litvinenko, RWTH Aachen, Germany, 
and Ronald Kriemann, MIS MPI Leipzig, Germany (www.hlibpro.com),  
January 2021

The theory is described here:

1) HLIBCov: Parallel hierarchical matrix approximation of large covariance matrices and likelihoods with applications in parameter identification
AlexanderLitvinenko, Ronald Kriemann, Marc G.Genton, Ying Sun, David E.Keyes,
MethodsX
Volume 7, 2020, 100600
https://www.sciencedirect.com/science/article/pii/S2215016119301761

2) Likelihood approximation with hierarchical matrices for large spatial datasets,
AlexanderLitvinenko, Ying Sun, Marc G. Genton, David E.Keyes, Computational Statistics & Data Analysis
Volume 137, September 2019, Pages 115-132,
https://www.sciencedirect.com/science/article/pii/S0167947319300374 

More details about "KAUST Competition on Spatial Statistics for Large Datasets" is here https://cemse.kaust.edu.sa/stsds/kaust-competition-spatial-statistics-large-datasets

With this package you can generate large Matern random fields, identify unknown parameters 
(variance, cov. length, smoothness, nugget) = (\sigma, \ell, \nu, \tau^2) simultaneously.

Also you can find procedures: HierarchicalMatrix-vector multiplication, Hierarchical Chiolesky, Computation of the joint Gaussian Loglikelihoood, matrix determinant, hierarchical matrix inverse, soution of a linear system, scalar product.

To install and run please follow these steps:

Install HLIBPro v2.9 from www.hlibpro.com  

Example 1:  
  
generate.cc is needed to generate "the truth"  
loglikelihood4.cc is needed to identify 4 unknown parameters of the Matern covariaance matrix

To try follow please the follwoing steps

Download and install HLIBPro from www.hlibpro.com. It shouldgo very fast (1 minute).


Copy generate.cc, prediction.cc and loglikelihood4.cc (and all other dataset files) to /hlibpro-2.9/examples/  

cd examples  

 vi SConscript  
 
 and modify it by adding 
 
 examples.append( cxxenv.Program( 'generate.cc' ) ) 
 
 examples.append( cxxenv.Program( 'prediction.cc' ) ) 
 
 examples.append( cxxenv.Program( 'loglikelihood4.cc' ) )  
 
 cd .. 
 scons 
 cd examples 
 and run 

./generate --data S1b_short_dataset1_training.csv 

As a result a new file 111gen_d.txt is generated.

Now run  

 ./loglikelihood4 111gen_d.txt  
 
 It will be generated this output:  
 
 ....  
 
  loglikelihood at  σ = 1.52958, ℓ = 0.102099, ν = 0.568078, tau = 6.06898e-06   is -607.012.  
  
  loglikelihood at  σ = 1.52958, ℓ = 0.102099, ν = 0.568078, tau = 6.06898e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52863, ℓ = 0.101924, ν = 0.568156, tau = 6.13021e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52863, ℓ = 0.101924, ν = 0.568156, tau = 6.13021e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52863, ℓ = 0.101924, ν = 0.568156, tau = 6.13021e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52942, ℓ = 0.102025, ν = 0.56812, tau = 6.07958e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52942, ℓ = 0.102025, ν = 0.56812, tau = 6.07958e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52942, ℓ = 0.102025, ν = 0.56812, tau = 6.07958e-06   is -607.012. 
  
  loglikelihood at  σ = 1.52942, ℓ = 0.102025, ν = 0.56812, tau = 6.07958e-06   is -607.012. 
  
converged to minimum. 

  loglikelihood at  σ = 1.52942, ℓ = 0.102025, ν = 0.56812, tau = 6.07958e-06   is -607.012. 
  
  


As a result four unknown parameters will be identified. They should be very close to parameters used in generate.cc, e.g., to:  

  double  sigma  = 2.0;  
  
  double  length = 0.1;  
  
  double  nu     = 0.5;  
  
  double  tau    = 0.0;  
  
  
  Example 2:
  To predict unknown values of a random field in some new locations run: 
  
   $ ./prediction --dataL dataset2_training.csv --dataT dataset2_testing.csv 
   
   or
   
   $ ./prediction --dataL dataset2_training.csv --dataT dataset2_testing.csv --ldl
   
   Here dataset2_training.csv contains the learning dataset (3 columsn, x,y,values) and 
   dataset2_testing.csv  (two columns, x,y ) contains the set of locations where the random field should be forecasted.
   
   $  ./prediction --dataL competition/S1a/dataset8_training.csv --dataT competition/S1b/dataset8_testing.csv --sigma=1.09429 --length=0.204309 --nu=0.60055 --tau=6.373e-9 --ldl


