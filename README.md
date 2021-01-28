# large_random_fields
Generate large Matern random fields, identify unknown parameters (variance, cov. length, smoothness, nugget) = (\sigma, \ell, \nu, \tau^2)\

Install HLIBPro v2.9 from www.hlibpro.com  
Example 1:  
  
generate.cc is needed to generate "the truth"  
loglikelihood4.cc is needed to identify 4 unknown parameters of the Matern   

Copy  

generate.cc and loglikelihood4.cc to /hlibpro-2.9/examples/  

cd examples  

 vi SConscript  
 
 and modify it by adding 
 
 examples.append( cxxenv.Program( 'generate.cc' ) ) 
 
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

