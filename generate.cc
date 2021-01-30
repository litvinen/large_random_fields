
// Run like this
//./generate --data competition/S1a/S1b_short_dataset1_training.csv 
// This code will forecast 4 unknown parameters of the Matern covariance matrix
// Developed by Alexander Litvinenko (RWTH Aachen) and Ronald Kriemann (MIS MPG Leipzig)
// Based on the HLIBPro library (v. 2.9) www.hlibpro.com
// No warranties.

#include <iostream>
#include <fstream>
#include <string>

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

#include <gsl/gsl_multimin.h>

#include "hlib.hh"

using namespace std;
using boost::format;
using namespace HLIB;
using namespace boost::program_options;
using  real_t    = HLIB::real;

enum {
    IDX_SIGMA  = 0,
    IDX_LENGTH = 1,
    IDX_NU     = 2,
    IDX_TAU    = 3
};


//Use a method described by Abramowitz and Stegun: 
double gaussrand_Stegun()
{
    static double U, V;
    static int phase = 0;
    double Z;

    if(phase == 0) {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
    } else
        Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);

    phase = 1 - phase;

    return Z;
}

//Use a method discussed in Knuth and due originally to Marsaglia:

double gaussrand_Knuth()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

// global options
int        nmin      = CFG::Cluster::nmin;
double     eps       = 1e-6;
double     fac_eps   = 1e-6;
double     shift     = 1e-7;
bool       use_ldl   = false;

//
// read dataset from file
//
void
read_data ( const std::string &       datafile,
            std::vector< T2Point > &  vertices,
            BLAS::Vector< double > &  Z_data )
{
    std::ifstream  in( datafile );
    
    if ( ! in ) // error
        exit( 1 );

    size_t  N_vtx = 0;
    
    #if 1

    std::string  line;
    
    std::getline( in, line );

    if ( line == "x,y" )
    {
        std::list< T2Point >  pos;
        std::list< double >   vals;

        while ( std::getline( in, line ) )
        {
            auto    parts = split( line, "," );
            double  x = atof( parts[0].c_str() );
            double  y = atof( parts[1].c_str() );
           
            
            pos.push_back( T2Point( x, y ) );
        }// while

        N_vtx = pos.size();

        std::cout << "learning dataset" << std::endl;
        std::cout << N_vtx << std::endl;
        
        vertices.resize( N_vtx );
        Z_data = BLAS::Vector< double >( N_vtx );

        int  i = 0;

        for ( auto  p : pos )
            vertices[ i++ ] = p;

        i = 0;
        
        for ( idx_t  i = 0; i < idx_t(N_vtx); ++i )
           Z_data( i++ ) =  gaussrand_Knuth(); //gaussrand_Stegun() ;
    }// if
    else
    {
        std::cout << "you should not be here, something is wrong with the input file" << std::endl;
        HERROR( ERR_NOT_IMPL, "", "" );
    }
    
    #else
    
    in >> N_vtx;

    std::cout << "reading " << N_vtx << " datapoints" << std::endl;
    
    vertices.resize( N_vtx );
    Z_data = BLAS::Vector< double >( N_vtx );
        
    for ( idx_t  i = 0; i < idx_t(N_vtx); ++i )
    {
        int     index = i;
        double  x, y, z;
        // double  v     = 0.0;
        
        in >> index >> x >> y >> z;
        // in >> index >> x >> y >> z >> v;

        vertices[ index ] = T2Point( x, y );
        //Z_data( index )   = v;
    }// for

    #endif
    
    //
    // for visualization of data, export 2D points with v value in csv file
    //

    // std::ofstream  out( "data.csv" );

    // out << "x,y,z,v" << std::endl;
    // out << "x,y,z" << std::endl;
    
    //for ( uint  i = 0; i < N_vtx; ++i )
    //   out << vertices[i].x() << "," << vertices[i].y() << ",0" << std::endl;
//      out << vertices[i].x() << "," << vertices[i].y() << ",0," << Z_data( i ) << std::endl;
}

//
// define PredictionProblem to forecast unknown values in new locations
//
struct GeneratingProblem
{
    std::vector< T2Point >                vertices;
    std::unique_ptr< TCoordinate >        coord;
    std::unique_ptr< TClusterTree >       ct;
    std::unique_ptr< TBlockClusterTree >  bct;
    std::unique_ptr< TVector >            Z;
    
    GeneratingProblem ( const std::string &  datafile )
    {
        init( datafile );
    }

    void
    init ( const std::string &  datafile )
    {
        BLAS::Vector< double >  Z_data;

        read_data( datafile, vertices, Z_data );
        std::cout << "the grid is successfully read" << std::endl;
        coord = std::make_unique< TCoordinate >( vertices );
        
        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder( & part_strat, nmin );
    
        ct = ct_builder.build( coord.get() );
        //print_vtk( & coord, "ct_coord" );
        //print_vtk( & coord_predict, "ct_coord_predict" );
        
    
        TStdGeomAdmCond    adm_cond( 2.0, use_min_diam );
        TBCBuilder         bct_builder;
   
        bct = bct_builder.build( ct.get(), ct.get(), & adm_cond );
        Z   = std::make_unique< TScalarVector >( *ct->root(), std::move( Z_data ) );
        ct->perm_e2i()->permute( Z.get() );
  }
    
    //BLAS::Vector< double > 
    std::unique_ptr< TVector > 
    eval ( const double  sigma,
           const double  length,
           const double  nu,
           const double tau )
    {
  
        TMaternCovCoeffFn< T2Point >  matern_coefffn( sigma, length, nu,  vertices );
        TPermCoeffFn< double >        coefffn( & matern_coefffn, ct->perm_i2e(), ct->perm_i2e() );
        
        TACAPlus< double >            aca( & coefffn );
        auto                          acc = fixed_prec( eps );
        TDenseMatBuilder< double >    h_builder( & coefffn, & aca );
        
        auto                          C        = h_builder.build( bct.get(), acc );
        TPSMatrixVis  mvis;
        
         //mvis.svd(true).print( C.get(), "myC" );
 
//        print_ps(bct->root(), "bct.eps");
//        print_ps(bct_predict->root(), "bct_predict.eps");
  
        //       mvis.svd(true).print( C_predict.get(), "myC_predict" );
  
        //if ( shift != 0.0 )
        //    add_identity( C.get(), tau*tau );
  
        auto                          fac_acc  = fixed_prec( fac_eps );
        auto                          C_fac    = C->copy();
        auto                          fac_opts = fac_options_t{ point_wise, CFG::Arith::storage_type, false };
    
        
        chol( C_fac.get(), fac_acc );
    

        //std::cout << "    |С|_F             = " << norm_F( C.get() ) << std::endl;
        //std::cout << "    |L|_F             = " << norm_F( C_fac.get() ) << std::endl;
        //std::cout << "    |С|_2             = " << norm_2( C.get() ) << std::endl;
        //std::cout << "    |L|_2             = " << norm_2( C_fac.get() ) << std::endl;
        //mvis.svd(true).print( C_fac.get(), "myL" );

        //std::unique_ptr< TFacInvMatrix >   C_inv;

     
      
        const size_t                  N     = vertices.size();
        
        
        
        auto                          Z_generated = C->row_vector();
    
      //  std::cout << "  size of C = " << C->rows() << "x"<< C->cols() << std::endl;
      //  std::cout << "  size of C_predict = " << C_predict->rows() <<"x"<<  C_predict->cols() << std::endl;
      //  std::cout << "  sol size = " << sol->size() << std::endl;
      //  std::cout << "  ||sol|| = " << sol->norm2() << std::endl;
      //  std::cout << "  ||Z_predict|| = " << Z_predict->norm2() << std::endl;
  
        //auto                          ZdotCZ = std::real( Z->dot( sol.get() ) );
        //Z_predict = C_predict * sol.get(); 
        mul_vec( real_t(1), C_fac.get(), Z.get(), real_t(0), Z_generated.get(), apply_normal );
        //C_predict->mul_vec( 1.0, sol.get(), 0.0, Z_predict.get(), apply_normal );
        
        ct->perm_i2e()->permute( Z_generated.get() );
        std::cout << "  ||Z_generated|| = " << Z_generated->norm2() << std::endl;
        
        //TMatlabVectorIO  vio;
 
        //vio.write( Z_predict,  "x.mat", "x" );

        FILE* f1;
        f1 = fopen("111gen_d.txt", "w");

        for ( size_t  i = 0; i < Z_generated->size(); i++ )
          fprintf(f1," %6.6e, %6.6e, %6.6e\n",   vertices[i].x(),  vertices[i].y(), Z_generated->entry(i));
        fclose(f1);
        return std::move( Z_generated );
    }
};

//
// wrapper from GSL to LogLikeliHoodProblem
//
/*double
  eval_logli ( const gsl_vector *  param,
  void *              data )
  {
  double sigma  = gsl_vector_get( param, IDX_SIGMA );
  double length = gsl_vector_get( param, IDX_LENGTH );
  double nu     = gsl_vector_get( param, IDX_NU );
  double tau    = gsl_vector_get( param, IDX_TAU );

  LogLikeliHoodProblem *  problem = static_cast< LogLikeliHoodProblem * >( data );

  return - problem->eval( sigma, length, nu, tau );
  }
*/
//
// optimization function using GSL
//

//
// main function
//
int
main ( int      argc,
       char **  argv )
{

    
    CFG::set_verbosity( 3 );
    INIT();
    
    //std::string  datafile = "datafile.txt";
    //std::string  datafile_predict = "datafile_predict.txt";
    std::string  datafile = "grid.txt";
    
    //
    // define command line options
    //

    options_description             all_opts;
    options_description             vis_opts( "usage: generatig [options] datafile\n  where options include" );
    options_description             hid_opts( "Hidden options" );
    positional_options_description  pos_opts;
    variables_map                   vm;

    // standard options
    vis_opts.add_options()
        ( "help,h",                       ": print this help text" )
        ( "threads,t",   value<int>(),    ": number of parallel threads" )
        ( "verbosity,v", value<int>(),    ": verbosity level" )
        ( "nmin",        value<int>(),    ": set minimal cluster size" )
        ( "eps,e",       value<double>(), ": set H accuracy" )
        ( "epslu",       value<double>(), ": set only H factorization accuracy" )
        ( "shift",       value<double>(), ": regularization parameter" )
        ( "ldl",                          ": use LDL factorization" )
        ;
    
    hid_opts.add_options()
        ( "data",        value<std::string>(), ": datafile " );

    // options for command line parsing
    all_opts.add( vis_opts ).add( hid_opts );

    // all "non-option" arguments should be "--data" arguments
    pos_opts.add( "data", -1 );

    //
    // parse command line options
    //

    try
    {
        store( command_line_parser( argc, argv ).options( all_opts ).positional( pos_opts ).run(), vm );
        notify( vm );
    }// try
    catch ( required_option &  e )
    {
        std::cout << e.get_option_name() << " requires an argument, try \"-h\"" << std::endl;
        exit( 1 );
    }// catch
    catch ( unknown_option &  e )
    {
        std::cout << e.what() << ", try \"-h\"" << std::endl;
        exit( 1 );
    }// catch

    //
    // eval command line options
    //

    if ( vm.count( "help") )
    {
        std::cout << vis_opts << std::endl;
        exit( 1 );
    }// if

    if ( vm.count( "nmin"      ) ) nmin     = vm["nmin"].as<int>();
    if ( vm.count( "eps"       ) ) eps      = vm["eps"].as<double>();
    if ( vm.count( "epslu"     ) ) fac_eps  = vm["epslu"].as<double>();
    if ( vm.count( "shift"     ) ) shift    = vm["shift"].as<double>();
    if ( vm.count( "threads"   ) ) CFG::set_nthreads( vm["threads"].as<int>() );
    if ( vm.count( "verbosity" ) ) CFG::set_verbosity( vm["verbosity"].as<int>() );
    if ( vm.count( "ldl"       ) ) use_ldl  = true;

    // default to general eps
    if ( fac_eps == -1 )
        fac_eps = eps;
    
    if ( vm.count( "data" ) )
        datafile = vm["data"].as<std::string>();
    else
    {
        std::cout << "usage: generating [options] datafile" << std::endl;
        exit( 1 );
    }// if


    double  sigma  = 2.0; //take these values from previous experiments (Part 1a)
    double  length = 0.1; 
    double  nu     = 0.5;
    double  tau    = 0.0;
    
    GeneratingProblem  problem( datafile );
    problem.eval( sigma, length, nu, tau);

    DONE();
}
