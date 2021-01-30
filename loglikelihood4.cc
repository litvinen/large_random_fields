// Run like this
// ./loglikelihood4 competition/S1a/dataset13_training.csv  --use_ldl
// This code will forecast 4 unknown parameters of the Matern covariance matrix
// Developed by Alexander Litvinenko (RWTH Aachen) and Ronald Kriemann (MIS MPG Leipzig)
// Based on the HLIBPro library (v. 2.9) www.hlibpro.com
// No warranties.
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>


#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <gsl/gsl_multimin.h>

#include "hlib.hh"

using namespace HLIB;
using namespace boost::program_options;

enum {
    IDX_SIGMA  = 0,
    IDX_LENGTH = 1,
    IDX_NU     = 2,
    IDX_TAU     = 3
};

// global options
int        nmin      = CFG::Cluster::nmin;
double     eps       = 1e-6;
double     fac_eps   = 1e-6;
double     shift     = 1e-8;
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

    if ( line == "x,y,values" )
    {
        std::list< T2Point >  pos;
        std::list< double >   vals;

        while ( std::getline( in, line ) )
        {
            auto    parts = split( line, "," );
            double  x = atof( parts[0].c_str() );
            double  y = atof( parts[1].c_str() );
            double  v = atof( parts[2].c_str() );
            
            pos.push_back( T2Point( x, y ) );
            vals.push_back( v );
        }// while

        N_vtx = pos.size();

        std::cout << N_vtx << std::endl;
        
        vertices.resize( N_vtx );
        Z_data = BLAS::Vector< double >( N_vtx );

        int  i = 0;

        for ( auto  p : pos )
            vertices[ i++ ] = p;

        i = 0;
        
        for ( auto  v : vals )
            Z_data( i++ ) = v;
    }// if
    else
        HERROR( ERR_NOT_IMPL, "", "" );
    
    #else
    
    in >> N_vtx;

    std::cout << "reading " << N_vtx << " datapoints" << std::endl;
    
    vertices.resize( N_vtx );
    Z_data = BLAS::Vector< double >( N_vtx );
        
    for ( idx_t  i = 0; i < idx_t(N_vtx); ++i )
    {
        int     index = i;
        double  x, y, z;
        double  v     = 0.0;

        in >> index >> x >> y >> z >> v;

        vertices[ index ] = T2Point( x, y );
        Z_data( index )   = v;
    }// for

    #endif
    
    //
    // for visualization of data, export 2D points with v value in csv file
    //

    std::ofstream  out( "data.csv" );

    out << "x,y,z,v" << std::endl;
    
    for ( uint  i = 0; i < N_vtx; ++i )
        out << vertices[i].x() << "," << vertices[i].y() << ",0," << Z_data( i ) << std::endl;
}

//
// define LogLikeliHood Problem to be evaluated at theta by minimization function
//
struct LogLikeliHoodProblem
{
    std::vector< T2Point >                vertices;
    std::unique_ptr< TCoordinate >        coord;
    std::unique_ptr< TClusterTree >       ct;
    std::unique_ptr< TBlockClusterTree >  bct;
    std::unique_ptr< TVector >            Z;

    LogLikeliHoodProblem ( const std::string &  datafile )
    {
        init( datafile );
    }

    void
    init ( const std::string &  datafile )
    {
        BLAS::Vector< double >  Z_data;

        read_data( datafile, vertices, Z_data );

        coord = std::make_unique< TCoordinate >( vertices );

        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder( & part_strat, nmin );
    
        ct = ct_builder.build( coord.get() );
    
        TStdGeomAdmCond    adm_cond( 2.0, use_min_diam );
        TBCBuilder         bct_builder;
    
        bct = bct_builder.build( ct.get(), ct.get(), & adm_cond );

        Z   = std::make_unique< TScalarVector >( *ct->root(), std::move( Z_data ) );

        ct->perm_e2i()->permute( Z.get() );
    }
    
  

    double
    eval ( const double  sigma,
           const double  length,
           const double  nu,
           const double tau )
    {
        TMaternCovCoeffFn< T2Point >  matern_coefffn( 2.0/pow(1.1,sigma), 1.0/pow(1.5,length), 1.0/pow(1.2, nu),  vertices );
       // TMaternCovCoeffFn< T2Point >  matern_coefffn( sigma, length, nu,  vertices );
        TPermCoeffFn< double >        coefffn( & matern_coefffn, ct->perm_i2e(), ct->perm_i2e() );

        TACAPlus< double >            aca( & coefffn );
        auto                          acc = fixed_prec( eps );
        TDenseMatBuilder< double >    h_builder( & coefffn, & aca );
    
        auto                          C        = h_builder.build( bct.get(), acc );

        if ( shift != 0.0 )
        {
           // std::cout << "Diagonal is added!!!!!!!!!!" << std::endl;
            add_identity( C.get(), 1/pow(2.0,tau) );
        }

        auto                          fac_acc  = fixed_prec( fac_eps );
        auto                          C_fac    = C->copy();
        auto                          fac_opts = fac_options_t{ point_wise, CFG::Arith::storage_type, false };
    
        if ( use_ldl )
            ldl( C_fac.get(), fac_acc, fac_opts );
        else
            chol( C_fac.get(), fac_acc );
    
        std::unique_ptr< TFacInvMatrix >   C_inv;

        if ( use_ldl )
            C_inv = std::make_unique< TLDLInvMatrix >( C_fac.get(), symmetric, point_wise );
        else
            C_inv = std::make_unique< TLLInvMatrix >( C_fac.get(), symmetric );
        
        const size_t                  N     = vertices.size();
        double                        log_det_C = 0.0;
    
        for ( idx_t  i = 0; i < idx_t(N); ++i )
        {
            if ( use_ldl ) log_det_C +=   std::log( C_fac->entry( i, i ) );
            else           log_det_C += 2*std::log( C_fac->entry( i, i ) ); // two factors L!
        }// if
        
        TStopCriterion                sstop( 250, 1e-16, 0.0 );
        TCG                           solver( sstop );
        auto                          sol = C->row_vector();
    
        solver.solve( C.get(), sol.get(), Z.get(), C_inv.get() );
        
        auto                          ZdotCZ = std::real( Z->dot( sol.get() ) );
        const double                  log2pi = std::log( 2.0 * Math::pi<double>() );
        auto                          LL     = -0.5 * ( N * log2pi + log_det_C + ZdotCZ );

        return LL;
    }
};

//
// wrapper from GSL to LogLikeliHoodProblem
//
double
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

//
// optimization function using GSL
//
double
maximize_likelihood ( double &                sigma,
                      double &                length,
                      double &                nu,
                      double &                tau,
                      LogLikeliHoodProblem &  problem )
{
    int        status   = 0;
    const int  max_iter = 500;

    const gsl_multimin_fminimizer_type *  T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *             s = NULL;
    gsl_vector *                          ss;
    gsl_vector *                          x;
    gsl_multimin_function                 minex_func;
    double                                size;

    x = gsl_vector_alloc( 4 );        // start value
    gsl_vector_set( x, IDX_SIGMA,  sigma );
    gsl_vector_set( x, IDX_LENGTH, length );
    gsl_vector_set( x, IDX_NU,     nu );
    gsl_vector_set( x, IDX_TAU,    tau );

    ss = gsl_vector_alloc( 4 );       // step sizes
    gsl_vector_set( ss, IDX_SIGMA,  1.0 );
    gsl_vector_set( ss, IDX_LENGTH, 1.0 );
    gsl_vector_set( ss, IDX_NU,     1.0 );
    gsl_vector_set( ss, IDX_TAU,    1.0 );
    //gsl_vector_set( ss, IDX_SIGMA,  0.05 );
    //gsl_vector_set( ss, IDX_LENGTH, 0.005 );
    //gsl_vector_set( ss, IDX_NU,     0.01 );
    //gsl_vector_set( ss, IDX_TAU,    1 );

    // Initialize method and iterate
    minex_func.n      = 4;
    minex_func.f      = & eval_logli;
    minex_func.params = & problem;

    s = gsl_multimin_fminimizer_alloc( T, 4 );
    gsl_multimin_fminimizer_set(s, & minex_func, x, ss );

    int     iter = 0;
    double  LL = 0;

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate( s );

        if ( status != 0 )
            break;

        size   = gsl_multimin_fminimizer_size( s );    // return eps for stopping criteria
        status = gsl_multimin_test_size( size, 1e-3 ); // This function tests the minimizer specific characteristic size 

        if ( status == GSL_SUCCESS )
        {
            std::cout << "converged to minimum" << std::endl;
            FILE* f4 = fopen( "111_S2b_results.txt", "a+");
           // fprintf(f4, "sigma, ell, nu, tau\n"); //2.0/pow(1.1,sigma), 1.0/pow(1.5,length), 1.0/pow(1.2, nu),
            fprintf(f4, "& %.6f & %.6f & %.6f & %.6f & %.6f \\\\ \\hline \n", 2.0/pow(1.1, sigma), 1.0/pow(1.5,length), 1.0/pow(1.2,nu), 1.0/pow(2.0, tau), LL );
            fclose(f4);
        }

        sigma  = gsl_vector_get( s->x, IDX_SIGMA );
        length = gsl_vector_get( s->x, IDX_LENGTH );
        nu     = gsl_vector_get( s->x, IDX_NU );
        tau    = gsl_vector_get( s->x, IDX_TAU );
        LL     = -s->fval;

        std::cout << "  loglikelihood at" //2.0/pow(1.1,sigma), 1.0/pow(1.5,length), 1.0/pow(1.2, nu),
                  << "  σ = " << 2.0/pow(1.1, sigma)
                  << ", ℓ = " << 1.0/pow(1.5, length)
                  << ", ν = " << 1.0/pow(1.2, nu)
                  << ", tau = " << 1.0/pow(2.0, tau)
                  << "   is " << LL
                  << std::endl;
    } while (( status == GSL_CONTINUE ) && ( iter < max_iter ));

    gsl_multimin_fminimizer_free( s );
    gsl_vector_free( ss );
    gsl_vector_free( x );

    return LL;
}

//
// main function
//
int
main ( int      argc,
       char **  argv )
{
    INIT();
    
    std::string  datafile = "datafile.txt";
    //Important ! sigma= 2.0/pow(1.1,sigma), length = 1.0/pow(1.5,length), 1.0/pow(1.2, nu),
    double  sigma  = 1.0; //take these values from previous experiments (Part 1a)
    double  length = 7.0; //means, cov length = 1/pow(1.5, 7)
    double  nu     = 4.0; //means, nu = 1/pow(1.2, 2)
    double  tau    = 20;
    
    //
    // define command line options
    //

    options_description             all_opts;
    options_description             vis_opts( "usage: loglikelihood [options] datafile\n  where options include" );
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
        ( "sigma",       value<double>(), ": sigma parameter" )
        ( "nu",          value<double>(), ": nu parameter" )
        ( "length",      value<double>(), ": length parameter" )
        ( "tau",         value<double>(), ": tau paramater" )
        ;
    
    hid_opts.add_options()
        ( "data",        value<std::string>(), ": datafile defining problem" )
        ;

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
    if ( vm.count( "sigma"     ) ) sigma    = vm["sigma"].as<double>();
    if ( vm.count( "nu"        ) ) nu       = vm["nu"].as<double>();
    if ( vm.count( "length"    ) ) length   = vm["length"].as<double>();
    if ( vm.count( "tau"       ) ) tau      = vm["tau"].as<double>();

    // default to general eps
    if ( fac_eps == -1 )
        fac_eps = eps;
    
    if ( vm.count( "data" ) )
        datafile = vm["data"].as<std::string>();
    else
    {
        std::cout << "usage: loglikelihood [options] datafile" << std::endl;
        exit( 1 );
    }// if


    std::cout << "initial parameters : " << sigma << ", " << nu << ", " << length << ", " << tau << std::endl;
    
    LogLikeliHoodProblem  problem( datafile );

    auto  LL = maximize_likelihood( sigma, length, nu, tau, problem );

    DONE();
    
}
