#include <iostream>
#include <fstream>
#include <string>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

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

// global options
int        nmin      = CFG::Cluster::nmin;
double     eps       = 1e-6;
double     fac_eps   = 1e-6;
double     shift     = 0.0;
bool       use_ldl   = false;

//
// read dataset from file
//
void
read_data ( const std::string &       datafile,
            std::vector< T2Point > &  vertices,
            BLAS::Vector< double > &  Z_data_ )
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

        std::cout << "learning dataset" << std::endl;
        std::cout << N_vtx << std::endl;
        
        vertices.resize( N_vtx );
        Z_data_ = BLAS::Vector< double >( N_vtx );

        int  i = 0;

        for ( auto  p : pos )
            vertices[ i++ ] = p;

        i = 0;
        int j=0;
        for ( auto  v : vals )
        {
          Z_data_( i++ ) = v;
        }
         
         
   }// if
    if ( line == "x,y,ignor" )
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
           // vals.push_back( v );
        }// while

        N_vtx = pos.size();

        std::cout << "learning dataset" << std::endl;
        std::cout << N_vtx << std::endl;
        
        vertices.resize( N_vtx );
        Z_data_ = BLAS::Vector< double >( N_vtx );

        int  i = 0;

        for ( auto  p : pos )
            vertices[ i++ ] = p;

        i = 0;
        
        for ( auto  v : vals )
           Z_data_( i++ ) = 0.0;
    }// if
   
    if ( line == "x,y" )
    {
        std::list< T2Point >  pos;

        while ( std::getline( in, line ) )
        {
            auto    parts = split( line, "," );
            double  x = atof( parts[0].c_str() );
            double  y = atof( parts[1].c_str() );
            
            pos.push_back( T2Point( x, y ) );
        }// while

        N_vtx = pos.size();
        std::cout << "testing dataset" << std::endl;
        std::cout << N_vtx << std::endl;
        
        vertices.resize( N_vtx );

        int  i = 0;

        for ( auto  p : pos )
            vertices[ i++ ] = p;
    }// if
    //else
    //{
    //    std::cout << "you should not be here" << std::endl;
    //    HERROR( ERR_NOT_IMPL, "", "" );
    //}
    
    #else
    
    in >> N_vtx;

    std::cout << "reading " << N_vtx << " datapoints" << std::endl;
    
    vertices.resize( N_vtx );
    Z_data_ = BLAS::Vector< double >( N_vtx );
        
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
struct PredictionProblem
{
    std::vector< T2Point >                vertices;
    std::vector< T2Point >                vertices_predict;
    std::unique_ptr< TCoordinate >        coord;
    std::unique_ptr< TCoordinate >        coord_predict;
    std::unique_ptr< TClusterTree >       ct;
    std::unique_ptr< TClusterTree >       ct_predict;
    std::unique_ptr< TBlockClusterTree >  bct;
    std::unique_ptr< TBlockClusterTree >  bct_predict;
    std::unique_ptr< TVector >            Z;
    std::unique_ptr< TVector >            Z_predict;

    PredictionProblem ( const std::string &  datafile, const std::string &  datafile_predict )
    {
        init( datafile, datafile_predict );
    }

    void
    init ( const std::string &  datafile, const std::string &  datafile_predict )
    {
        BLAS::Vector< double >  Z_data;
        BLAS::Vector< double >  Z_data_predict;

        read_data( datafile, vertices, Z_data );
        read_data( datafile_predict, vertices_predict, Z_data_predict );

        std::cout << "both datasets are succesfully read" << std::endl;
        coord = std::make_unique< TCoordinate >( vertices );
        coord_predict = std::make_unique< TCoordinate >( vertices_predict );
      
        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder( & part_strat, nmin );
    
        ct = ct_builder.build( coord.get() );
        
        //print_vtk( & coord, "ct_coord" );
        ct_predict = ct_builder.build( coord_predict.get() );
        //print_vtk( & coord_predict, "ct_coord_predict" );
        
    
        TStdGeomAdmCond    adm_cond( 2.0, use_min_diam );
        TBCBuilder         bct_builder;
   
        bct = bct_builder.build( ct.get(), ct.get(), & adm_cond );
        bct_predict = bct_builder.build( ct_predict.get(), ct.get(), & adm_cond );
  
        Z   = std::make_unique< TScalarVector >( *ct->root(), std::move( Z_data ) );
        Z_predict   = std::make_unique< TScalarVector >( *ct_predict->root(), std::move( Z_data_predict ) );
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
        
        TMaternCovCoeffFn< T2Point >  matern_coefffn_predict( sigma, length, nu, vertices_predict, vertices );
        TPermCoeffFn< double >        coefffn( & matern_coefffn, ct->perm_i2e(), ct->perm_i2e() );
        TPermCoeffFn< double >        coefffn_predict( & matern_coefffn_predict, ct_predict->perm_i2e(), ct->perm_i2e() );
  
        TACAPlus< double >            aca( & coefffn );
        TACAPlus< double >            aca_predict( & coefffn_predict );
        auto                          acc = fixed_prec( eps );
        TDenseMatBuilder< double >    h_builder( & coefffn, & aca );
        TDenseMatBuilder< double >    h_builder_predict( & coefffn_predict, & aca_predict );
    
        auto                          C        = h_builder.build( bct.get(), acc );
        //TPSMatrixVis  mvis;
        
        //mvis.svd(true).print( C.get(), "myC" ); Output matrix C as a .eps file 
        

       
        auto   C_predict= h_builder_predict.build( bct_predict.get(), unsymmetric, acc ); //A_21 - rectangular matrix
        
        if ( shift != 0.0 )
            add_identity( C.get(), tau*tau );
  
        auto                          fac_acc  = fixed_prec( fac_eps );
        auto                          C_fac    = C->copy();
        auto                          fac_opts = fac_options_t{ point_wise, CFG::Arith::storage_type, false };
    
        if ( use_ldl )
        {
            ldl( C_fac.get(), fac_acc, fac_opts );
            //std::cout <<  "LDL is used " << std::endl;
        }
        else
        {
            chol( C_fac.get(), fac_acc ); 
            //std::cout <<  "Cholesky is used " << std::endl;
        }
        //std::cout << "    |С|_F             = " << norm_F( C.get() ) << std::endl;
        //std::cout << "    |С|_2             = " << norm_2( C.get() ) << std::endl;
        //std::cout << "    |С_fac|_F             = " << norm_F( C_fac.get() ) << std::endl;
        //std::cout << "    |С_fac|_2             = " << norm_2( C_fac.get() ) << std::endl;
        
        std::unique_ptr< TFacInvMatrix >   C_inv;

        if ( use_ldl )
        {
            C_inv = std::make_unique< TLDLInvMatrix >( C_fac.get(), symmetric, point_wise );
        }
        else
        {
            C_inv = std::make_unique< TLLInvMatrix >( C_fac.get(), symmetric );
        }
        //std::cout << "    |С_inv|_2             = " << norm_2( C_inv.get() ) << std::endl;
        
        const size_t                  N     = vertices.size();
        
        TStopCriterion                sstop( 350, 1e-16, 0.0 );
        TCG                           solver( sstop );
        auto                          sol = C->row_vector();
        
        auto                          Z_predict = C_predict->row_vector();
        
        FILE* f1; 
  
        
        solver.solve( C.get(), sol.get(), Z.get(), C_inv.get() );
        
    
        //DBG::write( C_predict.get(), "A12.mat", "A12" );
        //std::cout << "  ||invC_Z|| = " << sol->norm2() << std::endl;
        mul_vec( real_t(1), C_predict.get(), sol.get(), real_t(0), Z_predict.get(), apply_normal );
        

        //std::cout << "  ||Z_predict|| = " << Z_predict->norm2() << std::endl;
        //C_predict->mul_vec( 1.0, sol.get(), 0.0, Z_predict.get(), apply_normal );
        

        ct_predict->perm_i2e()->permute( Z_predict.get() );
        //TMatlabVectorIO  vio;
        //vio.write( Z_predict,  "x.mat", "x" );

        f1 = fopen("111prediction.txt", "w");
        for ( size_t  i = 0; i < Z_predict->size(); i++ )
          fprintf(f1," %3.6e, %3.6e, %3.6e\n",   vertices_predict[i].x(),  vertices_predict[i].y(), Z_predict->entry(i));
        fclose(f1);
        return std::move( Z_predict );
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
    
    //CFG::set_verbosity( 3 );
    INIT();
    
    //std::string  datafile = "datafile.txt";
    //std::string  datafile_predict = "datafile_predict.txt";
    std::string  datafile = "LearningSet.txt";
    std::string  datafile_predict = "TestingSet.txt";
    
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
        ;
    
    hid_opts.add_options()
        ( "dataL",        value<std::string>(), ": datafile defining learning problem" )
        ( "dataT",        value<std::string>(), ": datafile defining testing problem" )
        ;

    // options for command line parsing
    all_opts.add( vis_opts ).add( hid_opts );

    // all "non-option" arguments should be "--data" arguments
    pos_opts.add( "dataL", -1 );
    pos_opts.add( "dataT", -1 );

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
    
    if ( vm.count( "dataL" ) )
        datafile = vm["dataL"].as<std::string>();
    else
    {
        std::cout << "usage: loglikelihood [options] datafile" << std::endl;
        exit( 1 );
    }// if

    if ( vm.count( "dataT" ) )
        datafile_predict = vm["dataT"].as<std::string>();
    else
    {
        std::cout << "usage: loglikelihood [options] datafile" << std::endl;
        exit( 1 );
    }// if

    double  sigma  = 0.538467; //take these values from previous experiments (Part 1a)
    double  length = 0.0106; 
    double  nu     = 2.4705;
    double  tau    = 1.5837e-7;
   
    PredictionProblem  problem( datafile, datafile_predict );
    problem.eval( sigma, length, nu, tau);

    DONE();
}
