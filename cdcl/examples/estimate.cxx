#include <vcl_fstream.h>
#include <vcl_string.h>
#include <vcl_iomanip.h>

#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <vnl/vnl_math.h>
#include <vnl/vnl_matlab_read.h>

#include <cdcl/cdcl_trans_similarity2d.h>
#include <cdcl/cdcl_trans_affine.h>
#include <cdcl/cdcl_trans_rigid3d.h>
#include <cdcl/cdcl_estimation_abs.h>
#include <cdcl/cdcl_estimation.h>
#include <cdcl/cdcl_estimation_ICP.h>
#include <cdcl/cdcl_estimation_symmetric.h>
#include <cdcl/cdcl_utils.h>
#include <cdcl/display/cdcl_display.h>

#define PIECEWISE 0

#if PIECEWISE
#include <cdcl/cdcl_region.h>
#include <cdcl/cdcl_estimation_piecewise.h>
#include <cdcl/display/cdcl_display_regions.h>
#endif

#include <vgui/vgui.h>
#include <vgui/vgui_window.h>
#include <vgui/vgui_adaptor.h>
#include <vgui/vgui_shell_tableau.h>
#include <vgui/vgui_viewer2D_tableau.h>
#include <vgui/vgui_clear_tableau.h>
#include <vgui/vgui_viewer3D_tableau.h>


//:
// \file
// \brief   Registration using CDC library.
// \author  Michal Sofka
// \date    Mar 2006


template < unsigned int dim, unsigned int dof >
void register_points( vnl_matrix< double > const                                   &  fixed,
                      vnl_matrix< double > const                                   &  moving,
                      vbl_smart_ptr< cdcl_trans< dim, dof > >                      &  transform,
                      unsigned int                                                    method,
                      void *                                                          caller = 0,
                      // cannot import type from cdcl_estimation.h because gcc complains
                      void   (*display_points) ( void *  caller,
                                                 vcl_vector< vbl_smart_ptr< cdcl_feature<dim> > > const &  moving,
                                                 vcl_vector< vbl_smart_ptr< cdcl_feature<dim> > > const &  fixed,
                                                 vcl_vector< vbl_smart_ptr< cdcl_match<dim> > >   const &  matches,
                                                 vbl_smart_ptr< cdcl_trans< dim, dof > >          const &  transform,
                                                 unsigned int                                              iteration ) = 0
#if PIECEWISE
                      // cannot import type from cdcl_estimation.h because gcc complains
                      , void   (*finalize_display) ( void *  caller, 
                                                     vcl_vector< vbl_smart_ptr< cdcl_region<dim, dof> > >  const &  regions,
                                                     unsigned int                                                   iteration ) = 0 )
#else
                                                 )
#endif
{
  vcl_vector< typename cdcl_feature< dim >::sptr >  fixed_fea, moving_fea;
  vnl_vector_fixed< double, dim >  center_fixed, center_moving;
  double  avg_rad_fixed, avg_rad_moving;

  // parse the raw data into feature set vectors
  cdcl_parse_raw_data( fixed, fixed_fea );
  cdcl_parse_raw_data( moving, moving_fea );

  if( dim == 2 ) {
    cdcl_normalize_data( fixed_fea, center_fixed, avg_rad_fixed );
    cdcl_normalize_data( moving_fea, center_moving, avg_rad_moving );

    // normalize initial transform
    transform->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

    // set up Covariance Driven Correspondence object
    typename cdcl_estimation_abs<dim, dof>::sptr  estimation;
    switch( method ) {
      case 0 :     {// SYMMETRIC, full covariance
                   typename cdcl_trans< dim, dof >::sptr  inverse_normalized = transform->inverse();
                   estimation = new cdcl_estimation_symmetric<dim, dof>( moving_fea, fixed_fea, transform, inverse_normalized );
                   break;}
      case 1 :     {// SYMMETRIC, transfer error covariance
                   typename cdcl_trans< dim, dof >::sptr  inverse_normalized = transform->inverse();
                   estimation = new cdcl_estimation_symmetric<dim, dof, cdcl_estimation_transfer< dim, dof> >( moving_fea, fixed_fea, transform, inverse_normalized );
                   break;}
      case 2 :     {// FORWARD, full covariance
                   estimation = new cdcl_estimation<dim, dof>( moving_fea, fixed_fea, transform );
                   break;}
      case 3 :     {// FORWARD, transfer error covariance
                   estimation = new cdcl_estimation_transfer<dim, dof>( moving_fea, fixed_fea, transform );
                   break;}
      case 4 :     {// FORWARD, ICP
                   estimation = new cdcl_estimation_ICP<dim, dof>( moving_fea, fixed_fea, transform );
                   break;}
      #if PIECEWISE
      case 5 :     {// PIECEWISE, transfer error covariance
                   typedef cdcl_estimation_transfer< dim, dof >  estimation_transfer_type;
                   estimation = new cdcl_estimation_piecewise<dim, dof, estimation_transfer_type>( moving_fea, fixed_fea, transform );
                   break;}
      #endif
    }

    // estimate
    if( caller ) {
      #if PIECEWISE
      if( finalize_display ) { // method 4: piecewise estimation
        //typedef cdcl_estimation_transfer< dim, dof >  estimation_transfer_type;
        //cdcl_estimation_piecewise< dim, dof, estimation_transfer_type >*  estimation_piecewise = dynamic_cast<cdcl_estimation_piecewise< dim, dof, estimation_transfer_type >*>( estimation.as_pointer() );
        //vbl_smart_ptr< cdcl_estimation_piecewise< dim, dof, estimation_transfer_type > >  estimation_piecewise_sptr( estimation_piecewise );
        //estimation_piecewise_sptr->run( caller, display_points, finalize_display );  // for piecewise estimation
      }
      else
      #endif
        estimation->run( caller, display_points );
    }
    else {
      estimation->run();
    }
  }
  else { // dim == 3
    cdcl_normalize_data( moving_fea, center_moving, avg_rad_moving ); // data need to be normalized for computing covariances inside subsample (or pass in centers)

    if( transform->is_type( cdcl_trans_rigid3d::type_id() ) ) {
      // RIGID
      // Same radius for fixed and moving
      avg_rad_fixed = avg_rad_moving;
      cdcl_normalize_data_known_radius( fixed_fea, center_fixed, avg_rad_fixed );
    }
    else {
      // AFFINE
      cdcl_normalize_data( fixed_fea, center_fixed, avg_rad_fixed );
    }

    // normalize initial transform
    transform->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

    vcl_vector< vcl_vector< typename cdcl_feature< dim >::sptr > >  fixed_mult, moving_mult;
    vcl_vector< double > dummy, fixed_spacing;
    subsample_data_fine_covariances( moving_fea, moving_mult, dummy );
    subsample_data_fine_covariances( fixed_fea, fixed_mult, fixed_spacing );

    // set up Covariance Driven Correspondence object
    typename cdcl_estimation_abs<dim, dof>::sptr  estimation;
    switch( method ) {
      //case 0 :     {// SYMMETRIC, full covariance
      //             typename cdcl_trans< dim, dof >::sptr  inverse_normalized = transform->inverse();
      //             estimation = new cdcl_estimation_symmetric<dim, dof>( moving_fea, fixed_fea, transform, inverse_normalized );
      //             break;}
      //case 1 :     {// SYMMETRIC, transfer error covariance
      //             typename cdcl_trans< dim, dof >::sptr  inverse_normalized = transform->inverse();
      //             estimation = new cdcl_estimation_symmetric<dim, dof, cdcl_estimation_transfer< dim, dof> >( moving_fea, fixed_fea, transform, inverse_normalized );
      //             break;}
      case 2 :     {// FORWARD, full covariance
                   estimation = new cdcl_estimation<dim, dof>( moving_mult, fixed_mult, transform, fixed_spacing );
                   break;}
      //case 3 :     {// FORWARD, transfer error covariance
      //             estimation = new cdcl_estimation_transfer<dim, dof>( moving_fea, fixed_fea, transform );
      //             break;}
      case 4 :     {// FORWARD, ICP
                   estimation = new cdcl_estimation_ICP<dim, dof>( moving_mult, fixed_mult, transform, fixed_spacing );
                   break;}
      //case 5 :     {// PIECEWISE, transfer error covariance
      //             estimation = new cdcl_estimation_piecewise<dim, dof>( moving_fea, fixed_fea, transform );
      //             break;}
      default:     vcl_cout << "Selected method: " << method << " not supported." << vcl_endl;
                   return;
    }

    if( caller )
      estimation->run( caller, display_points );
    else
      estimation->run();
  }

  // convert transformation back to unnormalized coordinates
  transform->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );

  vcl_cout << "Final transform: " << vcl_endl;
  transform->print( vcl_cout );
}


int main( int argc, char *argv[] )
{
#define ONLY_RIGID 1

#if !ONLY_RIGID
  vul_arg<vcl_string>           fixed_points             ( 0, "Fixed points" );
  vul_arg<vcl_string>           moving_points            ( 0, "Moving points" );
  vul_arg<bool>                 no_display               ( "-no_disp", "Do not use vgui display", false );  // only use "-no_disp" to set, no other arguments
  vul_arg<unsigned int>         model                    ( "-model", "Model: 2D: Similarity2D [0], Affine2D [1], 3D: Rigid3D [0], Affine3D [1]", 0 );
  #if PIECEWISE
  vul_arg<unsigned int>         method                   ( "-method", "Method: SYMMETRIC, full covariance [0]\nSYMMETRIC, transfer error covariance [1]\nFORWARD, full covariance [2]\nFORWARD, transfer error covariance [3]\nFORWARD, ICP [4]\nPIECEWISE, transfer error covariance [5]", 0 );
  #else
  vul_arg<unsigned int>         method                   ( "-method", "Method: SYMMETRIC, full covariance [0]\nSYMMETRIC, transfer error covariance [1]\nFORWARD, full covariance [2]\nFORWARD, transfer error covariance [3]\nFORWARD, ICP [4]", 0 );
  #endif
  vul_arg<vcl_vector<double> >  trans_par                ( "-trans", "Use initial transformation parameters" );
  vul_arg<vcl_string>           output_trans             ( "-output_trans", "Save final transformation into a given file", "" );
  vul_arg< const char* >        output_dir               ( "-out_dir", "Output directory name", "" );
#else
  vul_arg<vcl_string>           fixed_points             ( 0, "Fixed points" );
  vul_arg<vcl_string>           moving_points            ( 0, "Moving points" );
  vul_arg<vcl_vector<double> >  trans_par                ( "-trans", "Use initial transformation parameters" );
  vul_arg<vcl_string>           output_trans             ( "-output_trans", "Save final transformation into a given file", "" );
  vul_arg< const char* >        output_dir               ( "-out_dir", "Output directory name", "" );
#endif

  vul_arg_parse( argc, argv );

#if ONLY_RIGID
  vul_arg<bool>                 no_display               ( "", "", true );  // only use "-no_disp" to set, no other arguments
  vul_arg<unsigned int>         model                    ( "", "", 0 );
  vul_arg<unsigned int>         method                   ( "-method", "", 2 );
#endif

  if( !no_display() ) vgui::init(argc, argv);

  // read fixed points
 
  //vcl_string fixed_name = "Q";
  //// when saved with cmd:  save points Q -v4, complains about type, otherwise doesn't work at all
  //bool have_read = vnl_matlab_read_or_die( data_file, fixed, fixed_name.c_str() );
  
  // load ascii data that was saved as:
  // save Q.dat Q -ascii
  vcl_ifstream  fixed_data_file( fixed_points().c_str() );
  vnl_matrix< double >  fixed;
  bool have_read = fixed.read_ascii( fixed_data_file );
  fixed_data_file.close();

  if( !have_read ) {
    vcl_cerr << "Error: Can't read fixed data." << vcl_endl;
    return 1;
  }

  // read moving points
  // load ascii data that was saved in matlab as:
  // save P.dat P -ascii
  vcl_ifstream  moving_data_file( moving_points().c_str() );
  vnl_matrix< double >  moving;
  have_read = moving.read_ascii( moving_data_file );
  moving_data_file.close();

  if( !have_read ) {
    vcl_cerr << "Error: Can't read moving data." << vcl_endl;
    return 1;
  }

  if( fixed.cols() != moving.cols() ) {
    vcl_cerr << "Error: Moving and fixed data files do not have the same number of columns." << vcl_endl;
    return 1;
  }

  // fixed and moving file name
  vcl_string fixed_name = vul_file::strip_extension( vul_file::basename( fixed_points() ) );
  vcl_string moving_name = vul_file::strip_extension( vul_file::basename( moving_points() ) );

  vcl_string  cwd = vul_file::get_cwd();
  vul_file::change_directory( output_dir() );

  if( moving.columns() == 5 || moving.columns() == 2 ) {
    double scale, angle;
    vnl_double_2  translation;
    if( trans_par.set() ) { // parameters supplied from the command line
      if( trans_par().size() != 4 ) {
        vcl_cerr << "When supplying parameters, the total number must be 4." << vcl_endl;
        return 1;
      }
      scale = trans_par()[0];
      angle = trans_par()[1] * vnl_math::pi / 180.0;
      translation[0] = trans_par()[2];
      translation[1] = trans_par()[3];
    }
    else {
      scale = 1.0;
      angle = 0.0 * vnl_math::pi / 180.0;
      translation[0] = 10;
      translation[1] = 100;
    }

    // make directory with name indicating parameters
    vcl_ostringstream  dir;
    dir << fixed_name << "--" << moving_name << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << scale << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( angle*180/vnl_math::pi ) << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( translation[0] ) << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( translation[1] );

    vcl_string  new_dir = dir.str();

    vul_file::make_directory( new_dir.c_str() );
    vul_file::change_directory( new_dir.c_str() );


    // set up viewer
    vgui_clear_tableau_new     background_tableau;
    background_tableau->set_colour( 1.0, 1.0, 1.0 );
    vgui_easy2D_tableau_new    matches_easy2D( background_tableau );
    if( !no_display() ) {
      vgui_viewer2D_tableau_new  viewer( matches_easy2D );
      vgui_shell_tableau_new     shell( viewer );

      vgui_window* win = vgui::adapt( shell, 1024, 1200 );
      win->reposition( 0, 0 );
      win->show();
      viewer->center_image( 0, 0 );

    }

    // need a macro because trans defined in if statements (dim and dof are different)
    #define OUTPUT_TRANS \
    if( output_trans.set() ) { \
      vcl_ofstream  out_trans_file( output_trans().c_str() ); \
      out_trans_file << transform->get_A() << transform->get_translation() << vcl_endl; \
      out_trans_file.close(); \
    }

    // set up the initial transform
    if( model() == 0 ) {  // Similarity
      const unsigned int dim = cdcl_trans_similarity2d::dim;
      const unsigned int dof = cdcl_trans_similarity2d::dof;

      cdcl_trans<dim,dof>::sptr  transform = new cdcl_trans_similarity2d( scale, angle, translation );

      if( !no_display() ) {
        if( method() == 1 || method() == 3 ) { // SYMMETRIC, transfer error covariance  OR  FORWARD, transfer error covariance
          cdcl_display_callback2d_transfer< dim, dof >  callback( matches_easy2D );
          register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_callback2d_transfer< dim, dof >::display_points );
        }
        #if PIECEWISE
        else if( method() == 5 ) { // PIECEWISE, transfer error covariance
          cdcl_display_regions_callback2d< dim, dof >  callback( matches_easy2D );
          register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_regions_callback2d< dim, dof >::display_points, &cdcl_display_regions_callback2d< dim, dof >::finalize_display );
        }
        #endif
        else {
          cdcl_display_callback2d< dim, dof >  callback( matches_easy2D );
          register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_callback2d< dim, dof >::display_points );
        }
      }
      else {
        register_points( fixed, moving, transform, method() );
      }
      OUTPUT_TRANS;
    }
    else if( model() == 1 ) {  // Affine
      const unsigned int dim = cdcl_trans_affine<2>::dim;
      const unsigned int dof = cdcl_trans_affine<2>::dof;

      double sina = vcl_sin( angle );
      double cosa = vcl_cos( angle );

      vnl_matrix_fixed< double, dim, dim >  A;
      A( 0, 0 ) =  scale * cosa;
      A( 0, 1 ) = -scale * sina;
      A( 1, 0 ) =  scale * sina;
      A( 1, 1 ) =  scale * cosa;
      cdcl_trans<dim,dof>::sptr  transform = new cdcl_trans_affine<2>( A, translation, vnl_double_2( 0.0, 0.0 ) );

      if( !no_display() ) {
        if( method() == 1 || method() == 3 ) { // SYMMETRIC, transfer error covariance  OR  FORWARD, transfer error covariance
          cdcl_display_callback2d_transfer< dim, dof >  callback( matches_easy2D );
          register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_callback2d_transfer< dim, dof >::display_points );
        }
        #if PIECEWISE
        else if( method() == 5 ) { // PIECEWISE, transfer error covariance
          cdcl_display_regions_callback2d< dim, dof >  callback( matches_easy2D );
          register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_regions_callback2d< dim, dof >::display_points, &cdcl_display_regions_callback2d< dim, dof >::finalize_display );
        }
        #endif
        else {
          cdcl_display_callback2d< dim, dof >  callback( matches_easy2D );
          register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_callback2d< dim, dof >::display_points );
        }
      }
      else {
        register_points( fixed, moving, transform, method() );
      }
      OUTPUT_TRANS;
    }
    else {
      vcl_cerr << "Specified model is invalid ..." << vcl_endl;
      return 1;
    }

  }
  else if( moving.columns() == 9 || moving.columns() == 3 ) {
    const unsigned int dim = cdcl_trans_rigid3d::dim;
    const unsigned int dof = cdcl_trans_rigid3d::dof;

    double alpha, beta, gamma;
    vnl_vector_fixed< double, dim >  translation;

    if( trans_par.set() ) { // parameters supplied from the command line
      if( trans_par().size() != 6 ) {
        vcl_cerr << "When supplying parameters, the total number must be 6." << vcl_endl;
        return 1;
      }
      alpha = trans_par()[0] * vnl_math::pi / 180.0;
      beta  = trans_par()[1] * vnl_math::pi / 180.0;
      gamma = trans_par()[2] * vnl_math::pi / 180.0;

      translation[0] = trans_par()[3];
      translation[1] = trans_par()[4];
      translation[2] = trans_par()[5];
    }
    else {
      alpha =  5.0 * vnl_math::pi / 180.0;
      beta  = 10.0 * vnl_math::pi / 180.0;
      gamma = 45.0 * vnl_math::pi / 180.0;

      translation[0] = 0;//10;
      translation[1] = 0;// 5;
      translation[2] = 0;// 5;

      alpha = 0.0872665;
      beta  = 0.174533;
      gamma = 0.785398;
      translation[0] =  0.206294;
      translation[1] = -0.0229259;
      translation[2] =  0.00314579;
    }

    vcl_ostringstream  dir;
    dir << fixed_name << "--" << moving_name << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( alpha*180.0/vnl_math::pi ) << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( beta*180.0/vnl_math::pi ) << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( gamma*180.0/vnl_math::pi ) << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( 100.0*translation[0] )/100.0 << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( 100.0*translation[1] )/100.0 << "_"
        << vcl_setw( 2 ) << vcl_setfill( '0' ) << vnl_math_rnd( 100.0*translation[2] )/100.0;

    vcl_string  new_dir = dir.str();

    vul_file::make_directory( new_dir.c_str() );
    vul_file::change_directory( new_dir.c_str() );

    cdcl_trans<dim,dof>::sptr  transform_rigid = new cdcl_trans_rigid3d( alpha, beta, gamma, translation );

    vgui_easy3D_tableau_new    matches_easy3D;
    if( !no_display() ) {
      vgui_viewer3D_tableau_new  viewer( matches_easy3D );
      vgui_shell_tableau_new     shell( viewer );

      // set view parameters
      viewer->token.fov = 45.0f;
      viewer->token.scale = 0.02f;
      viewer->token.quat[0] = 0.30f;
      viewer->token.quat[1] = 0.28f;
      viewer->token.quat[2] = 0.20f;
      viewer->token.quat[3] = 0.89f;

      viewer->token.trans[0] = 0.0f;
      viewer->token.trans[1] = 0.0f;
      viewer->token.trans[2] = -10.0f;

      vgui_window* win = vgui::adapt( shell, 1200, 1600 );
      win->reposition( /*-120*/100, 0 );
      win->show();
    }

    if( model() == 0 ) {  // Rigid
      cdcl_trans<dim,dof>::sptr  transform = transform_rigid;
      if( !no_display() ) {
        cdcl_display_callback3d< dim, dof >  callback( matches_easy3D );
        register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_callback3d< dim, dof >::display_points );
      }
      else {
        register_points( fixed, moving, transform, method() );
      }
      OUTPUT_TRANS;
    }
    else if( model() == 1 ) {  // Affine
      const unsigned int dof = cdcl_trans_affine<3>::dof;
      cdcl_trans<dim,dof>::sptr  transform_affine = new cdcl_trans_affine<3>( transform_rigid->get_A(), transform_rigid->get_translation(), transform_rigid->center_moving_ );
      cdcl_trans<dim,dof>::sptr  transform = transform_affine;
      if( !no_display() ) {
        cdcl_display_callback3d< dim, dof >  callback( matches_easy3D );
        register_points( fixed, moving, transform, method(), (void*) &callback, &cdcl_display_callback3d< dim, dof >::display_points );
      }
      else {
        register_points( fixed, moving, transform, method() );
      }
      OUTPUT_TRANS;
    }

    vul_file::change_directory( cwd.c_str() );

  }
  else {
    vcl_cerr << "Error: Don't know what estimation to run for data with: " << moving.columns() << " columns." << vcl_endl;
    return 1;
  }

  return 0;
}
