#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>


#include <rrl/itkLocationRegistration.h>


//:
// \file
// \brief  Register and recognize volume location using itkLocationRegistration object.
// \author Michal Sofka
// \date   Oct 2007


int
main( int argc, char* argv[] )
{
  // e.g. location_registration.exe 1089/000002cropped 1091/000002cropped -trans 1089_000002cropped-1091_000002cropped.vtk
  // e.g. location_registration.exe 4748/2whole4747/2whole -matches  nodules -segment watershed
  //
  vul_arg< const char* >        fixedImageDir            ( 0, "Fixed Dicom volume directory" );
  vul_arg< const char* >        movingImageDir           ( 0, "Moving Dicom volume directory" );
  vul_arg< const char* >        transformFile            ( "-trans", "ITK transform file (transform used to warp moving points)", 0 );
  vul_arg< std::string >        prefix                   ( "-matches", "Prefix of the matches folder", "matches" );
  vul_arg< std::string >        segmentationSuffix       ( "-segment", "Suffix of the oversegmentation of fixed and moving Dicom volumes", "" );

  vul_arg_parse( argc, argv );

  vcl_cout << "Command: ";
  for( int i = 0; i < argc; ++i )
    vcl_cout << argv[i] << " ";
  vcl_cout << vcl_endl;


  typedef itk::LocationRegistration  LocationRegistrationType;
  LocationRegistrationType::Pointer  locreg = LocationRegistrationType::New();

  locreg->SetFixedImageDir( fixedImageDir() );
  locreg->SetMovingImageDir( movingImageDir() );
  if( transformFile.set() ) locreg->SetTransformFile( transformFile() );
  locreg->SetPrefix( prefix() );
  locreg->SetSegmentationSuffix( segmentationSuffix() );

  int exit_status = locreg->Run();

  return exit_status;
}
