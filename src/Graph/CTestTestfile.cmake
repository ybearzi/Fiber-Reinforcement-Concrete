# CMake generated Testfile for 
# Source directory: /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph
# Build directory: /usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(GraphTest "itkGraphTest")
ADD_TEST(GCTest_1 "itkBoykovGraphCutFilterTest" "2" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/r85slice.png" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_1.nii.gz")
ADD_TEST(GCTest_1_ImageCompare "ImageCompare" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_1.nii.gz" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/r85slice_out.nii.gz")
ADD_TEST(GCTest_2 "itkBoykovGraphCutFilterTest" "2" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/tools.jpg" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_2.nii.gz")
ADD_TEST(GCTest_2_ImageCompare "ImageCompare" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_2.nii.gz" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/tools_out.nii.gz")
ADD_TEST(GCTest_3 "itkBoykovGraphCutFilterTest" "2" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/peppers.jpg" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_3.nii.gz")
ADD_TEST(GCTest_3_ImageCompare "ImageCompare" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_3.nii.gz" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/peppers_out.nii.gz")
ADD_TEST(GCTest_4 "itkBoykovGraphCutFilterTest" "2" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/ct_scan.jpg" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_4.nii.gz")
ADD_TEST(GCTest_4_ImageCompare "ImageCompare" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_4.nii.gz" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/ct_scan_out.nii.gz")
ADD_TEST(GCTest_5 "itkBoykovGraphCutFilterTest" "3" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/lungs.nii.gz" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_5.nii.gz")
ADD_TEST(GCTest_5_ImageCompare "ImageCompare" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/test_5.nii.gz" "/usr/sci/projects/neuro/Pilot-Projects/FiberConcrete/src/Graph/lungs_out.nii.gz")
