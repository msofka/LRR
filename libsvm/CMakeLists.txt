
INCLUDE_DIRECTORIES( ${lrr_SOURCE_DIR} )

ADD_LIBRARY( libsvm
  svm.h  svm.cpp
  svm_interface.h  svm_interface.cpp
)

ADD_EXECUTABLE(svm-predict svm-predict.cpp)
TARGET_LINK_LIBRARIES(svm-predict libsvm)

ADD_EXECUTABLE(svm-train svm-train.cpp)
TARGET_LINK_LIBRARIES(svm-train libsvm)

ADD_EXECUTABLE(svm-scale svm-scale.cpp)
TARGET_LINK_LIBRARIES(svm-scale libsvm)

