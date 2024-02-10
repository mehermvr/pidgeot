set(EIGEN_BUILD_DOC
    OFF
    CACHE BOOL "Don't build Eigen docs")
set(EIGEN_BUILD_TESTING
    OFF
    CACHE BOOL "Don't build Eigen tests")
set(EIGEN_BUILD_PKGCONFIG
    OFF
    CACHE BOOL "Don't build Eigen pkg-config")
set(EIGEN_BUILD_BLAS
    OFF
    CACHE BOOL "Don't build blas module")
set(EIGEN_BUILD_LAPACK
    OFF
    CACHE BOOL "Don't build lapack module")

include(FetchContent)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 287c8017808fb4c7a651ba4c02363a773e2f0c46) # hash for v 3.4 commit
FetchContent_MakeAvailable(eigen)
