# set( libname "heedpp") set( lib_fullver ${PROJECT_VERSION}) set( lib_apiver 1 )
# add_library(${libname} SHARED "" )

# TODO Heed actually is interdependent from Garfield but just for the Random number
# production and in that just for the poisson. It can be split but the whole random number
# generation system should be reviewed. For now I just add the sources to the main
# garfield library. Needs the absolute path because it will be used from the garfield
# folder
# Because the sources are added to the sources of a target defined in another directory
# we need to provide the absolute path
target_sources(
  Garfield
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/heed++/code/BGMesh.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/ElElasticScat.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/EnTransfCS.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/EnTransfCS_BGM.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/EnergyMesh.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedCluster.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedCondElectron.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedDeltaElectron.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedDeltaElectronCS.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedMatterDef.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedParticle.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedParticle_BGM.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/HeedPhoton.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/PairProd.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/PhotoAbsCS.cpp
          ${CMAKE_CURRENT_LIST_DIR}/heed++/code/PhotoAbsCSLib.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/box.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/circumf.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/gparticle.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/mparticle.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/plane.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/polyline.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/straight.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/surface.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/trajestep.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/vec.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/geometry/volume.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/ioniz/bethe_bloch.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/ioniz/e_cont_enloss.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/DoubleAc.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/PolLeg.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/cubic.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/kinem.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/linexi2.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/lorgamma.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/math/parabola.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/matrix/inverse.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/matrix/multiply.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/matter/AtomDef.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/matter/GasDef.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/matter/MatterDef.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/matter/MoleculeDef.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/particle/eparticle.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/particle/particle_def.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/random/PointsRan.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/random/chisran.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/random/pois.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/random/rnorm.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/safetl/AbsArr.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/safetl/AbsPtr.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/stream/definp.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/stream/findmark.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/stream/prstream.cpp
          ${CMAKE_CURRENT_LIST_DIR}/wcpplib/util/FunNameStack.cpp
          ${CMAKE_CURRENT_LIST_DIR}/HeedFieldMap.cpp)

target_compile_definitions(Garfield PRIVATE GARFIELD_HEED_INTERFACE FUNNAMESTACK)
# garfield_set_all_default_properties(${libname}) # Sets the API version of the library
# set_property( TARGET ${libname} PROPERTY SOVERSION ${lib_apiver} ) set_property( TARGET
# ${libname} PROPERTY CXX_STANDARD 11) set_property( TARGET ${libname} PROPERTY
# CMAKE_CXX_STANDARD_REQUIRED ON) set_property( TARGET ${libname} PROPERTY CXX_EXTENSIONS
# OFF)

# target_include_directories(${libname} PUBLIC
# $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# install( TARGETS ${libname} EXPORT "${PROJECT_NAME}Targets" LIBRARY DESTINATION
# ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION
# ${CMAKE_INSTALL_BINDIR} )

# --- install the Heed database--------------------------------------------------
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/heed++/database DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Heed)
