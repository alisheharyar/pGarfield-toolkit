OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include/Garfield
LIBDIR = $(GARFIELD_HOME)/Library
HEEDDIR = $(GARFIELD_HOME)/Heed
MAGDIR = $(GARFIELD_HOME)/Magboltz
INCDIR_FLAG = $(GARFIELD_HOME)/Include
NEBDIR_FLAG = $(GARFIELD_HOME)/NeBem

HEADERS = $(wildcard $(INCDIR)/*.hh)

SOURCES = $(wildcard $(SRCDIR)/*.cc)

ifndef PARALLEL
	HEADERS := $(filter-out $(wildcard $(INCDIR)/RandomEngineMPI.hh),$(HEADERS))
	HEADERS := $(filter-out $(wildcard $(INCDIR)/RandomEngineMPIServer.hh),$(HEADERS))
	SOURCES := $(filter-out $(wildcard $(SRCDIR)/RandomEngineMPI.cc),$(SOURCES))
	SOURCES := $(filter-out $(wildcard $(SRCDIR)/RandomEngineMPIServer.cc),$(SOURCES))
endif

OBJECTS = $(subst $(SRCDIR),$(OBJDIR),$(SOURCES:.cc=.o))
OBJECTS += $(OBJDIR)/magboltz.o 
#OBJECTS += $(OBJDIR)/GarfieldDict.o

TARGETS = $(OBJECTS)
TARGETS += heed

# Fortran compiler
FC = gfortran

# C++ compiler if parallel mode is enabled
ifdef PARALLEL
CXX = mpicxx
endif

# Compilation flags
CFLAGS = -Wall -Wextra -pedantic -ansi -Wabi -Wno-long-long -Woverloaded-virtual \
	 `root-config --cflags` \
        -fpic -fno-common -c \
	-I$(INCDIR_FLAG) -I$(NEBDIR_FLAG) -I$(HEEDDIR) -DINS_CRETURN 

ifdef PARALLEL
CFLAGS += -DPARALLEL
endif

FFLAGS = -fpic -c

# Optimization flags
# CFLAGS += -Os
# FFLAGS += -Os
CFLAGS += -O2
FFLAGS += -O2

# Debug flags
# CFLAGS += -g
# FFLAGS += -g
# Profiling flag
 CFLAGS += -pg

print-%  : ; @echo $* = $($*)

# Linking flags
LDFLAGS = `root-config --glibs` `root-config --ldflags`-lGeom \
	-lgfortran -lm
LDFLAGS+= --evelibs` -lz
LDFLAGS+= --no-relax

all:	$(TARGETS)
	@echo Creating library libGarfield...
	@ar rcs $(LIBDIR)/libGarfield.a $(OBJECTS) \
	$(wildcard $(OBJDIR)/Heed/*.o)
	@ranlib $(LIBDIR)/libGarfield.a
	@touch $(OBJDIR)/last_updated_on
	@echo Finished.

.PHONY:	heed

installdirs : 
	@if [ ! -d $(GARFIELD_HOME)/Library/$(BFARCH) ] ; then mkdir -p $(GARFIELD_HOME)/Library/$(BFARCH); \
	    echo "   >>>> Create $(GARFIELD_HOME)/Library/$(BFARCH)"; fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); \
	    echo "   >>>> Create $(OBJDIR)"; \
        else echo " $(OBJDIR) already exists"; fi
	@if [ ! -d $(OBJDIR)/Heed/ ]; then mkdir -p $(OBJDIR)/Heed; \
	    echo "   >>>> Create $(OBJDIR)/Heed"; \
        else echo " $(OBJDIR)/Heed already exists"; fi
        
heed:	
	@echo Compiling Heed...
	@cd $(HEEDDIR); make; cd $(GARFIELD_HOME)
	@touch $(OBJDIR)/last_updated_on

clean:
	@echo Removing object files...
	@$(RM) $(OBJDIR)/*.o
	@echo Removing libraries...
	@$(RM) $(LIBDIR)/*.a
	@cd $(HEEDDIR); make clean
	@echo Removing dictionary...
	@$(RM) $(SRCDIR)/GarfieldDict.C 

$(OBJDIR)/AvalancheGrid.o: \
	$(SRCDIR)/AvalancheGrid.cc \
	$(INCDIR)/AvalancheGrid.hh \
	$(INCDIR)/FundamentalConstants.hh $(INCDIR)/GarfieldConstants.hh \
	$(INCDIR)/Random.hh \
	$(INCDIR)/Sensor.hh $(INCDIR)/Medium.hh $(INCDIR)/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/AvalancheMicroscopic.o: \
	$(SRCDIR)/AvalancheMicroscopic.cc \
	$(INCDIR)/AvalancheMicroscopic.hh \
	$(INCDIR)/FundamentalConstants.hh $(INCDIR)/GarfieldConstants.hh \
	$(INCDIR)/Random.hh \
	$(INCDIR)/Sensor.hh $(INCDIR)/Medium.hh $(INCDIR)/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/AvalancheMC.o: \
	$(SRCDIR)/AvalancheMC.cc $(INCDIR)/AvalancheMC.hh \
	$(INCDIR)/FundamentalConstants.hh $(INCDIR)/GarfieldConstants.hh \
	$(INCDIR)/Random.hh \
	$(INCDIR)/Sensor.hh $(INCDIR)/Medium.hh $(INCDIR)/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@      
$(OBJDIR)/DriftLineRKF.o: \
	$(SRCDIR)/DriftLineRKF.cc $(INCDIR)/DriftLineRKF.hh \
	$(INCDIR)/FundamentalConstants.hh \
	$(INCDIR)/Sensor.hh $(INCDIR)/Medium.hh $(INCDIR)/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
 
$(OBJDIR)/Track.o: \
	$(SRCDIR)/Track.cc $(INCDIR)/Track.hh \
	$(INCDIR)/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/TrackBichsel.o: \
	$(SRCDIR)/TrackBichsel.cc $(INCDIR)/TrackBichsel.hh \
	$(INCDIR)/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@       
$(OBJDIR)/TrackPAI.o: \
	$(SRCDIR)/TrackPAI.cc $(INCDIR)/TrackPAI.hh \
	$(INCDIR)/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/TrackSimple.o: \
	$(SRCDIR)/TrackSimple.cc $(INCDIR)/TrackSimple.hh \
	$(INCDIR)/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@  
$(OBJDIR)/TrackTrim.o: \
	$(SRCDIR)/TrackTrim.cc $(INCDIR)/TrackTrim.hh \
	$(INCDIR)/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@  
$(OBJDIR)/TrackSrim.o: \
	$(SRCDIR)/TrackSrim.cc $(INCDIR)/TrackSrim.hh \
  $(INCDIR)/FundamentalConstants.hh \
	$(INCDIR)/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@     
$(OBJDIR)/TrackHeed.o: \
	$(SRCDIR)/TrackHeed.cc $(INCDIR)/TrackHeed.hh \
	$(INCDIR)/Track.hh $(SRCDIR)/Track.cc \
	$(HEEDDIR)/HeedChamber.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/TrackElectron.o: \
	$(SRCDIR)/TrackElectron.cc $(INCDIR)/TrackElectron.hh \
	$(INCDIR)/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Component.o: \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh \
	$(INCDIR)/Medium.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentConstant.o: \
	$(SRCDIR)/ComponentConstant.cc $(INCDIR)/ComponentConstant.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentUser.o: \
	$(SRCDIR)/ComponentUser.cc $(INCDIR)/ComponentUser.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@       
$(OBJDIR)/ComponentAnalyticField.o: \
	$(SRCDIR)/ComponentAnalyticField.cc \
	$(INCDIR)/ComponentAnalyticField.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentNeBem2d.o: \
	$(SRCDIR)/ComponentNeBem2d.cc $(INCDIR)/ComponentNeBem2d.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/ComponentFieldMap.o: \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/ComponentFieldMap.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentAnsys121.o: \
	$(SRCDIR)/ComponentAnsys121.cc $(INCDIR)/ComponentAnsys121.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentAnsys123.o: \
	$(SRCDIR)/ComponentAnsys123.cc $(INCDIR)/ComponentAnsys123.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentCST.o: \
	$(SRCDIR)/ComponentCST.cc $(INCDIR)/ComponentCST.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentElmer.o: \
	$(SRCDIR)/ComponentElmer.cc $(INCDIR)/ComponentElmer.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentTcad2d.o: \
	$(SRCDIR)/ComponentTcad2d.cc $(INCDIR)/ComponentTcad2d.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/ComponentTcad3d.o: \
	$(SRCDIR)/ComponentTcad3d.cc $(INCDIR)/ComponentTcad3d.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/ComponentVoxel.o: \
	$(SRCDIR)/ComponentVoxel.cc $(INCDIR)/ComponentVoxel.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentTcadBase.o: \
	$(SRCDIR)/ComponentTcadBase.cc $(INCDIR)/ComponentTcadBase.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/ComponentParallelPlate.o: \
	$(SRCDIR)/ComponentParallelPlate.cc $(INCDIR)/ComponentParallelPlate.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@ 
$(OBJDIR)/ComponentNeBem3dMap.o: \
	$(SRCDIR)/ComponentNeBem3dMap.cc $(INCDIR)/ComponentNeBem3dMap.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@ 
$(OBJDIR)/ComponentNeBem3d.o: \
	$(SRCDIR)/ComponentNeBem3d.cc $(INCDIR)/ComponentNeBem3d.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@    
$(OBJDIR)/ComponentGrid.o: \
	$(SRCDIR)/ComponentGrid.cc $(INCDIR)/ComponentGrid.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@ 
$(OBJDIR)/ComponentElmer2d.o: \
	$(SRCDIR)/ComponentElmer2d.cc $(INCDIR)/ComponentElmer2d.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@  
$(OBJDIR)/ComponentComsol.o: \
	$(SRCDIR)/ComponentComsol.cc $(INCDIR)/ComponentComsol.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@       
$(OBJDIR)/ComponentUserMapBase.o: \
	$(SRCDIR)/ComponentUserMapBase.cc $(INCDIR)/ComponentUserMapBase.hh \
	$(SRCDIR)/Component.cc $(INCDIR)/Component.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@    
	
$(OBJDIR)/GeometrySimple.o: \
	$(SRCDIR)/GeometrySimple.cc $(INCDIR)/GeometrySimple.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/GeometryRoot.o: \
	$(SRCDIR)/GeometryRoot.cc $(INCDIR)/GeometryRoot.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@    

$(OBJDIR)/ViewFEMesh.o: \
	$(SRCDIR)/ViewFEMesh.cc $(INCDIR)/ViewFEMesh.hh \
	$(SRCDIR)/ViewDrift.cc $(INCDIR)/ViewDrift.hh \
	$(INCDIR)/ComponentFieldMap.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewField.o: \
	$(SRCDIR)/ViewField.cc $(INCDIR)/ViewField.hh \
	$(INCDIR)/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewDrift.o: \
	$(SRCDIR)/ViewDrift.cc $(INCDIR)/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewMedium.o: \
	$(SRCDIR)/ViewMedium.cc $(INCDIR)/ViewMedium.hh \
	$(INCDIR)/Medium.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewIsochrons.o: \
	$(SRCDIR)/ViewIsochrons.cc $(INCDIR)/ViewIsochrons.hh \
	$(INCDIR)/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewSignal.o: \
	$(SRCDIR)/ViewSignal.cc $(INCDIR)/ViewSignal.hh \
	$(INCDIR)/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewCell.o: \
	$(SRCDIR)/ViewCell.cc $(INCDIR)/ViewCell.hh \
	$(INCDIR)/ComponentAnalyticField.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewBase.o: \
	$(SRCDIR)/ViewBase.cc $(INCDIR)/ViewBase.hh \
	$(INCDIR)/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewGeometry.o: \
	$(SRCDIR)/ViewGeometry.cc $(INCDIR)/ViewGeometry.hh \
	$(INCDIR)/GeometrySimple.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Medium.o: \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh \
	$(INCDIR)/FundamentalConstants.hh \
	$(INCDIR)/Numerics.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumGas.o: \
	$(SRCDIR)/MediumGas.cc $(INCDIR)/MediumGas.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumMagboltz.o: \
	$(SRCDIR)/MediumMagboltz.cc $(INCDIR)/MediumMagboltz.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh $(SRCDIR)/OpticalData.cc \
	$(INCDIR)/FundamentalConstants.hh $(INCDIR)/Random.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/magboltz.o: \
	$(MAGDIR)/magboltz.f
	@echo $@
	@$(FC) $(FFLAGS) $< -o $@
$(OBJDIR)/MediumCdTe.o: \
	$(SRCDIR)/MediumCdTe.cc $(INCDIR)/MediumCdTe.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumGaAs.o: \
	$(SRCDIR)/MediumGaAs.cc $(INCDIR)/MediumGaAs.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumGaN.o: \
	$(SRCDIR)/MediumGaN.cc $(INCDIR)/MediumGaN.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumDiamond.o: \
	$(SRCDIR)/MediumDiamond.cc $(INCDIR)/MediumDiamond.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumSilicon.o: \
	$(SRCDIR)/MediumSilicon.cc $(INCDIR)/MediumSilicon.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Medium.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/OpticalData.o: \
	$(SRCDIR)/OpticalData.cc $(INCDIR)/OpticalData.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/SolidBox.o: \
	$(SRCDIR)/SolidBox.cc $(INCDIR)/SolidBox.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidTube.o: \
	$(SRCDIR)/SolidTube.cc $(INCDIR)/SolidTube.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidSphere.o: \
	$(SRCDIR)/SolidSphere.cc $(INCDIR)/SolidSphere.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidWire.o: \
	$(SRCDIR)/SolidWire.cc $(INCDIR)/SolidWire.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidRidge.o: \
	$(SRCDIR)/SolidRidge.cc $(INCDIR)/SolidRidge.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidHole.o: \
	$(SRCDIR)/SolidHole.cc $(INCDIR)/SolidHole.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidExtrusion.o: \
	$(SRCDIR)/SolidExtrusion.cc $(INCDIR)/SolidExtrusion.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/Solid.o: \
	$(SRCDIR)/Solid.cc $(INCDIR)/Solid.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/Shaper.o: \
	$(SRCDIR)/Shaper.cc $(INCDIR)/Shaper.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@


$(OBJDIR)/Random.o: \
	$(SRCDIR)/Random.cc $(INCDIR)/Random.hh
	$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/RandomEngineGSL.o: \
	$(SRCDIR)/RandomEngineGSL.cc $(INCDIR)/RandomEngineGSL.hh
	$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/RandomEngineRoot.o: \
	$(SRCDIR)/RandomEngineRoot.cc $(INCDIR)/RandomEngineRoot.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/PlottingEngineRoot.o: \
	$(SRCDIR)/PlottingEngineRoot.cc $(INCDIR)/PlottingEngineRoot.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@    
  
$(OBJDIR)/PlottingEngine.o: \
	$(SRCDIR)/PlottingEngine.cc $(INCDIR)/PlottingEngine.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@      

$(OBJDIR)/Numerics.o: \
	$(SRCDIR)/Numerics.cc $(INCDIR)/Numerics.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@  

$(OBJDIR)/Sensor.o: \
	$(SRCDIR)/Sensor.cc $(INCDIR)/Sensor.hh \
	$(INCDIR)/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Polygon.o: \
	$(SRCDIR)/Polygon.cc $(INCDIR)/Polygon.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/QuadTree.o: \
	$(SRCDIR)/QuadTree.cc $(INCDIR)/QuadTree.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/TetrahedralTree.o: \
	$(SRCDIR)/TetrahedralTree.cc $(INCDIR)/TetrahedralTree.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
 
$(OBJDIR)/KDTree.o: \
	$(SRCDIR)/KDTree.cc $(INCDIR)/KDTree.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
 
$(OBJDIR)/TGeoTet.o: \
	$(SRCDIR)/TGeoTet.cc $(INCDIR)/TGeoTet.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/tet_mesh.o: \
	$(SRCDIR)/tet_mesh.cc $(INCDIR)/tet_mesh.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

ifdef PARALLEL

$(OBJDIR)/RandomEngineMPI.o: \
	$(SRCDIR)/RandomEngineMPI.cc $(INCDIR)/RandomEngineMPI.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/RandomEngineMPIServer.o: \
	$(SRCDIR)/RandomEngineMPIServer.cc $(INCDIR)/RandomEngineMPIServer.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

endif

$(OBJDIR)/GarfieldDict.o: \
	$(SRCDIR)/GarfieldDict.C
	@echo $@
	@$(CXX) $(CFLAGS) -DDICT_SKIP_HEED $< -o $@

$(SRCDIR)/GarfieldDict.C: $(HEADERS) $(INCDIR)/LinkDef.h
	@echo Creating dictionary...
	@rootcint -f $@ -c $(CFLAGS) -p $^ 
