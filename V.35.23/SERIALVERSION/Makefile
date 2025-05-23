#################################################################################
# This makefile is for building the stand alone serial version of uquantchem	#
#         using the lapack-3.4.0 library that comes with the package.		#
#										#
# 	  	Petros Souvatzis: petros.souvatzis@fysik.uu.se 			#
# 	   	  Tue Apr 24 14:01:31 CEST 2012					#
#################################################################################

FC = gfortran

FFLAGS =  -ffree-form -ffree-line-length-none  

F77FLAGS = -ffixed-form

LAPACKPATH = /usr/lib64/

BLASPATH = /usr/lib64/

LIB = -L$(LAPACKPATH) -llapack -L$(BLASPATH) -lblas

OFLAG  = -O3 -ftree-vectorize 

SOURCE = datatypemodule.o exchcorrmodule.o random.o Afunc.o RHF.o Fnurec.o Fnu.o CISD.o getKv.o getK.o getJv.o getJ.o gammaf.o fj.o fac.o  \
	 gto.o gtop.o gtopp.o gtoppp.o eeints.o diaghHF.o diagh.o dfac.o checkinputfile.o chargedenssave.o exciteornot.o \
	 hforbitalsave.o binomfac.o basfunkval.o URHF.o primeeintegral.o potential.o overlap.o oneint.o normalize.o ROOTWEIGHT5.o \
	 leastsq.o  makedens.o kinetic.o ijkl.o homolumosave.o gettotalbasis.o rysquad.o rho.o readin.o ROOTWEIGHT4.o \
	 interpol.o gradbasfunkval.o laplacebasfunkval.o hforbitalval.o gradhforbitalval.o laplacehforbitalval.o det.o \
	 stoexponent.o slaterdet.o laplaceslaterdet.o gradslaterdet.o jastrowup.o jastrowdown.o jastrowud.o jastrow.o \
	 DIIS.o findclosestatom.o gradjastrow.o laplacejastrow.o guideforce.o laplacetrialfnk.o trialfnk.o EL.o dqmc.o \
	 preparedifusion.o ROOTWEIGHTMAX3.o readbasis.o primpotential.o primoverlap.o primkinetic.o EP2.o EP2so.o EP2pt5so.o EP3so.o \
         EP2plus3r.o EPP3so.o EPQ3so.o EPP3plus.o EPQ3plus.o EPL3so.o EPL3plusB.o EP2r.o EP3r.o EP2nD.o UEP2.o UEP2nD.o MP2.o uquantchem.o \
	 gradprimoverlap.o gradprimkinetic.o gradprimpotential.o gradprimeeintegral.o forces.o dAfunc.o relax.o linesearchmin.o \
	 massa.o moleculardynamics.o moleculardynamicssoft.o invert.o vxc.o vxcalt.o dvxc.o ngbasfunkval.o nggradbasfunkval.o nglaplacebasfunkval.o \
	 getvxcr.o atomicradii.o lebedev1.o lebedev2.o lebedev3.o lebedev4.o lebedev5.o chebgauss.o sofmu.o pvoronoi.o getvxc.o excdens.o excdr.o exc.o DFT.o quadcheck.o TRACE.o \
	 hessianbasfunk.o hessianrho.o dftforcedens.o getxcforce.o pnonorm.o tmu.o gradpvoronoi.o gradoverlap.o relaxf.o mulliken.o \
	 makedensT.o hessianrhonip.o dftforcedensnip.o getxcforcenip.o

uquantchem:  $(SOURCE) 
	$(FC) $(FFLAGS) $(SOURCE) $(LIB) -o uquantchem.s

clean: 
	-/bin/rm *.o *~ *.mod uquantchem.s

%.o: %.f
	$(FC) $(LIB) $(F77FLAGS) $(OFLAG) -c $<

%.o: %.f90
	$(FC) $(LIB) $(FFLAGS) $(OFLAG) -c $<

lebedev1.o: lebedev1.f90
	$(FC) $(LIB) $(FFLAGS) -O0 -c lebedev1.f90

lebedev2.o: lebedev2.f90
	$(FC) $(LIB) $(FFLAGS) -O0 -c lebedev2.f90

lebedev3.o: lebedev3.f90
	$(FC) $(LIB) $(FFLAGS) -O0 -c lebedev3.f90

lebedev4.o: lebedev4.f90
	$(FC) $(LIB) $(FFLAGS) -O0 -c lebedev4.f90

lebedev5.o: lebedev5.f90
	$(FC) $(LIB) $(FFLAGS) -O0 -c lebedev5.f90
