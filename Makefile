CC = gcc
CLINK = gcc
CFLAGS = -O2
CLIB= -L${HOME}/lib -lcutil -lm 
XDIR = ${HOME}/exec


OBJS23=	covar3.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o
covar3:	$(OBJS23) Makefile
	$(CLINK) $(CFLAGS) -o covar3 $(OBJS23) $(CLIB)
	mv covar3 $(XDIR)/$@

OBJS23a=	covar3_2halo.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o
covar3_2halo:	$(OBJS23a) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS23a) $(CLIB)
	mv $@ $(XDIR)/$@

OBJS24=	2dcorr.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dcorr:	$(OBJS24) Makefile
	$(CLINK) $(CFLAGS) -o 2dcorr $(OBJS24) $(CLIB)
	cp 2dcorr  $(XDIR)/$@

OBJS27=	cyl-2dcorr.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
cyl-2dcorr:	$(OBJS27) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS27) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS26=	los_pvd.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
los_pvd:	$(OBJS26) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS26) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS25=	check-distribution.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
check-distribution:	$(OBJS25) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS25) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS1=	model.o spline.o splint.o nrutil.o
model:	$(OBJS1) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1) $(CLIB)
	cp model $(HOME)/cosmo/exec/

OBJS2=	wrapper_2dcorr.o
wrapper_2dcorr:	$(OBJS2)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS2) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS3=	avg-2dcorr.o
avg-2dcorr:	$(OBJS3)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS3) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/$@

OBJS4=	model2.o spline.o splint.o nrutil.o
model2:	$(OBJS4) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS4) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS5=	model3.o spline.o splint.o nrutil.o
model3:	$(OBJS5) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS5) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS6=	model4.o spline.o splint.o nrutil.o qromo.o trapzd.o
model4:	$(OBJS6) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS6) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS7=	check-2dgauss.o velocity-profiles.o nrutil.o linklist.o rfind.o \
	rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
check-2dgauss:	$(OBJS7) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS7) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS77=	check-2dxgauss.o velocity-profiles.o nrutil.o linklist.o rfind.o \
	rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
check-2dxgauss:	$(OBJS77) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS77) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS8=	los_pdf.o spline.o splint.o nrutil.o
los_pdf:	$(OBJS8) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS8) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS9=	model5.o spline.o splint.o nrutil.o velocity_functions.o qromo.o trapzd.o integrate_los.o
model5:	$(OBJS9) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS9) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS10=	model6.o spline.o splint.o nrutil.o velocity_functions.o trapzd.o
model6:	$(OBJS10) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS10) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS11=	tan-histograms.o tan-velocity-profiles.o nrutil.o linklist.o rfind.o \
	rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
tan-histograms:	$(OBJS11) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS11) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS12=	2dx-corr.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dx-corr:	$(OBJS12) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS12) $(CLIB)
	cp $@  $(HOME)/cosmo/exec/

OBJS13=	wrapper_2dx-corr.o
wrapper_2dx-corr:	$(OBJS13)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS13) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS14=	avg-2dx-corr.o
avg-2dx-corr:	$(OBJS14)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS14) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS15=	model5x.o spline.o splint.o nrutil.o velocity_functions.o qromo.o trapzd.o
model5x:	$(OBJS15) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS15) $(CLIB)
	cp $@ $(HOME)/cosmo/exec/

OBJS16=	angular-2dcorr.o linklist3d.o rfind_3d.o rfind2_3d.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o
angular-2dcorr:	$(OBJS16) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS16) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS17=	avg-angular.o
avg-angular:	$(OBJS17)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS17) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS18=	xi-multipoles.o nrutil.o qtrap.o spline.o splint.o trapzd.o
xi-multipoles:	$(OBJS18)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS18) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS19= small_scale_measure.o nrutil.o splint.o spline.o least_squares.o
small_scale_measure:	$(OBJS19)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS19) $(CLIB)
	cp $@ $(XDIR)/$@

OBJS20=	2dcorr_jacknife.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dcorr_jacknife:	$(OBJS20) Makefile
	$(CLINK) $(CFLAGS) -o 2dcorr $(OBJS20) $(CLIB)
	cp 2dcorr  $(XDIR)/$@

OBJS21=	angular-2dcorr-jacknife.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o
angular-2dcorr-jacknife:	$(OBJS21) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS21) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS22=	wp.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
wp:	$(OBJS22) Makefile
	$(CLINK) $(CFLAGS) -o 2dcorr $(OBJS22) $(CLIB)
	cp 2dcorr  $(XDIR)/$@

OBJS124=	2dcorr_covar.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dcorr_covar:	$(OBJS124) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS124) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS125=	angular-2dcorr-covar.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-covar:	$(OBJS125) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS125) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS125a=	angular-2dcorr-covar-xcorr.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-covar-xcorr:	$(OBJS125a) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS125a) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS126=	wp_covar.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
wp_covar:	$(OBJS126) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS126) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS126b=	cluster_cic.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
cluster_cic:	$(OBJS126b) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS126b) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS126a=	wp_covarboot.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
wp_covarboot:	$(OBJS126a) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS126a) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS127=	wp_covar_random.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o getpos_random.o
wp_covar_random:	$(OBJS127) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS127) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS127=	wp_dispersion.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
wp_dispersion:	$(OBJS127) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS127) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS126c=	angular-2dcorr-SDSS.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-SDSS:	$(OBJS126c) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS126c) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS127= 2dcorr_mass.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dcorr_mass:	$(OBJS127) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS127) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS128= 2dcorr_2halo.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dcorr_2halo:	$(OBJS128) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS128) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS129= 2dcorr_1halo.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
2dcorr_1halo:	$(OBJS129) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS129) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS130=	wpXcorr.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
wpXcorr:	$(OBJS130) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS130) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS1301=	delta_sigma.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
delta_sigma:	$(OBJS1301) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1301) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS131= create_linear_plot.o spline.o splint.o splie2.o splin2.o nrutil.o
create_linear_plot:	$(OBJS131) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS131) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS132=	wp_test.o linklist.o rfind.o rfind2.o getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o
wp_test:	$(OBJS132) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS132) $(CLIB)
	cp $@  $(XDIR)/$@


OBJS133=	angular-2dcorr-covar-1halo.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-covar-1halo:	$(OBJS133) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS133) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS134=	angular-2dcorr-covar-2halo.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-covar-2halo:	$(OBJS134) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS134) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS135=	angular-2dcorr-linear.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-linear:	$(OBJS135) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS135) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS136=	angular-2dcorr-linear-mu.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-linear-mu:	$(OBJS136) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS136) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS137=	angular-2dcorr-linear-mu-old.o linklist3d.o rfind_3d.o rfind2_3d.o \
	getpos_ascii.o getpos_tipsy.o   \
	getpos_fastfood.o random_selection.o brute_force.o redshift_distortions.o \
	ran1.o gasdev.o i3tensor_2.o spline.o splint.o qtrap.o trapzd.o
angular-2dcorr-linear-mu-old:	$(OBJS137) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS137) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS138= cluster_assignment.o gasdev.o ran1.o
cluster_assignment:	$(OBJS138) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS138) $(CLIB)
	cp $@  $(XDIR)/$@

OBJS139= cluster_MN.o gasdev.o ran1.o
cluster_MN:	$(OBJS139) Makefile
	$(CLINK) $(CFLAGS) -o $@ $(OBJS139) $(CLIB)
	cp $@  $(XDIR)/$@


clean:
	rm -f *.o
