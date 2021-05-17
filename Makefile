
#cidArr = b2kidston b12unifit1 b13resi
#zyorg = $(patsubst %, ensemble-wise/zyorg_%.nc, $(cidArr))
#zyorg = $(patsubst %,ensemble-wise/zyorg_%.nc,$(1))
zyorg = ensemble-wise/zyorg_$(cid).nc
cy300 = ensemble-wise/cy.cospectrum.$(cid).300.nc
lininst = ensemble-wise/lininst_$(cid).nc
zyorg2 = fms-output/$(cid)ens$(ii)/zyorg_$(cid)ens$(ii).nc

 # makefile snakefile n bren
 # makefile software dependency? https://medium.com/@ionrock/makefile-madness-and-software-dependencies-7c38962793c9
 # makefile force replot

plot: fig1 fig2 fig3 fig4
.PHONY: plot
fig%: fig%.pdf
	@echo $< is ready
fig2: fig3and2.pdf
	@echo $< is ready
fig3: fig3and2.pdf
	@echo $< is ready
fig1.pdf: yp.diff.panel.pdf
	cp -a $< $@
fig3and2.pdf: cy.cospectrum.panel.300.pdf
	cp -a $< $@
fig4.pdf: yl.lengthscale.pub300.pdf
	cp -a $< $@
yp.diff.panel.pdf: yp.diff.panel.ncl gsn_csm.ncl $(foreach cid,b1ctrl b2kidston b12unifit1 b13resi,$(zyorg))
	ncl -Q $<
cy.cospectrum.panel.300.pdf: cy.cospectrum.panel.ncl gsn_csm.ncl $(foreach cid,b2kidston b12unifit1 b13resi,$(cy300))
	ncl -Q $<
yl.lengthscale.pub300.pdf: yl.lengthscale.ncl gsn_csm.ncl $(foreach cid,b1ctrl,$(zyorg) $(lininst)) $(foreach cid,b12unifit1 b15unifp15 b16unifm15 b14unifm3 b2kidston b13resi,$(cy300) $(zyorg) $(lininst))
	ncl -Q $<
yp.uncertainty.pdf: yp.uncertainty.ncl gsn_csm.ncl $(foreach cid,b1ctrl b12unifit1,$(foreach ii,$(shell seq 20),$(zyorg2)))
	@echo prerequisites not complete
#	ncl -Q $<
#ensemble-wise/zyorg_%.nc:
#	@echo recipe not yet ready $@, use crontab-ensemble-fms-monitor.sh and nces
#ensemble-wise/cy.cospectrum.%.300.nc: cy.cospectrum.ncl
#	@echo rule not yet ready $@, use crontab-ensemble-fms-monitor.sh
#	ncl -Q $<
ensemble-wise/lininst_%.nc: lininst_pchan.m lininst.m $(foreach cid,%,$(zyorg))
	export CASENAME=$*; \matlab -nodesktop -nodisplay -nosplash -batch "lininst_pchan"
gridfrac.nc: | gridfrac.m
	@echo recipe not yet ready $@
%.ncl: ${projroot}/%.ncl
	rsync -a $< $(@D)/
#ensemble-wise/%.nc: ${projroot}/ensemble-wise/%.nc
#	rsync -a $< $(@D)/
%.nc:
	mkdir -p $(@D)/
	rsync -a ${projroot}/$*.nc $(@D)/
#test: $(call zyorg,b1ctrl b2kidston b12unifit1 b13resi)
#	@echo $^
#	@echo $(call cy300,b2kidston b12unifit1 b13resi)
#	@echo $(call zyorg,b1ctrl b2kidston b12unifit1 b13resi)
#	@echo $(patsubst %,$(zyorg) $(lininst),b1ctrl b2kidston b12unifit1 b13resi)
#	echo $(foreach cid,b1ctrl b2kidston b12unifit1 b13resi,ensemble-wise/zyorg_$(cid).nc ensemble-wise/lininst_$(cid).nc)
#	echo $(zyorg)
#	@echo $(foreach cid,b1ctrl b2kidston b12unifit1 b13resi,$(zyorg) $(lininst))
#	@echo $(foreach cid,b1ctrl,$(zyorg) $(lininst)) $(foreach cid,b12unifit1 b15unifp15 b16unifm15 b14unifm3,$(cy300) $(zyorg) $(lininst))
ensemble-wise/test_%.nc: lininst_pchan.m lininst.m $(foreach cid,%,$(zyorg))
	@echo prerequisite list $^
#test:
#	@echo $(foreach cid,b1ctrl b12unifit1,$(foreach ii,$(shell seq 20),$(zyorg2)))

