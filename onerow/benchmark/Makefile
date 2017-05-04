miplib2010-top := \
	instances/gmu-35-40.pre.done \
	instances/eil33-2.pre.done \
	instances/neos-1337307.pre.done \
	instances/opm2-z7-s2.pre.done \
	instances/mik-250-1-100-1.pre.done \
	instances/neos-686190.pre.done \
	instances/mine-90-10.pre.done \
	instances/cov1075.pre.done \
	instances/mine-166-5.pre.done \
	instances/n3div36.pre.done

miplib2010 := \
	instances/30n20b8.pre.done \
	instances/aflow40b.pre.done \
	instances/air04.pre.done \
	instances/app1-2.pre.done \
	instances/bab5.pre.done \
	instances/beasleyC3.pre.done \
	instances/biella1.pre.done \
	instances/bienst2.pre.done \
	instances/binkar10_1.pre.done \
	instances/bley_xl1.pre.done \
	instances/core2536-691.pre.done \
	instances/cov1075.pre.done \
	instances/csched010.pre.done \
	instances/danoint.pre.done \
	instances/dfn-gwin-UUM.pre.done \
	instances/eil33-2.pre.done \
	instances/eilB101.pre.done \
	instances/enlight13.pre.done \
	instances/ex9.pre.done \
	instances/glass4.pre.done \
	instances/gmu-35-40.pre.done \
	instances/iis-100-0-cov.pre.done \
	instances/iis-bupa-cov.pre.done \
	instances/iis-pima-cov.pre.done \
	instances/lectsched-4-obj.pre.done \
	instances/macrophage.pre.done \
	instances/map18.pre.done \
	instances/map20.pre.done \
	instances/mcsched.pre.done \
	instances/mik-250-1-100-1.pre.done \
	instances/mine-166-5.pre.done \
	instances/mine-90-10.pre.done \
	instances/msc98-ip.pre.done \
	instances/mspp16.pre.done \
	instances/mzzv11.pre.done \
	instances/n3div36.pre.done \
	instances/n3seq24.pre.done \
	instances/n4-3.pre.done \
	instances/neos-1109824.pre.done \
	instances/neos-1337307.pre.done \
	instances/neos-1396125.pre.done \
	instances/neos13.pre.done \
	instances/neos-1601936.pre.done \
	instances/neos18.pre.done \
	instances/neos-476283.pre.done \
	instances/neos-686190.pre.done \
	instances/neos-916792.pre.done \
	instances/neos-934278.pre.done \
	instances/net12.pre.done \
	instances/netdiversion.pre.done \
	instances/newdano.pre.done \
	instances/noswot.pre.done \
	instances/ns1208400.pre.done \
	instances/ns1688347.pre.done \
	instances/ns1758913.pre.done \
	instances/ns1830653.pre.done \
	instances/opm2-z7-s2.pre.done \
	instances/pg5_34.pre.done \
	instances/pigeon-10.pre.done \
	instances/pw-myciel4.pre.done \
	instances/qiu.pre.done \
	instances/rail507.pre.done \
	instances/ran16x16.pre.done \
	instances/reblock67.pre.done \
	instances/rmatr100-p10.pre.done \
	instances/rmatr100-p5.pre.done \
	instances/rmine6.pre.done \
	instances/rocII-4-11.pre.done \
	instances/rococoC10-001000.pre.done \
	instances/roll3000.pre.done \
	instances/sp98ic.pre.done \
	instances/sp98ir.pre.done \
	instances/tanglegram1.pre.done \
	instances/tanglegram2.pre.done \
	instances/timtab1.pre.done \
	instances/triptim1.pre.done \
	instances/unitcal_7.pre.done \
	instances/vpphard.pre.done \
	instances/zib54-UUE.pre.done

mir := $(patsubst instances/%, mir/%, $(patsubst %.done,%.yaml,$(miplib3-all-instances)))
bases := $(patsubst instances/%, bases/%, $(patsubst %.done,%.bas,$(miplib3-all-instances)))
solutions := $(patsubst instances/%, solutions/%, $(patsubst %.done,%.x,$(miplib3-all-instances)))

benchmark := ../../build/onerow/benchmark/onerow-benchmark.run

.SECONDARY:

miplib2010: $(miplib2010)
miplib2010-top: $(miplib2010-top)

instances/%.done: $(benchmark) out/%.yaml
	@touch $@

out/%.yaml: $(benchmark) instances/%.mps.gz
	@$(benchmark) -mw -f instances/$*.mps.gz -b bases/$*.bas -s out/$*.yaml 2>&1 | tee out/$*.log

clean:
	rm -f out/* instances/*.done
