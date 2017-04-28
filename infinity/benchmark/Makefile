miplib3-small := \
	instances/bell5.pre.done \
	instances/blend2.pre.done \
	instances/dcmulti.pre.done \
	instances/egout.pre.done \
	instances/flugpl.pre.done \
	instances/gen.pre.done \
	instances/gt2.pre.done \
	instances/lseu.pre.done \
	instances/mas74.pre.done \
	instances/modglob.pre.done \
	instances/p0033.pre.done \
	instances/p0548.pre.done \
	instances/pp08a.pre.done \
	instances/pp08aCUTS.pre.done \
	instances/set1ch.pre.done \
	instances/vpm1.pre.done \
	instances/vpm2.pre.done \

miplib3-all-instances := \
	instances/air03.pre.done \
	instances/bell3a.pre.done \
	instances/bell5.pre.done \
	instances/blend2.pre.done \
	instances/cap6000.pre.done \
	instances/danoint.pre.done \
	instances/dcmulti.pre.done \
	instances/egout.pre.done \
	instances/fiber.pre.done \
	instances/fixnet6.pre.done \
	instances/flugpl.pre.done \
	instances/gen.pre.done \
	instances/gesa2_o.pre.done \
	instances/gesa2.pre.done \
	instances/gesa3_o.pre.done \
	instances/gt2.pre.done \
	instances/harp2.pre.done \
	instances/khb05250.pre.done \
	instances/l152lav.pre.done \
	instances/lseu.pre.done \
	instances/markshare1.pre.done \
	instances/markshare2.pre.done \
	instances/mas74.pre.done \
	instances/misc03.pre.done \
	instances/misc06.pre.done \
	instances/misc07.pre.done \
	instances/mitre.pre.done \
	instances/mkc.pre.done \
	instances/mod008.pre.done \
	instances/mod010.pre.done \
	instances/mod011.pre.done \
	instances/modglob.pre.done \
	instances/nw04.pre.done \
	instances/p0033.pre.done \
	instances/p0201.pre.done \
	instances/p0282.pre.done \
	instances/p2756.pre.done \
	instances/pk1.pre.done \
	instances/pp08aCUTS.pre.done \
	instances/pp08a.pre.done \
	instances/qiu.pre.done \
	instances/qnet1_o.pre.done \
	instances/rentacar.pre.done \
	instances/rgn.pre.done \
	instances/rout.pre.done \
	instances/set1ch.pre.done \
	instances/stein27.pre.done \
	instances/stein45.pre.done \
	instances/swath.pre.done \
	instances/vpm1.pre.done \
	instances/vpm2.pre.done \
	instances/10teams.pre.done \
	instances/gesa3.pre.done \
	instances/p0548.pre.done \
	instances/qnet1.pre.done \
	instances/seymour.pre.done
#	instances/fast0507.pre.done \
#	instances/air04.pre.done \
#	instances/air05.pre.done \

mir := $(patsubst instances/%, mir/%, $(patsubst %.done,%.yaml,$(miplib3-all-instances)))
bases := $(patsubst instances/%, bases/%, $(patsubst %.done,%.bas,$(miplib3-all-instances)))
solutions := $(patsubst instances/%, solutions/%, $(patsubst %.done,%.x,$(miplib3-all-instances)))

2row-cont := $(patsubst instances/%, 2row-cont/%, $(patsubst %.done,%.yaml, $(miplib3-all-instances)))
2row-lift := $(patsubst instances/%, 2row-lift/%, $(patsubst %.done,%.yaml, $(miplib3-all-instances)))
3row-cont := $(patsubst instances/%, 3row-cont/%, $(patsubst %.done,%.yaml, $(miplib3-all-instances)))
3row-lift := $(patsubst instances/%, 3row-lift/%, $(patsubst %.done,%.yaml, $(miplib3-all-instances)))
5row-cont := $(patsubst instances/%, 5row-cont/%, $(patsubst %.done,%.yaml, $(miplib3-all-instances)))
5row-lift := $(patsubst instances/%, 5row-lift/%, $(patsubst %.done,%.yaml, $(miplib3-all-instances)))

.SECONDARY:

multirow := ../../build/infinity/benchmark/infinity-benchmark.run

miplib3-all-instances: $(miplib3-all-instances)
miplib3-small: $(miplib3-small)

miplib3-2row-cont: $(2row-cont)
miplib3-2row-lift: $(2row-lift)
miplib3-3row-cont: $(3row-cont)
miplib3-3row-lift: $(3row-lift)
miplib3-5row-cont: $(5row-cont)
miplib3-5row-lift: $(5row-lift)

bases: $(bases)

solutions: $(solutions)

mir: $(mir)

instances/%.done: $(multirow) 2row-cont/%.yaml 2row-lift/%.yaml 3row-cont/%.yaml
	@touch $@

2row-cont/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas
	@echo "   2row-cont  $*..."
	@$(multirow) --mir --greedy --problem instances/$*.mps.gz --basis bases/$*.bas --log 2row-cont/$*.log --stats $@ --solution solutions/$*.x || touch 2row-cont/$*.err.$$?
	@echo "   2row-cont  $* [done]"

2row-lift/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas
	@echo "   2row-lift  $*..."
	@$(multirow) --mir --greedy --problem instances/$*.mps.gz --basis bases/$*.bas --log 2row-lift/$*.log --stats $@ --solution solutions/$*.x --lift || touch 2row-lift/$*.err.$$?
	@echo "   2row-lift  $* [done]"

3row-cont/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas
	@echo "   3row-cont  $*..."
	@$(multirow) --mir --greedy --problem instances/$*.mps.gz --basis bases/$*.bas --log 3row-cont/$*.log --stats $@ --solution solutions/$*.x --rows 3 || touch 3row-cont/$*.err.$$?
	@echo "   3row-cont  $* [done]"

3row-lift/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas
	@echo "   3row-lift  $*..."
	@$(multirow) --mir --greedy --problem instances/$*.mps.gz --basis bases/$*.bas --log 3row-lift/$*.log --stats $@ --solution solutions/$*.x --rows 3 --lift || touch 3row-lift/$*.err.$$?
	@echo "   3row-lift  $* [done]"

5row-cont/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas
	@echo "   5row-cont  $*..."
	@$(multirow) --mir --greedy --problem instances/$*.mps.gz --basis bases/$*.bas --log 5row-cont/$*.log --stats $@ --solution solutions/$*.x --rows 5 || touch 5row-cont/$*.err.$$?
	@echo "   5row-cont  $* [done]"

5row-lift/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas
	@echo "   5row-lift  $*..."
	@$(multirow) --mir --greedy --problem instances/$*.mps.gz --basis bases/$*.bas --log 5row-lift/$*.log --stats $@ --solution solutions/$*.x --rows 5 --lift || touch 5row-lift/$*.err.$$?
	@echo "   5row-lift  $* [done]"


mir/%.yaml: $(multirow) instances/%.mps.gz bases/%.bas solutions/%.x
	@echo "      mir  $*"
	@$(multirow) --mir --problem instances/$*.mps.gz --basis bases/$*.bas --log mir/$*.log --stats $@ --solution solutions/$*.x  > /dev/null

bases/%.bas:
	@echo "     base  $*"
	@$(multirow) --problem instances/$*.mps.gz --write-basis $@ --log $@.log > /dev/null

solutions/%.x:
	@echo " solution  $*"
	@$(multirow) --keep-integrality --problem instances/$*.mps.gz --write-solution $@ --log $@.log > /dev/null

clean:
	rm -fv [2-9]row-lift/*
	rm -fv [2-9]row-cont/*
	rm -fv instances/*.done

table:
	scripts/table.rb 2row-cont/ pi-2 > tables/2row-cont.txt
	scripts/table.rb 2row-lift/ pi-2 > tables/2row-lift.txt
	scripts/table.rb 3row-cont/ pi-3 > tables/3row-cont.txt
	scripts/table.rb 3row-lift/ pi-3 > tables/3row-lift.txt
	scripts/table.rb 5row-cont/ pi-3 > tables/5row-cont.txt
	scripts/table.rb 5row-lift/ pi-3 > tables/5row-lift.txt
	scripts/table-csv.rb 2row-cont/ pi-2 > csv/2row-cont.csv
	scripts/table-csv.rb 3row-cont/ pi-3 > csv/3row-cont.csv
	scripts/table-csv.rb 2row-lift/ pi-2 > csv/2row-lift.csv
	scripts/table-csv.rb 3row-lift/ pi-3 > csv/3row-lift.csv
