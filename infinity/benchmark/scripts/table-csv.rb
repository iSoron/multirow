#!/usr/bin/env ruby
require 'yaml'

def gmean(x, shift)
    sum=0.0
    x.each {|v| sum += Math.log(v+shift)}
    sum /= x.size
    Math.exp(sum)-shift
end

def mean(x)
    sum=0.0
    x.each {|v| sum += v}
    sum /= x.size
end

data = YAML::load(File.open("data/data.yaml"))

n = 0

exact_cuts = []
exact_gaps = []
exact_times = []
greedy_cuts = []
greedy_gaps = []
greedy_times = []
mir_cuts = []
mir_gaps = []
percut_times = []
improv = []

missing = []
no_improvement = []
mir_differs = []
count = 1

print("instance,")
print("gmir,")
print("ginfty,")
print("gexact,")
print("improv,")
print("ttotal,")
print("tpercut,")
print("ncuts")
print("\n")


Dir[ARGV[0] + "*yaml"].sort.each do |f|
    yaml = YAML::load(File.open(f))

    instance = yaml['input-filename'].gsub("instances/", "").gsub(".pre.mps.gz", "")
    d = data[instance]

    obj_lp = yaml['obj-value'][0]
    obj_ip = d['ip']
    gap_pi = d[ARGV[1]]
    obj_mir = yaml['obj-value'][1]
    obj_greedy = yaml['obj-value'][2]

    gap_mir = (obj_mir - obj_lp) / (obj_ip - obj_lp) * 100
    gap_greedy = (obj_greedy - obj_lp) / (obj_ip - obj_lp) * 100

    #gap_pi = [gap_pi, gap_greedy].max

    improv_greedy = (gap_greedy - gap_mir) / (gap_pi - gap_mir) * 100
    improv_greedy = 0 if improv_greedy.abs < 0.01

    if (gap_mir - d['gmi']).abs > 1
        mir_differs.push([instance, gap_mir, d['gmi']])
    end

    begin
        time_percut = yaml['time-per-cut'][2] * 1000
    rescue NoMethodError
        time_percut = 0.0
    end
    percut_times.push(time_percut)

    mir_gaps.push(gap_mir)
    greedy_gaps.push(gap_greedy)
    exact_gaps.push(gap_pi)
    improv.push(improv_greedy)

    print("%s," % instance.gsub("_", "\\_"))
    print("%.1f," % (gap_mir))
    print("%.1f," % (gap_greedy))
    print("%.1f," % (gap_pi))

    if (gap_pi - gap_mir).abs > 0.01
        print("%.1f," % (improv_greedy))
    else
        print("%s," % "---")
    end

    print("%d:%02d," % [yaml['user-cpu-time'][2] / 60, yaml['user-cpu-time'][2] % 60])
    print("%d ms," % (time_percut))
    print("%d" % (yaml['added-cuts'][2]))
    print("\n")

    count += 1
end

mean_exact = mean(exact_gaps)
mean_greedy = mean(greedy_gaps)
mean_mir = mean(mir_gaps)

#print("\n%22s" % "")
#print("%8.2f%% " % mean_mir)
#print("%8.2f%% " % mean_greedy)
#print("%8.2f%% " % mean_exact)
