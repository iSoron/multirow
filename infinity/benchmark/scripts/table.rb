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

print("%4s " % "")
print("%16s " % "instance")
print("%9s " % "g-mir")
print("%9s " % "g-infty")
print("%9s " % "g-exact")
print("%9s " % "improv")
print("%9s " % "t-total")
print("%9s " % "t-per-cut")
print("%6s " % "n-cuts")
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

    print("%4d " % count)
    print("%16s " % instance)
    print("%8.2f%% " % (gap_mir))
    print("%8.2f%% " % (gap_greedy))
    print("%8.2f%% " % (gap_pi))

    if (gap_pi - gap_mir).abs > 0.01
        print("%8.2f%% " % (improv_greedy))
    else
        print("%8s  " % "---")
    end

    print("%6d:%02d " % [yaml['user-cpu-time'][2] / 60, yaml['user-cpu-time'][2] % 60])
    print("%6d ms " % (time_percut))
    print("%6d " % (yaml['added-cuts'][2]))
    print("\n")

    count += 1
end

mean_exact = mean(exact_gaps)
mean_greedy = mean(greedy_gaps)
mean_mir = mean(mir_gaps)

#ncuts_exact = gmean(exact_cuts, 10)
#ncuts_greedy = gmean(greedy_cuts, 10)
#ncuts_mir = gmean(mir_cuts, 10)
#
#time_exact = gmean(exact_times, 5)
#time_greedy = gmean(greedy_times, 5)
#time_percut = gmean(percut_times, 5)
#
print("\n%22s" % "")
print("%8.2f%% " % mean_mir)
print("%8.2f%% " % mean_greedy)
print("%8.2f%% " % mean_exact)
#print("%8.1f%% " % mean(improv))
#print("  %8d " % ncuts_mir)
#print("%8d " % ncuts_greedy)
#print("%8d   " % ncuts_exact)
#print("%4d:%02d " % [time_greedy / 60, time_greedy % 60])
#print("%4d:%02d " % [time_exact / 60, time_exact % 60])
#print("  %4.0f ms" % time_percut)
print("\n")

#print("%30s" % "")
#print("%5.1f%% " % ((mean_greedy / mean_mir - 1) * 100))
#print("%5.1f%% " % ((mean_exact / mean_mir - 1) * 100))
#print("\n")

print("%34s" % "")
print("%6.2f%% " % ((mean_greedy - mean_mir) / (mean_exact - mean_mir) * 100))
print("\n")

if missing.size > 0
    print("Missing data:\n")
    missing.each do |m|
        print("  %s\n" % m)
    end
end

#if mir_differs.size > 0
#    print("MIR differs:\n")
#    mir_differs.each do |m|
#        print("%4d %16s %6.2f %6.2f\n" % [count, m[0], m[1], m[2]])
#        count += 1
#    end
#    print("\n")
#end

print("Disabled:\n")
Dir["instances/off/*"].each do |f|
    print("%4d %16s\n" % [count, f.gsub("instances/off/", "").gsub(".pre.mps.gz", "")])
    count += 1
end
print("\n")


print("Error:\n")
Dir[ARGV[0] + "*.err*"].sort.each do |f|
    instance = f.gsub("greedy/", "").gsub(".pre.err.", " ").split(" ")
    print("%4d %-30s %4d\n" % [count, instance[0], instance[1]])
    count += 1
end
