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

mip_values = YAML::load(File.open("mip-value.yaml"))

n = 0

times = []
mir_gaps = []
greedy_gaps = []
i=1

ARGV.each do |f|
    yaml = YAML::load(File.open(f))

    instance = yaml['input-filename'].gsub("instances/", "").gsub(".pre.mps.gz", "")
    time_greedy = yaml['user-cpu-time'][2]
    times.push(time_greedy)

    mip_value = mip_values[instance]
    lp_value = yaml['obj-value'][0]
    mir_value = yaml['obj-value'][1]
    greedy_value = yaml['obj-value'][2]

    next if lp_value == mip_value

    mir_gap = (mir_value - lp_value) / (mip_value - lp_value)
    greedy_gap = (greedy_value - lp_value) / (mip_value - lp_value)

    mir_ncuts = yaml['added-cuts'][1]
    greedy_ncuts = yaml['added-cuts'][2]

    mir_gaps.push(mir_gap)
    greedy_gaps.push(greedy_gap)

    print("%20s " % instance)
    print("%14.2f " % lp_value)
    print("%14.2f " % mip_value)
    print("| %14.2f " % mir_value)
    print("%14.2f " % greedy_value)
    print("| %5.1f%% " % (mir_gap * 100))
    print("%5.1f%% " % (greedy_gap * 100))
    print("| %8d " % mir_ncuts)
    print("%8d " % greedy_ncuts)
    print("| %4d:%02d " % [yaml['user-cpu-time'][2] / 60, yaml['user-cpu-time'][2] % 60])
    print("\n")

    i += 1
end

greedy_gap = mean(greedy_gaps) * 100
mir_gap = mean(mir_gaps) * 100

time_greedy = gmean(times, 10)

print("\n%85s" % "")
print("%5.1f%% " % mir_gap)
print("%5.1f%% " % greedy_gap)
print("%21s %4d:%02d" % ["", time_greedy / 60, time_greedy % 60])

print("\n%92s" % "")
print("%5.1f%% " % ((greedy_gap / mir_gap - 1) * 100))
print("\n")
